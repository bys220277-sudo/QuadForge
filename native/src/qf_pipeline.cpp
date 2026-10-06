#include "qf_pipeline.h"

#include "qf_extract.h"
#include "qf_features.h"
#include "qf_field.h"
#include "qf_param.h"
#include "qf_remesh.h"
#include "qf_quadpost.h"

#include <queue>
#include <unordered_map>

namespace qf {

static void store(Result& res, const char* name, const std::vector<double>& v) { res.darr[name] = v; }
static void store(Result& res, const char* name, const std::vector<int>& v) { res.iarr[name] = v; }
template <class T>
static std::vector<int> as_int(const std::vector<T>& v) { return std::vector<int>(v.begin(), v.end()); }

static void log_index_stats(Ctx& ctx, const char* label, const std::vector<int>& idx) {
    int n1 = 0, hi = 0;
    for (int v : idx) { if (v != 0) ++n1; if (std::abs(v) >= 2) ++hi; }
    ctx.logf("поле (%s): сингулярных вершин %d, с |индексом|>=2: %d", label, n1, hi);
}


// Ограничения поля по фичерям: жёсткие (soft_w<=0) или мягкие с весом soft_w.
static void build_field_input(const TriMesh& M, const std::vector<uint8_t>& feat, double soft_w, FieldInput& fi,
                              std::vector<int>* feat_faces) {
    const int nf = M.nF();
    std::vector<C2> fixed;
    field_constraints_from_features(M, feat, fixed);
    if (feat_faces) {
        feat_faces->assign(nf, 0);
        for (int f = 0; f < nf; ++f) (*feat_faces)[f] = fixed[f].abs2() > 0.5 ? 1 : 0;
    }
    if (soft_w > 0) {
        fi.fixed.clear();
        fi.soft_z.assign(nf, C2(1, 0));
        fi.soft_w.assign(nf, 0.0);
        for (int f = 0; f < nf; ++f)
            if (fixed[f].abs2() > 0.5) { fi.soft_z[f] = fixed[f]; fi.soft_w[f] = soft_w; }
    } else {
        fi.fixed = fixed;
    }
}

// Мягкий приор: выравнивание по главным кривизнам там, где анизотропия кривизны надёжна.
static void add_curvature_prior(const TriMesh& M, FieldInput& fi, double strength, double L, const QFParams& p) {
    const int nF = M.nF();
    CurvatureData cd;
    compute_face_curvature(M, cd, p.i_reserved[2] > 0 ? p.i_reserved[2] : 4);
    if (fi.soft_w.empty()) { fi.soft_z.assign(nF, C2(1, 0)); fi.soft_w.assign(nF, 0.0); }
    const double gate0 = p.d_reserved[1] > 0 ? p.d_reserved[1] : 0.3;
    for (int f = 0; f < nF; ++f) {
        if (fi.soft_w[f] > 0 || (!fi.fixed.empty() && fi.fixed[f].abs2() > 0.5)) continue;
        double a = std::fabs(cd.k1[f] - cd.k2[f]);
        double sm = std::fabs(cd.k1[f]) + std::fabs(cd.k2[f]) + 1e-12;
        double aniso = a / sm;
        double gate = std::min(1.0, std::max(0.0, (a * L - gate0) / (2.0 * gate0)));
        fi.soft_z[f] = polar(4.0 * cd.theta[f]);
        fi.soft_w[f] = 2.0 * strength * aniso * aniso * gate;
    }
}

// Поле размеров (длина ребра квада) на грань исходника: адаптивность по кривизне + карта плотности.
// Детализация поверхности в масштабе квада: угловой разброс нормалей (рад) в окрестности
// радиуса ~r вокруг каждой грани. Складки, морщины, мелкие формы и тонкие части дают большой
// разброс, гладкие участки — малый. Не зависит от плотности триангуляции исходника.
static std::vector<double> surface_detail(const TriMesh& S, double r) {
    // Ячейки сетки размера r; в каждой — 6 корзин по преобладающему направлению нормали,
    // чтобы противоположные стороны тонких частей (уши, пальцы, пластины) не смешивались.
    const int nf = S.nF();
    struct Bin { V3 n, c; double a = 0; };
    struct Cell { Bin b[6]; };
    std::unordered_map<u64, Cell> grid;
    grid.reserve(nf / 4 + 16);
    auto key = [&](int x, int y, int z) {
        return ((u64)(uint32_t)(x + (1 << 20)) << 42) | ((u64)(uint32_t)(y + (1 << 20)) << 21) | (u64)(uint32_t)(z + (1 << 20));
    };
    auto dir_bin = [](const V3& n) {
        const double ax = std::fabs(n.x), ay = std::fabs(n.y), az = std::fabs(n.z);
        if (ax >= ay && ax >= az) return n.x >= 0 ? 0 : 1;
        if (ay >= az) return n.y >= 0 ? 2 : 3;
        return n.z >= 0 ? 4 : 5;
    };
    std::vector<V3> cen(nf);
    std::vector<int> cx(nf), cy(nf), cz(nf);
    for (int f = 0; f < nf; ++f) {
        cen[f] = (S.P[S.F[3 * f]] + S.P[S.F[3 * f + 1]] + S.P[S.F[3 * f + 2]]) * (1.0 / 3.0);
        cx[f] = (int)std::floor(cen[f].x / r);
        cy[f] = (int)std::floor(cen[f].y / r);
        cz[f] = (int)std::floor(cen[f].z / r);
        Bin& c = grid[key(cx[f], cy[f], cz[f])].b[dir_bin(S.FN[f])];
        const double a = S.FA[f];
        c.n += S.FN[f] * a;
        c.c += cen[f] * a;
        c.a += a;
    }
    for (auto& kv : grid)
        for (Bin& c : kv.second.b)
            if (c.a > 0) c.c = c.c / c.a;
    std::vector<double> det(nf, 0.0);
    const double inv2s2 = 1.0 / (2.0 * (0.6 * r) * (0.6 * r));
    parallel_each(nf, [&](i64 f) {
        const V3& nf0 = S.FN[f];
        V3 m;
        double w = 0;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = grid.find(key(cx[f] + dx, cy[f] + dy, cz[f] + dz));
                    if (it == grid.end()) continue;
                    for (const Bin& c : it->second.b) {
                        if (c.a <= 0) continue;
                        // другая сторона тонкой части: нормаль корзины смотрит против нормали грани
                        if (dot(c.n, nf0) < -0.2 * c.a) continue;
                        const double g = std::exp(-norm2(c.c - cen[f]) * inv2s2);
                        m += c.n * g;
                        w += c.a * g;
                    }
                }
        const double len = w > 0 ? norm(m) / w : 1.0;
        det[f] = std::sqrt(std::max(0.0, 2.0 * (1.0 - std::min(1.0, len))));
    }, 1024);
    return det;
}

// Ограничение градиента размера: h(f) <= h(g) + gamma * |c_f - c_g| (Дейкстра по смежности граней).
static void limit_size_gradient(const TriMesh& S, std::vector<double>& h, double gamma) {
    const int nf = S.nF();
    std::vector<V3> cen(nf);
    for (int f = 0; f < nf; ++f) cen[f] = (S.P[S.F[3 * f]] + S.P[S.F[3 * f + 1]] + S.P[S.F[3 * f + 2]]) * (1.0 / 3.0);
    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    for (int f = 0; f < nf; ++f) pq.push({h[f], f});
    while (!pq.empty()) {
        auto [hv, f] = pq.top();
        pq.pop();
        if (hv > h[f]) continue;
        for (int i = 0; i < 3; ++i) {
            const int g = S.FF[3 * f + i];
            if (g < 0) continue;
            const double cand = hv + gamma * norm(cen[g] - cen[f]);
            if (cand < h[g]) { h[g] = cand; pq.push({cand, g}); }
        }
    }
}

static std::vector<double> compute_source_sizes(const TriMesh& S, const QFParams& p, const std::vector<double>& vdens,
                                                int target_count, double L, Ctx& ctx) {
    const int nf = S.nF();
    std::vector<double> h(nf, L);
    const bool use_adapt = p.adaptive > 1e-3;
    const bool use_dens = !vdens.empty() && p.density_strength > 1e-3;
    if (!use_adapt && !use_dens) return h;
    std::vector<double> lh(nf, 0.0);  // log-множитель
    if (use_adapt) {
        // детализация на двух масштабах: квад и полтора квада
        std::vector<double> d1 = surface_detail(S, 1.0 * L), d2 = surface_detail(S, 1.6 * L);
        const double gain = 6.0 * p.adaptive;
        for (int f = 0; f < nf; ++f) {
            const double d = 0.5 * (d1[f] + d2[f]);
            lh[f] += -std::log(1.0 + gain * d);
        }
    }
    if (use_dens) {
        for (int f = 0; f < nf; ++f) {
            double d = (vdens[S.F[3 * f]] + vdens[S.F[3 * f + 1]] + vdens[S.F[3 * f + 2]]) / 3.0;
            d = std::max(0.05, 1.0 + (d - 1.0) * p.density_strength);
            lh[f] += -0.5 * std::log(d);
        }
    }
    // предварительная нормировка к целевому числу и ограничение перепада размеров
    for (double& v : lh) v = std::max(std::log(0.15), std::min(std::log(4.0), v));
    double cnt = 0;
    for (int f = 0; f < nf; ++f) {
        h[f] = L * std::exp(lh[f]);
        cnt += S.FA[f] / (h[f] * h[f]);
    }
    if (target_count > 0 && cnt > 0) {
        const double sc = std::sqrt(cnt / target_count);
        for (double& v : h) v *= sc;
    }
    // плавный переход плотности: размер растёт не быстрее gamma на единицу длины
    const double gamma = p.d_reserved[11] > 0 ? p.d_reserved[11] : 0.3;
    limit_size_gradient(S, h, gamma);
    cnt = 0;
    for (int f = 0; f < nf; ++f) cnt += S.FA[f] / (h[f] * h[f]);
    if (target_count > 0 && cnt > 0) {
        const double sc = std::sqrt(cnt / target_count);
        for (double& v : h) v *= sc;
    }
    double mn = 1e300, mx = 0;
    for (double v : h) { mn = std::min(mn, v); mx = std::max(mx, v); }
    ctx.logf("поле размеров: от %.4g до %.4g (равномерное %.4g), перепад x%.2f", mn, mx, L, mx / std::max(mn, 1e-300));
    return h;
}

// Мягкое выравнивание поля вдоль направляющих (guides): грани в пределах ~L от отрезка.
static void add_guide_prior(const TriMesh& M, FieldInput& fi, const SegBvh& guides, double strength, double radius) {
    if (guides.bvh.empty() || strength <= 0) return;
    const int nF = M.nF();
    if (fi.soft_w.empty()) { fi.soft_z.assign(nF, C2(1, 0)); fi.soft_w.assign(nF, 0.0); }
    parallel_each(nF, [&](i64 f) {
        if (!fi.fixed.empty() && fi.fixed[f].abs2() > 0.5) return;
        V3 c = (M.P[M.F[3 * f]] + M.P[M.F[3 * f + 1]] + M.P[M.F[3 * f + 2]]) * (1.0 / 3.0), q;
        double d2;
        int sgi = guides.nearest(c, q, d2);
        if (sgi < 0) return;
        double d = std::sqrt(d2);
        if (d > radius) return;
        V3 dir = guides.B[sgi] - guides.A[sgi];
        const V3& n = M.FN[f];
        dir = dir - n * dot(n, dir);
        if (norm(dir) < 1e-12) return;
        double w = 6.0 * strength * (1.0 - d / radius);
        if (w > fi.soft_w[f]) {
            fi.soft_z[f] = polar(4.0 * angle_in_face(M, (int)f, dir));
            fi.soft_w[f] = w;
        }
    });
}

void run_pipeline(const QFParams& p, const QFInput& in, Ctx& ctx, Result& res) {
    const double t0 = now_seconds();
    ctx.progress(0.01f, 1, "Подготовка сетки");
    const int nv = in.nv, nf = in.nf;
    if (!in.V || !in.F || nv < 3 || nf < 1) {
        res.error = "пустая сетка";
        return;
    }
    std::vector<V3> P(nv);
    for (int i = 0; i < nv; ++i) P[i] = V3(in.V[3 * i], in.V[3 * i + 1], in.V[3 * i + 2]);
    std::vector<int> Fi(in.F, in.F + 3 * (size_t)nf);
    TriMesh S;  // исходная (очищенная) сетка
    std::vector<int> vert_origin;
    std::vector<int> weld_rep;
    {
        CleanResult cr = clean_triangle_soup(P, Fi);
        ctx.logf("очистка: сварено %d, T-стыков %d, вырожд. %d, дубликатов %d, неманифольдных %d, расщеплено вершин %d",
                 cr.welded_vertices, cr.tjunction_splits, cr.removed_degenerate, cr.removed_duplicate,
                 cr.removed_nonmanifold, cr.split_vertices);
        weld_rep = std::move(cr.weld_rep);
        S.P = std::move(P);
        S.F = std::move(Fi);
        vert_origin = std::move(cr.vert_origin);
    }
    std::string err;
    if (!S.build(&err)) {
        res.error = "Сетка непригодна: " + err;
        return;
    }
    const double area = S.total_area();
    const double L = p.target_edge_length > 0 ? p.target_edge_length
                                               : std::sqrt(area / std::max(1, p.target_faces));
    const int target_count = p.target_edge_length > 0 ? (int)(area / (p.target_edge_length * p.target_edge_length))
                                                      : p.target_faces;
    ctx.logf("исходник: %d вершин, %d граней, площадь %.4g, целевое ребро %.4g", S.nV(), S.nF(), area, L);

    FeatureOptions fo;
    fo.sharp_angle = p.sharp_angle;
    fo.detect_sharp = p.detect_sharp != 0;
    fo.boundary = p.preserve_boundary != 0;
    // пары острых рёбер пользователя — в нумерации после сварки
    std::vector<int> user_edges;
    if (in.sharp_edges && in.n_sharp > 0) {
        user_edges.assign(in.sharp_edges, in.sharp_edges + 2 * (size_t)in.n_sharp);
        for (auto& v : user_edges)
            if (v >= 0 && v < (int)weld_rep.size()) v = weld_rep[v];
    }
    fo.user_edges = user_edges.empty() ? nullptr : user_edges.data();
    fo.n_user = (int)(user_edges.size() / 2);
    fo.origin = &vert_origin;
    std::vector<uint8_t> featS;
    detect_features(S, fo, featS);
    {
        int nd = drop_close_parallel_features(S, featS, L, 0.5);
        if (nd) ctx.logf("сняты близкие параллельные фичери: %d рёбер", nd);
        int ns = drop_short_feature_chains(S, featS, L, 1.0);
        if (ns) ctx.logf("сняты короткие оторванные фичери: %d рёбер", ns);
    }
    size_t n_feat_he = 0;
    for (auto b : featS) n_feat_he += b;
    ctx.logf("фичери исходника: %zu полурёбер", n_feat_he);

    // плотность и размеры
    std::vector<double> vdens;
    if (in.vdensity) {
        vdens.resize(S.nV());
        for (int v = 0; v < S.nV(); ++v) vdens[v] = std::max(0.0, (double)in.vdensity[vert_origin[v]]);
    }
    std::vector<double> hS = compute_source_sizes(S, p, vdens, target_count, L, ctx);
    const bool adaptive_sizes = p.adaptive > 1e-3 || (!vdens.empty() && p.density_strength > 1e-3);
    TriBvh srcBvh;
    srcBvh.build(S.P, S.F);

    // направляющие
    SegBvh guides;
    if (in.guide_edges && in.n_guide > 0) {
        std::vector<V3> A, B;
        for (int k = 0; k < in.n_guide; ++k) {
            int a = in.guide_edges[2 * k], b = in.guide_edges[2 * k + 1];
            if (a < 0 || b < 0 || a >= nv || b >= nv) continue;
            A.push_back(V3(in.V[3 * a], in.V[3 * a + 1], in.V[3 * a + 2]));
            B.push_back(V3(in.V[3 * b], in.V[3 * b + 1], in.V[3 * b + 2]));
        }
        if (!A.empty()) guides.build(A, B);
        ctx.logf("направляющие: %zu отрезков", A.size());
    }

    // --- рабочая сетка
    TriMesh M;
    std::vector<uint8_t> feat;
    double kscale = p.working_scale > 0 ? p.working_scale : (p.d_reserved[0] > 0 ? p.d_reserved[0] : 0.0);
    const bool skip_remesh = p.working_scale < 0;
    if (kscale <= 0) kscale = n_feat_he > 0 ? 3.5 : 3.0;
    if (skip_remesh) {
        M = S;
        feat = featS;
        kscale = 0;
    } else {
        ctx.progress(0.05f, 2, "Изотропный ремешинг");
        RemeshOptions ro;
        ro.target_len = L / kscale;
        std::vector<double> srcsz;
        if (adaptive_sizes) {
            srcsz.resize(hS.size());
            for (size_t f = 0; f < hS.size(); ++f) srcsz[f] = hS[f] / kscale;
            ro.src_size = &srcsz;
        }
        ro.iterations = p.i_reserved[0] > 0 ? p.i_reserved[0] : 6;
        RemeshStats rs;
        if (!isotropic_remesh(S, featS, ro, M, feat, rs, ctx)) {
            res.error = "Ремешинг не удался";
            return;
        }
    }
    if (p.i_reserved[1] == 1) {  // отладка: вернуть рабочую сетку
        res.V.resize(M.P.size() * 3);
        for (size_t i = 0; i < M.P.size(); ++i) {
            res.V[3 * i] = (float)M.P[i].x; res.V[3 * i + 1] = (float)M.P[i].y; res.V[3 * i + 2] = (float)M.P[i].z;
        }
        for (int f = 0; f < M.nF(); ++f) {
            res.face_sizes.push_back(3);
            for (int k = 0; k < 3; ++k) res.indices.push_back(M.F[3 * f + k]);
        }
        res.ok = true;
        return;
    }
    const int nF = M.nF();
    ctx.logf("рабочая сетка: %d вершин, %d граней (%.2f с)", M.nV(), nF, now_seconds() - t0);

    ctx.progress(0.15f, 3, "Поле направлений");
    FieldGeometry G;
    compute_field_geometry(M, G);
    const double soft_feat_w = p.d_reserved[6];
    std::vector<int> feat_faces;
    std::vector<C2> z;
    std::vector<int> fidx;
    // simple = поле без жёстких фичеров и без выравнивания по кривизне: запасной путь, когда
    // фичери/шум дают тысячи особых точек и параметризация не помещается в память.
    auto build_field = [&](bool simple) {
    const std::vector<uint8_t> nofeatM(simple ? feat.size() : 0, 0);
    FieldInput fi;
    build_field_input(M, simple ? nofeatM : feat, soft_feat_w, fi, &feat_faces);
    if (p.curvature_align > 0 && !simple) add_curvature_prior(M, fi, p.curvature_align, L, p);
    add_guide_prior(M, fi, guides, p.guide_strength, 1.2 * L);
    fi.seed = (unsigned)p.seed;
    if (p.d_reserved[3] > 0) fi.eps = p.d_reserved[3];
    if (p.d_reserved[4] > 0) fi.cg_tol = p.d_reserved[4];
    if (p.i_reserved[7] > 0) fi.cg_max_iters = p.i_reserved[7];
    if (p.i_reserved[6] > 0) fi.outer_iters = p.i_reserved[6];
    const double kfield = p.field_scale > 0 ? p.field_scale : (p.d_reserved[5] > 0 ? p.d_reserved[5] : 1.0);
    bool coarse_done = false;
    if (kfield > 0 && kscale > 0 && L * kfield > 1.5 * (L / kscale)) {
        // поле считаем на грубой сетке (масштаб квада): исключает субмасштабные сингулярности
        TriMesh Cm;
        std::vector<uint8_t> featC;
        RemeshOptions rc;
        rc.target_len = L * kfield;
        rc.iterations = 6;
        std::vector<double> srcszC;
        if (adaptive_sizes) {
            srcszC.resize(hS.size());
            for (size_t f = 0; f < hS.size(); ++f) srcszC[f] = hS[f] * kfield;
            rc.src_size = &srcszC;
        }
        RemeshStats rsc;
        const std::vector<uint8_t> nofeatS(simple ? featS.size() : 0, 0);
        if (isotropic_remesh(S, simple ? nofeatS : featS, rc, Cm, featC, rsc, ctx)) {
            FieldGeometry Gc;
            compute_field_geometry(Cm, Gc);
            FieldInput fic;
            if (simple) std::fill(featC.begin(), featC.end(), 0);
            build_field_input(Cm, featC, soft_feat_w, fic, nullptr);
            if (p.curvature_align > 0 && !simple) add_curvature_prior(Cm, fic, p.curvature_align, L, p);
            add_guide_prior(Cm, fic, guides, p.guide_strength, 1.2 * L);
            fic.seed = fi.seed; fic.eps = fi.eps; fic.cg_tol = fi.cg_tol; fic.cg_max_iters = fi.cg_max_iters;
            fic.outer_iters = fi.outer_iters;
            std::vector<C2> zc;
            compute_cross_field(Cm, Gc, fic, zc, ctx);
            polish_field(Cm, Gc, fic, zc, 200, ctx);
            std::vector<int> ci;
            compute_field_index(Cm, zc, ci);
            log_index_stats(ctx, "грубая сетка", ci);
            prolong_field(Cm, zc, M, z);
            for (int f = 0; f < nF; ++f)
                if (!fi.fixed.empty() && fi.fixed[f].abs2() > 0.5) z[f] = fi.fixed[f];
            coarse_done = true;
        }
    }
    if (!coarse_done) {
        compute_cross_field(M, G, fi, z, ctx);
        compute_field_index(M, z, fidx);
        log_index_stats(ctx, "линейное", fidx);
    }
    polish_field(M, G, fi, z, p.i_reserved[5] >= 0 ? (p.i_reserved[5] > 0 ? p.i_reserved[5] : (coarse_done ? 30 : 300)) : 0, ctx, coarse_done ? 0.0 : 0.08);
    compute_field_index(M, z, fidx);
    log_index_stats(ctx, "после релаксации", fidx);
    };
    build_field(false);
    ctx.logf("поле готово (%.2f с)", now_seconds() - t0);
    if (p.i_reserved[1] == 2) {  // отладка: остановиться после поля
        std::vector<double> d;
        for (auto& q : M.P) { d.push_back(q.x); d.push_back(q.y); d.push_back(q.z); }
        store(res, "dbg_P", d);
        store(res, "dbg_F", M.F);
        store(res, "dbg_fieldindex", fidx);
        store(res, "dbg_fixed", feat_faces);
        store(res, "dbg_feat", as_int(feat));
        std::vector<double> zr(nF * 2);
        for (int f = 0; f < nF; ++f) { zr[2 * f] = z[f].re; zr[2 * f + 1] = z[f].im; }
        store(res, "dbg_z", zr);
        res.ok = true;
        return;
    }

    ctx.progress(0.30f, 4, "Параметризация");
    std::vector<double> size(nF, L);
    if (adaptive_sizes) {
        parallel_each(nF, [&](i64 f) {
            V3 c = (M.P[M.F[3 * f]] + M.P[M.F[3 * f + 1]] + M.P[M.F[3 * f + 2]]) * (1.0 / 3.0), q;
            double d2;
            int sf = srcBvh.nearest(c, q, d2);
            if (sf >= 0) size[f] = hS[sf];
        });
    }
    ParamInput pi;
    pi.M = &M;
    pi.G = &G;
    pi.z = &z;
    pi.size = &size;
    pi.feature = &feat;
    pi.verbose = p.verbose != 0;
    if (p.i_reserved[4] > 0) pi.cg_max_iters = p.i_reserved[4];
    if (p.d_reserved[2] > 0) pi.cg_tol = p.d_reserved[2];
    pi.skip_rounding = p.i_reserved[9] == 1;
    pi.direct = p.i_reserved[10] != 1;
    if (p.i_reserved[11] > 0) pi.dense_row = p.i_reserved[11];
    if (p.d_reserved[10] != 0) pi.sliver_clamp = std::max(0.0, p.d_reserved[10]);
    pi.conformal = p.d_reserved[7];
    pi.integrable = p.d_reserved[9];
    std::vector<uint8_t> fixed_faces_u8(nF, 0);
    for (int f = 0; f < nF; ++f) fixed_faces_u8[f] = feat_faces[f] ? 1 : 0;
    pi.fixed_faces = &fixed_faces_u8;
    if (p.d_reserved[8] > 0) pi.conformal_clamp = p.d_reserved[8];
    if (pi.skip_rounding) pi.repair_rounds = 0;

    const int max_iters = 3 + 2 * std::max(0, std::min(3, p.quality));
    int grow_left = std::max(1, p.quality + 1);
    int calib_left = 2;
    ParamResult pr, best_pr;
    ExtractResult xr, best_xr;
    double best_score = 1e300;
    std::vector<double> best_size = size;
    int prev_flipped = -1;
    bool last_grow = false;
    bool field_simple = false;
    for (int iter = 0; iter < max_iters; ++iter) {
        ctx.check_cancel();
        ctx.progress(0.30f + 0.08f * iter, 4, "Параметризация");
        pr = ParamResult();
        compute_param(pi, pr, ctx);
        if (!pr.ok && pr.too_complex && !field_simple) {
            field_simple = true;
            ctx.logf("поле слишком сложное (%d особых точек): пересчёт без жёстких фичеров и кривизны", pr.nsingular);
            ctx.progress(0.15f, 3, "Поле направлений (упрощённое)");
            build_field(true);
            std::fill(fixed_faces_u8.begin(), fixed_faces_u8.end(), 0);
            pi.feature = nullptr;
            continue;
        }
        if (!pr.ok) {
            if (pr.too_complex)
                res.error = "Поле направлений слишком сложное (" + std::to_string(pr.nsingular) +
                            " особых точек): отключите «Помеченные острые»/«Складки» или уменьшите число квадов";
            break;
        }
        ctx.progress(0.7f, 5, "Извлечение квадов");
        xr = ExtractResult();
        extract_quads(M, pr.uv, &feat, xr, ctx);
        const double count = (double)xr.quads.size();
        const double ratio = target_count > 0 && xr.ok ? count / target_count : 1.0;
        double score = pr.flipped * 10.0 + xr.nonquad * 10.0 + xr.dead_ends * 10.0 + std::fabs(ratio - 1.0) * 5.0;
        if (!xr.ok) score += 1e6;
        if (score < best_score) { best_score = score; best_pr = pr; best_xr = xr; best_size = size; }
        if (iter + 1 >= max_iters) break;
        const int flip_allow = std::max(2, (int)(2e-4 * nF));
        // укрупнение зон с перевёрнутыми гранями — пока оно заметно помогает
        const bool improving = prev_flipped < 0 || pr.flipped < 0.8 * prev_flipped;
        prev_flipped = pr.flipped;
        if (pr.flipped > flip_allow && grow_left > 0 && improving) {
            --grow_left;
            std::vector<double> mult(nF, 1.0);
            std::vector<uint8_t> vmark(M.nV(), 0);
            for (int f : pr.bad_faces)
                for (int c = 0; c < 3; ++c) vmark[M.F[3 * f + c]] = 1;
            for (int dilate = 0; dilate < 3; ++dilate) {
                std::vector<uint8_t> nv2 = vmark;
                for (int f = 0; f < nF; ++f) {
                    bool any = false;
                    for (int c = 0; c < 3; ++c) any = any || vmark[M.F[3 * f + c]];
                    if (any) for (int c = 0; c < 3; ++c) nv2[M.F[3 * f + c]] = 1;
                }
                vmark.swap(nv2);
            }
            for (int f = 0; f < nF; ++f) {
                bool any = false;
                for (int c = 0; c < 3; ++c) any = any || vmark[M.F[3 * f + c]];
                if (any) mult[f] = 1.4;
            }
            for (int sm = 0; sm < 6; ++sm) {  // сглаживание градиента размеров
                std::vector<double> m2 = mult;
                for (int f = 0; f < nF; ++f) {
                    double a = mult[f]; int n = 1;
                    for (int i = 0; i < 3; ++i) { int g = M.FF[3 * f + i]; if (g >= 0) { a += mult[g]; ++n; } }
                    m2[f] = a / n;
                }
                mult.swap(m2);
            }
            for (int f = 0; f < nF; ++f) size[f] *= mult[f];
            ctx.logf("итерация %d: %d перевёрнутых — укрупняем зоны сетки (осталось попыток %d)", iter, pr.flipped, grow_left);
            last_grow = true;
            continue;
        }
        // калибровка числа квадов: число ~ 1/L^2. Если укрупнение не помогло — от лучшего варианта.
        int ref_count = (int)count;
        bool ref_ok = xr.ok;
        if (last_grow && score > best_score) {
            size = best_size;
            ref_count = (int)best_xr.quads.size();
            ref_ok = best_xr.ok;
        }
        grow_left = 0;
        last_grow = false;
        const double ref_ratio = target_count > 0 ? ref_count / (double)target_count : 1.0;
        if (target_count > 0 && calib_left > 0 && ref_ok && std::fabs(ref_ratio - 1.0) > 0.06) {
            const double scale = std::sqrt(std::max(0.2, std::min(5.0, ref_ratio)));
            --calib_left;
            for (double& sz : size) sz *= scale;
            ctx.logf("калибровка: получено %d квадов при цели %d -> масштаб размера x%.3f", ref_count,
                     target_count, scale);
            continue;
        }
        break;
    }
    if (!best_xr.ok && pr.ok && xr.ok) { best_pr = pr; best_xr = xr; }
    pr = best_pr;
    xr = best_xr;
    if (xr.ok) {
        ctx.progress(0.85f, 6, "Ремонт и релаксация");
        PolyMesh pm;
        pm.V = xr.V;
        pm.vfeat = xr.vfeature;
        pm.vfeat.resize(pm.V.size(), 0);
        for (auto& q : xr.quads) pm.faces.push_back({q[0], q[1], q[2], q[3]});
        for (auto& o : xr.others) pm.faces.push_back(o);
        RepairStats rst;
        repair_quad_mesh(pm, rst, ctx);
        const int ndoublets = remove_doublets(pm);
        if (ndoublets) ctx.logf("ремонт: удалено дублетов %d", ndoublets);
        RelaxOptions ro;
        ro.iterations = p.relax_iterations > 0 ? p.relax_iterations : (p.relax_iterations < 0 ? 0 : 3);
        relax_and_project(pm, S, featS, L, ro, ctx);
        QuadQuality qq = evaluate_quad_mesh(pm);
        ctx.logf("итог: %d граней (%d квадов, %d треуг., %d n-угольн.), %d вершин, нерегулярных %d, "
                 "углы %.1f..%.1f (ср. откл. %.1f°), CV рёбер %.3f, граничных рёбер %d, неманифольдных %d",
                 qq.faces, qq.quads, qq.tris, qq.ngons, qq.verts, qq.irregular, qq.min_angle, qq.max_angle,
                 qq.angle_dev_mean, qq.edge_cv, qq.boundary_edges, qq.nonmanifold_edges);
        res.V.resize(pm.V.size() * 3);
        for (size_t i = 0; i < pm.V.size(); ++i) {
            res.V[3 * i] = (float)pm.V[i].x;
            res.V[3 * i + 1] = (float)pm.V[i].y;
            res.V[3 * i + 2] = (float)pm.V[i].z;
        }
        for (auto& f : pm.faces) {
            res.face_sizes.push_back((int)f.size());
            for (int k : f) res.indices.push_back(k);
        }
        res.darr["quality"] = {(double)qq.faces, (double)qq.quads, (double)qq.tris, (double)qq.ngons,
                               (double)qq.verts, (double)qq.irregular, qq.min_angle, qq.max_angle,
                               qq.angle_dev_mean, qq.edge_cv, (double)qq.boundary_edges,
                               (double)qq.nonmanifold_edges, (double)pr.flipped, (double)rst.removed_faces};
    }

    if (p.debug) {
        std::vector<double> d;
        for (auto& q : M.P) { d.push_back(q.x); d.push_back(q.y); d.push_back(q.z); }
        store(res, "dbg_P", d);
        store(res, "dbg_F", M.F);
        store(res, "dbg_uv", pr.uv);
        store(res, "dbg_rep", pr.rep);
        if (!pr.phi.empty()) store(res, "dbg_phi", pr.phi);
        if (!pr.foreign.empty()) store(res, "dbg_foreign", pr.foreign);
        store(res, "dbg_cut", as_int(pr.cut));
        store(res, "dbg_vindex", as_int(pr.vindex));
        store(res, "dbg_feat", as_int(feat));
        store(res, "dbg_fieldindex", fidx);
        std::vector<double> zr(nF * 2);
        for (int f = 0; f < nF; ++f) { zr[2 * f] = z[f].re; zr[2 * f + 1] = z[f].im; }
        store(res, "dbg_z", zr);
        store(res, "dbg_fixed", feat_faces);
    }
    ctx.logf("всего %.2f с", now_seconds() - t0);
    res.ok = pr.ok && xr.ok && !res.face_sizes.empty();
}

}  // namespace qf
