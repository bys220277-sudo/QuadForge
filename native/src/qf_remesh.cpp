#include "qf_remesh.h"

#include <unordered_map>
#include <unordered_set>

namespace qf {

namespace {

struct Tri {
    int v[3];
    int nb[3];
    int nbe[3];
    uint8_t ef[3];
};

struct Dyn {
    std::vector<V3> P;
    std::vector<int> vface;       // -1: удалена
    std::vector<double> vh;       // целевая длина в вершине
    std::vector<Tri> T;
    std::vector<uint8_t> tdead;
    int alive_faces = 0;

    int corner_of(int f, int v) const {
        const Tri& t = T[f];
        return t.v[0] == v ? 0 : (t.v[1] == v ? 1 : 2);
    }
    void link(int f, int i, int g, int j) {
        T[f].nb[i] = g; T[f].nbe[i] = j;
        T[g].nb[j] = f; T[g].nbe[j] = i;
    }
    // связать ребро i грани f с гранью g (ищем обратное ребро); g<0 -> граница
    bool link_v(int f, int i, int g) {
        if (g < 0) { T[f].nb[i] = -1; T[f].nbe[i] = -1; return true; }
        int a = T[f].v[i], b = T[f].v[(i + 1) % 3];
        for (int j = 0; j < 3; ++j)
            if (T[g].v[j] == b && T[g].v[(j + 1) % 3] == a) { link(f, i, g, j); return true; }
        return false;
    }
    int add_vertex(const V3& p, double h) {
        P.push_back(p);
        vface.push_back(-1);
        vh.push_back(h);
        return (int)P.size() - 1;
    }
    int add_tri() {
        Tri t;
        std::memset(&t, 0, sizeof(t));
        T.push_back(t);
        tdead.push_back(0);
        ++alive_faces;
        return (int)T.size() - 1;
    }
    void kill(int f) {
        if (!tdead[f]) { tdead[f] = 1; --alive_faces; }
    }

    // Грани вокруг v (порядок как в TriMesh веерах). Возвращает true если замкнуто.
    bool ring(int v, std::vector<std::pair<int, int>>& out) const {
        out.clear();
        int f0 = vface[v];
        int c0 = corner_of(f0, v);
        int f = f0, c = c0;
        bool closed = false;
        for (int guard = 0; guard < 1000000; ++guard) {
            int g = T[f].nb[(c + 2) % 3];
            if (g < 0) break;
            int j = T[f].nbe[(c + 2) % 3];
            f = g; c = j;
            if (f == f0) { closed = true; break; }
        }
        int sf = f, sc = c;
        for (int guard = 0; guard < 1000000; ++guard) {
            out.push_back({f, c});
            int g = T[f].nb[c];
            if (g < 0) break;
            int j = T[f].nbe[c];
            int nf = g, nc = (j + 1) % 3;
            if (nf == sf && nc == sc) break;
            f = nf; c = nc;
        }
        return closed;
    }

    void neighbors(int v, std::vector<int>& out, bool* closed_out = nullptr) const {
        thread_local std::vector<std::pair<int, int>> r;
        bool closed = ring(v, r);
        out.clear();
        for (auto& fc : r) out.push_back(T[fc.first].v[(fc.second + 1) % 3]);
        if (!closed) out.push_back(T[r[0].first].v[(r[0].second + 2) % 3]);
        if (closed_out) *closed_out = closed;
    }

    // тип вершины: 0 поверхность, 1 цепочка фичера, 2 угол/узел (фиксирован)
    int vtype(int v, double corner_angle) const {
        thread_local std::vector<std::pair<int, int>> r;
        bool closed = ring(v, r);
        int nfeat = 0;
        int nbv[3] = {-1, -1, -1};
        for (auto& fc : r)
            if (T[fc.first].ef[fc.second]) {
                if (nfeat < 3) nbv[nfeat] = T[fc.first].v[(fc.second + 1) % 3];
                ++nfeat;
            }
        if (!closed && T[r[0].first].ef[(r[0].second + 2) % 3]) {
            if (nfeat < 3) nbv[nfeat] = T[r[0].first].v[(r[0].second + 2) % 3];
            ++nfeat;
        }
        if (nfeat == 0) return 0;
        if (nfeat != 2) return 2;
        V3 a = P[nbv[0]] - P[v], b = P[nbv[1]] - P[v];
        double la = norm(a), lb = norm(b);
        if (la < 1e-300 || lb < 1e-300) return 2;
        double cs = dot(a, b) / (la * lb);
        double ang = std::acos(std::max(-1.0, std::min(1.0, cs)));
        if (kPi - ang > corner_angle) return 2;
        return 1;
    }

    // Отладочная проверка согласованности. Возвращает описание первой ошибки или "".
    std::string validate() const {
        char buf[200];
        for (size_t f = 0; f < T.size(); ++f) {
            if (tdead[f]) continue;
            for (int i = 0; i < 3; ++i) {
                int g = T[f].nb[i];
                if (g == -1) continue;
                if (g < 0 || g >= (int)T.size() || tdead[g]) {
                    std::snprintf(buf, sizeof(buf), "face %zu edge %d -> dead/invalid neighbour %d", f, i, g);
                    return buf;
                }
                int j = T[f].nbe[i];
                if (j < 0 || j > 2 || T[g].nb[j] != (int)f || T[g].nbe[j] != i) {
                    std::snprintf(buf, sizeof(buf), "face %zu edge %d <-> %d asym (j=%d)", f, i, g, j);
                    return buf;
                }
                if (T[g].v[j] != T[f].v[(i + 1) % 3] || T[g].v[(j + 1) % 3] != T[f].v[i]) {
                    std::snprintf(buf, sizeof(buf), "face %zu edge %d <-> %d vertex mismatch", f, i, g);
                    return buf;
                }
            }
            for (int i = 0; i < 3; ++i) {
                int vf = vface[T[f].v[i]];
                if (vf < 0 || tdead[vf]) {
                    std::snprintf(buf, sizeof(buf), "vertex %d has dead vface (face %zu)", T[f].v[i], f);
                    return buf;
                }
                const Tri& tv = T[vf];
                if (tv.v[0] != T[f].v[i] && tv.v[1] != T[f].v[i] && tv.v[2] != T[f].v[i]) {
                    std::snprintf(buf, sizeof(buf), "vertex %d vface %d does not contain it (face %zu)", T[f].v[i], vf, f);
                    return buf;
                }
            }
        }
        return "";
    }

    V3 face_normal(int f) const {
        const Tri& t = T[f];
        return cross(P[t.v[1]] - P[t.v[0]], P[t.v[2]] - P[t.v[0]]);
    }

    double edge_len(int f, int i) const {
        const Tri& t = T[f];
        return norm(P[t.v[(i + 1) % 3]] - P[t.v[i]]);
    }
};

// ---------------------------------------------------------------------------
int split_edge(Dyn& D, int f, int i, const V3& pos, double h) {
    const Tri F0 = D.T[f];
    int a = F0.v[i], b = F0.v[(i + 1) % 3], c = F0.v[(i + 2) % 3];
    int g = F0.nb[i], j = F0.nbe[i];
    uint8_t e_ab = F0.ef[i], e_bc = F0.ef[(i + 1) % 3], e_ca = F0.ef[(i + 2) % 3];
    int n_bc = F0.nb[(i + 1) % 3], n_ca = F0.nb[(i + 2) % 3];
    int m = D.add_vertex(pos, h);
    int fn = D.add_tri();
    // f=(a,m,c), fn=(m,b,c)
    {
        Tri& t = D.T[f];
        t.v[0] = a; t.v[1] = m; t.v[2] = c;
        t.ef[0] = e_ab; t.ef[1] = 0; t.ef[2] = e_ca;
        Tri& u = D.T[fn];
        u.v[0] = m; u.v[1] = b; u.v[2] = c;
        u.ef[0] = e_ab; u.ef[1] = e_bc; u.ef[2] = 0;
    }
    D.link(f, 1, fn, 2);
    D.link_v(f, 2, n_ca);
    D.link_v(fn, 1, n_bc);
    if (g >= 0) {
        const Tri G0 = D.T[g];
        int d = G0.v[(j + 2) % 3];
        int n_ad = G0.nb[(j + 1) % 3], n_db = G0.nb[(j + 2) % 3];
        uint8_t e_ad = G0.ef[(j + 1) % 3], e_db = G0.ef[(j + 2) % 3];
        int gn = D.add_tri();
        Tri& t = D.T[g];
        t.v[0] = b; t.v[1] = m; t.v[2] = d;
        t.ef[0] = e_ab; t.ef[1] = 0; t.ef[2] = e_db;
        Tri& u = D.T[gn];
        u.v[0] = m; u.v[1] = a; u.v[2] = d;
        u.ef[0] = e_ab; u.ef[1] = e_ad; u.ef[2] = 0;
        D.link(f, 0, gn, 0);
        D.link(fn, 0, g, 0);
        D.link(g, 1, gn, 2);
        D.link_v(g, 2, n_db);
        D.link_v(gn, 1, n_ad);
        D.vface[d] = g;
        D.vface[a] = f;
        D.vface[b] = fn;
        D.vface[m] = f;
        D.vface[c] = f;
    } else {
        D.T[f].nb[0] = -1; D.T[f].nbe[0] = -1;
        D.T[fn].nb[0] = -1; D.T[fn].nbe[0] = -1;
        D.vface[a] = f;
        D.vface[b] = fn;
        D.vface[m] = f;
        D.vface[c] = f;
    }
    return m;
}

bool flip_edge(Dyn& D, int f, int i) {
    const Tri F0 = D.T[f];
    int g = F0.nb[i], j = F0.nbe[i];
    if (g < 0) return false;
    const Tri G0 = D.T[g];
    int a = F0.v[i], b = F0.v[(i + 1) % 3], c = F0.v[(i + 2) % 3];
    int d = G0.v[(j + 2) % 3];
    int n1 = F0.nb[(i + 1) % 3], n2 = F0.nb[(i + 2) % 3];
    int n3 = G0.nb[(j + 1) % 3], n4 = G0.nb[(j + 2) % 3];
    if (n1 == f || n1 == g || n2 == f || n2 == g || n3 == f || n3 == g || n4 == f || n4 == g) return false;
    uint8_t ef1 = F0.ef[(i + 1) % 3], ef2 = F0.ef[(i + 2) % 3];
    uint8_t ef3 = G0.ef[(j + 1) % 3], ef4 = G0.ef[(j + 2) % 3];
    Tri& t = D.T[f];
    t.v[0] = a; t.v[1] = d; t.v[2] = c;
    t.ef[0] = ef3; t.ef[1] = 0; t.ef[2] = ef2;
    Tri& u = D.T[g];
    u.v[0] = d; u.v[1] = b; u.v[2] = c;
    u.ef[0] = ef4; u.ef[1] = ef1; u.ef[2] = 0;
    D.link(f, 1, g, 2);
    D.link_v(f, 0, n3);
    D.link_v(f, 2, n2);
    D.link_v(g, 0, n4);
    D.link_v(g, 1, n1);
    D.vface[a] = f; D.vface[d] = f; D.vface[c] = f; D.vface[b] = g;
    return true;
}

// Проверка и выполнение коллапса a=T[f].v[i] -> b=T[f].v[i+1].
bool try_collapse(Dyn& D, int f, int i, double max_len, double corner_angle, bool do_it) {
    const Tri F0 = D.T[f];
    int a = F0.v[i], b = F0.v[(i + 1) % 3], c = F0.v[(i + 2) % 3];
    int g = F0.nb[i], j = F0.nbe[i];
    bool edge_feat = F0.ef[i] != 0;
    int ta = D.vtype(a, corner_angle);
    if (edge_feat) { if (ta != 1) return false; }
    else {
        if (ta != 0) return false;
    }
    thread_local std::vector<std::pair<int, int>> ra;
    thread_local std::vector<int> na, nb_;
    bool closed_a;
    D.neighbors(a, na, &closed_a);
    bool closed_b;
    D.neighbors(b, nb_, &closed_b);
    int d = -1;
    if (g >= 0) {
        d = D.T[g].v[(j + 2) % 3];
        // оба на границе, а ребро внутреннее — нельзя
        if (!closed_a && !closed_b) return false;
    }
    // link condition
    {
        std::vector<int> sa = na, sb = nb_;
        std::sort(sa.begin(), sa.end());
        std::sort(sb.begin(), sb.end());
        int common = 0;
        size_t x = 0, y = 0;
        while (x < sa.size() && y < sb.size()) {
            if (sa[x] == sb[y]) { ++common; ++x; ++y; }
            else if (sa[x] < sb[y]) ++x; else ++y;
        }
        int need = g >= 0 ? 2 : 1;
        if (common != need) return false;
    }
    // валентность не должна падать до < 3
    if (closed_a && na.size() <= 3) return false;
    D.ring(a, ra);
    const V3 pb = D.P[b];
    for (auto& fc : ra) {
        int h = fc.first;
        if (h == f || h == g) continue;
        const Tri& th = D.T[h];
        V3 n0 = D.face_normal(h);
        V3 q[3] = {D.P[th.v[0]], D.P[th.v[1]], D.P[th.v[2]]};
        q[fc.second] = pb;
        V3 n1 = cross(q[1] - q[0], q[2] - q[0]);
        double l0 = norm(n0), l1 = norm(n1);
        if (l1 < 1e-14 * std::max(l0, 1e-300) || l1 < 1e-30) return false;
        if (dot(n0, n1) < 0.25 * l0 * l1) return false;
    }
    for (int w : na) {
        if (w == b) continue;
        if (norm(D.P[w] - pb) > max_len) return false;
    }
    if (!do_it) return true;

    // --- выполнение
    int n_bc = F0.nb[(i + 1) % 3], n_ca = F0.nb[(i + 2) % 3];
    uint8_t ef_bc = F0.ef[(i + 1) % 3], ef_ca = F0.ef[(i + 2) % 3];
    int n_ad = -2, n_db = -2;
    uint8_t ef_ad = 0, ef_db = 0;
    if (g >= 0) {
        const Tri& G0 = D.T[g];
        n_ad = G0.nb[(j + 1) % 3]; n_db = G0.nb[(j + 2) % 3];
        ef_ad = G0.ef[(j + 1) % 3]; ef_db = G0.ef[(j + 2) % 3];
    }
    // замена a -> b в остальных гранях
    for (auto& fc : ra) {
        if (fc.first == f || fc.first == g) continue;
        D.T[fc.first].v[fc.second] = b;
    }
    D.kill(f);
    if (g >= 0) D.kill(g);
    // перелинковка за f: рёбра (b,c) и (c,a)
    {
        uint8_t fl = ef_bc | ef_ca;
        if (n_bc >= 0 && n_ca >= 0) {
            // найти рёбра в соседях
            int jb = -1, jc = -1;
            for (int k = 0; k < 3; ++k) {
                if (D.T[n_bc].nb[k] == f) jb = k;
                if (D.T[n_ca].nb[k] == f) jc = k;
            }
            if (jb < 0 || jc < 0) return true;  // не должно происходить
            D.link(n_bc, jb, n_ca, jc);
            D.T[n_bc].ef[jb] = fl;
            D.T[n_ca].ef[jc] = fl;
        } else if (n_bc >= 0) {
            for (int k = 0; k < 3; ++k) if (D.T[n_bc].nb[k] == f) { D.T[n_bc].nb[k] = -1; D.T[n_bc].nbe[k] = -1; D.T[n_bc].ef[k] = 1; }
        } else if (n_ca >= 0) {
            for (int k = 0; k < 3; ++k) if (D.T[n_ca].nb[k] == f) { D.T[n_ca].nb[k] = -1; D.T[n_ca].nbe[k] = -1; D.T[n_ca].ef[k] = 1; }
        }
    }
    if (g >= 0) {
        uint8_t fl = ef_ad | ef_db;
        if (n_ad >= 0 && n_db >= 0) {
            int ja = -1, jd = -1;
            for (int k = 0; k < 3; ++k) {
                if (D.T[n_ad].nb[k] == g) ja = k;
                if (D.T[n_db].nb[k] == g) jd = k;
            }
            if (ja < 0 || jd < 0) return true;
            D.link(n_ad, ja, n_db, jd);
            D.T[n_ad].ef[ja] = fl;
            D.T[n_db].ef[jd] = fl;
        } else if (n_ad >= 0) {
            for (int k = 0; k < 3; ++k) if (D.T[n_ad].nb[k] == g) { D.T[n_ad].nb[k] = -1; D.T[n_ad].nbe[k] = -1; D.T[n_ad].ef[k] = 1; }
        } else if (n_db >= 0) {
            for (int k = 0; k < 3; ++k) if (D.T[n_db].nb[k] == g) { D.T[n_db].nb[k] = -1; D.T[n_db].nbe[k] = -1; D.T[n_db].ef[k] = 1; }
        }
    }
    D.vface[a] = -1;
    for (auto& fc : ra) {
        if (fc.first == f || fc.first == g) continue;
        for (int k = 0; k < 3; ++k) D.vface[D.T[fc.first].v[k]] = fc.first;
    }
    (void)c; (void)d;
    return true;
}


// Геометрическая допустимость флипа ребра (f,i): внутреннее, не фичер, нет дубликата ребра c-d,
// нормали новых граней близки к старым, площади положительны.
// Если новое ребро c-d длиннее целевого размера — строгий порог kFlipCosLong (~14°): при 45°
// каскад флипов Делоне на длинных узких треугольниках CAD строил «перепонки» сквозь деталь.
// Для коротких рёбер оставляем мягкий kFlipCos: строгий везде ломал мелкие отверстия (flange, gear).
static constexpr double kFlipCos = 0.7;
static constexpr double kFlipCosLong = 0.97;
bool flip_geometry_ok(const Dyn& D, int f, int i) {
    const Tri& t = D.T[f];
    int g = t.nb[i];
    if (g < 0 || t.ef[i]) return false;
    int j = t.nbe[i];
    int a = t.v[i], b = t.v[(i + 1) % 3], c = t.v[(i + 2) % 3];
    int d = D.T[g].v[(j + 2) % 3];
    if (a == d || b == c || c == d) return false;
    {
        std::vector<int> nc;
        D.neighbors(c, nc);
        if (std::find(nc.begin(), nc.end(), d) != nc.end()) return false;
    }
    V3 n_old_f = D.face_normal(f), n_old_g = D.face_normal(g);
    V3 n_new_f = cross(D.P[d] - D.P[a], D.P[c] - D.P[a]);
    V3 n_new_g = cross(D.P[b] - D.P[d], D.P[c] - D.P[d]);
    double lf = norm(n_new_f), lg = norm(n_new_g);
    double lof = norm(n_old_f), log_ = norm(n_old_g);
    if (lf < 1e-12 * (lof + log_) || lg < 1e-12 * (lof + log_)) return false;
    const double fc = norm(D.P[c] - D.P[d]) > std::min(D.vh[c], D.vh[d]) ? kFlipCosLong : kFlipCos;
    if (dot(n_new_f, n_old_f) < fc * lf * lof || dot(n_new_f, n_old_g) < fc * lf * log_) return false;
    if (dot(n_new_g, n_old_f) < fc * lg * lof || dot(n_new_g, n_old_g) < fc * lg * log_) return false;
    return true;
}

// Критерий Делоне (сумма противолежащих углов > pi)
bool delaunay_violated(const Dyn& D, int f, int i) {
    const Tri& t = D.T[f];
    int g = t.nb[i];
    if (g < 0) return false;
    int j = t.nbe[i];
    const V3& pa = D.P[t.v[i]];
    const V3& pb = D.P[t.v[(i + 1) % 3]];
    const V3& pc = D.P[t.v[(i + 2) % 3]];
    const V3& pd = D.P[D.T[g].v[(j + 2) % 3]];
    auto cot_at = [](const V3& p, const V3& q, const V3& r) {  // угол при p между q и r
        V3 u = q - p, w = r - p;
        double cr = norm(cross(u, w));
        return cr > 1e-300 ? dot(u, w) / cr : 0.0;
    };
    return cot_at(pc, pa, pb) + cot_at(pd, pa, pb) < -1e-9;
}

// Каскад флипов Делоне вокруг ребра после вставки вершины.
void delaunay_fix(Dyn& D, std::vector<std::pair<int, int>>& stack) {
    int guard = 0;
    while (!stack.empty() && guard++ < 200000) {
        auto [f, i] = stack.back();
        stack.pop_back();
        if (D.tdead[f]) continue;
        if (!delaunay_violated(D, f, i) || !flip_geometry_ok(D, f, i)) continue;
        const Tri F0 = D.T[f];
        int j = F0.nbe[i];
        int g = F0.nb[i];
        (void)j;
        if (!flip_edge(D, f, i)) continue;
        // после флипа f=(a,d,c): внешние рёбра e0 (a,d) и e2 (c,a); g=(d,b,c): e0 (d,b), e1 (b,c)
        stack.push_back({f, 0});
        stack.push_back({f, 2});
        stack.push_back({g, 0});
        stack.push_back({g, 1});
    }
}

}  // namespace

bool isotropic_remesh(const TriMesh& src, const std::vector<uint8_t>& feat, const RemeshOptions& opt,
                      TriMesh& out, std::vector<uint8_t>& out_feat, RemeshStats& stats, Ctx& ctx) {
    double t0 = now_seconds();
    Dyn D;
    const int nf0 = src.nF(), nv0 = src.nV();
    stats.faces_in = nf0;
    D.P = src.P;
    D.vface.assign(nv0, -1);
    D.vh.assign(nv0, opt.target_len);
    D.T.resize(nf0);
    D.tdead.assign(nf0, 0);
    D.alive_faces = nf0;
    for (int f = 0; f < nf0; ++f) {
        Tri& t = D.T[f];
        for (int i = 0; i < 3; ++i) {
            t.v[i] = src.F[3 * f + i];
            t.nb[i] = src.FF[3 * f + i];
            t.nbe[i] = src.FE[3 * f + i];
            t.ef[i] = (src.FF[3 * f + i] < 0) ? 1 : (feat[3 * f + i] ? 1 : 0);
            D.vface[t.v[i]] = f;
        }
    }
    // BVH: поверхность и фичер-отрезки исходника
    TriBvh surf;
    surf.build(src.P, src.F);
    SegBvh fsegs;
    {
        std::vector<V3> A, B;
        for (int f = 0; f < nf0; ++f)
            for (int i = 0; i < 3; ++i) {
                int g = src.FF[3 * f + i];
                bool is = (g < 0) || feat[3 * f + i];
                if (!is) continue;
                if (g >= 0 && g < f) continue;
                A.push_back(src.P[src.F[3 * f + i]]);
                B.push_back(src.P[src.F[3 * f + (i + 1) % 3]]);
            }
        if (!A.empty()) fsegs.build(A, B);
    }
    auto size_at = [&](int src_face) {
        return opt.src_size ? (*opt.src_size)[src_face] : opt.target_len;
    };
    // начальные размеры в вершинах
    if (opt.src_size) {
        for (int v = 0; v < nv0; ++v) {
            V3 q;
            double d2;
            int tri = surf.nearest(D.P[v], q, d2);
            if (tri >= 0) D.vh[v] = size_at(tri);
        }
    }
    const double corner_angle = opt.corner_angle;
    const bool dbg = std::getenv("QF_REMESH_DEBUG") != nullptr;
    const bool dbg_fold = std::getenv("QF_REMESH_FOLDS") != nullptr;
    auto count_folds = [&](const char* where, int it) {
        if (!dbg_fold) return;
        int nfold = 0, nsharp = 0;
        for (size_t f = 0; f < D.T.size(); ++f) {
            if (D.tdead[f]) continue;
            for (int i = 0; i < 3; ++i) {
                int g = D.T[f].nb[i];
                if (g < 0 || (size_t)g < f || D.tdead[g] || D.T[f].ef[i]) continue;
                V3 a = normalized(D.face_normal((int)f), V3(0, 0, 0)), b = normalized(D.face_normal(g), V3(0, 0, 0));
                double c = dot(a, b);
                if (c < 0) ++nfold;
                else if (c < 0.5) ++nsharp;
            }
        }
        ctx.logf("  folds %s it %d: обратных %d, острых %d", where, it, nfold, nsharp);
    };
    if (dbg) {
        std::string e = D.validate();
        if (!e.empty()) { ctx.logf("ОШИБКА на входе: %s", e.c_str()); return false; }
    }

    auto edge_h = [&](int a, int b) { return 0.5 * (D.vh[a] + D.vh[b]); };

    for (int it = 0; it < opt.iterations; ++it) {
        ctx.check_cancel();
        // ------------------------------------------------------ расщепление
        for (int pass = 0; pass < 4; ++pass) {
            int nsplit = 0;
            const int nfaces_now = (int)D.T.size();
            for (int f = 0; f < nfaces_now; ++f) {
                if (D.tdead[f]) continue;
                for (int i = 0; i < 3; ++i) {
                    if (D.tdead[f]) break;
                    const Tri& t = D.T[f];
                    int a = t.v[i], b = t.v[(i + 1) % 3];
                    int g = t.nb[i];
                    if (g >= 0 && g < f) continue;  // обрабатываем ребро один раз (но g может быть новым)
                    double len = norm(D.P[b] - D.P[a]);
                    double h = edge_h(a, b);
                    if (len > (4.0 / 3.0) * h) {
                        if (D.alive_faces > opt.max_faces) continue;
                        V3 mid = (D.P[a] + D.P[b]) * 0.5;
                        int m = split_edge(D, f, i, mid, h);
                        ++nsplit;
                        {
                            // рёбра, противолежащие m в окружающих гранях
                            thread_local std::vector<std::pair<int, int>> st;
                            thread_local std::vector<std::pair<int, int>> rr;
                            st.clear();
                            D.ring(m, rr);
                            for (auto& fc : rr) st.push_back({fc.first, (fc.second + 1) % 3});
                            delaunay_fix(D, st);
                        }
                        if (dbg) { std::string e = D.validate(); if (!e.empty()) { ctx.logf("ОШИБКА после split: %s", e.c_str()); return false; } }
                    }
                }
            }
            stats.splits += nsplit;
            if (ctx.verbose) ctx.logf("  split pass %d: %d (граней %d)", pass, nsplit, D.alive_faces);
            if (nsplit == 0) break;
        }
        count_folds("split", it);
        // ------------------------------------------------------ схлопывание
        for (int pass = 0; pass < 8; ++pass) {
            int ncol = 0;
            const int nfaces_now = (int)D.T.size();
            for (int f = 0; f < nfaces_now; ++f) {
                if (D.tdead[f]) continue;
                for (int i = 0; i < 3; ++i) {
                    if (D.tdead[f]) break;
                    const Tri& t = D.T[f];
                    int a = t.v[i], b = t.v[(i + 1) % 3];
                    double len = norm(D.P[b] - D.P[a]);
                    double h = edge_h(a, b);
                    if (len < 0.8 * h) {
                        double maxl = (4.0 / 3.0) * h;
                        if (try_collapse(D, f, i, maxl, corner_angle, true)) {
                            ++ncol;
                            if (dbg) { std::string e = D.validate(); if (!e.empty()) { ctx.logf("ОШИБКА после collapse: %s", e.c_str()); return false; } }
                            break;
                        }
                        // пробуем обратное направление
                        const Tri& t2 = D.T[f];
                        int g = t2.nb[i];
                        if (g >= 0 && !D.tdead[g]) {
                            int j = t2.nbe[i];
                            if (try_collapse(D, g, j, maxl, corner_angle, true)) {
                                ++ncol;
                                if (dbg) { std::string e = D.validate(); if (!e.empty()) { ctx.logf("ОШИБКА после collapse(rev): %s", e.c_str()); return false; } }
                                break;
                            }
                        }
                    }
                }
            }
            stats.collapses += ncol;
            if (ctx.verbose) ctx.logf("  collapse pass %d: %d (граней %d)", pass, ncol, D.alive_faces);
            if (ncol == 0) break;
        }
        count_folds("collapse", it);
        // ------------------------------------------------------ флипы
        {
            int nflip = 0;
            // 1) Делоне
            for (int sweep = 0; sweep < 3; ++sweep) {
                int n = 0;
                const int nfaces_now = (int)D.T.size();
                for (int f = 0; f < nfaces_now; ++f) {
                    if (D.tdead[f]) continue;
                    for (int i = 0; i < 3; ++i) {
                        if (D.tdead[f]) break;
                        const Tri& t = D.T[f];
                        int g = t.nb[i];
                        if (g < 0 || g < f) continue;
                        if (delaunay_violated(D, f, i) && flip_geometry_ok(D, f, i) && flip_edge(D, f, i)) ++n;
                    }
                }
                nflip += n;
                if (n == 0) break;
            }
            // 2) выравнивание валентностей
            const int nfaces_now = (int)D.T.size();
            for (int f = 0; f < nfaces_now; ++f) {
                if (D.tdead[f]) continue;
                for (int i = 0; i < 3; ++i) {
                    if (D.tdead[f]) break;
                    const Tri& t = D.T[f];
                    int g = t.nb[i];
                    if (g < 0 || g < f || t.ef[i]) continue;
                    int j = t.nbe[i];
                    int a = t.v[i], b = t.v[(i + 1) % 3], c = t.v[(i + 2) % 3];
                    int d = D.T[g].v[(j + 2) % 3];
                    auto val_target = [&](int v) {
                        bool closed;
                        thread_local std::vector<int> nbv;
                        D.neighbors(v, nbv, &closed);
                        int tv = closed ? 6 : 4;
                        return std::make_pair((int)nbv.size(), tv);
                    };
                    auto pa = val_target(a), pb = val_target(b), pc = val_target(c), pd = val_target(d);
                    int before = (pa.first - pa.second) * (pa.first - pa.second) +
                                 (pb.first - pb.second) * (pb.first - pb.second) +
                                 (pc.first - pc.second) * (pc.first - pc.second) +
                                 (pd.first - pd.second) * (pd.first - pd.second);
                    int after = (pa.first - 1 - pa.second) * (pa.first - 1 - pa.second) +
                                (pb.first - 1 - pb.second) * (pb.first - 1 - pb.second) +
                                (pc.first + 1 - pc.second) * (pc.first + 1 - pc.second) +
                                (pd.first + 1 - pd.second) * (pd.first + 1 - pd.second);
                    if (after >= before) continue;
                    if (!flip_geometry_ok(D, f, i)) continue;
                    double hh = edge_h(a, b);
                    if (norm(D.P[c] - D.P[d]) > 2.0 * hh) continue;
                    if (flip_edge(D, f, i)) {
                        ++nflip;
                        if (dbg) { std::string e = D.validate(); if (!e.empty()) { ctx.logf("ОШИБКА после flip: %s", e.c_str()); return false; } }
                    }
                }
            }
            stats.flips += nflip;
        }
        count_folds("flip", it);
        // ------------------------------------------------------ сглаживание + проекция
        {
            const int nv = (int)D.P.size();
            std::vector<V3> newP(nv);
            std::vector<uint8_t> vt(nv, 255);
            std::vector<int> valid;
            for (int v = 0; v < nv; ++v)
                if (D.vface[v] >= 0 && !D.tdead[D.vface[v]]) valid.push_back(v);
            parallel_for((i64)valid.size(), [&](i64 b, i64 e) {
                std::vector<std::pair<int, int>> r;
                for (i64 k = b; k < e; ++k) {
                    int v = valid[k];
                    newP[v] = D.P[v];
                    int type = D.vtype(v, corner_angle);
                    vt[v] = (uint8_t)type;
                    if (type == 2) continue;
                    bool closed = D.ring(v, r);
                    V3 q(0, 0, 0);
                    double wsum = 0;
                    V3 nrm(0, 0, 0);
                    if (type == 0) {
                        for (auto& fc : r) {
                            const Tri& t = D.T[fc.first];
                            int w = t.v[(fc.second + 1) % 3];
                            q += D.P[w];
                            wsum += 1;
                            nrm += D.face_normal(fc.first);
                        }
                        if (!closed) {
                            q += D.P[D.T[r[0].first].v[(r[0].second + 2) % 3]];
                            wsum += 1;
                        }
                    } else {
                        for (auto& fc : r) {
                            const Tri& t = D.T[fc.first];
                            if (t.ef[fc.second]) { q += D.P[t.v[(fc.second + 1) % 3]]; wsum += 1; }
                        }
                        if (!closed && D.T[r[0].first].ef[(r[0].second + 2) % 3]) {
                            q += D.P[D.T[r[0].first].v[(r[0].second + 2) % 3]];
                            wsum += 1;
                        }
                    }
                    if (wsum < 1) continue;
                    q = q / wsum;
                    V3 delta = (q - D.P[v]) * opt.smooth_lambda;
                    if (type == 0) {
                        V3 n = normalized(nrm, V3(0, 0, 1));
                        delta = delta - n * dot(n, delta);
                    }
                    newP[v] = D.P[v] + delta;
                }
            }, 256);
            // проекция (только на грани исходника с согласованной нормалью — защита тонких частей)
            std::vector<V3> vnorm(D.P.size());
            parallel_for((i64)valid.size(), [&](i64 b, i64 e) {
                std::vector<std::pair<int, int>> r;
                for (i64 k = b; k < e; ++k) {
                    int v = valid[k];
                    D.ring(v, r);
                    V3 n;
                    for (auto& fc : r) n += D.face_normal(fc.first);
                    vnorm[v] = normalized(n, V3(0, 0, 0));
                }
            }, 256);
            parallel_for((i64)valid.size(), [&](i64 b, i64 e) {
                for (i64 k = b; k < e; ++k) {
                    int v = valid[k];
                    if (vt[v] == 2 || vt[v] == 255) continue;
                    V3 qpt;
                    double d2;
                    if (vt[v] == 0) {
                        const V3 nv = vnorm[v];
                        int tri = -1;
                        if (norm2(nv) > 0.5)
                            tri = surf.nearest_if(newP[v], qpt, d2, [&](int f) { return dot(src.FN[f], nv) > 0.0; });
                        if (tri < 0) tri = surf.nearest(newP[v], qpt, d2);
                        if (tri >= 0) { newP[v] = qpt; if (opt.src_size) D.vh[v] = size_at(tri); }
                    } else if (!fsegs.bvh.empty()) {
                        int s = fsegs.nearest(newP[v], qpt, d2);
                        if (s >= 0) newP[v] = qpt;
                    }
                }
            }, 256);
            // откат перемещений, переворачивающих грани. Вершины двигаются одновременно, поэтому
            // проверяем грани с учётом новых позиций всех трёх вершин и откатываем итеративно,
            // пока ни одна грань не перевёрнута относительно своего положения до сглаживания.
            std::vector<uint8_t> revert(D.P.size(), 0);
            {
                const int nT = (int)D.T.size();
                std::vector<uint8_t> bad(nT, 0);
                for (int rit = 0; rit < 20; ++rit) {
                    parallel_for(nT, [&](i64 b, i64 e) {
                        for (i64 f = b; f < e; ++f) {
                            bad[f] = 0;
                            if (D.tdead[f]) continue;
                            const Tri& t = D.T[f];
                            bool moved = false;
                            V3 q[3];
                            for (int c = 0; c < 3; ++c) {
                                const int v = t.v[c];
                                const bool mv = !revert[v] && vt[v] != 255;
                                moved = moved || mv;
                                q[c] = mv ? newP[v] : D.P[v];
                            }
                            if (!moved) continue;
                            V3 n0 = D.face_normal((int)f);
                            V3 n1 = cross(q[1] - q[0], q[2] - q[0]);
                            if (dot(n0, n1) <= 0.1 * norm(n0) * norm(n1)) bad[f] = 1;
                        }
                    }, 1024);
                    int nnew = 0;
                    for (int f = 0; f < nT; ++f)
                        if (bad[f])
                            for (int c = 0; c < 3; ++c) {
                                const int v = D.T[f].v[c];
                                if (!revert[v] && vt[v] != 255) { revert[v] = 1; ++nnew; }
                            }
                    if (nnew == 0) break;
                }
            }
            for (int v : valid) if (!revert[v] && vt[v] != 255) D.P[v] = newP[v];
        }
        count_folds("smooth", it);
        if (ctx.verbose)
            ctx.logf("remesh: итерация %d: граней %d (разбиений %d, схлопываний %d, флипов %d) (%.2f с)", it,
                     D.alive_faces, stats.splits, stats.collapses, stats.flips, now_seconds() - t0);
    }

    // ------------------------------------------------------------ выгрузка
    std::vector<int> remap(D.P.size(), -1);
    out = TriMesh();
    out_feat.clear();
    std::unordered_set<u64> fkeys;
    for (size_t f = 0; f < D.T.size(); ++f) {
        if (D.tdead[f]) continue;
        for (int i = 0; i < 3; ++i) {
            int v = D.T[f].v[i];
            if (remap[v] < 0) { remap[v] = (int)out.P.size(); out.P.push_back(D.P[v]); }
            out.F.push_back(remap[v]);
        }
        for (int i = 0; i < 3; ++i)
            if (D.T[f].ef[i]) {
                int a = remap[D.T[f].v[i]], b = remap[D.T[f].v[(i + 1) % 3]];
                if (a > b) std::swap(a, b);
                fkeys.insert(((u64)(uint32_t)a << 32) | (uint32_t)b);
            }
    }
    stats.faces_out = (int)(out.F.size() / 3);
    std::string err;
    if (!out.build(&err)) {
        ctx.logf("remesh: результат не прошёл проверку топологии: %s", err.c_str());
        return false;
    }
    out_feat.assign(3 * (size_t)out.nF(), 0);
    for (int f = 0; f < out.nF(); ++f)
        for (int i = 0; i < 3; ++i) {
            int a = out.F[3 * f + i], b = out.F[3 * f + (i + 1) % 3];
            if (a > b) std::swap(a, b);
            if (out.FF[3 * f + i] < 0 || fkeys.count(((u64)(uint32_t)a << 32) | (uint32_t)b)) out_feat[3 * f + i] = 1;
        }
    ctx.logf("remesh: %d -> %d граней за %.2f с (разбиений %d, схлопываний %d, флипов %d)", stats.faces_in,
             stats.faces_out, now_seconds() - t0, stats.splits, stats.collapses, stats.flips);
    return true;
}

}  // namespace qf
