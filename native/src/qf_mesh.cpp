#include "qf_mesh.h"

#include <map>
#include <unordered_map>
#include <unordered_set>

namespace qf {

namespace {
inline u64 edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((u64)(uint32_t)a << 32) | (uint32_t)b;
}
}  // namespace

void TriMesh::compute_face_geometry() {
    const int nf = nF();
    FN.assign(nf, V3());
    FA.assign(nf, 0.0);
    parallel_each(nf, [&](i64 f) {
        const V3& a = P[F[3 * f]];
        const V3& b = P[F[3 * f + 1]];
        const V3& c = P[F[3 * f + 2]];
        V3 n = cross(b - a, c - a);
        double len = norm(n);
        FA[f] = 0.5 * len;
        FN[f] = len > 1e-300 ? n / len : V3(0, 0, 1);
    });
}

bool TriMesh::build_adjacency(std::string* err) {
    const int nf = nF();
    FF.assign(3 * (size_t)nf, -1);
    FE.assign(3 * (size_t)nf, -1);
    struct HE { u64 key; int slot; };
    std::vector<HE> hes;
    hes.reserve(3 * (size_t)nf);
    for (int f = 0; f < nf; ++f)
        for (int i = 0; i < 3; ++i)
            hes.push_back({edge_key(F[3 * f + i], F[3 * f + (i + 1) % 3]), 3 * f + i});
    std::sort(hes.begin(), hes.end(), [](const HE& a, const HE& b) {
        return a.key != b.key ? a.key < b.key : a.slot < b.slot;
    });
    size_t i = 0;
    while (i < hes.size()) {
        size_t j = i;
        while (j < hes.size() && hes[j].key == hes[i].key) ++j;
        size_t cnt = j - i;
        if (cnt > 2) {
            if (err) *err = "неманифольдное ребро (больше двух граней)";
            return false;
        }
        if (cnt == 2) {
            int s0 = hes[i].slot, s1 = hes[i + 1].slot;
            FF[s0] = s1 / 3;
            FE[s0] = s1 % 3;
            FF[s1] = s0 / 3;
            FE[s1] = s0 % 3;
        }
        i = j;
    }
    return true;
}

bool TriMesh::orient_consistently(std::string* err) {
    const int nf = nF();
    std::vector<uint8_t> visited(nf, 0), flip(nf, 0);
    std::vector<int> stack;
    for (int s = 0; s < nf; ++s) {
        if (visited[s]) continue;
        visited[s] = 1;
        stack.push_back(s);
        while (!stack.empty()) {
            int f = stack.back();
            stack.pop_back();
            for (int i = 0; i < 3; ++i) {
                int g = FF[3 * f + i];
                if (g < 0) continue;
                // ребро f: (a,b). Согласованное соседство: в g ребро идёт b->a.
                int a = F[3 * f + i], b = F[3 * f + (i + 1) % 3];
                int j = FE[3 * f + i];
                int ga = F[3 * g + j], gb = F[3 * g + (j + 1) % 3];
                bool fa = flip[f], fg = flip[g];
                // эффективное направление ребра с учётом флипов
                int ea = fa ? b : a, eb = fa ? a : b;
                int na = fg ? gb : ga, nb = fg ? ga : gb;
                bool consistent = (ea == nb && eb == na);
                if (!visited[g]) {
                    visited[g] = 1;
                    if (!consistent) flip[g] ^= 1;
                    stack.push_back(g);
                } else if (!consistent) {
                    if (err) *err = "поверхность неориентируема";
                    return false;
                }
            }
        }
    }
    bool any = false;
    for (int f = 0; f < nf; ++f)
        if (flip[f]) {
            std::swap(F[3 * f + 1], F[3 * f + 2]);
            any = true;
        }
    if (any) return build_adjacency(err);
    return true;
}

bool TriMesh::build_fans(std::string* err) {
    const int nv = nV(), nf = nF();
    // Все углы по вершинам.
    std::vector<int> cnt(nv + 1, 0);
    for (int k = 0; k < 3 * nf; ++k) ++cnt[F[k] + 1];
    for (int v = 0; v < nv; ++v) cnt[v + 1] += cnt[v];
    std::vector<int> corners(3 * (size_t)nf);
    {
        std::vector<int> fill(cnt.begin(), cnt.end() - 1);
        for (int k = 0; k < 3 * nf; ++k) corners[fill[F[k]]++] = k;  // k = 3f+i
    }
    vfan_off.assign(nv + 1, 0);
    vfan_face.assign(3 * (size_t)nf, 0);
    vfan_corner.assign(3 * (size_t)nf, 0);
    vbound.assign(nv, 0);
    std::atomic<int> bad{0};
    // Порядок обхода: next(f,i) = пересечь ребро i (от v), попасть в соседа.
    parallel_each(nv, [&](i64 v) {
        int b = cnt[v], e = cnt[v + 1];
        int n = e - b;
        if (n == 0) return;
        // найти начало: идём назад (prev) пока не граница/цикл
        int s = corners[b];
        int f = s / 3, i = s % 3;
        int start_f = f, start_i = i;
        bool closed = false;
        for (int step = 0; step <= n; ++step) {
            int ep = (i + 2) % 3;  // ребро, оканчивающееся в v
            int g = FF[3 * f + ep];
            if (g < 0) { start_f = f; start_i = i; closed = false; break; }
            int j = FE[3 * f + ep];
            f = g; i = j;  // в g это ребро идёт от v, угол v = j
            if (f == s / 3 && i == s % 3) { closed = true; break; }
            if (step == n) { bad = 1; return; }
        }
        if (closed) { start_f = s / 3; start_i = s % 3; }
        // обход вперёд
        int out = cnt[v];
        f = start_f; i = start_i;
        int walked = 0;
        for (;;) {
            vfan_face[out + walked] = f;
            vfan_corner[out + walked] = (uint8_t)i;
            ++walked;
            if (walked > n) { bad = 1; return; }
            int g = FF[3 * f + i];
            if (g < 0) break;
            int j = FE[3 * f + i];
            int ni = (j + 1) % 3;
            f = g; i = ni;
            if (f == start_f && i == start_i) break;
        }
        if (walked != n) { bad = 1; return; }
        vbound[v] = closed ? 0 : 1;
    }, 512);
    if (bad) {
        if (err) *err = "неманифольдная вершина";
        return false;
    }
    vfan_off = cnt;
    return true;
}

double TriMesh::total_area() const {
    double s = 0;
    for (double a : FA) s += a;
    return s;
}

bool TriMesh::build(std::string* err) {
    if (!build_adjacency(err)) return false;
    if (!orient_consistently(err)) return false;
    compute_face_geometry();
    if (!build_fans(err)) return false;
    return true;
}

// ---------------------------------------------------------------------------
// Очистка
// ---------------------------------------------------------------------------
// Сварка вершин, совпадающих в пределах tol (пространственный хеш). rep[i] — представитель.
static int weld_vertices(const std::vector<V3>& P, double tol, std::vector<int>& rep) {
    const int n = (int)P.size();
    rep.resize(n);
    for (int i = 0; i < n; ++i) rep[i] = i;
    if (tol <= 0) return 0;
    auto cell = [&](double v) { return (int64_t)std::floor(v / tol); };
    auto key = [](int64_t x, int64_t y, int64_t z) {
        return ((u64)(x & 0x1FFFFF) << 42) | ((u64)(y & 0x1FFFFF) << 21) | (u64)(z & 0x1FFFFF);
    };
    std::unordered_map<u64, std::vector<int>> grid;
    grid.reserve(n * 2);
    const double t2 = tol * tol;
    int nw = 0;
    for (int i = 0; i < n; ++i) {
        int64_t X = cell(P[i].x), Y = cell(P[i].y), Z = cell(P[i].z);
        int found = -1;
        for (int dx = -1; dx <= 1 && found < 0; ++dx)
            for (int dy = -1; dy <= 1 && found < 0; ++dy)
                for (int dz = -1; dz <= 1 && found < 0; ++dz) {
                    auto it = grid.find(key(X + dx, Y + dy, Z + dz));
                    if (it == grid.end()) continue;
                    for (int j : it->second)
                        if (norm2(P[j] - P[i]) <= t2) { found = j; break; }
                }
        if (found >= 0) { rep[i] = found; ++nw; }
        else grid[key(X, Y, Z)].push_back(i);
    }
    return nw;
}

// Сшивка T-стыков: граничная вершина, лежащая на внутренности чужого граничного ребра,
// разбивает грань этого ребра на две. Повторяем, пока есть что сшивать.
static int stitch_tjunctions(const std::vector<V3>& P, std::vector<int>& F, double tol) {
    int total = 0;
    for (int round = 0; round < 8; ++round) {
        const int nf = (int)(F.size() / 3);
        std::unordered_map<u64, int> cnt;
        cnt.reserve(nf * 3);
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) ++cnt[edge_key(F[3 * f + i], F[3 * f + (i + 1) % 3])];
        struct BE { int f, i; };
        std::vector<BE> be;
        std::vector<int> bv;
        std::unordered_set<int> bvs;
        double sumlen = 0;
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) {
                int a = F[3 * f + i], b = F[3 * f + (i + 1) % 3];
                if (cnt[edge_key(a, b)] != 1) continue;
                be.push_back({f, i});
                sumlen += norm(P[b] - P[a]);
                if (bvs.insert(a).second) bv.push_back(a);
                if (bvs.insert(b).second) bv.push_back(b);
            }
        if (be.empty()) break;
        const double cs = std::max(sumlen / be.size(), tol * 4);
        auto cell = [&](double v) { return (int64_t)std::floor(v / cs); };
        auto key = [](int64_t x, int64_t y, int64_t z) {
            return ((u64)(x & 0x1FFFFF) << 42) | ((u64)(y & 0x1FFFFF) << 21) | (u64)(z & 0x1FFFFF);
        };
        std::unordered_map<u64, std::vector<int>> grid;
        for (int v : bv) grid[key(cell(P[v].x), cell(P[v].y), cell(P[v].z))].push_back(v);
        // для каждой грани — не больше одного разбиения за раунд
        std::vector<int> split_at(3 * (size_t)nf, -1);
        std::vector<uint8_t> fbusy(nf, 0);
        int nsplit = 0;
        for (const BE& e : be) {
            if (fbusy[e.f]) continue;
            int a = F[3 * e.f + e.i], b = F[3 * e.f + (e.i + 1) % 3];
            V3 pa = P[a], d = P[b] - P[a];
            double l2 = norm2(d);
            if (l2 <= tol * tol) continue;
            int best = -1;
            double bt = 2;
            int64_t x0 = cell(std::min(pa.x, P[b].x) - tol), x1 = cell(std::max(pa.x, P[b].x) + tol);
            int64_t y0 = cell(std::min(pa.y, P[b].y) - tol), y1 = cell(std::max(pa.y, P[b].y) + tol);
            int64_t z0 = cell(std::min(pa.z, P[b].z) - tol), z1 = cell(std::max(pa.z, P[b].z) + tol);
            if ((x1 - x0 + 1) * (y1 - y0 + 1) * (z1 - z0 + 1) > 4096) continue;
            for (int64_t x = x0; x <= x1; ++x)
                for (int64_t y = y0; y <= y1; ++y)
                    for (int64_t z = z0; z <= z1; ++z) {
                        auto it = grid.find(key(x, y, z));
                        if (it == grid.end()) continue;
                        for (int w : it->second) {
                            if (w == a || w == b) continue;
                            double t = dot(P[w] - pa, d) / l2;
                            if (t <= 1e-6 || t >= 1 - 1e-6) continue;
                            if (norm2(pa + d * t - P[w]) > tol * tol) continue;
                            if (t < bt) { bt = t; best = w; }  // ближайшая к a — остальные в след. раундах
                        }
                    }
            if (best < 0) continue;
            split_at[3 * e.f + e.i] = best;
            fbusy[e.f] = 1;
            ++nsplit;
        }
        if (!nsplit) break;
        for (int f = 0; f < nf; ++f) {
            if (!fbusy[f]) continue;
            for (int i = 0; i < 3; ++i) {
                int w = split_at[3 * f + i];
                if (w < 0) continue;
                int a = F[3 * f + i], b = F[3 * f + (i + 1) % 3], c = F[3 * f + (i + 2) % 3];
                F[3 * f] = a; F[3 * f + 1] = w; F[3 * f + 2] = c;
                F.push_back(w); F.push_back(b); F.push_back(c);
                break;
            }
        }
        total += nsplit;
    }
    return total;
}

CleanResult clean_triangle_soup(std::vector<V3>& P, std::vector<int>& F) {
    CleanResult res;
    const int nv0 = (int)P.size();
    int nf = (int)(F.size() / 3);

    // 1. вырожденные и дубликаты
    std::vector<int> keep;
    keep.reserve(nf);
    std::unordered_set<u64> seen;
    seen.reserve(nf * 2);
    // bbox для относительных допусков
    V3 mn(1e300, 1e300, 1e300), mx(-1e300, -1e300, -1e300);
    for (auto& p : P)
        for (int k = 0; k < 3; ++k) { mn[k] = std::min(mn[k], p[k]); mx[k] = std::max(mx[k], p[k]); }
    double diag = norm(mx - mn);
    double amin = diag * diag * 1e-18;

    // 0. сварка совпадающих вершин (суп STL, расщеплённые нормали, дубли после булевых)
    //    и сшивка T-стыков (независимо триангулированные грани CAD)
    {
        const double wtol = 1e-6 * diag;
        res.welded_vertices = weld_vertices(P, wtol, res.weld_rep);
        if (res.welded_vertices)
            for (auto& v : F)
                if (v >= 0 && v < nv0) v = res.weld_rep[v];
        res.tjunction_splits = stitch_tjunctions(P, F, 1e-5 * diag);
        nf = (int)(F.size() / 3);
    }
    std::vector<int> F2;
    F2.reserve(F.size());
    for (int f = 0; f < nf; ++f) {
        int a = F[3 * f], b = F[3 * f + 1], c = F[3 * f + 2];
        if (a < 0 || b < 0 || c < 0 || a >= nv0 || b >= nv0 || c >= nv0 || a == b || b == c || a == c) {
            ++res.removed_degenerate;
            continue;
        }
        double area = 0.5 * norm(cross(P[b] - P[a], P[c] - P[a]));
        if (area <= amin) { ++res.removed_degenerate; continue; }
        int s[3] = {a, b, c};
        std::sort(s, s + 3);
        u64 key = ((u64)s[0] * 1000003ULL) ^ ((u64)s[1] * 998244353ULL << 1) ^ ((u64)s[2] * 12582917ULL << 2);
        // ключ приблизительный; точная проверка ниже
        (void)key;
        F2.push_back(a); F2.push_back(b); F2.push_back(c);
    }
    nf = (int)(F2.size() / 3);
    {
        std::map<std::array<int, 3>, int> uniq;
        std::vector<int> F3;
        F3.reserve(F2.size());
        for (int f = 0; f < nf; ++f) {
            std::array<int, 3> k = {F2[3 * f], F2[3 * f + 1], F2[3 * f + 2]};
            std::sort(k.begin(), k.end());
            if (uniq.count(k)) { ++res.removed_duplicate; continue; }
            uniq[k] = f;
            F3.push_back(F2[3 * f]); F3.push_back(F2[3 * f + 1]); F3.push_back(F2[3 * f + 2]);
        }
        F2.swap(F3);
        nf = (int)(F2.size() / 3);
    }

    // 2. неманифольдные рёбра: оставляем две грани с наибольшей площадью
    {
        std::unordered_map<u64, std::vector<int>> edges;
        edges.reserve(nf * 2);
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) edges[edge_key(F2[3 * f + i], F2[3 * f + (i + 1) % 3])].push_back(f);
        std::vector<uint8_t> dead(nf, 0);
        for (auto& kv : edges) {
            auto& fl = kv.second;
            if (fl.size() <= 2) continue;
            std::sort(fl.begin(), fl.end(), [&](int x, int y) {
                double ax = norm(cross(P[F2[3 * x + 1]] - P[F2[3 * x]], P[F2[3 * x + 2]] - P[F2[3 * x]]));
                double ay = norm(cross(P[F2[3 * y + 1]] - P[F2[3 * y]], P[F2[3 * y + 2]] - P[F2[3 * y]]));
                return ax > ay;
            });
            for (size_t k = 2; k < fl.size(); ++k)
                if (!dead[fl[k]]) { dead[fl[k]] = 1; ++res.removed_nonmanifold; }
        }
        if (res.removed_nonmanifold) {
            std::vector<int> F3;
            for (int f = 0; f < nf; ++f)
                if (!dead[f]) { F3.push_back(F2[3 * f]); F3.push_back(F2[3 * f + 1]); F3.push_back(F2[3 * f + 2]); }
            F2.swap(F3);
            nf = (int)(F2.size() / 3);
        }
    }

    // 3. неманифольдные вершины: дублируем вершину для каждой компоненты веера
    std::vector<V3> P2 = P;
    std::vector<int> origin(P.size());
    for (size_t i = 0; i < P.size(); ++i) origin[i] = (int)i;
    {
        // ребро -> два соседа
        std::unordered_map<u64, std::array<int, 2>> emap;
        emap.reserve(nf * 2);
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) {
                u64 k = edge_key(F2[3 * f + i], F2[3 * f + (i + 1) % 3]);
                auto it = emap.find(k);
                if (it == emap.end()) emap[k] = {f, -1};
                else it->second[1] = f;
            }
        const std::vector<int> F0 = F2;  // топология до переименований (emap построен по ней)
        // грани по вершинам
        std::vector<std::vector<int>> vf(P.size());
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) vf[F2[3 * f + i]].push_back(f);
        for (int v = 0; v < (int)P.size(); ++v) {
            auto& fl = vf[v];
            if (fl.size() < 2) continue;
            // компоненты связности граней вокруг v через рёбра, содержащие v
            std::unordered_map<int, int> comp;
            int ncomp = 0;
            for (int s : fl) {
                if (comp.count(s)) continue;
                std::vector<int> st = {s};
                comp[s] = ncomp;
                while (!st.empty()) {
                    int f = st.back();
                    st.pop_back();
                    for (int i = 0; i < 3; ++i) {
                        int a = F0[3 * f + i], b = F0[3 * f + (i + 1) % 3];
                        if (a != v && b != v) continue;
                        const auto& nb = emap.at(edge_key(a, b));
                        int g = nb[0] == f ? nb[1] : nb[0];
                        if (g >= 0 && !comp.count(g)) { comp[g] = ncomp; st.push_back(g); }
                    }
                }
                ++ncomp;
            }
            if (ncomp <= 1) continue;
            std::vector<int> newv(ncomp, v);
            for (int c = 1; c < ncomp; ++c) {
                newv[c] = (int)P2.size();
                P2.push_back(P[v]);
                origin.push_back(origin[v]);
                ++res.split_vertices;
            }
            for (int f : fl) {
                int c = comp[f];
                if (c == 0) continue;
                for (int i = 0; i < 3; ++i)
                    if (F2[3 * f + i] == v) F2[3 * f + i] = newv[c];
            }
        }
    }

    // 4. удалить неиспользуемые вершины
    std::vector<int> remap(P2.size(), -1);
    std::vector<V3> P3;
    std::vector<int> orig3;
    for (int f = 0; f < nf; ++f)
        for (int i = 0; i < 3; ++i) {
            int v = F2[3 * f + i];
            if (remap[v] < 0) {
                remap[v] = (int)P3.size();
                P3.push_back(P2[v]);
                orig3.push_back(origin[v]);
            }
            F2[3 * f + i] = remap[v];
        }
    P.swap(P3);
    F.swap(F2);
    res.vert_origin.swap(orig3);
    return res;
}

// ---------------------------------------------------------------------------
// Нормали вершин
// ---------------------------------------------------------------------------
void compute_vertex_normals(const TriMesh& M, std::vector<V3>& VN) {
    VN.assign(M.nV(), V3());
    parallel_each(M.nV(), [&](i64 v) {
        V3 n;
        M.for_fan((int)v, [&](int f, int c) {
            const V3& p = M.P[M.F[3 * f + c]];
            V3 a = M.P[M.F[3 * f + (c + 1) % 3]] - p;
            V3 b = M.P[M.F[3 * f + (c + 2) % 3]] - p;
            double la = norm(a), lb = norm(b);
            if (la < 1e-300 || lb < 1e-300) return;
            double ang = std::acos(std::max(-1.0, std::min(1.0, dot(a, b) / (la * lb))));
            n += M.FN[f] * ang;
        });
        VN[v] = normalized(n, V3(0, 0, 1));
    }, 512);
}

// ---------------------------------------------------------------------------
// Кривизна по граням (Rusinkiewicz 2004), сглаживание по соседям.
// ---------------------------------------------------------------------------
void compute_face_curvature(const TriMesh& M, CurvatureData& out, int smooth_iters) {
    const int nf = M.nF();
    std::vector<V3> VN;
    compute_vertex_normals(M, VN);
    // тензор (a b; b c) в базисе грани -> 3x3 в мировых координатах (6 чисел)
    struct T6 { double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0; };
    std::vector<T6> T(nf), T2(nf);
    parallel_each(nf, [&](i64 f) {
        V3 e1, e2;
        M.face_frame((int)f, e1, e2);
        const int* v = &M.F[3 * f];
        // рёбра и разности нормалей (в базисе грани)
        double A[6][3];
        double rhs[6];
        for (int i = 0; i < 3; ++i) {
            V3 e = M.P[v[(i + 2) % 3]] - M.P[v[(i + 1) % 3]];
            V3 dn = VN[v[(i + 2) % 3]] - VN[v[(i + 1) % 3]];
            double eu = dot(e, e1), ev = dot(e, e2);
            double nu = dot(dn, e1), nv = dot(dn, e2);
            // II = [[a b],[b c]]: II*e = dn
            A[2 * i][0] = eu; A[2 * i][1] = ev; A[2 * i][2] = 0; rhs[2 * i] = nu;
            A[2 * i + 1][0] = 0; A[2 * i + 1][1] = eu; A[2 * i + 1][2] = ev; rhs[2 * i + 1] = nv;
        }
        double ata[3][3] = {{0}}, atb[3] = {0};
        for (int r = 0; r < 6; ++r)
            for (int i = 0; i < 3; ++i) {
                atb[i] += A[r][i] * rhs[r];
                for (int j = 0; j < 3; ++j) ata[i][j] += A[r][i] * A[r][j];
            }
        // решить 3x3 (Крамер)
        auto det3 = [](double m[3][3]) {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        };
        double d = det3(ata);
        double a = 0, b = 0, c = 0;
        if (std::fabs(d) > 1e-300) {
            double m[3][3];
            for (int col = 0; col < 3; ++col) {
                std::memcpy(m, ata, sizeof(m));
                for (int r = 0; r < 3; ++r) m[r][col] = atb[r];
                double s = det3(m) / d;
                if (col == 0) a = s; else if (col == 1) b = s; else c = s;
            }
        }
        T6 t;
        t.xx = a * e1.x * e1.x + 2 * b * e1.x * e2.x + c * e2.x * e2.x;
        t.yy = a * e1.y * e1.y + 2 * b * e1.y * e2.y + c * e2.y * e2.y;
        t.zz = a * e1.z * e1.z + 2 * b * e1.z * e2.z + c * e2.z * e2.z;
        t.xy = a * e1.x * e1.y + b * (e1.x * e2.y + e1.y * e2.x) + c * e2.x * e2.y;
        t.xz = a * e1.x * e1.z + b * (e1.x * e2.z + e1.z * e2.x) + c * e2.x * e2.z;
        t.yz = a * e1.y * e1.z + b * (e1.y * e2.z + e1.z * e2.y) + c * e2.y * e2.z;
        T[f] = t;
    });
    // сглаживание: усреднение с соседями (по площади)
    for (int it = 0; it < smooth_iters; ++it) {
        parallel_each(nf, [&](i64 f) {
            T6 s;
            double wsum = 0;
            auto add = [&](const T6& t, double w) {
                s.xx += w * t.xx; s.xy += w * t.xy; s.xz += w * t.xz;
                s.yy += w * t.yy; s.yz += w * t.yz; s.zz += w * t.zz;
                wsum += w;
            };
            add(T[f], M.FA[f]);
            for (int i = 0; i < 3; ++i) {
                int g = M.FF[3 * f + i];
                if (g >= 0) add(T[g], M.FA[g]);
            }
            double inv = 1.0 / wsum;
            s.xx *= inv; s.xy *= inv; s.xz *= inv; s.yy *= inv; s.yz *= inv; s.zz *= inv;
            T2[f] = s;
        });
        T.swap(T2);
    }
    out.k1.assign(nf, 0);
    out.k2.assign(nf, 0);
    out.theta.assign(nf, 0);
    parallel_each(nf, [&](i64 f) {
        V3 e1, e2;
        M.face_frame((int)f, e1, e2);
        const T6& t = T[f];
        auto apply = [&](const V3& v) {
            return V3(t.xx * v.x + t.xy * v.y + t.xz * v.z,
                      t.xy * v.x + t.yy * v.y + t.yz * v.z,
                      t.xz * v.x + t.yz * v.y + t.zz * v.z);
        };
        double a = dot(e1, apply(e1)), b = dot(e1, apply(e2)), c = dot(e2, apply(e2));
        double mean = 0.5 * (a + c);
        double diff = std::sqrt(0.25 * (a - c) * (a - c) + b * b);
        out.k1[f] = mean + diff;
        out.k2[f] = mean - diff;
        out.theta[f] = 0.5 * std::atan2(2 * b, a - c);  // направление k1 (макс. кривизна)
    });
}

}  // namespace qf
