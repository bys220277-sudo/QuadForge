#include "qf_features.h"

#include <unordered_map>
#include <unordered_set>

namespace qf {

static inline u64 ekey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((u64)(uint32_t)a << 32) | (uint32_t)b;
}

void detect_features(const TriMesh& M, const FeatureOptions& opt, std::vector<uint8_t>& feat) {
    const int nf = M.nF();
    feat.assign(3 * (size_t)nf, 0);
    std::unordered_set<u64> user;
    for (int k = 0; k < opt.n_user; ++k) user.insert(ekey(opt.user_edges[2 * k], opt.user_edges[2 * k + 1]));
    const double cs = std::cos(opt.sharp_angle);
    parallel_each(nf, [&](i64 f) {
        for (int i = 0; i < 3; ++i) {
            int g = M.FF[3 * f + i];
            bool is = false;
            if (g < 0) {
                is = opt.boundary;
            } else {
                if (opt.detect_sharp && opt.sharp_angle > 0 && dot(M.FN[f], M.FN[g]) < cs) is = true;
                if (!is && !user.empty()) {
                    int a = M.F[3 * f + i], b = M.F[3 * f + (i + 1) % 3];
                    if (opt.origin) { a = (*opt.origin)[a]; b = (*opt.origin)[b]; }
                    if (user.count(ekey(a, b))) is = true;
                }
            }
            feat[3 * f + i] = is ? 1 : 0;
        }
    });
}

void field_constraints_from_features(const TriMesh& M, const std::vector<uint8_t>& feat,
                                     std::vector<C2>& fixed) {
    const int nf = M.nF();
    fixed.assign(nf, C2(0, 0));
    parallel_each(nf, [&](i64 f) {
        double best = -1;
        int bi = -1;
        for (int i = 0; i < 3; ++i) {
            if (!feat[3 * f + i]) continue;
            double len = norm(M.P[M.F[3 * f + (i + 1) % 3]] - M.P[M.F[3 * f + i]]);
            if (len > best) { best = len; bi = i; }
        }
        if (bi < 0) return;
        V3 e = M.P[M.F[3 * f + (bi + 1) % 3]] - M.P[M.F[3 * f + bi]];
        fixed[f] = polar(4.0 * angle_in_face(M, (int)f, e));
    });
}

int drop_close_parallel_features(const TriMesh& M, std::vector<uint8_t>& feat, double L, double fac) {
    struct E { int f, i, a, b; V3 m, d; };
    std::vector<E> es;
    for (int f = 0; f < M.nF(); ++f)
        for (int i = 0; i < 3; ++i) {
            if (!feat[3 * f + i] || M.FF[3 * f + i] < 0) continue;  // границы не трогаем
            int g = M.FF[3 * f + i];
            if (g < f) continue;
            int a = M.F[3 * f + i], b = M.F[3 * f + (i + 1) % 3];
            V3 d = M.P[b] - M.P[a];
            es.push_back({f, i, a, b, (M.P[a] + M.P[b]) * 0.5, normalized(d)});
        }
    const double cell = fac * L;
    std::unordered_map<u64, std::vector<int>> grid;
    if (es.empty()) return 0;
    auto key = [&](int x, int y, int z) {
        return ((u64)(uint32_t)(x + (1 << 20)) << 42) | ((u64)(uint32_t)(y + (1 << 20)) << 21) | (u64)(uint32_t)(z + (1 << 20));
    };
    auto cx = [&](double v) { return (int)std::floor(v / cell); };
    for (int k = 0; k < (int)es.size(); ++k) grid[key(cx(es[k].m.x), cx(es[k].m.y), cx(es[k].m.z))].push_back(k);
    std::vector<uint8_t> drop(es.size(), 0);
    {
        // цепочки фичеров: компоненты рёбер через вершины фичер-степени 2
        std::unordered_map<int, std::vector<int>> vedges;
        for (int k = 0; k < (int)es.size(); ++k) { vedges[es[k].a].push_back(k); vedges[es[k].b].push_back(k); }
        std::vector<int> chain(es.size(), -1);
        int nch = 0;
        for (int k0 = 0; k0 < (int)es.size(); ++k0) {
            if (chain[k0] >= 0) continue;
            std::vector<int> st = {k0};
            chain[k0] = nch;
            while (!st.empty()) {
                int k = st.back(); st.pop_back();
                for (int v : {es[k].a, es[k].b}) {
                    auto& ve = vedges[v];
                    if (ve.size() != 2) continue;  // угол/узел — граница цепочки
                    for (int o : ve) if (chain[o] < 0) { chain[o] = nch; st.push_back(o); }
                }
            }
            ++nch;
        }
        std::vector<double> clen(nch, 0.0);
        for (int k = 0; k < (int)es.size(); ++k) clen[chain[k]] += norm(M.P[es[k].b] - M.P[es[k].a]);
        // партнёры: для каждого ребра — цепочки близких параллельных рёбер
        std::vector<std::vector<int>> partners(es.size());
        for (int k = 0; k < (int)es.size(); ++k) {
            const E& e = es[k];
            int X = cx(e.m.x), Y = cx(e.m.y), Z = cx(e.m.z);
            for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) for (int dz = -1; dz <= 1; ++dz) {
                auto it = grid.find(key(X + dx, Y + dy, Z + dz));
                if (it == grid.end()) continue;
                for (int o : it->second) {
                    if (chain[o] == chain[k]) continue;
                    const E& q = es[o];
                    if (std::fabs(dot(e.d, q.d)) < 0.9) continue;
                    V3 dv = q.m - e.m; double pd = norm(dv - e.d * dot(dv, e.d));
                    if (pd > 0.02 * L && pd < fac * L && std::fabs(dot(dv, e.d)) < 0.5 * L) partners[k].push_back(chain[o]);
                }
            }
        }
        std::vector<std::vector<int>> cedges(nch);
        for (int k = 0; k < (int)es.size(); ++k) cedges[chain[k]].push_back(k);
        std::vector<int> order(nch);
        for (int c = 0; c < nch; ++c) order[c] = c;
        std::sort(order.begin(), order.end(), [&](int a, int b) { return clen[a] > clen[b]; });
        std::vector<uint8_t> kept(nch, 0), dropc(nch, 0);
        for (int c : order) {
            double close = 0, tot = 0;
            for (int k : cedges[c]) {
                double l = norm(M.P[es[k].b] - M.P[es[k].a]);
                tot += l;
                for (int pc : partners[k]) if (kept[pc]) { close += l; break; }
            }
            if (tot > 0 && close > 0.3 * tot) dropc[c] = 1; else kept[c] = 1;
        }
        for (int k = 0; k < (int)es.size(); ++k) drop[k] = dropc[chain[k]];
    }
    int n = 0;
    for (int k = 0; k < (int)es.size(); ++k)
        if (drop[k]) {
            const E& e = es[k];
            int g = M.FF[3 * e.f + e.i], j = M.FE[3 * e.f + e.i];
            feat[3 * e.f + e.i] = 0;
            feat[3 * g + j] = 0;
            ++n;
        }
    return n;
}

int drop_short_feature_chains(const TriMesh& M, std::vector<uint8_t>& feat, double L, double min_len) {
    struct E { int f, i, a, b; };
    std::vector<E> es;
    std::unordered_map<int, int> deg;  // фичер-степень вершины (включая границы)
    for (int f = 0; f < M.nF(); ++f)
        for (int i = 0; i < 3; ++i) {
            if (!feat[3 * f + i]) continue;
            int g = M.FF[3 * f + i];
            if (g >= 0 && g < f) continue;
            int a = M.F[3 * f + i], b = M.F[3 * f + (i + 1) % 3];
            ++deg[a]; ++deg[b];
            if (g >= 0) es.push_back({f, i, a, b});
        }
    if (es.empty()) return 0;
    std::unordered_map<int, std::vector<int>> vedges;
    for (int k = 0; k < (int)es.size(); ++k) { vedges[es[k].a].push_back(k); vedges[es[k].b].push_back(k); }
    std::vector<int> chain(es.size(), -1);
    std::vector<double> clen;
    std::vector<uint8_t> cfree;  // все концы цепочки свободны (степень 1)
    for (int k0 = 0; k0 < (int)es.size(); ++k0) {
        if (chain[k0] >= 0) continue;
        const int c = (int)clen.size();
        clen.push_back(0.0);
        cfree.push_back(1);
        std::vector<int> st = {k0};
        chain[k0] = c;
        while (!st.empty()) {
            int k = st.back(); st.pop_back();
            clen[c] += norm(M.P[es[k].b] - M.P[es[k].a]);
            for (int v : {es[k].a, es[k].b}) {
                const int d = deg[v];
                if (d == 2) {
                    for (int o : vedges[v]) if (chain[o] < 0) { chain[o] = c; st.push_back(o); }
                    if ((int)vedges[v].size() < 2) cfree[c] = 0;  // продолжается граничным ребром
                } else if (d != 1) {
                    cfree[c] = 0;  // узел/угол — цепочка часть сети
                }
            }
        }
    }
    int n = 0;
    for (int k = 0; k < (int)es.size(); ++k) {
        const int c = chain[k];
        if (!cfree[c] || clen[c] >= min_len * L) continue;
        const E& e = es[k];
        int g = M.FF[3 * e.f + e.i], j = M.FE[3 * e.f + e.i];
        feat[3 * e.f + e.i] = 0;
        feat[3 * g + j] = 0;
        ++n;
    }
    return n;
}

}  // namespace qf
