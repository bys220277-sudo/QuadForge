#include "qf_extract.h"

#include <unordered_map>

namespace qf {

namespace {

constexpr double kEps = 1e-8;       // допуск классификации по значению iso-координаты
constexpr double kLatEps = 1e-7;    // допуск «целая координата»

struct NodeStore {
    std::vector<V3> pos;
    std::vector<V3> nrm;
    std::vector<uint8_t> lattice;
    std::vector<uint8_t> feat;  // на фичер-ребре/границе
    std::vector<uint8_t> bnd;   // на границе меша
    int add(const V3& p, const V3& n, bool lat, bool ft, bool bd) {
        pos.push_back(p);
        nrm.push_back(n);
        lattice.push_back(lat ? 1 : 0);
        feat.push_back(ft ? 1 : 0);
        bnd.push_back(bd ? 1 : 0);
        return (int)pos.size() - 1;
    }
};

struct Builder {
    const TriMesh& M;
    const std::vector<double>& uv;
    const std::vector<uint8_t>* feature;
    std::vector<V3> VN;
    NodeStore nodes;
    std::vector<int> corner_node;
    std::vector<uint8_t> vfeat;
    std::unordered_map<u64, int> edge_map;
    std::vector<double> node_t;  // для узлов на рёбрах (проверка близости)
    std::unordered_map<u64, int> interior_map;
    std::unordered_map<int, std::vector<std::pair<double, int>>> slot_nodes;
    std::vector<std::pair<int, int>> edges;
    struct Collinear { int slot, va, vb; };
    std::vector<Collinear> collinear;
    std::unordered_map<int, uint8_t> collinear_seen;
    int skipped = 0;

    Builder(const TriMesh& m, const std::vector<double>& u, const std::vector<uint8_t>* ft)
        : M(m), uv(u), feature(ft) {
        compute_vertex_normals(M, VN);
        corner_node.assign(M.nV(), -1);
        vfeat.assign(M.nV(), 0);
        if (feature)
            for (int f = 0; f < M.nF(); ++f)
                for (int i = 0; i < 3; ++i)
                    if ((*feature)[3 * f + i]) {
                        vfeat[M.F[3 * f + i]] = 1;
                        vfeat[M.F[3 * f + (i + 1) % 3]] = 1;
                    }
    }

    int slot_of(int f, int i) const {
        int g = M.FF[3 * f + i];
        if (g < 0) return 3 * f + i;
        return std::min(3 * f + i, 3 * g + M.FE[3 * f + i]);
    }
    bool is_feature_edge(int f, int i) const {
        return feature && (*feature)[3 * f + i];
    }

    int corner(int v, bool lat) {
        int& id = corner_node[v];
        if (id < 0) {
            id = nodes.add(M.P[v], VN[v], lat, vfeat[v] != 0, M.vbound[v] != 0);
        } else if (lat) {
            nodes.lattice[id] = 1;
        }
        return id;
    }

    // узел на ребре (f,i) в параметре t от меньшей вершины к большей
    int edge_node(int f, int i, double t, bool lat) {
        int a = M.F[3 * f + i], b = M.F[3 * f + (i + 1) % 3];
        int lo = std::min(a, b), hi = std::max(a, b);
        int slot = slot_of(f, i);
        int64_t tq = (int64_t)(t * 1e9);
        for (int64_t d = -1; d <= 1; ++d) {
            u64 key = ((u64)slot << 32) | (u64)(uint32_t)(tq + d);
            auto it = edge_map.find(key);
            if (it != edge_map.end() && std::fabs(node_t[it->second] - t) < 3e-9) {
                if (lat) nodes.lattice[it->second] = 1;
                return it->second;
            }
        }
        V3 p = M.P[lo] * (1.0 - t) + M.P[hi] * t;
        int g = M.FF[3 * f + i];
        V3 n = M.FN[f];
        if (g >= 0) n = normalized(M.FN[f] + M.FN[g], M.FN[f]);
        bool ft = g < 0 || is_feature_edge(f, i);
        int id = nodes.add(p, n, lat, ft, g < 0);
        node_t.resize(nodes.pos.size(), 0.0);
        node_t[id] = t;
        edge_map[((u64)slot << 32) | (u64)(uint32_t)tq] = id;
        slot_nodes[slot].push_back({t, id});
        return id;
    }

    int interior_node(int f, int iu, int iv, const V3& p) {
        u64 key = ((u64)f << 32) | ((u64)(uint16_t)(iu + 32768) << 16) | (u64)(uint16_t)(iv + 32768);
        auto it = interior_map.find(key);
        if (it != interior_map.end()) return it->second;
        int id = nodes.add(p, M.FN[f], true, false, false);
        interior_map[key] = id;
        return id;
    }

    void add_edge(int a, int b) {
        if (a != b) edges.push_back({std::min(a, b), std::max(a, b)});
    }

    void process_triangle(int f) {
        const double* q = &uv[6 * (size_t)f];
        const int* v = &M.F[3 * (size_t)f];
        for (int axis = 0; axis < 2; ++axis) {
            double w0 = q[axis], w1 = q[2 + axis], w2 = q[4 + axis];
            double wmin = std::min(w0, std::min(w1, w2)), wmax = std::max(w0, std::max(w1, w2));
            if (wmax - wmin > 5000.0) { ++skipped; return; }
            long cmin = (long)std::ceil(wmin - kEps), cmax = (long)std::floor(wmax + kEps);
            for (long c = cmin; c <= cmax; ++c) {
                double w[3] = {w0 - c, w1 - c, w2 - c};
                int s[3];
                int nz = 0;
                for (int i = 0; i < 3; ++i) {
                    s[i] = std::fabs(w[i]) < kEps ? 0 : (w[i] > 0 ? 1 : -1);
                    nz += s[i] == 0;
                }
                if (nz == 3) continue;
                if (nz == 2) {
                    int za = -1, zb = -1;
                    for (int i = 0; i < 3; ++i)
                        if (s[i] == 0) { if (za < 0) za = i; else zb = i; }
                    int ei = (zb == (za + 1) % 3) ? za : zb;  // ребро между za и zb
                    int slot = slot_of(f, ei);
                    if (!collinear_seen.count(slot)) {
                        collinear_seen[slot] = 1;
                        collinear.push_back({slot, v[ei], v[(ei + 1) % 3]});
                        // убедимся, что узлы-углы существуют и помечены решёточными при целой другой координате
                    }
                    // углы как узлы (решёточные ли — определит другая ось)
                    for (int k : {za, zb}) {
                        double other = q[2 * k + (1 - axis)];
                        corner(v[k], std::fabs(other - std::round(other)) < kLatEps);
                    }
                    continue;
                }
                // конечные точки сегмента
                struct Pt { int node; double o; V3 pos; };
                Pt pts[2];
                int np = 0;
                for (int i = 0; i < 3 && np < 3; ++i) {
                    if (s[i] == 0) {
                        double other = q[2 * i + (1 - axis)];
                        bool lat = std::fabs(other - std::round(other)) < kLatEps;
                        if (np < 2) pts[np] = {corner(v[i], lat), other, M.P[v[i]]};
                        ++np;
                    }
                }
                for (int i = 0; i < 3; ++i) {
                    int j = (i + 1) % 3;
                    if (s[i] * s[j] == -1) {
                        // ребро i: (v[i], v[j]); параметр от меньшей вершины
                        int lo = std::min(v[i], v[j]);
                        double wl = (lo == v[i]) ? w[i] : w[j];
                        double wh = (lo == v[i]) ? w[j] : w[i];
                        double t = wl / (wl - wh);
                        double ol = q[2 * (lo == v[i] ? i : j) + (1 - axis)];
                        double oh = q[2 * (lo == v[i] ? j : i) + (1 - axis)];
                        double other = ol * (1.0 - t) + oh * t;
                        bool lat = std::fabs(other - std::round(other)) < kLatEps;
                        int id = edge_node(f, i, t, lat);
                        V3 p = nodes.pos[id];
                        if (np < 2) pts[np] = {id, other, p};
                        ++np;
                    }
                }
                if (np != 2) continue;
                // внутренние решёточные точки
                double oa = pts[0].o, ob = pts[1].o;
                double lo_o = std::min(oa, ob), hi_o = std::max(oa, ob);
                long dmin = (long)std::floor(lo_o + kLatEps) + 1, dmax = (long)std::ceil(hi_o - kLatEps) - 1;
                int prev = pts[0].node;
                if (dmin <= dmax) {
                    std::vector<long> ds;
                    for (long d = dmin; d <= dmax; ++d) ds.push_back(d);
                    if (ob < oa) std::reverse(ds.begin(), ds.end());
                    for (long d : ds) {
                        double s_par = (d - oa) / (ob - oa);
                        V3 p = pts[0].pos * (1.0 - s_par) + pts[1].pos * s_par;
                        int iu = axis == 0 ? (int)c : (int)d;
                        int iv = axis == 0 ? (int)d : (int)c;
                        int id = interior_node(f, iu, iv, p);
                        add_edge(prev, id);
                        prev = id;
                    }
                }
                add_edge(prev, pts[1].node);
            }
        }
    }
};

}  // namespace

bool extract_quads(const TriMesh& M, const std::vector<double>& uv, const std::vector<uint8_t>* feature,
                   ExtractResult& out, Ctx& ctx) {
    double t0 = now_seconds();
    out = ExtractResult();
    Builder B(M, uv, feature);
    const int nf = M.nF();
    for (int f = 0; f < nf; ++f) {
        B.process_triangle(f);
        if ((f & 8191) == 8191) ctx.check_cancel();
    }
    out.skipped_triangles = B.skipped;
    // коллинеарные сегменты: цепочка вдоль ребра с разбиением по узлам на ребре
    for (const auto& cl : B.collinear) {
        int lo = std::min(cl.va, cl.vb), hi = std::max(cl.va, cl.vb);
        int nlo = B.corner(lo, false), nhi = B.corner(hi, false);
        auto it = B.slot_nodes.find(cl.slot);
        std::vector<std::pair<double, int>> list;
        if (it != B.slot_nodes.end()) list = it->second;
        std::sort(list.begin(), list.end());
        int prev = nlo;
        for (auto& pr : list) {
            B.add_edge(prev, pr.second);
            prev = pr.second;
        }
        B.add_edge(prev, nhi);
    }
    // уникальные рёбра
    std::sort(B.edges.begin(), B.edges.end());
    B.edges.erase(std::unique(B.edges.begin(), B.edges.end()), B.edges.end());
    const int nn = (int)B.nodes.pos.size();
    std::vector<std::vector<int>> adj(nn);
    for (auto& e : B.edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    ctx.logf("extract: %d узлов, %zu сегментов (%.2f с)", nn, B.edges.size(), now_seconds() - t0);

    // тупики и узлы степени != 2 без решёточного флага считаем вершинами
    std::vector<uint8_t> is_vert(nn, 0);
    for (int i = 0; i < nn; ++i) {
        int d = (int)adj[i].size();
        if (d == 0) continue;
        if (B.nodes.lattice[i]) is_vert[i] = 1;
        else if (d != 2) { is_vert[i] = 1; if (d == 1) ++out.dead_ends; }
    }
    // цепочки
    struct Chain { int a, b; int first_a, first_b; bool bnd; };
    std::vector<Chain> chains;
    {
        std::vector<std::pair<int, int>> used;  // (node, neighbor) направленные использованные
        std::unordered_map<u64, uint8_t> seen;
        auto dkey = [](int a, int b) { return ((u64)(uint32_t)a << 32) | (uint32_t)b; };
        for (int a = 0; a < nn; ++a) {
            if (!is_vert[a]) continue;
            for (int nb : adj[a]) {
                if (seen.count(dkey(a, nb))) continue;
                int prev = a, cur = nb;
                bool bnd = B.nodes.bnd[a] && B.nodes.bnd[nb];
                int last_prev = prev;
                while (!is_vert[cur]) {
                    // степень 2: идём дальше
                    int nxt = adj[cur][0] == prev ? adj[cur][1] : adj[cur][0];
                    bnd = bnd && B.nodes.bnd[cur];
                    prev = cur;
                    cur = nxt;
                    last_prev = prev;
                }
                seen[dkey(a, nb)] = 1;
                seen[dkey(cur, last_prev)] = 1;
                if (cur == a && last_prev == nb) continue;  // вырожденная петля из одного ребра
                chains.push_back({a, cur, nb, last_prev, bnd});
            }
        }
    }
    ctx.logf("extract: %zu цепочек между вершинами", chains.size());

    // система вращения
    const int nc = (int)chains.size();
    std::vector<std::vector<int>> rot(nn);  // исходящие полурёбра (2c: a->b, 2c+1: b->a)
    for (int c = 0; c < nc; ++c) {
        rot[chains[c].a].push_back(2 * c);
        rot[chains[c].b].push_back(2 * c + 1);
    }
    auto head = [&](int h) { return (h & 1) ? chains[h >> 1].a : chains[h >> 1].b; };
    auto tail = [&](int h) { return (h & 1) ? chains[h >> 1].b : chains[h >> 1].a; };
    auto first_nb = [&](int h) { return (h & 1) ? chains[h >> 1].first_b : chains[h >> 1].first_a; };
    std::vector<int> pos_in_rot(2 * (size_t)nc, -1);
    for (int v = 0; v < nn; ++v) {
        auto& r = rot[v];
        if (r.empty()) continue;
        V3 n = B.nodes.nrm[v];
        V3 r1 = std::fabs(n.x) < 0.9 ? cross(n, V3(1, 0, 0)) : cross(n, V3(0, 1, 0));
        r1 = normalized(r1);
        V3 r2 = cross(n, r1);
        std::vector<std::pair<double, int>> ang;
        for (int h : r) {
            V3 d = B.nodes.pos[first_nb(h)] - B.nodes.pos[v];
            ang.push_back({std::atan2(dot(d, r2), dot(d, r1)), h});
        }
        std::sort(ang.begin(), ang.end());
        for (size_t i = 0; i < ang.size(); ++i) {
            r[i] = ang[i].second;
            pos_in_rot[ang[i].second] = (int)i;
        }
    }
    // трассировка граней
    std::vector<uint8_t> visited(2 * (size_t)nc, 0);
    std::vector<int> vmap(nn, -1);
    auto vid = [&](int node) {
        if (vmap[node] < 0) {
            vmap[node] = (int)out.V.size();
            out.V.push_back(B.nodes.pos[node]);
            out.vfeature.push_back(B.nodes.feat[node]);
        }
        return vmap[node];
    };
    for (int h0 = 0; h0 < 2 * nc; ++h0) {
        if (visited[h0]) continue;
        std::vector<int> cyc;
        int h = h0;
        bool all_bnd = true;
        int guard = 0;
        while (!visited[h]) {
            visited[h] = 1;
            cyc.push_back(tail(h));
            all_bnd = all_bnd && chains[h >> 1].bnd;
            int b = head(h);
            int tw = h ^ 1;
            int idx = pos_in_rot[tw];
            int deg = (int)rot[b].size();
            h = rot[b][(idx - 1 + deg) % deg];
            if (++guard > 100000) break;
        }
        if (cyc.size() == 4) {
            std::array<int, 4> q;
            for (int k = 0; k < 4; ++k) q[k] = cyc[k];
            // повторяющиеся вершины -> дефект
            bool dup = q[0] == q[1] || q[0] == q[2] || q[0] == q[3] || q[1] == q[2] || q[1] == q[3] || q[2] == q[3];
            if (dup) { ++out.nonquad; continue; }
            out.quads.push_back({vid(q[0]), vid(q[1]), vid(q[2]), vid(q[3])});
        } else if (all_bnd && cyc.size() > 4) {
            ++out.boundary_loops;
        } else {
            ++out.nonquad;
            std::vector<int> poly;
            for (int n : cyc) poly.push_back(vid(n));
            out.others.push_back(poly);
        }
    }
    ctx.logf("extract: %zu квадов, %d не-квадов, %d внешних циклов, %d тупиков (%.2f с)",
             out.quads.size(), out.nonquad, out.boundary_loops, out.dead_ends, now_seconds() - t0);
    out.ok = !out.quads.empty();
    return out.ok;
}

}  // namespace qf
