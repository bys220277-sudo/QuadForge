#include "qf_quadpost.h"

#include <unordered_map>
#include <unordered_set>

namespace qf {

namespace {

inline u64 ekey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((u64)(uint32_t)a << 32) | (uint32_t)b;
}

struct Filler {
    PolyMesh& pm;
    RepairStats& st;
    std::vector<std::vector<int>> out;

    int add_vertex(const V3& p) {
        pm.V.push_back(p);
        pm.vfeat.push_back(0);
        return (int)pm.V.size() - 1;
    }

    void quad(int a, int b, int c, int d) {
        out.push_back({a, b, c, d});
        ++st.added_quads;
    }

    void fill_even(std::vector<int> loop) {
        const int n = (int)loop.size();
        if (n < 4) return;
        if (n == 4) { quad(loop[0], loop[1], loop[2], loop[3]); return; }
        if (n <= 8) {
            V3 c;
            for (int v : loop) c += pm.V[v];
            c = c / n;
            int cv = add_vertex(c);
            for (int i = 0; i < n; i += 2) quad(cv, loop[i], loop[(i + 1) % n], loop[(i + 2) % n]);
            return;
        }
        // делим по самой короткой «противоположной» хорде
        int best_i = 0;
        double best_d = 1e300;
        for (int i = 0; i < n / 2; ++i) {
            double d = norm(pm.V[loop[i]] - pm.V[loop[(i + n / 2) % n]]);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        int i = best_i, j = (best_i + n / 2) % n;
        std::vector<int> A, B;
        for (int k = i;; k = (k + 1) % n) { A.push_back(loop[k]); if (k == j) break; }
        for (int k = j;; k = (k + 1) % n) { B.push_back(loop[k]); if (k == i) break; }
        if ((n / 2) % 2 == 0) {
            // нужна промежуточная вершина на хорде
            int m = add_vertex((pm.V[loop[i]] + pm.V[loop[j]]) * 0.5);
            A.push_back(m);
            B.push_back(m);
        }
        fill_even(A);
        fill_even(B);
    }

    void fill(std::vector<int> loop) {
        int n = (int)loop.size();
        if (n < 3) return;
        if (n == 3) {
            out.push_back({loop[0], loop[1], loop[2]});
            ++st.added_tris;
            return;
        }
        if (n % 2 == 1) {
            // отрезаем угол с наименьшим углом треугольником
            int bi = 0;
            double best = 1e300;
            for (int i = 0; i < n; ++i) {
                const V3& p = pm.V[loop[i]];
                V3 a = pm.V[loop[(i + n - 1) % n]] - p, b = pm.V[loop[(i + 1) % n]] - p;
                double cs = dot(a, b) / std::max(norm(a) * norm(b), 1e-300);
                double ang = std::acos(std::max(-1.0, std::min(1.0, cs)));
                if (ang < best) { best = ang; bi = i; }
            }
            out.push_back({loop[(bi + n - 1) % n], loop[bi], loop[(bi + 1) % n]});
            ++st.added_tris;
            loop.erase(loop.begin() + bi);
            --n;
        }
        fill_even(loop);
    }
};

}  // namespace

void compact_mesh(PolyMesh& pm) {
    std::vector<int> remap(pm.V.size(), -1);
    std::vector<V3> V;
    std::vector<uint8_t> vf;
    for (auto& f : pm.faces)
        for (int& v : f) {
            if (remap[v] < 0) {
                remap[v] = (int)V.size();
                V.push_back(pm.V[v]);
                vf.push_back(v < (int)pm.vfeat.size() ? pm.vfeat[v] : 0);
            }
            v = remap[v];
        }
    pm.V.swap(V);
    pm.vfeat.swap(vf);
}

void repair_quad_mesh(PolyMesh& pm, RepairStats& st, Ctx& ctx) {
    const int nf = (int)pm.faces.size();
    std::vector<uint8_t> bad(nf, 0);
    // исходные граничные рёбра (до удаления)
    std::unordered_map<u64, int> ecount;
    std::unordered_map<u64, int> dircount;  // (a,b) направленное
    ecount.reserve(nf * 4);
    for (int f = 0; f < nf; ++f) {
        auto& F = pm.faces[f];
        int n = (int)F.size();
        if (n != 4) bad[f] = 1;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (F[i] == F[j]) bad[f] = 1;
        for (int i = 0; i < n; ++i) {
            int a = F[i], b = F[(i + 1) % n];
            ++ecount[ekey(a, b)];
            ++dircount[((u64)(uint32_t)a << 32) | (uint32_t)b];
        }
    }
    std::unordered_set<u64> orig_boundary;
    for (auto& kv : ecount)
        if (kv.second == 1) orig_boundary.insert(kv.first);
    for (int f = 0; f < nf; ++f) {
        auto& F = pm.faces[f];
        int n = (int)F.size();
        for (int i = 0; i < n; ++i) {
            int a = F[i], b = F[(i + 1) % n];
            if (ecount[ekey(a, b)] > 2) bad[f] = 1;
            if (dircount[((u64)(uint32_t)a << 32) | (uint32_t)b] > 1) bad[f] = 1;
        }
    }
    // Нерегулярные «бабочки»: вершина, где грани образуют несколько вееров — удалим грани с плохим
    // окружением грубо: оставляем как есть (обрабатывается разбиением петель ниже).
    std::vector<std::vector<int>> keep;
    keep.reserve(nf);
    for (int f = 0; f < nf; ++f) {
        if (bad[f]) { ++st.removed_faces; continue; }
        keep.push_back(pm.faces[f]);
    }
    if (st.removed_faces == 0) return;
    // граничные полурёбра после удаления
    std::unordered_map<u64, int> ec2;
    for (auto& F : keep)
        for (size_t i = 0; i < F.size(); ++i) ++ec2[ekey(F[i], F[(i + 1) % F.size()])];
    // исходящие полурёбра дыр: для грани с ребром a->b, граничным, дыра обходит b->a
    std::unordered_map<int, std::vector<int>> outgoing;  // вершина -> концы
    std::unordered_map<u64, uint8_t> used;
    int nhe = 0;
    for (auto& F : keep)
        for (size_t i = 0; i < F.size(); ++i) {
            int a = F[i], b = F[(i + 1) % F.size()];
            if (ec2[ekey(a, b)] == 1) { outgoing[b].push_back(a); ++nhe; }
        }
    Filler fl{pm, st, {}};
    for (auto& kv : outgoing) {
        for (int start_to : kv.second) {
            int s = kv.first;
            u64 k0 = ((u64)(uint32_t)s << 32) | (uint32_t)start_to;
            if (used.count(k0)) continue;
            // обход петли
            std::vector<int> loop;
            int cur = s, nxt = start_to;
            bool ok = true;
            for (int guard = 0; guard < 100000; ++guard) {
                used[((u64)(uint32_t)cur << 32) | (uint32_t)nxt] = 1;
                loop.push_back(cur);
                if (nxt == s) break;
                auto it = outgoing.find(nxt);
                if (it == outgoing.end()) { ok = false; break; }
                int cand = -1;
                for (int t : it->second)
                    if (!used.count(((u64)(uint32_t)nxt << 32) | (uint32_t)t)) { cand = t; break; }
                if (cand < 0) { ok = false; break; }
                cur = nxt;
                nxt = cand;
            }
            if (!ok || loop.size() < 3) { ++st.holes_left; continue; }
            // исходная граница? не закрываем
            bool touches_orig = false;
            for (size_t i = 0; i < loop.size(); ++i)
                if (orig_boundary.count(ekey(loop[i], loop[(i + 1) % loop.size()]))) { touches_orig = true; break; }
            if (touches_orig) { ++st.holes_left; continue; }
            // разбиение несimple-петель по повторяющимся вершинам
            std::vector<std::vector<int>> parts;
            std::vector<int> stack;
            std::unordered_map<int, int> pos;
            for (int v : loop) {
                auto it = pos.find(v);
                if (it != pos.end()) {
                    const int at = it->second;  // копия: erase ниже инвалидирует итератор
                    std::vector<int> sub(stack.begin() + at, stack.end());
                    parts.push_back(sub);
                    for (size_t q = at; q < stack.size(); ++q) pos.erase(stack[q]);
                    stack.resize(at);
                }
                pos[v] = (int)stack.size();
                stack.push_back(v);
            }
            if (!stack.empty()) parts.push_back(stack);
            for (auto& part : parts) {
                if (part.size() >= 3) { fl.fill(part); ++st.holes_filled; }
            }
        }
    }
    keep.insert(keep.end(), fl.out.begin(), fl.out.end());
    pm.faces.swap(keep);
    compact_mesh(pm);
    ctx.logf("ремонт: удалено %d граней, заполнено дыр %d (+%d квадов, +%d треуг.), оставлено дыр %d",
             st.removed_faces, st.holes_filled, st.added_quads, st.added_tris, st.holes_left);
}

int remove_doublets(PolyMesh& pm) {
    int total = 0;
    for (int pass = 0; pass < 8; ++pass) {
        const int nv = (int)pm.V.size();
        std::vector<std::vector<int>> vf(nv);
        std::unordered_map<u64, int> ecount;
        for (int f = 0; f < (int)pm.faces.size(); ++f) {
            const auto& F = pm.faces[f];
            for (size_t i = 0; i < F.size(); ++i) {
                vf[F[i]].push_back(f);
                ++ecount[ekey(F[i], F[(i + 1) % F.size()])];
            }
        }
        std::vector<uint8_t> dead(pm.faces.size(), 0), touched(pm.faces.size(), 0);
        int removed = 0;
        for (int v = 0; v < nv; ++v) {
            if (vf[v].size() != 2) continue;
            if (v < (int)pm.vfeat.size() && pm.vfeat[v]) continue;
            const int f1 = vf[v][0], f2 = vf[v][1];
            if (f1 == f2 || touched[f1] || touched[f2]) continue;
            const auto& A = pm.faces[f1];
            const auto& B = pm.faces[f2];
            if (A.size() != 4 || B.size() != 4) continue;
            int ia = (int)(std::find(A.begin(), A.end(), v) - A.begin());
            int ib = (int)(std::find(B.begin(), B.end(), v) - B.begin());
            // A = (v, a, x, b), B должна быть (v, b, y, a)
            const int a = A[(ia + 1) % 4], x = A[(ia + 2) % 4], b = A[(ia + 3) % 4];
            if (B[(ib + 1) % 4] != b || B[(ib + 3) % 4] != a) continue;
            const int y = B[(ib + 2) % 4];
            if (ecount[ekey(v, a)] != 2 || ecount[ekey(v, b)] != 2) continue;
            if (a == b || x == y || x == a || x == b || y == a || y == b) continue;
            pm.faces[f1] = {a, x, b, y};
            dead[f2] = 1;
            touched[f1] = touched[f2] = 1;
            ++removed;
        }
        if (removed == 0) break;
        std::vector<std::vector<int>> keep;
        keep.reserve(pm.faces.size());
        for (size_t f = 0; f < pm.faces.size(); ++f)
            if (!dead[f]) keep.push_back(std::move(pm.faces[f]));
        pm.faces.swap(keep);
        total += removed;
    }
    if (total) compact_mesh(pm);
    return total;
}

void relax_and_project(PolyMesh& pm, const TriMesh& src, const std::vector<uint8_t>& src_feat, double L,
                       const RelaxOptions& opt, Ctx& ctx) {
    const int nv = (int)pm.V.size();
    if (nv == 0) return;
    // BVH исходника и фичер-отрезков
    TriBvh surf;
    surf.build(src.P, src.F);
    SegBvh fseg;
    {
        std::vector<V3> A, B;
        for (int f = 0; f < src.nF(); ++f)
            for (int i = 0; i < 3; ++i) {
                int g = src.FF[3 * f + i];
                bool is = g < 0 || src_feat[3 * f + i];
                if (!is || (g >= 0 && g < f)) continue;
                A.push_back(src.P[src.F[3 * f + i]]);
                B.push_back(src.P[src.F[3 * f + (i + 1) % 3]]);
            }
        if (!A.empty()) fseg.build(A, B);
    }
    // рёбра и смежность
    std::unordered_map<u64, int> ecount;
    for (auto& F : pm.faces)
        for (size_t i = 0; i < F.size(); ++i) ++ecount[ekey(F[i], F[(i + 1) % F.size()])];
    std::vector<std::vector<int>> nbr(nv);
    std::vector<int> fdeg(nv, 0);
    std::vector<std::vector<int>> chain(nv);  // соседи по фичер-рёбрам
    for (auto& kv : ecount) {
        int a = (int)(kv.first >> 32), b = (int)(kv.first & 0xffffffffu);
        nbr[a].push_back(b);
        nbr[b].push_back(a);
        bool fe = kv.second == 1;  // граница квад-сетки
        if (!fe && pm.vfeat[a] && pm.vfeat[b] && !fseg.bvh.empty()) {
            V3 mid = (pm.V[a] + pm.V[b]) * 0.5, q;
            double d2;
            if (fseg.nearest(mid, q, d2) >= 0 && d2 < (0.2 * L) * (0.2 * L)) fe = true;
        }
        if (fe) { chain[a].push_back(b); chain[b].push_back(a); }
    }
    std::vector<uint8_t> type(nv, 0);
    for (int v = 0; v < nv; ++v) {
        size_t c = chain[v].size();
        type[v] = c == 0 ? 0 : (c == 2 ? 1 : 2);
        if (c == 2) {  // острый поворот цепочки — угол
            V3 a = pm.V[chain[v][0]] - pm.V[v], b = pm.V[chain[v][1]] - pm.V[v];
            double cs = dot(a, b) / std::max(norm(a) * norm(b), 1e-300);
            if (cs > -0.7) type[v] = 2;
        }
    }
    // нормали вершин по граням
    auto vertex_normals = [&](std::vector<V3>& N) {
        N.assign(nv, V3());
        for (auto& F : pm.faces) {
            V3 n;
            for (size_t i = 0; i < F.size(); ++i) n += cross(pm.V[F[i]], pm.V[F[(i + 1) % F.size()]]);
            for (int v : F) N[v] += n;
        }
        for (auto& n : N) n = normalized(n);
    };
    auto project = [&](int v, const V3& p) -> V3 {
        V3 q;
        double d2;
        if (type[v] >= 1 && !fseg.bvh.empty()) {
            if (fseg.nearest(p, q, d2) >= 0) return q;
        }
        if (surf.nearest(p, q, d2) >= 0) return q;
        return p;
    };
    std::vector<V3> N, Q(nv);
    for (int it = 0; it < opt.iterations; ++it) {
        ctx.check_cancel();
        vertex_normals(N);
        parallel_each(nv, [&](i64 v) {
            Q[v] = pm.V[v];
            if (type[v] == 2) return;
            const auto& nb = type[v] == 1 ? chain[v] : nbr[v];
            if (nb.empty()) return;
            V3 c;
            for (int u : nb) c += pm.V[u];
            c = c / (double)nb.size();
            V3 d = (c - pm.V[v]) * opt.lambda;
            if (type[v] == 0) d = d - N[v] * dot(N[v], d);
            Q[v] = pm.V[v] + d;
        }, 256);
        if (opt.project) {
            parallel_each(nv, [&](i64 v) { pm.V[v] = type[v] == 2 ? pm.V[v] : project((int)v, Q[v]); }, 256);
        } else {
            pm.V.swap(Q);
        }
    }
    if (opt.iterations == 0 && opt.project) {
        parallel_each(nv, [&](i64 v) { pm.V[v] = project((int)v, pm.V[v]); }, 256);
    } else if (opt.project) {
        // финальная проекция углов на поверхность
        parallel_each(nv, [&](i64 v) {
            if (type[v] == 2) {
                V3 q;
                double d2;
                if (surf.nearest(pm.V[v], q, d2) >= 0) pm.V[v] = q;
            }
        }, 256);
    }
}

QuadQuality evaluate_quad_mesh(const PolyMesh& pm) {
    QuadQuality q;
    q.faces = (int)pm.faces.size();
    std::unordered_map<u64, int> ec;
    std::vector<int> val(pm.V.size(), 0);
    double amin = 360, amax = 0, adev = 0;
    int nang = 0;
    double lsum = 0, l2sum = 0;
    int ln = 0;
    for (auto& F : pm.faces) {
        int n = (int)F.size();
        if (n == 4) ++q.quads; else if (n == 3) ++q.tris; else ++q.ngons;
        for (int i = 0; i < n; ++i) {
            ++ec[ekey(F[i], F[(i + 1) % n])];
            const V3& p = pm.V[F[i]];
            V3 a = pm.V[F[(i + n - 1) % n]] - p, b = pm.V[F[(i + 1) % n]] - p;
            double cs = dot(a, b) / std::max(norm(a) * norm(b), 1e-300);
            double ang = std::acos(std::max(-1.0, std::min(1.0, cs))) * 180.0 / kPi;
            amin = std::min(amin, ang);
            amax = std::max(amax, ang);
            if (n == 4) { adev += std::fabs(ang - 90.0); ++nang; }
        }
    }
    for (auto& kv : ec) {
        int a = (int)(kv.first >> 32), b = (int)(kv.first & 0xffffffffu);
        ++val[a];
        ++val[b];
        if (kv.second == 1) ++q.boundary_edges;
        if (kv.second > 2) ++q.nonmanifold_edges;
        double l = norm(pm.V[a] - pm.V[b]);
        lsum += l;
        l2sum += l * l;
        ++ln;
    }
    q.verts = 0;
    for (size_t v = 0; v < pm.V.size(); ++v) {
        if (val[v] == 0) continue;
        ++q.verts;
        bool on_b = false;  // упрощённо: граничные вершины считаем регулярными при валентности 3
        (void)on_b;
        if (val[v] != 4) ++q.irregular;
    }
    q.min_angle = amin;
    q.max_angle = amax;
    q.angle_dev_mean = nang ? adev / nang : 0;
    if (ln) {
        double m = lsum / ln;
        double var = std::max(0.0, l2sum / ln - m * m);
        q.edge_cv = m > 0 ? std::sqrt(var) / m : 0;
    }
    return q;
}

}  // namespace qf
