#include "qf_sparse.h"

namespace qf {

// ---------------------------------------------------------------------------
// Вложенное рассечение по координатам
// ---------------------------------------------------------------------------
std::vector<int> nested_dissection_order(const SymCsr& A, const std::vector<double>& coords, int leaf) {
    const int n = A.n;
    std::vector<int> order;
    order.reserve(n);
    std::vector<int> mark(n, -1);   // идентификатор текущего подмножества
    std::vector<int> side(n, 0);
    struct Task { std::vector<int> set; };
    // рекурсия через явный стек: (множество, фаза)
    struct Frame { std::vector<int> set; std::vector<int> sep; int stage; std::vector<int> a, b; };
    std::vector<Frame> st;
    st.push_back({std::vector<int>(), {}, 0, {}, {}});
    st.back().set.resize(n);
    for (int i = 0; i < n; ++i) st.back().set[i] = i;
    int stamp = 0;
    while (!st.empty()) {
        Frame& fr = st.back();
        if (fr.stage == 0) {
            if ((int)fr.set.size() <= leaf) {
                for (int v : fr.set) order.push_back(v);
                st.pop_back();
                continue;
            }
            // ось наибольшей протяжённости
            double mn[3] = {1e300, 1e300, 1e300}, mx[3] = {-1e300, -1e300, -1e300};
            for (int v : fr.set)
                for (int k = 0; k < 3; ++k) {
                    mn[k] = std::min(mn[k], coords[3 * v + k]);
                    mx[k] = std::max(mx[k], coords[3 * v + k]);
                }
            int ax = 0;
            for (int k = 1; k < 3; ++k) if (mx[k] - mn[k] > mx[ax] - mn[ax]) ax = k;
            std::vector<int> s = fr.set;
            size_t mid = s.size() / 2;
            std::nth_element(s.begin(), s.begin() + mid, s.end(),
                             [&](int x, int y) { return coords[3 * x + ax] < coords[3 * y + ax]; });
            ++stamp;
            for (size_t q = 0; q < s.size(); ++q) { mark[s[q]] = stamp; side[s[q]] = q < mid ? 0 : 1; }
            // разделитель: вершины стороны 0, смежные со стороной 1
            std::vector<int> a, b, sep;
            for (int v : s) {
                if (side[v] == 1) { b.push_back(v); continue; }
                bool touch = false;
                for (i64 p = A.off[v]; p < A.off[v + 1]; ++p) {
                    int u = A.col[p];
                    if (mark[u] == stamp && side[u] == 1) { touch = true; break; }
                }
                if (touch) sep.push_back(v); else a.push_back(v);
            }
            if (a.empty() || b.empty()) {  // вырожденное разбиение
                for (int v : fr.set) order.push_back(v);
                st.pop_back();
                continue;
            }
            fr.sep = std::move(sep);
            fr.a = std::move(a);
            fr.b = std::move(b);
            fr.stage = 1;
            std::vector<int> aa = fr.a;
            st.push_back({std::move(aa), {}, 0, {}, {}});
        } else if (fr.stage == 1) {
            fr.stage = 2;
            std::vector<int> bb = fr.b;
            st.push_back({std::move(bb), {}, 0, {}, {}});
        } else {
            for (int v : fr.sep) order.push_back(v);
            st.pop_back();
        }
    }
    return order;
}

// ---------------------------------------------------------------------------
// Супернодальный мультифронтальный Холецкий
// ---------------------------------------------------------------------------
namespace {

// Частичное разложение фронта F (m x m, по столбцам, нижний треугольник): первые k столбцов
// раскладываются, хвост F[k:,k:] получает дополнение Шура.
bool partial_factor(double* F, int m, int k) {
    constexpr int B = 48;
    for (int kb = 0; kb < k; kb += B) {
        const int ke = std::min(k, kb + B);
        for (int t = kb; t < ke; ++t) {
            double* ct = F + (size_t)t * m;
            double d = ct[t];
            if (!(d > 0) || !std::isfinite(d)) return false;
            d = std::sqrt(d);
            ct[t] = d;
            const double inv = 1.0 / d;
            for (int i = t + 1; i < m; ++i) ct[i] *= inv;
            for (int t2 = t + 1; t2 < ke; ++t2) {
                const double a = ct[t2];
                if (a == 0.0) continue;
                double* c2 = F + (size_t)t2 * m;
                for (int i = t2; i < m; ++i) c2[i] -= a * ct[i];
            }
        }
        if (ke >= m) continue;
        auto update = [&](i64 jb, i64 je) {
            for (i64 j = jb; j < je; ++j) {
                double* cj = F + (size_t)j * m;
                for (int t = kb; t < ke; ++t) {
                    const double* ct = F + (size_t)t * m;
                    const double a = ct[j];
                    if (a == 0.0) continue;
                    for (int i = (int)j; i < m; ++i) cj[i] -= a * ct[i];
                }
            }
        };
        const i64 rest = m - ke;
        if ((double)rest * rest * (ke - kb) > 2e6) parallel_for(rest, [&](i64 b, i64 e) { update(ke + b, ke + e); }, 4);
        else update(ke, m);
    }
    return true;
}

}  // namespace

bool SparseCholesky::factor(const SymCsr& A, const std::vector<int>& perm, double reg, i64 max_nnz) {
    n = A.n;
    ok_ = false;
    ns = 0;
    nnz_ = 0;
    std::vector<int> p0 = perm;
    if ((int)p0.size() != n) {
        p0.resize(n);
        for (int i = 0; i < n; ++i) p0[i] = i;
    }
    std::vector<int> ip0(n);
    for (int i = 0; i < n; ++i) ip0[p0[i]] = i;
    // дерево исключения в порядке p0 (верхняя часть столбца k = соседи с меньшим номером)
    auto etree = [&](const std::vector<int>& pm, const std::vector<int>& ipm, std::vector<int>& parent) {
        parent.assign(n, -1);
        std::vector<int> ancestor(n, -1);
        for (int k = 0; k < n; ++k) {
            const int ko = pm[k];
            for (i64 p = A.off[ko]; p < A.off[ko + 1]; ++p) {
                int i = ipm[A.col[p]];
                while (i != -1 && i < k) {
                    const int inext = ancestor[i];
                    ancestor[i] = k;
                    if (inext == -1) parent[i] = k;
                    i = inext;
                }
            }
        }
    };
    std::vector<int> parent0;
    etree(p0, ip0, parent0);
    // постпорядок
    std::vector<int> head(n, -1), next(n, -1), post(n);
    for (int j = n - 1; j >= 0; --j)
        if (parent0[j] != -1) { next[j] = head[parent0[j]]; head[parent0[j]] = j; }
    {
        int k = 0;
        std::vector<int> stack;
        for (int j = 0; j < n; ++j) {
            if (parent0[j] != -1) continue;
            stack.push_back(j);
            while (!stack.empty()) {
                const int p = stack.back();
                const int c = head[p];
                if (c == -1) { stack.pop_back(); post[k++] = p; }
                else { head[p] = next[c]; stack.push_back(c); }
            }
        }
    }
    perm_.resize(n);
    for (int k = 0; k < n; ++k) perm_[k] = p0[post[k]];
    iperm_.assign(n, 0);
    for (int i = 0; i < n; ++i) iperm_[perm_[i]] = i;
    std::vector<int> parent;
    etree(perm_, iperm_, parent);
    // число элементов в столбцах L (ereach)
    std::vector<i64> cc(n, 1);
    {
        std::vector<int> w(n, -1);
        for (int k = 0; k < n; ++k) {
            w[k] = k;
            const int ko = perm_[k];
            for (i64 p = A.off[ko]; p < A.off[ko + 1]; ++p) {
                int i = iperm_[A.col[p]];
                if (i >= k) continue;
                for (; w[i] != k; i = parent[i]) { w[i] = k; ++cc[i]; }
            }
        }
    }
    // супернодам: фундаментальные + ослабленная амальгамация цепочек
    std::vector<int> nchild(n, 0);
    for (int j = 0; j < n; ++j) if (parent[j] != -1) ++nchild[parent[j]];
    sn_c0.clear();
    sn_c0.push_back(0);
    {
        int sz = 1;
        for (int j = 1; j < n; ++j) {
            bool merge = false;
            if (parent[j - 1] == j) {
                if (cc[j - 1] == cc[j] + 1 && nchild[j] == 1) merge = true;
                else if (sz < 16 && (cc[j] + 1 - cc[j - 1]) * sz <= 8 * sz + 16) merge = true;
            }
            if (merge && sz < 128) ++sz;
            else { sn_c0.push_back(j); sz = 1; }
        }
    }
    ns = (int)sn_c0.size();
    sn_c0.push_back(n);
    std::vector<int> snode_of(n);
    for (int s = 0; s < ns; ++s)
        for (int j = sn_c0[s]; j < sn_c0[s + 1]; ++j) snode_of[j] = s;
    std::vector<int> sn_parent(ns, -1);
    std::vector<std::vector<int>> sn_children(ns);
    for (int s = 0; s < ns; ++s) {
        const int pj = parent[sn_c0[s + 1] - 1];
        if (pj != -1) { sn_parent[s] = snode_of[pj]; sn_children[snode_of[pj]].push_back(s); }
    }
    // структуры строк
    sn_roff.assign(ns + 1, 0);
    sn_rows.clear();
    {
        std::vector<int> mark(n, -1), extra;
        for (int s = 0; s < ns; ++s) {
            const int c0 = sn_c0[s], c1 = sn_c0[s + 1];
            extra.clear();
            for (int j = c0; j < c1; ++j) mark[j] = s;
            for (int j = c0; j < c1; ++j) {
                const int jo = perm_[j];
                for (i64 p = A.off[jo]; p < A.off[jo + 1]; ++p) {
                    const int i = iperm_[A.col[p]];
                    if (i >= c1 && mark[i] != s) { mark[i] = s; extra.push_back(i); }
                }
            }
            for (int c : sn_children[s]) {
                const int kc = sn_c0[c + 1] - sn_c0[c];
                for (i64 q = sn_roff[c] + kc; q < sn_roff[c + 1]; ++q) {
                    const int i = sn_rows[q];
                    if (i >= c1 && mark[i] != s) { mark[i] = s; extra.push_back(i); }
                }
            }
            std::sort(extra.begin(), extra.end());
            for (int j = c0; j < c1; ++j) sn_rows.push_back(j);
            sn_rows.insert(sn_rows.end(), extra.begin(), extra.end());
            sn_roff[s + 1] = (i64)sn_rows.size();
        }
    }
    sn_xoff.assign(ns + 1, 0);
    for (int s = 0; s < ns; ++s) {
        const i64 m = sn_roff[s + 1] - sn_roff[s], k = sn_c0[s + 1] - sn_c0[s];
        sn_xoff[s + 1] = sn_xoff[s] + m * k;
    }
    nnz_ = sn_xoff[ns];
    if (max_nnz > 0 && nnz_ > max_nnz) return false;
    Lx.assign(nnz_, 0.0);
    // численный мультифронтальный проход (супернод в постпорядке: дети раньше родителей)
    std::vector<int> pos(n, -1);
    std::vector<double> F;
    std::vector<double> ustack;
    struct Upd { i64 off; int mu; int sn; };
    std::vector<Upd> ust;
    std::vector<int> idx;
    for (int s = 0; s < ns; ++s) {
        const int c0 = sn_c0[s], c1 = sn_c0[s + 1], k = c1 - c0;
        const int m = (int)(sn_roff[s + 1] - sn_roff[s]);
        const int* R = &sn_rows[sn_roff[s]];
        for (int a = 0; a < m; ++a) pos[R[a]] = a;
        F.assign((size_t)m * m, 0.0);
        for (int j = c0; j < c1; ++j) {
            const int t = j - c0;
            const int jo = perm_[j];
            double* col = &F[(size_t)t * m];
            for (i64 p = A.off[jo]; p < A.off[jo + 1]; ++p) {
                const int i = iperm_[A.col[p]];
                if (i >= j) col[pos[i]] += A.val[p];
            }
            col[t] += reg;
        }
        // сложение обновлений детей (вершина стека)
        const int nc = (int)sn_children[s].size();
        for (int q = 0; q < nc; ++q) {
            const Upd u = ust.back();
            ust.pop_back();
            const int kc = sn_c0[u.sn + 1] - sn_c0[u.sn];
            const int* Rc = &sn_rows[sn_roff[u.sn] + kc];
            idx.resize(u.mu);
            for (int a = 0; a < u.mu; ++a) idx[a] = pos[Rc[a]];
            const double* U = &ustack[u.off];
            for (int b = 0; b < u.mu; ++b) {
                double* col = &F[(size_t)idx[b] * m];
                const double* ub = U + (size_t)b * u.mu;
                for (int a = b; a < u.mu; ++a) col[idx[a]] += ub[a];
            }
            ustack.resize(u.off);
        }
        if (!partial_factor(F.data(), m, k)) return false;
        double* X = &Lx[sn_xoff[s]];
        std::memcpy(X, F.data(), sizeof(double) * (size_t)m * k);
        if (m > k) {
            const int mu = m - k;
            const i64 off = (i64)ustack.size();
            ustack.resize(off + (i64)mu * mu);
            double* U = &ustack[off];
            for (int b = 0; b < mu; ++b)
                std::memcpy(U + (size_t)b * mu + b, &F[(size_t)(k + b) * m + k + b], sizeof(double) * (mu - b));
            ust.push_back({off, mu, s});
        }
    }
    ok_ = true;
    return true;
}

void SparseCholesky::solve_inplace(double* bx, double* y) const {
    for (int k = 0; k < n; ++k) y[k] = bx[perm_[k]];
    for (int s = 0; s < ns; ++s) {
        const int c0 = sn_c0[s], k = sn_c0[s + 1] - c0;
        const int m = (int)(sn_roff[s + 1] - sn_roff[s]);
        const int* R = &sn_rows[sn_roff[s]];
        const double* X = &Lx[sn_xoff[s]];
        for (int t = 0; t < k; ++t) {
            const double* ct = X + (size_t)t * m;
            const double yt = y[c0 + t] / ct[t];
            y[c0 + t] = yt;
            for (int i = t + 1; i < m; ++i) y[R[i]] -= ct[i] * yt;
        }
    }
    for (int s = ns - 1; s >= 0; --s) {
        const int c0 = sn_c0[s], k = sn_c0[s + 1] - c0;
        const int m = (int)(sn_roff[s + 1] - sn_roff[s]);
        const int* R = &sn_rows[sn_roff[s]];
        const double* X = &Lx[sn_xoff[s]];
        for (int t = k - 1; t >= 0; --t) {
            const double* ct = X + (size_t)t * m;
            double acc = y[c0 + t];
            for (int i = t + 1; i < m; ++i) acc -= ct[i] * y[R[i]];
            y[c0 + t] = acc / ct[t];
        }
    }
    for (int k = 0; k < n; ++k) bx[perm_[k]] = y[k];
}

void SparseCholesky::solve(const std::vector<double>& b, std::vector<double>& x) const {
    x = b;
    std::vector<double> y(n);
    solve_inplace(x.data(), y.data());
}

// ---------------------------------------------------------------------------
// Плотный Холецкий
// ---------------------------------------------------------------------------
bool dense_cholesky(double* A, int n) {
    constexpr int B = 64;
    for (int kb = 0; kb < n; kb += B) {
        const int ke = std::min(n, kb + B);
        // диагональный блок
        for (int k = kb; k < ke; ++k) {
            double* rk = A + (size_t)k * n;
            double d = rk[k];
            for (int p = kb; p < k; ++p) d -= rk[p] * rk[p];
            if (!(d > 0) || !std::isfinite(d)) return false;
            d = std::sqrt(d);
            rk[k] = d;
            for (int i = k + 1; i < ke; ++i) {
                double* ri = A + (size_t)i * n;
                double v = ri[k];
                for (int p = kb; p < k; ++p) v -= ri[p] * rk[p];
                ri[k] = v / d;
            }
        }
        if (ke >= n) break;
        // панель: строки ниже блока, L_ik = A_ik L_kk^{-T}
        parallel_for(n - ke, [&](i64 b, i64 e) {
            for (i64 q = b; q < e; ++q) {
                double* ri = A + (size_t)(ke + q) * n;
                for (int k = kb; k < ke; ++k) {
                    const double* rk = A + (size_t)k * n;
                    double v = ri[k];
                    for (int p = kb; p < k; ++p) v -= ri[p] * rk[p];
                    ri[k] = v / rk[k];
                }
            }
        }, 16);
        // обновление хвоста: A_ij -= L_i,blk · L_j,blk (j <= i)
        parallel_for(n - ke, [&](i64 b, i64 e) {
            for (i64 q = b; q < e; ++q) {
                const int i = ke + (int)q;
                double* ri = A + (size_t)i * n;
                for (int j = ke; j <= i; ++j) {
                    const double* rj = A + (size_t)j * n;
                    double v = 0;
                    for (int p = kb; p < ke; ++p) v += ri[p] * rj[p];
                    ri[j] -= v;
                }
            }
        }, 8);
    }
    return true;
}

void dense_cholesky_solve(const double* L, int n, double* b) {
    for (int i = 0; i < n; ++i) {
        const double* ri = L + (size_t)i * n;
        double v = b[i];
        for (int p = 0; p < i; ++p) v -= ri[p] * b[p];
        b[i] = v / ri[i];
    }
    for (int i = n - 1; i >= 0; --i) {
        double v = b[i] / L[(size_t)i * n + i];
        b[i] = v;
        for (int p = 0; p < i; ++p) b[p] -= L[(size_t)i * n + p] * v;
    }
}

}  // namespace qf
