// MIQ: комбинг поля, граф разреза, сектора, линейные связи (элиминация),
// квадратичная система, CG, жадное округление целочисленных переменных.
#include "qf_param.h"
#include "qf_sparse.h"

#include <cassert>
#include <unordered_map>

namespace qf {

namespace {

struct Term {
    int var;
    double coef;
};

// ---------------------------------------------------------------------------
// Элиминатор линейных однородных связей Σ a_j y_j = 0.
// ---------------------------------------------------------------------------
class Eliminator {
public:
    Eliminator(int n, std::vector<uint8_t> prio)
        : dep_(n, 0), expr_(n), prio_(std::move(prio)), stamp_(n, ~(u64)0) {}

    // Полное разложение c*var по текущим свободным переменным (добавляется в out).
    // Попутно сжимает пути (выражения зависимых переменных переписываются в развёрнутом виде).
    void expand(int var, double c, std::vector<Term>& out, std::vector<Term>& stack) {
        if (!dep_[var]) { out.push_back({var, c}); return; }
        compress(var);
        for (const Term& u : expr_[var]) out.push_back({u.var, c * u.coef});
        (void)stack;
    }

    // Переписывает expr_[var] так, чтобы оно ссылалось только на свободные переменные.
    void compress(int var) {
        // итеративный обход в глубину по зависимостям
        std::vector<std::pair<int, size_t>> st;
        st.push_back({var, 0});
        while (!st.empty()) {
            auto& top = st.back();
            int v = top.first;
            auto& ex = expr_[v];
            bool pushed = false;
            while (top.second < ex.size()) {
                int u = ex[top.second].var;
                if (dep_[u] && !clean(u)) { st.push_back({u, 0}); pushed = true; break; }
                ++top.second;
            }
            if (pushed) continue;
            // все ссылки v либо свободны, либо уже «чистые» — разворачиваем
            if (!clean(v)) {
                std::vector<Term> acc;
                acc.reserve(ex.size());
                for (const Term& t : ex) {
                    if (!dep_[t.var]) acc.push_back(t);
                    else for (const Term& w : expr_[t.var]) acc.push_back({w.var, t.coef * w.coef});
                }
                merge(acc);
                ex.swap(acc);
            }
            stamp_[v] = epoch_;
            st.pop_back();
        }
    }

    static void merge(std::vector<Term>& v) {
        std::sort(v.begin(), v.end(), [](const Term& a, const Term& b) { return a.var < b.var; });
        size_t w = 0;
        for (size_t i = 0; i < v.size();) {
            size_t j = i;
            double s = 0;
            while (j < v.size() && v[j].var == v[i].var) s += v[j++].coef;
            if (std::fabs(s) > 1e-10) v[w++] = {v[i].var, s};
            i = j;
        }
        v.resize(w);
    }

    void resolved(int var, std::vector<Term>& out) const {
        out.clear();
        if (!dep_[var]) { out.push_back({var, 1.0}); return; }
        out = expr_[var];
    }

    // Возвращает pivot или -1 (избыточное уравнение).
    int add_equation(const std::vector<Term>& eq) {
        std::vector<Term> acc, stack;
        for (const Term& t : eq) expand(t.var, t.coef, acc, stack);
        merge(acc);
        if (acc.empty()) return -1;
        int best = -1;
        uint32_t bestKey = 0;
        for (size_t i = 0; i < acc.size(); ++i) {
            bool unit = std::fabs(std::fabs(acc[i].coef) - 1.0) < 1e-9;
            uint32_t key = ((uint32_t)prio_[acc[i].var] << 28) | ((unit ? 0u : 1u) << 27) |
                           (uint32_t)((0x7FFFFFF) - (acc[i].var & 0x7FFFFFF));
            if (best < 0 || key < bestKey) { best = (int)i; bestKey = key; }
        }
        const Term p = acc[best];
        std::vector<Term> e;
        e.reserve(acc.size() - 1);
        for (size_t i = 0; i < acc.size(); ++i)
            if ((int)i != best) e.push_back({acc[i].var, -acc[i].coef / p.coef});
        dep_[p.var] = 1;
        expr_[p.var] = std::move(e);
        ++neliminated_;
        ++epoch_;  // все ранее «чистые» выражения могли устареть (p.var стала зависимой)
        stamp_[p.var] = epoch_;  // e ссылается только на свободные
        order_.push_back(p.var);
        return p.var;
    }

    // После всех уравнений: развернуть все выражения до свободных переменных.
    void finalize() {
        ++epoch_;
        for (auto it = order_.rbegin(); it != order_.rend(); ++it) {
            int v = *it;
            std::vector<Term> acc;
            for (const Term& t : expr_[v]) {
                if (!dep_[t.var]) acc.push_back(t);
                else for (const Term& w : expr_[t.var]) acc.push_back({w.var, t.coef * w.coef});
            }
            merge(acc);
            expr_[v].swap(acc);
            stamp_[v] = epoch_;
        }
    }

    bool is_free(int v) const { return !dep_[v]; }
    int neliminated() const { return neliminated_; }

private:
    bool clean(int v) const { return stamp_[v] == epoch_; }
    std::vector<uint8_t> dep_;
    std::vector<std::vector<Term>> expr_;
    std::vector<uint8_t> prio_;
    std::vector<int> order_;
    std::vector<u64> stamp_;
    u64 epoch_ = 0;
    int neliminated_ = 0;
};

// ---------------------------------------------------------------------------
// CSR
// ---------------------------------------------------------------------------
struct Csr {
    int n = 0;
    std::vector<i64> off;
    std::vector<int> col;
    std::vector<double> val;
};

void spmv(const Csr& A, const std::vector<double>& x, std::vector<double>& y) {
    parallel_for(A.n, [&](i64 b, i64 e) {
        for (i64 i = b; i < e; ++i) {
            double s = 0;
            for (i64 k = A.off[i]; k < A.off[i + 1]; ++k) s += A.val[k] * x[A.col[k]];
            y[i] = s;
        }
    }, 1024);
}

double dotv(const std::vector<double>& a, const std::vector<double>& b) {
    return parallel_sum((i64)a.size(), [&](i64 s, i64 e) {
        double r = 0;
        for (i64 i = s; i < e; ++i) r += a[i] * b[i];
        return r;
    });
}

// Неполное разложение Холецкого IC(0) для маскированной системы (фиксированные строки — единичные).
struct IC0 {
    int n = 0;
    std::vector<int> perm, iperm;  // perm[new] = old, iperm[old] = new (пусто — тождественная)
    mutable std::vector<double> tmp_r, tmp_z;
    std::vector<i64> off;
    std::vector<int> col;
    std::vector<double> val;
    std::vector<double> d;
    bool ok = false;
    double shift = 0.0;
    // Сдвинутое IC(0) (Manteuffel): при пробое увеличиваем сдвиг диагонали.
    void factor(const Csr& H0, const std::vector<uint8_t>& fixed0) {
        static const double shifts[] = {0.0, 1e-3, 1e-2, 0.05, 0.2, 1.0};
        ok = false;
        // плотные строки (глобальные переменные) — в конец упорядочения
        const int N = H0.n;
        perm.resize(N);
        iperm.resize(N);
        {
            std::vector<std::pair<i64, int>> key(N);
            for (int i = 0; i < N; ++i) {
                i64 len = H0.off[i + 1] - H0.off[i];
                key[i] = {len > 64 ? len : 0, i};
            }
            std::stable_sort(key.begin(), key.end(), [](const std::pair<i64, int>& a, const std::pair<i64, int>& b) {
                return a.first < b.first;
            });
            for (int k = 0; k < N; ++k) { perm[k] = key[k].second; iperm[key[k].second] = k; }
        }
        Csr H;
        H.n = N;
        H.off.assign(N + 1, 0);
        for (int k = 0; k < N; ++k) H.off[k + 1] = H.off[k] + (H0.off[perm[k] + 1] - H0.off[perm[k]]);
        H.col.resize(H.off[N]);
        H.val.resize(H.off[N]);
        for (int k = 0; k < N; ++k) {
            int i = perm[k];
            std::vector<std::pair<int, double>> row;
            for (i64 q = H0.off[i]; q < H0.off[i + 1]; ++q) row.push_back({iperm[H0.col[q]], H0.val[q]});
            std::sort(row.begin(), row.end());
            i64 o = H.off[k];
            for (auto& e : row) { H.col[o] = e.first; H.val[o] = e.second; ++o; }
        }
        std::vector<uint8_t> fixed(N);
        for (int k = 0; k < N; ++k) fixed[k] = fixed0[perm[k]];
        for (double a : shifts) {
            if (try_factor(H, fixed, a)) { ok = true; shift = a; return; }
        }
    }
    bool try_factor(const Csr& H, const std::vector<uint8_t>& fixed, double alpha) {
        n = H.n;
        off.assign(n + 1, 0);
        for (int i = 0; i < n; ++i) {
            i64 c = 0;
            if (!fixed[i])
                for (i64 k = H.off[i]; k < H.off[i + 1]; ++k)
                    if (H.col[k] < i && !fixed[H.col[k]]) ++c;
            off[i + 1] = off[i] + c;
        }
        col.resize(off[n]);
        val.resize(off[n]);
        d.assign(n, 1.0);
        for (int i = 0; i < n; ++i) {
            if (fixed[i]) continue;
            i64 o = off[i];
            double aii = 0;
            for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) {
                int c = H.col[k];
                if (c == i) aii = H.val[k];
                else if (c < i && !fixed[c]) { col[o] = c; val[o] = H.val[k]; ++o; }
            }
            if (!(aii > 0)) {  // пустая строка: изолируем
                for (i64 a = off[i]; a < off[i + 1]; ++a) val[a] = 0.0;
                d[i] = 1.0;
                continue;
            }
            for (i64 a = off[i]; a < off[i + 1]; ++a) {
                int k = col[a];
                double sacc = 0;
                i64 p1 = off[i], p2 = off[k];
                while (p1 < a && p2 < off[k + 1]) {
                    int c1 = col[p1], c2 = col[p2];
                    if (c1 == c2) { sacc += val[p1] * val[p2]; ++p1; ++p2; }
                    else if (c1 < c2) ++p1; else ++p2;
                }
                val[a] = (val[a] - sacc) / d[k];
            }
            double s2 = aii * (1.0 + alpha);
            for (i64 a = off[i]; a < off[i + 1]; ++a) s2 -= val[a] * val[a];
            if (!(s2 > 1e-8 * aii) || !std::isfinite(s2)) return false;
            d[i] = std::sqrt(s2);
        }
        return true;
    }
    void apply(const std::vector<double>& r_in, std::vector<double>& z_out) const {
        tmp_r.resize(n);
        tmp_z.resize(n);
        for (int k = 0; k < n; ++k) tmp_r[k] = r_in[perm[k]];
        std::vector<double>& z = tmp_z;
        const std::vector<double>& r = tmp_r;
        // L y = r
        for (int i = 0; i < n; ++i) {
            double sacc = r[i];
            for (i64 a = off[i]; a < off[i + 1]; ++a) sacc -= val[a] * z[col[a]];
            z[i] = sacc / d[i];
        }
        // L^T z = y
        for (int i = n - 1; i >= 0; --i) {
            z[i] /= d[i];
            const double zi = z[i];
            for (i64 a = off[i]; a < off[i + 1]; ++a) z[col[a]] -= val[a] * zi;
        }
        for (int k = 0; k < n; ++k) z_out[perm[k]] = z[k];
    }
};

// Решает H x = g при фиксированных x[i] (mask). x — начальное приближение (на фиксированных —
// значения). Предобусловливатель: IC(0) (если ic != nullptr) или Якоби.
int cg_solve(const Csr& H, const std::vector<double>& g, const std::vector<uint8_t>& fixed,
             const std::vector<double>& invdiag, std::vector<double>& x, int max_iters, double tol,
             Ctx& ctx, double* res_out = nullptr, const IC0* ic = nullptr) {
    const int n = H.n;
    std::vector<double> r(n), z(n), p(n), Ap(n);
    parallel_for(n, [&](i64 b, i64 e) {
        for (i64 i = b; i < e; ++i) {
            double s = 0;
            for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) s += H.val[k] * x[H.col[k]];
            r[i] = fixed[i] ? 0.0 : g[i] - s;
            z[i] = r[i] * invdiag[i];
        }
    }, 1024);
    if (ic) {
        ic->apply(r, z);
        for (int i = 0; i < n; ++i) if (fixed[i]) z[i] = 0.0;
    }
    p = z;
    double rz = dotv(r, z);
    double gnorm = std::sqrt(parallel_sum(n, [&](i64 s, i64 e) {
        double v = 0;
        for (i64 i = s; i < e; ++i) if (!fixed[i]) v += g[i] * g[i];
        return v;
    })) + 1e-30;
    int it = 0;
    double rn2 = dotv(r, r);
    const i64 grain = 1024;
    const i64 nchunks = (n + grain - 1) / grain;
    std::vector<double> part1(nchunks), part2(nchunks), part3(nchunks);
    for (; it < max_iters; ++it) {
        if (std::sqrt(rn2) < tol * gnorm) break;
        parallel_for(nchunks, [&](i64 cb, i64 ce) {
            for (i64 c = cb; c < ce; ++c) {
                i64 b = c * grain, e = std::min<i64>(n, b + grain);
                double acc = 0;
                for (i64 i = b; i < e; ++i) {
                    double s = 0;
                    for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) s += H.val[k] * p[H.col[k]];
                    if (fixed[i]) s = 0.0;
                    Ap[i] = s;
                    acc += p[i] * s;
                }
                part1[c] = acc;
            }
        }, 1);
        double pAp = 0;
        for (double v : part1) pAp += v;
        if (!(pAp > 0)) {
            if (it == 0) ctx.logf("cg: остановка на первой итерации, pAp=%g rz=%g |r|=%g", pAp, rz, std::sqrt(rn2));
            break;
        }
        const double alpha = rz / pAp;
        double rz2 = 0;
        rn2 = 0;
        if (ic) {
            parallel_for(nchunks, [&](i64 cb, i64 ce) {
                for (i64 c = cb; c < ce; ++c) {
                    i64 b = c * grain, e = std::min<i64>(n, b + grain);
                    double a3 = 0;
                    for (i64 i = b; i < e; ++i) {
                        x[i] += alpha * p[i];
                        r[i] -= alpha * Ap[i];
                        a3 += r[i] * r[i];
                    }
                    part3[c] = a3;
                }
            }, 1);
            ic->apply(r, z);
            for (i64 c = 0; c < nchunks; ++c) rn2 += part3[c];
            for (int i = 0; i < n; ++i) { if (fixed[i]) z[i] = 0.0; rz2 += r[i] * z[i]; }
        } else {
            parallel_for(nchunks, [&](i64 cb, i64 ce) {
                for (i64 c = cb; c < ce; ++c) {
                    i64 b = c * grain, e = std::min<i64>(n, b + grain);
                    double a2 = 0, a3 = 0;
                    for (i64 i = b; i < e; ++i) {
                        x[i] += alpha * p[i];
                        r[i] -= alpha * Ap[i];
                        z[i] = r[i] * invdiag[i];
                        a2 += r[i] * z[i];
                        a3 += r[i] * r[i];
                    }
                    part2[c] = a2;
                    part3[c] = a3;
                }
            }, 1);
            for (i64 c = 0; c < nchunks; ++c) { rz2 += part2[c]; rn2 += part3[c]; }
        }
        const double beta = rz2 / rz;
        rz = rz2;
        parallel_for(n, [&](i64 b, i64 e) {
            for (i64 i = b; i < e; ++i) p[i] = z[i] + beta * p[i];
        }, 4096);
        if ((it & 63) == 63) ctx.check_cancel();
    }
    if (res_out) *res_out = std::sqrt(rn2) / gnorm;
    return it;
}

// ---------------------------------------------------------------------------
// Прямой решатель: переменные делятся на G (целочисленные и с плотными строками) и L (остальные).
// H_LL факторизуется разреженным Холецким (вложенное рассечение), по G строится плотное
// дополнение Шура S = H_GG - H_GL H_LL^{-1} H_LG. Округление меняет только набор
// фиксированных переменных G, поэтому каждый раунд — плотный Холецкий подматрицы S.
// ---------------------------------------------------------------------------
struct SchurSolver {
    int n = 0, nl = 0, ng = 0;
    std::vector<int> loc, lvars, gvars;
    std::vector<uint8_t> ing;
    SparseCholesky chol;
    std::vector<i64> lg_off;     // H_LG по строкам L: (локальный индекс G, значение)
    std::vector<int> lg_col;
    std::vector<double> lg_val;
    std::vector<double> S, s, gl;
    std::vector<double> dense, rhs, xl, work;
    bool ok = false;

    bool setup(const Csr& H, const std::vector<double>& g, const std::vector<uint8_t>& in_g,
               const std::vector<double>& coords, i64 max_fill, Ctx& ctx) {
        const double t0 = now_seconds();
        ok = false;
        n = H.n;
        ing = in_g;
        loc.assign(n, -1);
        lvars.clear();
        gvars.clear();
        nl = ng = 0;
        for (int i = 0; i < n; ++i) {
            if (ing[i]) { loc[i] = ng++; gvars.push_back(i); }
            else { loc[i] = nl++; lvars.push_back(i); }
        }
        SymCsr A;
        A.n = nl;
        A.off.assign(nl + 1, 0);
        lg_off.assign(nl + 1, 0);
        for (int l = 0; l < nl; ++l) {
            const int i = lvars[l];
            i64 a = 0, b = 0;
            for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) (ing[H.col[k]] ? b : a) += 1;
            A.off[l + 1] = A.off[l] + a;
            lg_off[l + 1] = lg_off[l] + b;
        }
        A.col.resize(A.off[nl]);
        A.val.resize(A.off[nl]);
        lg_col.resize(lg_off[nl]);
        lg_val.resize(lg_off[nl]);
        double diag_sum = 0;
        for (int l = 0; l < nl; ++l) {
            const int i = lvars[l];
            i64 a = A.off[l], b = lg_off[l];
            for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) {
                const int c = H.col[k];
                if (ing[c]) { lg_col[b] = loc[c]; lg_val[b] = H.val[k]; ++b; }
                else {
                    A.col[a] = loc[c];
                    A.val[a] = H.val[k];
                    if (c == i) diag_sum += H.val[k];
                    ++a;
                }
            }
        }
        std::vector<double> cl(3 * (size_t)nl);
        for (int l = 0; l < nl; ++l)
            for (int q = 0; q < 3; ++q) cl[3 * (size_t)l + q] = coords[3 * (size_t)lvars[l] + q];
        std::vector<int> perm = nested_dissection_order(A, cl, 64);
        const double t1 = now_seconds();
        // малая регуляризация: H вырождена по сдвигам компонент без целочисленных связей
        double reg = nl > 0 ? 1e-10 * diag_sum / nl : 0.0;
        bool fok = nl == 0;
        for (int attempt = 0; attempt < 3 && !fok; ++attempt, reg *= 1e3) fok = chol.factor(A, perm, reg, max_fill);
        if (!fok) {
            ctx.logf("param: прямой решатель: разложение H_LL не удалось (L=%d, заполнение %lld)", nl,
                     (long long)chol.nnz());
            return false;
        }
        const double t2 = now_seconds();
        // столбцы H_LG по G
        std::vector<i64> gc_off(ng + 1, 0);
        for (i64 p = 0; p < lg_off[nl]; ++p) ++gc_off[lg_col[p] + 1];
        for (int a = 0; a < ng; ++a) gc_off[a + 1] += gc_off[a];
        std::vector<int> gc_row(gc_off[ng]);
        std::vector<double> gc_val(gc_off[ng]);
        {
            std::vector<i64> fill(gc_off.begin(), gc_off.end() - 1);
            for (int l = 0; l < nl; ++l)
                for (i64 p = lg_off[l]; p < lg_off[l + 1]; ++p) {
                    i64 o = fill[lg_col[p]]++;
                    gc_row[o] = l;
                    gc_val[o] = lg_val[p];
                }
        }
        S.assign((size_t)ng * ng, 0.0);
        for (int a = 0; a < ng; ++a) {
            const int i = gvars[a];
            for (i64 k = H.off[i]; k < H.off[i + 1]; ++k)
                if (ing[H.col[k]]) S[(size_t)a * ng + loc[H.col[k]]] += H.val[k];
        }
        std::atomic<bool> cancelled{false};
        parallel_for(ng, [&](i64 b, i64 e) {
            std::vector<double> z(nl), w(nl);
            for (i64 a = b; a < e; ++a) {
                if (cancelled.load(std::memory_order_relaxed)) return;
                std::fill(z.begin(), z.end(), 0.0);
                for (i64 p = gc_off[a]; p < gc_off[a + 1]; ++p) z[gc_row[p]] = gc_val[p];
                chol.solve_inplace(z.data(), w.data());
                double* Sa = &S[(size_t)a * ng];
                for (int l = 0; l < nl; ++l) {
                    const double zl = z[l];
                    if (zl == 0.0) continue;
                    for (i64 p = lg_off[l]; p < lg_off[l + 1]; ++p) Sa[lg_col[p]] -= lg_val[p] * zl;
                }
                if ((a & 15) == 0 && ctx.st && ctx.st->cancel) cancelled = true;
            }
        }, 1);
        ctx.check_cancel();
        for (int a = 0; a < ng; ++a)
            for (int c = 0; c < a; ++c) {
                const double v = 0.5 * (S[(size_t)a * ng + c] + S[(size_t)c * ng + a]);
                S[(size_t)a * ng + c] = S[(size_t)c * ng + a] = v;
            }
        gl.assign(nl, 0.0);
        for (int l = 0; l < nl; ++l) gl[l] = g[lvars[l]];
        s.assign(ng, 0.0);
        for (int a = 0; a < ng; ++a) s[a] = g[gvars[a]];
        if (nl > 0) {
            std::vector<double> z0 = gl, w(nl);
            chol.solve_inplace(z0.data(), w.data());
            for (int l = 0; l < nl; ++l)
                for (i64 p = lg_off[l]; p < lg_off[l + 1]; ++p) s[lg_col[p]] -= lg_val[p] * z0[l];
        }
        ctx.logf("param: прямой решатель: L=%d (заполнение %lld, упорядочение %.2f с, разложение %.2f с), "
                 "G=%d (Шур %.2f с)", nl, (long long)chol.nnz(), t1 - t0, t2 - t1, ng, now_seconds() - t2);
        ok = true;
        return true;
    }

    // Ремонт: H изменилась локально (веса граней), все переменные G зафиксированы — PCG по L
    // с прежним разложением H_LL в роли предобусловливателя. Возвращает число итераций или -1.
    int refine(const Csr& H, const std::vector<double>& g, const std::vector<uint8_t>& fixed,
               std::vector<double>& x, int max_it, double tol) {
        for (int a = 0; a < ng; ++a)
            if (!fixed[gvars[a]]) return -1;
        if (nl == 0) return 0;
        std::vector<double> r(nl), z(nl), p(nl), Ap(nl), w(nl), xl0(nl);
        // r = g_L - H_L* x (по всем столбцам, G — зафиксированные значения)
        parallel_for(nl, [&](i64 b, i64 e) {
            for (i64 l = b; l < e; ++l) {
                const int i = lvars[l];
                double s = 0;
                for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) s += H.val[k] * x[H.col[k]];
                r[l] = g[i] - s;
                xl0[l] = x[i];
            }
        }, 512);
        double gn = 0;
        for (int l = 0; l < nl; ++l) gn += g[lvars[l]] * g[lvars[l]];
        gn = std::sqrt(std::max(gn, 1e-300));
        auto apply_LL = [&](const std::vector<double>& v, std::vector<double>& out) {
            parallel_for(nl, [&](i64 b, i64 e) {
                for (i64 l = b; l < e; ++l) {
                    const int i = lvars[l];
                    double s = 0;
                    for (i64 k = H.off[i]; k < H.off[i + 1]; ++k)
                        if (!ing[H.col[k]]) s += H.val[k] * v[loc[H.col[k]]];
                    out[l] = s;
                }
            }, 512);
        };
        z = r;
        chol.solve_inplace(z.data(), w.data());
        p = z;
        double rz = 0;
        for (int l = 0; l < nl; ++l) rz += r[l] * z[l];
        std::vector<double> xl = xl0;
        for (int it = 0; it < max_it; ++it) {
            double rn = 0;
            for (int l = 0; l < nl; ++l) rn += r[l] * r[l];
            if (std::sqrt(rn) <= tol * gn) {
                for (int l = 0; l < nl; ++l) x[lvars[l]] = xl[l];
                return it;
            }
            apply_LL(p, Ap);
            double pAp = 0;
            for (int l = 0; l < nl; ++l) pAp += p[l] * Ap[l];
            if (!(pAp > 0)) break;
            const double alpha = rz / pAp;
            for (int l = 0; l < nl; ++l) { xl[l] += alpha * p[l]; r[l] -= alpha * Ap[l]; }
            z = r;
            chol.solve_inplace(z.data(), w.data());
            double rz2 = 0;
            for (int l = 0; l < nl; ++l) rz2 += r[l] * z[l];
            const double beta = rz2 / rz;
            rz = rz2;
            for (int l = 0; l < nl; ++l) p[l] = z[l] + beta * p[l];
        }
        return -1;
    }

    // Минимизация при фиксированных переменных (только из G); x — на фиксированных значения.
    bool solve(const std::vector<uint8_t>& fixed, std::vector<double>& x) {
        std::vector<int> P;
        P.reserve(ng);
        for (int a = 0; a < ng; ++a)
            if (!fixed[gvars[a]]) P.push_back(a);
        const int m = (int)P.size();
        if (m > 0) {
            rhs.assign(m, 0.0);
            parallel_each(m, [&](i64 q) {
                const double* Sa = &S[(size_t)P[q] * ng];
                double v = s[P[q]];
                for (int c = 0; c < ng; ++c)
                    if (fixed[gvars[c]]) v -= Sa[c] * x[gvars[c]];
                rhs[q] = v;
            }, 64);
            double dmax = 0;
            for (int q = 0; q < m; ++q) dmax = std::max(dmax, S[(size_t)P[q] * ng + P[q]]);
            // S вырождена по сдвигу параметризации, пока целые не зафиксированы: без регуляризации
            // крошечный ведущий элемент раздувает калибровочную компоненту и губит точность x_L.
            double reg = 1e-10 * std::max(dmax, 1e-300);
            bool fok = false;
            for (int attempt = 0; attempt < 4 && !fok; ++attempt) {
                dense.resize((size_t)m * m);
                parallel_each(m, [&](i64 q) {
                    const double* Sa = &S[(size_t)P[q] * ng];
                    double* row = &dense[(size_t)q * m];
                    for (int r = 0; r <= q; ++r) row[r] = Sa[P[r]];
                    row[q] += reg;
                }, 64);
                fok = dense_cholesky(dense.data(), m);
                reg *= 1e3;
            }
            if (!fok) return false;
            dense_cholesky_solve(dense.data(), m, rhs.data());
            for (int q = 0; q < m; ++q) x[gvars[P[q]]] = rhs[q];
        }
        if (nl > 0) {
            xl.resize(nl);
            work.resize(nl);
            parallel_each(nl, [&](i64 l) {
                double v = gl[l];
                for (i64 p = lg_off[l]; p < lg_off[l + 1]; ++p) v -= lg_val[p] * x[gvars[lg_col[p]]];
                xl[l] = v;
            }, 1024);
            chol.solve_inplace(xl.data(), work.data());
            for (int l = 0; l < nl; ++l) x[lvars[l]] = xl[l];
        }
        return true;
    }
};

inline void rot_entries(int k, double& r00, double& r01, double& r10, double& r11) {
    switch (k & 3) {
        case 0: r00 = 1; r01 = 0; r10 = 0; r11 = 1; break;
        case 1: r00 = 0; r01 = 1; r10 = -1; r11 = 0; break;
        case 2: r00 = -1; r01 = 0; r10 = 0; r11 = -1; break;
        default: r00 = 0; r01 = -1; r10 = 1; r11 = 0; break;
    }
}

}  // namespace

bool compute_param(const ParamInput& in, ParamResult& out, Ctx& ctx) {
    const TriMesh& M = *in.M;
    const FieldGeometry& G = *in.G;
    const int nf = M.nF(), nv = M.nV();
    out = ParamResult();
    double t0 = now_seconds();

    // ---------------------------------------------------------------- комбинг
    std::vector<double> theta0(nf);
    for (int f = 0; f < nf; ++f) theta0[f] = 0.25 * (*in.z)[f].arg();
    std::vector<double>& rep = out.rep;
    rep.assign(nf, 0.0);
    std::vector<uint8_t> tree(3 * (size_t)nf, 0);
    {
        std::vector<uint8_t> vis(nf, 0);
        std::vector<int> queue;
        queue.reserve(nf);
        for (int s = 0; s < nf; ++s) {
            if (vis[s]) continue;
            vis[s] = 1;
            rep[s] = theta0[s];
            size_t qh = queue.size();
            queue.push_back(s);
            while (qh < queue.size()) {
                int f = queue[qh++];
                for (int i = 0; i < 3; ++i) {
                    int g = M.FF[3 * f + i];
                    if (g < 0 || vis[g]) continue;
                    int j = M.FE[3 * f + i];
                    double pred = rep[f] + G.transport[3 * g + j];
                    double k0 = std::round((pred - theta0[g]) / kHalfPi);
                    rep[g] = theta0[g] + k0 * kHalfPi;
                    vis[g] = 1;
                    tree[3 * f + i] = 1;
                    tree[3 * g + j] = 1;
                    queue.push_back(g);
                }
            }
        }
    }
    out.kjump.assign(3 * (size_t)nf, 0);
    parallel_each(nf, [&](i64 f) {
        for (int i = 0; i < 3; ++i) {
            int g = M.FF[3 * f + i];
            if (g < 0 || g < f) continue;
            int j = M.FE[3 * f + i];
            double d = (rep[g] - rep[f] - G.transport[3 * g + j]) / kHalfPi;
            int k = ((int)std::llround(d)) & 3;
            out.kjump[3 * f + i] = (uint8_t)k;
            out.kjump[3 * g + j] = (uint8_t)((4 - k) & 3);
        }
    });

    // ------------------------------------------------- интегрируемость поля
    if (in.integrable > 1e-9) {
        // min Σ w (α_g - α_f - δ_e)^2, α = 0 на зафиксированных гранях
        std::vector<double> alpha(nf, 0.0), diag(nf, 0.0), rhs(nf, 0.0);
        std::vector<std::array<double, 3>> wv(nf), dv(nf);
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) {
                wv[f][i] = 0;
                dv[f][i] = 0;
                int g = M.FF[3 * f + i];
                if (g < 0) continue;
                int j = M.FE[3 * f + i];
                double delta = rep[g] - rep[f] - G.transport[3 * g + j];
                delta -= std::round(delta / kHalfPi) * kHalfPi;
                double w = G.weight[3 * f + i];
                wv[f][i] = w;
                dv[f][i] = delta;  // поворот f->g
                diag[f] += w;
                rhs[f] -= w * delta;  // d/dα_f (α_g-α_f-δ)^2 -> -(α_g-α_f-δ)
            }
        auto isfix = [&](int f) { return in.fixed_faces && (*in.fixed_faces)[f]; };
        auto apply = [&](const std::vector<double>& x, std::vector<double>& y) {
            parallel_each(nf, [&](i64 f) {
                if (isfix((int)f)) { y[f] = x[f]; return; }
                double sacc = (diag[f] + 1e-9) * x[f];
                for (int i = 0; i < 3; ++i) {
                    int g = M.FF[3 * f + i];
                    if (g < 0 || isfix(g)) continue;
                    sacc -= wv[f][i] * x[g];
                }
                y[f] = sacc;
            }, 1024);
        };
        std::vector<double> b(nf);
        for (int f = 0; f < nf; ++f) b[f] = isfix(f) ? 0.0 : rhs[f];
        std::vector<double> r(b), z(nf), pd(nf), Ap(nf);
        for (int f = 0; f < nf; ++f) z[f] = r[f] / (diag[f] + 1e-9);
        pd = z;
        double rz = dotv(r, z), r0 = std::sqrt(dotv(r, r)) + 1e-30;
        int it = 0;
        for (; it < 3000; ++it) {
            apply(pd, Ap);
            double pAp = dotv(pd, Ap);
            if (!(pAp > 0)) break;
            double al = rz / pAp;
            for (int f = 0; f < nf; ++f) { alpha[f] += al * pd[f]; r[f] -= al * Ap[f]; z[f] = r[f] / (diag[f] + 1e-9); }
            if (std::sqrt(dotv(r, r)) < 1e-7 * r0) break;
            double rz2 = dotv(r, z);
            double be = rz2 / rz;
            rz = rz2;
            for (int f = 0; f < nf; ++f) pd[f] = z[f] + be * pd[f];
        }
        double amax = 0;
        for (int f = 0; f < nf; ++f) {
            rep[f] -= in.integrable * alpha[f];
            amax = std::max(amax, std::fabs(alpha[f]));
        }
        ctx.logf("param: интегрируемость: поворот поля до %.1f°, CG %d", amax * 180.0 / kPi, it);
    }

    // ------------------------------------------------- конформная поправка размеров
    std::vector<double> size_eff(in.size->begin(), in.size->end());
    if (std::fabs(in.conformal) > 1e-12) {
        // φ_b - φ_a = δ_e * ℓ/ℓ*, вес ℓ*/ℓ (сопряжённая гармоника угла поля)
        std::vector<V3> C(nf);
        for (int f = 0; f < nf; ++f)
            C[f] = (M.P[M.F[3 * f]] + M.P[M.F[3 * f + 1]] + M.P[M.F[3 * f + 2]]) * (1.0 / 3.0);
        struct Eq { int a, b; double w, rhs; };
        std::vector<Eq> eqs;
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) {
                int g = M.FF[3 * f + i];
                if (g < 0 || g < f) continue;
                int j = M.FE[3 * f + i];
                double k = out.kjump[3 * f + i];
                double delta = rep[g] - rep[f] - G.transport[3 * g + j];
                delta -= std::round(delta / kHalfPi) * kHalfPi;  // остаток после периодического скачка
                (void)k;
                int a = M.F[3 * f + i], b = M.F[3 * f + (i + 1) % 3];
                double l = norm(M.P[b] - M.P[a]);
                double ls = std::max(norm(C[g] - C[f]), 1e-12 * l);
                eqs.push_back({a, b, ls / std::max(l, 1e-300), delta * l / ls});
            }
        // нормальные уравнения: L φ = r (CG Якоби)
        std::vector<double> diag(nv, 0.0), rhs(nv, 0.0), phi(nv, 0.0);
        std::vector<std::vector<std::pair<int, double>>> adj(nv);
        for (auto& e : eqs) {
            diag[e.a] += e.w; diag[e.b] += e.w;
            adj[e.a].push_back({e.b, e.w}); adj[e.b].push_back({e.a, e.w});
            rhs[e.b] += e.w * e.rhs; rhs[e.a] -= e.w * e.rhs;
        }
        auto apply = [&](const std::vector<double>& x, std::vector<double>& y) {
            parallel_each(nv, [&](i64 v) {
                double s2 = (diag[v] + 1e-9) * x[v];
                for (auto& q : adj[v]) s2 -= q.second * x[q.first];
                y[v] = s2;
            }, 1024);
        };
        std::vector<double> r(rhs), z(nv), pdir(nv), Ap(nv);
        for (int v = 0; v < nv; ++v) z[v] = r[v] / (diag[v] + 1e-9);
        pdir = z;
        double rz = dotv(r, z), r0 = std::sqrt(dotv(r, r)) + 1e-30;
        for (int it = 0; it < 2000; ++it) {
            apply(pdir, Ap);
            double pAp = dotv(pdir, Ap);
            if (!(pAp > 0)) break;
            double al = rz / pAp;
            for (int v = 0; v < nv; ++v) { phi[v] += al * pdir[v]; r[v] -= al * Ap[v]; z[v] = r[v] / (diag[v] + 1e-9); }
            double rz2 = dotv(r, z);
            if (std::sqrt(dotv(r, r)) < 1e-6 * r0) break;
            double be = rz2 / rz;
            rz = rz2;
            for (int v = 0; v < nv; ++v) pdir[v] = z[v] + be * pdir[v];
        }
        // среднее по площади -> 0, ограничение диапазона
        double acc = 0, wa = 0;
        std::vector<double> phif(nf);
        for (int f = 0; f < nf; ++f) {
            phif[f] = (phi[M.F[3 * f]] + phi[M.F[3 * f + 1]] + phi[M.F[3 * f + 2]]) / 3.0;
            acc += phif[f] * M.FA[f];
            wa += M.FA[f];
        }
        acc /= std::max(wa, 1e-300);
        double cnt0 = 0, cnt1 = 0, pmin = 1e300, pmax = -1e300;
        for (int f = 0; f < nf; ++f) {
            double ph = in.conformal * (phif[f] - acc);
            ph = std::max(-in.conformal_clamp, std::min(in.conformal_clamp, ph));
            pmin = std::min(pmin, ph);
            pmax = std::max(pmax, ph);
            size_eff[f] = (*in.size)[f] * std::exp(ph);
            cnt0 += M.FA[f] / ((*in.size)[f] * (*in.size)[f]);
            cnt1 += M.FA[f] / (size_eff[f] * size_eff[f]);
        }
        double sc = std::sqrt(cnt1 / std::max(cnt0, 1e-300));
        for (double& v : size_eff) v *= sc;
        out.phi.resize(nf);
        for (int f = 0; f < nf; ++f) out.phi[f] = phif[f] - acc;
        ctx.logf("param: конформная поправка размеров: log-масштаб от %.2f до %.2f", pmin, pmax);
    }

    // ------------------------------------------------- сингулярности по вершинам
    out.vindex.assign(nv, 0);
    std::vector<uint8_t> singular(nv, 0);
    parallel_each(nv, [&](i64 v) {
        if (M.vbound[v]) { out.vindex[v] = -1; return; }
        int h = 0;
        M.for_fan((int)v, [&](int f, int c) { h += out.kjump[3 * f + c]; });
        h &= 3;
        out.vindex[v] = (int8_t)h;
        singular[v] = h != 0;
    }, 512);
    for (int v = 0; v < nv; ++v) out.nsingular += singular[v];

    // ------------------------------------------------------------ граф разреза
    std::vector<uint8_t>& cut = out.cut;
    cut.assign(3 * (size_t)nf, 0);
    std::vector<int> cutdeg(nv, 0);
    for (int f = 0; f < nf; ++f)
        for (int i = 0; i < 3; ++i) {
            int g = M.FF[3 * f + i];
            if (g < 0 || tree[3 * f + i] || g < f) continue;
            int j = M.FE[3 * f + i];
            cut[3 * f + i] = 1;
            cut[3 * g + j] = 1;
            ++cutdeg[M.F[3 * f + i]];
            ++cutdeg[M.F[3 * f + (i + 1) % 3]];
        }
    {
        // обрезка листьев в несингулярных внутренних вершинах
        std::vector<int> stack;
        auto is_prunable = [&](int v) { return !M.vbound[v] && !singular[v] && cutdeg[v] == 1; };
        for (int v = 0; v < nv; ++v)
            if (is_prunable(v)) stack.push_back(v);
        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            if (!is_prunable(v)) continue;
            // найти единственное разрезанное ребро в веере v
            int ef = -1, ei = -1;
            M.for_fan(v, [&](int f, int c) {
                if (cut[3 * f + c]) { ef = f; ei = c; }
            });
            if (ef < 0) { cutdeg[v] = 0; continue; }
            int g = M.FF[3 * ef + ei], j = M.FE[3 * ef + ei];
            int w = M.F[3 * ef + (ei + 1) % 3];
            cut[3 * ef + ei] = 0;
            cut[3 * g + j] = 0;
            --cutdeg[v];
            --cutdeg[w];
            if (is_prunable(w)) stack.push_back(w);
        }
    }
    // ------------------------------------------------------------ сектора
    std::vector<int> corner_copy(3 * (size_t)nf, -1);
    int ncopies = 0;
    std::vector<uint8_t> copy_singular;  // копия принадлежит сингулярной вершине
    std::vector<int> copy_vertex;
    for (int v = 0; v < nv; ++v) {
        int b = M.vfan_off[v], e = M.vfan_off[v + 1];
        int m = e - b;
        if (m == 0) continue;
        auto cutflag = [&](int t) {
            int f = M.vfan_face[b + t], c = M.vfan_corner[b + t];
            return (int)cut[3 * f + c];
        };
        auto new_copy = [&]() {
            copy_singular.push_back(singular[v]);
            copy_vertex.push_back(v);
            return ncopies++;
        };
        if (!M.vbound[v]) {
            int ncut = 0, first = -1;
            for (int t = 0; t < m; ++t)
                if (cutflag(t)) { ++ncut; if (first < 0) first = t; }
            if (ncut == 0) {
                int c = new_copy();
                for (int t = 0; t < m; ++t) corner_copy[3 * M.vfan_face[b + t] + M.vfan_corner[b + t]] = c;
            } else {
                int t0s = (first + 1) % m;
                int c = new_copy();
                for (int s = 0; s < m; ++s) {
                    int t = (t0s + s) % m;
                    corner_copy[3 * M.vfan_face[b + t] + M.vfan_corner[b + t]] = c;
                    if (cutflag(t) && s < m - 1) c = new_copy();
                }
            }
        } else {
            int c = new_copy();
            for (int t = 0; t < m; ++t) {
                corner_copy[3 * M.vfan_face[b + t] + M.vfan_corner[b + t]] = c;
                if (t < m - 1 && cutflag(t)) c = new_copy();
            }
        }
    }
    // индексы разрезанных рёбер
    std::vector<int> cut_id(3 * (size_t)nf, -1);
    int ncut_edges = 0;
    for (int f = 0; f < nf; ++f)
        for (int i = 0; i < 3; ++i) {
            int g = M.FF[3 * f + i];
            if (!cut[3 * f + i] || g < f) continue;
            cut_id[3 * f + i] = ncut_edges;
            cut_id[3 * g + M.FE[3 * f + i]] = ncut_edges;
            ++ncut_edges;
        }
    out.ncut_edges = ncut_edges;
    const int nslots = 2 * ncopies + 2 * ncut_edges;
    ctx.logf("param: комбинг+разрез: %d сингулярностей, %d рёбер разреза, %d копий вершин (%.2f с)",
             out.nsingular, ncut_edges, ncopies, now_seconds() - t0);

    // ------------------------------------------- приоритеты и «обязательно целые»
    std::vector<uint8_t> prio(nslots, 0), required(nslots, 0);
    for (int c = 0; c < ncopies; ++c)
        if (copy_singular[c]) {
            required[2 * c] = required[2 * c + 1] = 1;
            prio[2 * c] = prio[2 * c + 1] = 2;
        }
    for (int e = 0; e < ncut_edges; ++e) {
        prio[2 * ncopies + 2 * e] = prio[2 * ncopies + 2 * e + 1] = 1;
        required[2 * ncopies + 2 * e] = required[2 * ncopies + 2 * e + 1] = 1;
    }
    // фичери: какая координата постоянна вдоль ребра
    struct FeatEq { int slot_a, slot_b; };
    std::vector<FeatEq> feat_eqs;
    if (in.feature) {
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) {
                if (!(*in.feature)[3 * f + i]) continue;
                V3 e = M.P[M.F[3 * f + (i + 1) % 3]] - M.P[M.F[3 * f + i]];
                double alpha = angle_in_face(M, f, e);
                double d = (alpha - rep[f]) / kHalfPi;
                double m = std::round(d);
                if (std::fabs(d - m) * kHalfPi > 0.35) continue;  // не выровнено
                int mi = ((int)m) & 3;
                int axis = (mi & 1) ? 0 : 1;  // e ∥ d (u-направление) => v постоянна
                int ca = corner_copy[3 * f + i], cb = corner_copy[3 * f + (i + 1) % 3];
                int sa = 2 * ca + axis, sb = 2 * cb + axis;
                feat_eqs.push_back({sa, sb});
                required[sa] = required[sb] = 1;
                if (prio[sa] < 2) prio[sa] = 2;
                if (prio[sb] < 2) prio[sb] = 2;
            }
    }

    // ------------------------------------------------------------ уравнения
    Eliminator elim(nslots, prio);
    {
        std::vector<Term> eq;
        for (int f = 0; f < nf; ++f)
            for (int i = 0; i < 3; ++i) {
                int g = M.FF[3 * f + i];
                if (!cut[3 * f + i] || g < f) continue;
                int j = M.FE[3 * f + i];
                int k = out.kjump[3 * f + i];
                double r00, r01, r10, r11;
                rot_entries(k, r00, r01, r10, r11);
                int ce = cut_id[3 * f + i];
                int ts = 2 * ncopies + 2 * ce;
                int fa = corner_copy[3 * f + i], fb = corner_copy[3 * f + (i + 1) % 3];
                int ga = corner_copy[3 * g + (j + 1) % 3], gb = corner_copy[3 * g + j];
                for (int p = 0; p < 2; ++p) {
                    int fc = p == 0 ? fa : fb;
                    int gc = p == 0 ? ga : gb;
                    // u
                    eq.clear();
                    eq.push_back({2 * gc, 1.0});
                    eq.push_back({2 * fc, -r00});
                    eq.push_back({2 * fc + 1, -r01});
                    eq.push_back({ts, -1.0});
                    elim.add_equation(eq);
                    // v
                    eq.clear();
                    eq.push_back({2 * gc + 1, 1.0});
                    eq.push_back({2 * fc, -r10});
                    eq.push_back({2 * fc + 1, -r11});
                    eq.push_back({ts + 1, -1.0});
                    elim.add_equation(eq);
                }
            }
        for (const FeatEq& fe : feat_eqs) {
            if (fe.slot_a == fe.slot_b) continue;
            eq.clear();
            eq.push_back({fe.slot_a, 1.0});
            eq.push_back({fe.slot_b, -1.0});
            elim.add_equation(eq);
        }
    }
    elim.finalize();
    // свободные переменные
    std::vector<int> free_idx(nslots, -1);
    int nfree = 0;
    for (int s = 0; s < nslots; ++s)
        if (elim.is_free(s)) free_idx[s] = nfree++;
    out.n_free = nfree;
    // разложения всех слотов
    std::vector<i64> sl_off(nslots + 1, 0);
    std::vector<int> sl_var;
    std::vector<double> sl_coef;
    {
        std::vector<std::vector<Term>> tmp(nslots);
        parallel_each(nslots, [&](i64 s) { elim.resolved((int)s, tmp[s]); }, 512);
        for (int s = 0; s < nslots; ++s) sl_off[s + 1] = sl_off[s] + (i64)tmp[s].size();
        sl_var.resize(sl_off[nslots]);
        sl_coef.resize(sl_off[nslots]);
        for (int s = 0; s < nslots; ++s) {
            i64 o = sl_off[s];
            for (const Term& t : tmp[s]) {
                sl_var[o] = free_idx[t.var];
                sl_coef[o] = t.coef;
                ++o;
            }
        }
    }
    // отладка: «чужие» ссылки — слот копии вершины v зависит от свободной переменной-копии другой вершины
    if (in.verbose) {
        std::vector<int> free_owner(nfree, -1);  // вершина-владелец свободной переменной (копийный слот)
        for (int sl = 0; sl < 2 * ncopies; ++sl)
            if (free_idx[sl] >= 0) free_owner[free_idx[sl]] = copy_vertex[sl / 2];
        int foreign = 0, foreign_regular = 0;
        std::vector<int> fv_count(nv, 0);
        for (int sl = 0; sl < 2 * ncopies; ++sl) {
            int v = copy_vertex[sl / 2];
            for (i64 k = sl_off[sl]; k < sl_off[sl + 1]; ++k) {
                int fvv = sl_var[k];
                if (free_owner[fvv] >= 0 && free_owner[fvv] != v) {
                    ++foreign;
                    if (!singular[v] && !singular[free_owner[fvv]]) ++foreign_regular;
                    ++fv_count[v];
                    break;
                }
            }
        }
        out.foreign.assign(fv_count.begin(), fv_count.end());
        ctx.logf("param: отладка: слотов с чужими ссылками %d (между регулярными вершинами %d)", foreign, foreign_regular);
    }
    // целочисленные свободные переменные: замыкание по required
    std::vector<uint8_t> is_int(nfree, 0);
    for (int s = 0; s < nslots; ++s)
        if (required[s])
            for (i64 k = sl_off[s]; k < sl_off[s + 1]; ++k) is_int[sl_var[k]] = 1;
    int nint = 0;
    for (int i = 0; i < nfree; ++i) nint += is_int[i];
    out.n_int = nint;
    ctx.logf("param: %d переменных, %d свободных, %d целочисленных (%.2f с)", nslots, nfree, nint,
             now_seconds() - t0);
    ctx.check_cancel();

    // ------------------------------------------------- локальные системы по граням
    std::vector<i64> fv_off(nf + 1, 0);
    std::vector<int> fv_var;      // локальные переменные грани
    std::vector<i64> fk_off(nf + 1, 0);
    std::vector<double> fk;       // локальные матрицы (|S|^2)
    std::vector<double> fg;       // локальные правые части (|S|)
    std::vector<double> face_c0(nf, 0.0);  // константа энергии
    std::vector<double> wmul(nf, 1.0);     // множители весов граней (ремонт перевёрнутых)

    // Локальная система грани: переменные (динамически), K = Σ w a a^T, g = Σ w τ a.
    struct FaceLocal {
        int nvar = 0;
        std::vector<int> vars;
        std::vector<double> K, g, rows;  // rows: 4 x nvar
        double c0 = 0;
    };
    auto build_face = [&](int f, FaceLocal& L) {
        V3 e1, e2;
        M.face_frame(f, e1, e2);
        const int* v = &M.F[3 * f];
        double px[3], py[3];
        for (int c = 0; c < 3; ++c) {
            V3 d = M.P[v[c]] - M.P[v[0]];
            px[c] = dot(d, e1);
            py[c] = dot(d, e2);
        }
        double area2 = (px[1] - px[0]) * (py[2] - py[0]) - (px[2] - px[0]) * (py[1] - py[0]);
        // вырожденные грани (нулевая площадь / неверная нормаль) не участвуют в энергии
        const bool degenerate = !(area2 > 1e-14 * (norm2(M.P[v[1]] - M.P[v[0]]) + norm2(M.P[v[2]] - M.P[v[0]])));
        double A = degenerate ? 0.0 : 0.5 * area2 * wmul[f];
        // Иглы (почти нулевая площадь) дают жёсткость ~ e^2/площадь и делают систему плохо
        // обусловленной: оператор градиента для них считаем по ограниченной снизу площади.
        const double es = norm2(M.P[v[1]] - M.P[v[0]]) + norm2(M.P[v[2]] - M.P[v[1]]) + norm2(M.P[v[0]] - M.P[v[2]]);
        double invA2 = degenerate ? 0.0 : 1.0 / std::max(area2, in.sliver_clamp * es);
        double gN[3][2];
        for (int c = 0; c < 3; ++c) {
            int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
            gN[c][0] = (py[c1] - py[c2]) * invA2;
            gN[c][1] = (px[c2] - px[c1]) * invA2;
        }
        double sc = 1.0 / size_eff[f];
        double th = rep[f];
        double tx[2][2] = {{sc * std::cos(th), sc * std::sin(th)}, {-sc * std::sin(th), sc * std::cos(th)}};
        L.vars.clear();
        auto var_pos = [&](int var) {
            for (int i = 0; i < (int)L.vars.size(); ++i)
                if (L.vars[i] == var) return i;
            L.vars.push_back(var);
            return (int)L.vars.size() - 1;
        };
        for (int c = 0; c < 3; ++c) {
            int copy = corner_copy[3 * f + c];
            for (int a = 0; a < 2; ++a) {
                int slot = 2 * copy + a;
                for (i64 k = sl_off[slot]; k < sl_off[slot + 1]; ++k) var_pos(sl_var[k]);
            }
        }
        const int n = (int)L.vars.size();
        L.nvar = n;
        L.rows.assign(4 * (size_t)n, 0.0);
        double rhs[4];
        for (int a = 0; a < 2; ++a)
            for (int q = 0; q < 2; ++q) {
                int r = 2 * a + q;
                rhs[r] = tx[a][q];
                for (int c = 0; c < 3; ++c) {
                    int copy = corner_copy[3 * f + c];
                    int slot = 2 * copy + a;
                    for (i64 k = sl_off[slot]; k < sl_off[slot + 1]; ++k)
                        L.rows[(size_t)r * n + var_pos(sl_var[k])] += gN[c][q] * sl_coef[k];
                }
            }
        L.K.assign((size_t)n * n, 0.0);
        L.g.assign(n, 0.0);
        L.c0 = 0;
        for (int r = 0; r < 4; ++r) {
            const double* row = &L.rows[(size_t)r * n];
            for (int i = 0; i < n; ++i) {
                double wi = A * row[i];
                if (wi == 0) continue;
                L.g[i] += wi * rhs[r];
                for (int j = 0; j < n; ++j) L.K[(size_t)i * n + j] += wi * row[j];
            }
            L.c0 += A * rhs[r] * rhs[r];
        }
    };
    {
        std::vector<int> nvar(nf);
        parallel_for(nf, [&](i64 b, i64 e) {
            std::vector<int> vars;
            for (i64 f = b; f < e; ++f) {
                vars.clear();
                for (int c = 0; c < 3; ++c) {
                    int copy = corner_copy[3 * f + c];
                    for (int a = 0; a < 2; ++a) {
                        int slot = 2 * copy + a;
                        for (i64 k = sl_off[slot]; k < sl_off[slot + 1]; ++k) vars.push_back(sl_var[k]);
                    }
                }
                std::sort(vars.begin(), vars.end());
                nvar[f] = (int)(std::unique(vars.begin(), vars.end()) - vars.begin());
            }
        }, 512);
        {
            int mx = 0, over = 0;
            for (int f = 0; f < nf; ++f) { mx = std::max(mx, nvar[f]); if (nvar[f] >= 40) ++over; }
            i64 maxterms = 0;
            for (int sl = 0; sl < nslots; ++sl) maxterms = std::max(maxterms, sl_off[sl + 1] - sl_off[sl]);
            ctx.logf("param: локальных переменных на грань max %d (переполнений %d), max термов в слоте %lld", mx, over, (long long)maxterms);
            double bytes = 0;
            for (int f = 0; f < nf; ++f) bytes += 8.0 * nvar[f] * (nvar[f] + 2.0) + 8.0 * nvar[f];
            if (bytes > (double)in.max_local_bytes) {
                ctx.logf("param: система слишком сложна: %d особых точек, нужно ~%.1f ГБ — прерываю", out.nsingular,
                         bytes / 1e9);
                out.too_complex = true;
                return false;
            }
        }
        for (int f = 0; f < nf; ++f) {
            fv_off[f + 1] = fv_off[f] + nvar[f];
            fk_off[f + 1] = fk_off[f] + (i64)nvar[f] * nvar[f];
        }
        fv_var.resize(fv_off[nf]);
        fk.resize(fk_off[nf]);
        fg.resize(fv_off[nf]);
        parallel_for(nf, [&](i64 b, i64 e) {
          FaceLocal L;
          for (i64 f = b; f < e; ++f) {
            build_face((int)f, L);
            for (int i = 0; i < L.nvar; ++i) {
                fv_var[fv_off[f] + i] = L.vars[i];
                fg[fv_off[f] + i] = L.g[i];
            }
            for (int i = 0; i < L.nvar * L.nvar; ++i) fk[fk_off[f] + i] = L.K[i];
            face_c0[f] = L.c0;
          }
        }, 256);
    }
    ctx.check_cancel();

    // ------------------------------------------------- сборка H (CSR) и g
    Csr H;
    H.n = nfree;
    std::vector<double> gvec(nfree, 0.0);
    std::vector<double> invdiag(nfree, 1.0);
    std::vector<i64> vf_off(nfree + 1, 0);
    std::vector<int> vf_face, vf_pos;
    {
        for (i64 k = 0; k < fv_off[nf]; ++k) ++vf_off[fv_var[k] + 1];
        for (int i = 0; i < nfree; ++i) vf_off[i + 1] += vf_off[i];
        vf_face.resize(vf_off[nfree]);
        vf_pos.resize(vf_off[nfree]);
        std::vector<i64> fill(vf_off.begin(), vf_off.end() - 1);
        for (int f = 0; f < nf; ++f)
            for (i64 k = fv_off[f]; k < fv_off[f + 1]; ++k) {
                i64 o = fill[fv_var[k]]++;
                vf_face[o] = f;
                vf_pos[o] = (int)(k - fv_off[f]);
            }
    }
    auto assemble_H = [&]() {
        std::vector<std::vector<std::pair<int, double>>> rows(nfree);
        parallel_for(nfree, [&](i64 b, i64 e) {
            std::unordered_map<int, double> mp;
            for (i64 i = b; i < e; ++i) {
                mp.clear();
                double gi = 0;
                for (i64 k = vf_off[i]; k < vf_off[i + 1]; ++k) {
                    int f = vf_face[k], pos = vf_pos[k];
                    int n = (int)(fv_off[f + 1] - fv_off[f]);
                    const double* Kf = &fk[fk_off[f]];
                    const int* vars = &fv_var[fv_off[f]];
                    gi += fg[fv_off[f] + pos];
                    for (int j = 0; j < n; ++j) mp[vars[j]] += Kf[pos * n + j];
                }
                gvec[i] = gi;
                auto& r = rows[i];
                r.reserve(mp.size());
                for (auto& kv : mp) r.push_back(kv);
                std::sort(r.begin(), r.end());
            }
        }, 128);
        H.off.assign(nfree + 1, 0);
        for (int i = 0; i < nfree; ++i) H.off[i + 1] = H.off[i] + (i64)rows[i].size();
        H.col.resize(H.off[nfree]);
        H.val.resize(H.off[nfree]);
        parallel_each(nfree, [&](i64 i) {
            i64 o = H.off[i];
            for (auto& kv : rows[i]) { H.col[o] = kv.first; H.val[o] = kv.second; ++o; }
        }, 1024);
        for (int i = 0; i < nfree; ++i) {
            double d = 0;
            for (i64 k = H.off[i]; k < H.off[i + 1]; ++k)
                if (H.col[k] == i) d = H.val[k];
            invdiag[i] = d > 1e-14 ? 1.0 / d : 1.0;
        }
    };
    assemble_H();
    ctx.logf("param: матрица %d x %d, nnz=%lld (%.2f с)", nfree, nfree, (long long)H.off[nfree],
             now_seconds() - t0);

    // ------------------------------------------------- решение + округление
    std::vector<double> x(nfree, 0.0);
    std::vector<uint8_t> fixed(nfree, 0);
    double res = 0;
    IC0 ic;
    // прямой решатель: G — целые, t-слоты и плотные строки; координаты переменных для упорядочения
    SchurSolver schur;
    bool use_direct = in.direct;
    bool schur_dirty = true;
    bool schur_refine = false;  // H изменилась локально: сначала PCG со старым разложением
    std::vector<uint8_t> in_g(nfree, 0);
    std::vector<double> var_xyz(3 * (size_t)nfree, 0.0);
    if (use_direct) {
        for (int sl = 0; sl < nslots; ++sl) {
            const int i = free_idx[sl];
            if (i < 0) continue;
            if (sl < 2 * ncopies) {
                const V3& p = M.P[copy_vertex[sl / 2]];
                var_xyz[3 * (size_t)i] = p.x;
                var_xyz[3 * (size_t)i + 1] = p.y;
                var_xyz[3 * (size_t)i + 2] = p.z;
            } else {
                in_g[i] = 1;
            }
        }
        int ng = 0;
        for (int i = 0; i < nfree; ++i) {
            if (is_int[i] || H.off[i + 1] - H.off[i] > in.dense_row) in_g[i] = 1;
            ng += in_g[i];
        }
        if (ng > in.direct_max_g) {
            ctx.logf("param: плотная часть %d > %d — используется CG", ng, in.direct_max_g);
            use_direct = false;
        }
    }
    auto solve = [&](double tol) {
        if (use_direct) {
            if (schur_refine && !schur_dirty && schur.ok) {
                schur_refine = false;
                const double tr = now_seconds();
                const int itr = schur.refine(H, gvec, fixed, x, 200, 1e-10);
                if (itr >= 0) {
                    ctx.logf("param: прямой решатель: уточнение PCG со старым разложением, %d итераций (%.2f с)", itr,
                             now_seconds() - tr);
                    res = 0;
                    return itr;
                }
                schur_dirty = true;
            }
            schur_refine = false;
            if (schur_dirty) {
                schur_dirty = false;
                if (!schur.setup(H, gvec, in_g, var_xyz, in.direct_max_fill, ctx)) use_direct = false;
            }
            if (use_direct && schur.solve(fixed, x)) {
                res = 0;
                if (in.verbose) {
                    double rn = 0, gn = 0, rl = 0, rg = 0, xm = 0;
                    for (int i = 0; i < nfree; ++i) {
                        xm = std::max(xm, std::fabs(x[i]));
                        if (fixed[i]) continue;
                        double a = 0;
                        for (i64 k = H.off[i]; k < H.off[i + 1]; ++k) a += H.val[k] * x[H.col[k]];
                        rn += (a - gvec[i]) * (a - gvec[i]);
                        (in_g[i] ? rg : rl) += (a - gvec[i]) * (a - gvec[i]);
                        gn += gvec[i] * gvec[i];
                    }
                    res = std::sqrt(rn / std::max(gn, 1e-300));
                    ctx.logf("param: прямой решатель: невязка %.2e (L %.2e, G %.2e), max|x| %.3g", res,
                             std::sqrt(rl / std::max(gn, 1e-300)), std::sqrt(rg / std::max(gn, 1e-300)), xm);
                }
                return 0;
            }
            if (use_direct) ctx.logf("param: прямой решатель: плотная часть не разложилась — CG");
            use_direct = false;
        }
        if (in.use_ic) ic.factor(H, fixed);
        const bool use = in.use_ic && ic.ok;
        return cg_solve(H, gvec, fixed, invdiag, x, in.cg_max_iters, tol, ctx, &res, use ? &ic : nullptr);
    };
    {
        int bad_g = 0, bad_h = 0;
        for (double v : gvec) if (!std::isfinite(v)) ++bad_g;
        for (double v : H.val) if (!std::isfinite(v)) ++bad_h;
        int bad_rep = 0;
        for (double v : rep) if (!std::isfinite(v)) ++bad_rep;
        if (bad_g || bad_h || bad_rep)
            ctx.logf("param: ВНИМАНИЕ нечисловые значения: g=%d H=%d rep=%d", bad_g, bad_h, bad_rep);
    }
    int it = solve(in.cg_tol);
    if (use_direct) ctx.logf("param: релаксированное решение: прямое (%.2f с)", now_seconds() - t0);
    else ctx.logf("param: релаксированное решение: CG %d итераций, невязка %.2e (%.2f с)", it, res,
                  now_seconds() - t0);
    if (in.verbose) {
        // статистика релаксированного решения
        int flipped = 0;
        for (int f = 0; f < nf; ++f) {
            double q[6];
            for (int c = 0; c < 3; ++c) {
                int copy = corner_copy[3 * f + c];
                for (int a = 0; a < 2; ++a) {
                    int slot = 2 * copy + a;
                    double sacc = 0;
                    for (i64 k = sl_off[slot]; k < sl_off[slot + 1]; ++k) sacc += sl_coef[k] * x[sl_var[k]];
                    q[2 * c + a] = sacc;
                }
            }
            double a2 = (q[2] - q[0]) * (q[5] - q[1]) - (q[4] - q[0]) * (q[3] - q[1]);
            if (a2 <= 0) ++flipped;
        }
        ctx.logf("param: релаксированное решение: перевёрнуто %d из %d", flipped, nf);
    }
    {
        std::vector<int> pending;
        for (int i = 0; i < nfree; ++i)
            if (is_int[i] && !in.skip_rounding) pending.push_back(i);
        int round = 0;
        std::vector<std::pair<double, int>> cand;
        while (!pending.empty()) {
            ctx.check_cancel();
            cand.clear();
            for (int i : pending) cand.push_back({std::fabs(x[i] - std::round(x[i])), i});
            std::sort(cand.begin(), cand.end());
            size_t take = 0;
            while (take < cand.size() && cand[take].first < in.round_tol) ++take;
            size_t minb = std::max<size_t>(1, (size_t)std::ceil(cand.size() * in.round_frac));
            if (round >= in.max_rounds - 1) take = cand.size();
            else if (take < minb) take = minb;
            for (size_t q = 0; q < take; ++q) {
                int i = cand[q].second;
                x[i] = std::round(x[i]);
                fixed[i] = 1;
            }
            std::vector<int> rest;
            for (size_t q = take; q < cand.size(); ++q) rest.push_back(cand[q].second);
            pending.swap(rest);
            it = solve(pending.empty() ? in.cg_tol : in.cg_tol_round);
            ++round;
            if (in.verbose)
                ctx.logf("param: округление, раунд %d: зафиксировано %zu, осталось %zu, CG %d (%.2e)",
                         round, take, pending.size(), it, res);
        }
        ctx.logf("param: округление завершено за %d раундов (%.2f с)", round, now_seconds() - t0);
    }

    // ------------------------------------------------- вычисление uv и статистики
    out.uv.assign(6 * (size_t)nf, 0.0);
    auto eval_uv = [&]() {
        parallel_each(nf, [&](i64 f) {
            for (int c = 0; c < 3; ++c) {
                int copy = corner_copy[3 * f + c];
                for (int a = 0; a < 2; ++a) {
                    int slot = 2 * copy + a;
                    double s = 0;
                    for (i64 k = sl_off[slot]; k < sl_off[slot + 1]; ++k) s += sl_coef[k] * x[sl_var[k]];
                    out.uv[6 * f + 2 * c + a] = s;
                }
            }
        });
        std::vector<int> bad;
        double minr = 1e300;
        for (int f = 0; f < nf; ++f) {
            const double* q = &out.uv[6 * f];
            double a2 = (q[2] - q[0]) * (q[5] - q[1]) - (q[4] - q[0]) * (q[3] - q[1]);
            double s = 1.0 / size_eff[f];
            double ref = M.FA[f] * s * s * 2.0;
            if (a2 <= 0) bad.push_back(f);
            if (ref > 0) minr = std::min(minr, a2 / ref);
        }
        out.flipped = (int)bad.size();
        out.min_area_ratio = minr;
        return bad;
    };
    std::vector<int> bad = eval_uv();
    ctx.logf("param: перевёрнутых треугольников: %d из %d, min отношение площади %.3f", out.flipped, nf,
             out.min_area_ratio);
    std::vector<double> x_best = x;
    int best_bad = out.flipped;
    // Ремонт: повышаем вес окрестности перевёрнутых граней и пересчитываем (целые зафиксированы).
    for (int rep_round = 0; rep_round < in.repair_rounds && !bad.empty(); ++rep_round) {
        ctx.check_cancel();
        std::vector<uint8_t> vmark(nv, 0), fmark(nf, 0);
        for (int f : bad)
            for (int c = 0; c < 3; ++c) vmark[M.F[3 * f + c]] = 1;
        for (int depth = 0; depth < 3; ++depth) {
            std::vector<uint8_t> nvmark = vmark;
            for (int f = 0; f < nf; ++f) {
                bool any = false;
                for (int c = 0; c < 3; ++c) any = any || vmark[M.F[3 * f + c]];
                if (!any) continue;
                fmark[f] = 1;
                for (int c = 0; c < 3; ++c) nvmark[M.F[3 * f + c]] = 1;
            }
            vmark.swap(nvmark);
        }
        const double mult = std::pow(6.0, rep_round + 1);
        std::vector<int> changed;
        for (int f = 0; f < nf; ++f)
            if (fmark[f]) { wmul[f] = mult; changed.push_back(f); }
        parallel_each((i64)changed.size(), [&](i64 k) {
            int f = changed[k];
            thread_local FaceLocal L;
            build_face(f, L);
            for (int i = 0; i < L.nvar; ++i) {
                fv_var[fv_off[f] + i] = L.vars[i];
                fg[fv_off[f] + i] = L.g[i];
            }
            for (int i = 0; i < L.nvar * L.nvar; ++i) fk[fk_off[f] + i] = L.K[i];
        }, 64);
        assemble_H();
        schur_refine = true;
        solve(in.cg_tol);
        bad = eval_uv();
        ctx.logf("param: ремонт %d (вес x%.0f на %zu граней): перевёрнуто %d, min отношение %.3f", rep_round + 1,
                 mult, changed.size(), out.flipped, out.min_area_ratio);
        if (out.flipped < best_bad) {
            best_bad = out.flipped;
            x_best = x;
        } else {
            x = x_best;
            bad = eval_uv();
            break;
        }
    }
    out.bad_faces = bad;
    ctx.logf("param: готово за %.2f с", now_seconds() - t0);
    out.ok = true;
    return true;
}

}  // namespace qf
