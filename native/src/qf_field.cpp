#include "qf_field.h"

#include <random>

#include "qf_bvh.h"

namespace qf {

void compute_field_geometry(const TriMesh& M, FieldGeometry& G) {
    const int nf = M.nF();
    G.transport.assign(3 * (size_t)nf, 0.0);
    G.weight.assign(3 * (size_t)nf, 0.0);
    // центроиды для весов
    std::vector<V3> C(nf);
    parallel_each(nf, [&](i64 f) {
        C[f] = (M.P[M.F[3 * f]] + M.P[M.F[3 * f + 1]] + M.P[M.F[3 * f + 2]]) * (1.0 / 3.0);
    });
    parallel_each(nf, [&](i64 f) {
        for (int i = 0; i < 3; ++i) {
            int g = M.FF[3 * f + i];
            if (g < 0) { G.transport[3 * f + i] = 0; G.weight[3 * f + i] = 0; continue; }
            int j = M.FE[3 * f + i];
            V3 a = M.P[M.F[3 * f + i]], b = M.P[M.F[3 * f + (i + 1) % 3]];
            V3 e = b - a;
            double phi_f = angle_in_face(M, (int)f, e);
            double phi_g = angle_in_face(M, g, e);
            (void)j;
            G.transport[3 * f + i] = phi_f - phi_g;  // g -> f
            double len = norm(e);
            double dist = norm(C[f] - C[g]);
            G.weight[3 * f + i] = dist > 1e-300 ? len / dist : 1.0;
        }
    });
}

void field_angles_from_z(const std::vector<C2>& z, std::vector<double>& theta) {
    theta.resize(z.size());
    for (size_t f = 0; f < z.size(); ++f) theta[f] = 0.25 * z[f].arg();
}

namespace {

struct Op {
    const TriMesh* M;
    const FieldGeometry* G;
    std::vector<C2> c;       // 3*nF: коэффициент e^{i4 r_gf}
    std::vector<double> diag;  // диагональ (для предобусловливателя)
    std::vector<uint8_t> fixed;
    std::vector<double> lam;  // мягкий вес + eps*A на грань

    void apply(const std::vector<C2>& x, std::vector<C2>& y) const {
        const int nf = M->nF();
        parallel_each(nf, [&](i64 f) {
            if (fixed[f]) { y[f] = x[f]; return; }
            C2 s = x[f] * lam[f];
            for (int i = 0; i < 3; ++i) {
                int g = M->FF[3 * f + i];
                if (g < 0) continue;
                double w = G->weight[3 * f + i];
                if (fixed[g]) {
                    s += x[f] * w;  // вклад g идёт в правую часть
                } else {
                    s += (x[f] - c[3 * f + i] * x[g]) * w;
                }
            }
            y[f] = s;
        });
    }
};

double cdot_re(const std::vector<C2>& a, const std::vector<C2>& b) {
    return parallel_sum((i64)a.size(), [&](i64 s, i64 e) {
        double r = 0;
        for (i64 i = s; i < e; ++i) r += a[i].re * b[i].re + a[i].im * b[i].im;
        return r;
    });
}

}  // namespace

void compute_cross_field(const TriMesh& M, const FieldGeometry& G, const FieldInput& in,
                         std::vector<C2>& z, Ctx& ctx) {
    const int nf = M.nF();
    Op op;
    op.M = &M;
    op.G = &G;
    op.c.resize(3 * (size_t)nf);
    for (size_t k = 0; k < op.c.size(); ++k) op.c[k] = polar(4.0 * G.transport[k]);
    op.fixed.assign(nf, 0);
    bool any_fixed = !in.fixed.empty();
    for (int f = 0; f < nf; ++f)
        if (any_fixed && in.fixed[f].abs2() > 0.5) op.fixed[f] = 1;

    // средняя диагональ
    double avg = 0;
    int cnt = 0;
    for (int f = 0; f < nf; ++f) {
        double s = 0;
        for (int i = 0; i < 3; ++i) s += G.weight[3 * f + i];
        avg += s;
        ++cnt;
    }
    avg = cnt ? avg / cnt : 3.0;
    const double eps = in.eps * avg;

    op.lam.assign(nf, eps);
    const bool has_soft = !in.soft_w.empty();
    for (int f = 0; f < nf; ++f)
        if (has_soft) op.lam[f] += in.soft_w[f];
    op.diag.assign(nf, 1.0);
    for (int f = 0; f < nf; ++f) {
        double s = op.lam[f];
        for (int i = 0; i < 3; ++i)
            if (M.FF[3 * f + i] >= 0) s += G.weight[3 * f + i];
        op.diag[f] = op.fixed[f] ? 1.0 : s;
    }

    // начальное приближение
    z.assign(nf, C2(1, 0));
    {
        std::mt19937 rng(in.seed);
        std::uniform_real_distribution<double> U(0.0, 2.0 * kPi);
        for (int f = 0; f < nf; ++f) {
            if (op.fixed[f]) z[f] = in.fixed[f];
            else if (has_soft && in.soft_w[f] > 0.05) z[f] = in.soft_z[f];
            else z[f] = polar(U(rng));
        }
    }

    std::vector<C2> b(nf), r(nf), p(nf), Ap(nf), x(nf), zp(nf);
    for (int outer = 0; outer < in.outer_iters; ++outer) {
        ctx.check_cancel();
        // нормализуем текущее приближение (proximal point)
        for (int f = 0; f < nf; ++f) {
            double a = z[f].abs();
            zp[f] = a > 1e-12 ? z[f] * (1.0 / a) : C2(1, 0);
            if (op.fixed[f]) zp[f] = in.fixed[f];
        }
        // правая часть
        parallel_each(nf, [&](i64 f) {
            if (op.fixed[f]) { b[f] = in.fixed[f]; return; }
            C2 s = zp[f] * eps;
            if (has_soft) s += in.soft_z[f] * in.soft_w[f];
            for (int i = 0; i < 3; ++i) {
                int g = M.FF[3 * f + i];
                if (g < 0) continue;
                if (op.fixed[g]) s += op.c[3 * f + i] * in.fixed[g] * G.weight[3 * f + i];
            }
            b[f] = s;
        });
        // CG с предобусловливателем Якоби; x0 = zp
        x = zp;
        op.apply(x, Ap);
        for (int f = 0; f < nf; ++f) r[f] = b[f] - Ap[f];
        std::vector<C2> zres(nf);
        for (int f = 0; f < nf; ++f) zres[f] = r[f] * (1.0 / op.diag[f]);
        p = zres;
        double rz = cdot_re(r, zres);
        double bnorm = std::sqrt(cdot_re(b, b)) + 1e-30;
        int it = 0;
        for (; it < in.cg_max_iters; ++it) {
            if (std::sqrt(cdot_re(r, r)) < in.cg_tol * bnorm) break;
            op.apply(p, Ap);
            double pAp = cdot_re(p, Ap);
            if (pAp <= 0) break;
            double alpha = rz / pAp;
            parallel_each(nf, [&](i64 f) {
                x[f] += p[f] * alpha;
                r[f] -= Ap[f] * alpha;
            });
            for (int f = 0; f < nf; ++f) zres[f] = r[f] * (1.0 / op.diag[f]);
            double rz2 = cdot_re(r, zres);
            double beta = rz2 / rz;
            rz = rz2;
            parallel_each(nf, [&](i64 f) { p[f] = zres[f] + p[f] * beta; });
            if ((it & 31) == 31) ctx.check_cancel();
        }
        z = x;
        if (ctx.verbose) ctx.logf("поле: внешняя итерация %d, CG %d итераций", outer, it);
    }
    // финальная нормализация
    for (int f = 0; f < nf; ++f) {
        if (op.fixed[f]) { z[f] = in.fixed[f]; continue; }
        double a = z[f].abs();
        z[f] = a > 1e-12 ? z[f] * (1.0 / a) : C2(1, 0);
    }
}


void polish_field(const TriMesh& M, const FieldGeometry& G, const FieldInput& in,
                  std::vector<C2>& z, int iterations, Ctx& ctx, double noise) {
    const int nf = M.nF();
    std::vector<C2> c(3 * (size_t)nf);
    for (size_t k = 0; k < c.size(); ++k) c[k] = polar(4.0 * G.transport[k]);
    std::vector<uint8_t> fixed(nf, 0);
    const bool any_fixed = !in.fixed.empty();
    for (int f = 0; f < nf; ++f)
        if (any_fixed && in.fixed[f].abs2() > 0.5) fixed[f] = 1;
    const bool has_soft = !in.soft_w.empty();
    // небольшое возмущение для нарушения симметрии
    std::mt19937 rng(in.seed * 2654435761u + 17u);
    std::normal_distribution<double> N(0.0, std::max(noise, 1e-12));
    for (int f = 0; f < nf && noise > 0; ++f) {
        if (fixed[f]) continue;
        z[f] = z[f] + C2(N(rng), N(rng));
        double a = z[f].abs();
        z[f] = a > 1e-12 ? z[f] * (1.0 / a) : C2(1, 0);
    }
    std::vector<C2> zn(nf);
    for (int it = 0; it < iterations; ++it) {
        parallel_each(nf, [&](i64 f) {
            if (fixed[f]) { zn[f] = in.fixed[f]; return; }
            C2 s(0, 0);
            double wsum = 0;
            for (int i = 0; i < 3; ++i) {
                int g = M.FF[3 * f + i];
                if (g < 0) continue;
                double w = G.weight[3 * f + i];
                const C2& zg = fixed[g] ? in.fixed[g] : z[g];
                s += c[3 * f + i] * zg * w;
                wsum += w;
            }
            if (has_soft) s += in.soft_z[f] * in.soft_w[f];
            s += z[f] * (0.5 * wsum);  // демпфирование
            double a = s.abs();
            zn[f] = a > 1e-12 ? s * (1.0 / a) : z[f];
        });
        z.swap(zn);
        if ((it & 15) == 15) ctx.check_cancel();
    }
}

void compute_field_index(const TriMesh& M, const std::vector<C2>& z, std::vector<int>& index) {
    const int nv = M.nV();
    index.assign(nv, 0);
    parallel_each(nv, [&](i64 v) {
        if (M.vbound[v]) return;
        const int b = M.vfan_off[v], e = M.vfan_off[v + 1];
        const int m = e - b;
        if (m < 3) return;
        auto wrap = [](double x) {  // в [-pi/4, pi/4)
            double q = x / kHalfPi;
            return (q - std::floor(q + 0.5)) * kHalfPi;
        };
        double Lam = 0, S = 0, psi_lift = 0, psi0 = 0;
        for (int t = 0; t <= m; ++t) {
            int tt = t % m;
            int f = M.vfan_face[b + tt], c = M.vfan_corner[b + tt];
            const V3& pv = M.P[M.F[3 * f + c]];
            V3 o = M.P[M.F[3 * f + (c + 1) % 3]] - pv;
            V3 q = M.P[M.F[3 * f + (c + 2) % 3]] - pv;
            double lo = norm(o), lq = norm(q);
            double alpha = std::acos(std::max(-1.0, std::min(1.0, dot(o, q) / std::max(lo * lq, 1e-300))));
            if (t >= 1) Lam -= alpha;
            if (t < m) S += alpha;
            double rho = 0.25 * z[f].arg() - angle_in_face(M, f, o);
            double psi = Lam + rho;
            if (t == 0) { psi0 = psi; psi_lift = psi; }
            else psi_lift = psi_lift + wrap(psi - psi_lift);
            if (t == m) {
                double n = std::round((psi_lift - psi) / kHalfPi);
                int nat = (int)std::lround(S / kHalfPi);
                index[v] = (int)n - nat;
            }
        }
        (void)psi0;
    }, 256);
}


void prolong_field(const TriMesh& C, const std::vector<C2>& zc, const TriMesh& W, std::vector<C2>& zw) {
    TriBvh bvh;
    bvh.build(C.P, C.F);
    const int nf = W.nF();
    zw.assign(nf, C2(1, 0));
    parallel_each(nf, [&](i64 f) {
        V3 c = (W.P[W.F[3 * f]] + W.P[W.F[3 * f + 1]] + W.P[W.F[3 * f + 2]]) * (1.0 / 3.0);
        V3 q;
        double d2;
        const V3& nw = W.FN[f];
        int fc = bvh.nearest_if(c, q, d2, [&](int g) { return dot(C.FN[g], nw) > 0.7; });
        if (fc < 0) fc = bvh.nearest(c, q, d2);
        if (fc < 0) return;
        V3 e1, e2;
        C.face_frame(fc, e1, e2);
        double th = 0.25 * zc[fc].arg();
        V3 d = e1 * std::cos(th) + e2 * std::sin(th);
        // проекция на касательную плоскость грани W
        const V3& n = W.FN[f];
        V3 dp = d - n * dot(n, d);
        if (norm(dp) < 1e-9) { zw[f] = C2(1, 0); return; }
        zw[f] = polar(4.0 * angle_in_face(W, (int)f, dp));
    });
}

}  // namespace qf
