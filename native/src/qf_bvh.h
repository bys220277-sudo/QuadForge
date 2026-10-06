// BVH для треугольников и отрезков: ближайшая точка (CPU, double).
#pragma once

#include "qf_common.h"

namespace qf {

// Ближайшая точка на треугольнике (Ericson). Возвращает барицентрики (u,v,w) через ref-параметры.
inline V3 closest_point_triangle(const V3& p, const V3& a, const V3& b, const V3& c, double& bu, double& bv,
                                 double& bw) {
    V3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { bu = 1; bv = 0; bw = 0; return a; }
    V3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { bu = 0; bv = 1; bw = 0; return b; }
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        bu = 1 - v; bv = v; bw = 0;
        return a + ab * v;
    }
    V3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { bu = 0; bv = 0; bw = 1; return c; }
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        bu = 1 - w; bv = 0; bw = w;
        return a + ac * w;
    }
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        bu = 0; bv = 1 - w; bw = w;
        return b + (c - b) * w;
    }
    double denom = va + vb + vc;
    if (std::fabs(denom) < 1e-300) { bu = 1; bv = 0; bw = 0; return a; }
    double inv = 1.0 / denom;
    double v = vb * inv, w = vc * inv;
    bu = 1 - v - w; bv = v; bw = w;
    return a + ab * v + ac * w;
}

inline V3 closest_point_segment(const V3& p, const V3& a, const V3& b, double& t) {
    V3 ab = b - a;
    double l2 = dot(ab, ab);
    t = l2 > 1e-300 ? dot(p - a, ab) / l2 : 0.0;
    t = std::max(0.0, std::min(1.0, t));
    return a + ab * t;
}

struct BvhNode {
    double bmin[3], bmax[3];
    int left, right;   // -1 для листа
    int start, count;  // диапазон в order для листа
};

// Универсальное BVH над примитивами, заданными боксами.
class Bvh {
public:
    // boxes: 6 double на примитив (min xyz, max xyz)
    void build(const std::vector<double>& boxes, int leaf_size = 4);
    std::vector<BvhNode> nodes;
    std::vector<int> order;
    bool empty() const { return nodes.empty(); }
};

// Облако треугольников + запросы.
struct TriBvh {
    const std::vector<V3>* P = nullptr;
    const std::vector<int>* F = nullptr;  // 3*nF
    Bvh bvh;
    void build(const std::vector<V3>& P_, const std::vector<int>& F_);
    // Ближайшая точка; возвращает индекс треугольника и барицентрики.
    int nearest(const V3& p, V3& out, double& d2, double* bary = nullptr) const;
    // Ближайшая среди треугольников, для которых accept(f) истинно.
    int nearest_if(const V3& p, V3& out, double& d2, const std::function<bool(int)>& accept) const;
};

struct SegBvh {
    std::vector<V3> A, B;
    Bvh bvh;
    void build(const std::vector<V3>& a, const std::vector<V3>& b);
    int nearest(const V3& p, V3& out, double& d2) const;
};

}  // namespace qf
