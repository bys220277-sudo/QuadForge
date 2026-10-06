#include "qf_bvh.h"

namespace qf {

namespace {

struct BuildCtx {
    const std::vector<double>* boxes;
    std::vector<double> centers;  // 3 на примитив
    std::vector<int>* order;
    std::vector<BvhNode>* nodes;
    int leaf;
};

int build_rec(BuildCtx& c, int start, int count) {
    BvhNode node;
    for (int k = 0; k < 3; ++k) { node.bmin[k] = 1e300; node.bmax[k] = -1e300; }
    for (int i = start; i < start + count; ++i) {
        int p = (*c.order)[i];
        for (int k = 0; k < 3; ++k) {
            node.bmin[k] = std::min(node.bmin[k], (*c.boxes)[6 * p + k]);
            node.bmax[k] = std::max(node.bmax[k], (*c.boxes)[6 * p + 3 + k]);
        }
    }
    int id = (int)c.nodes->size();
    c.nodes->push_back(node);
    if (count <= c.leaf) {
        (*c.nodes)[id].left = (*c.nodes)[id].right = -1;
        (*c.nodes)[id].start = start;
        (*c.nodes)[id].count = count;
        return id;
    }
    // ось по максимальной протяжённости центров
    double cmin[3] = {1e300, 1e300, 1e300}, cmax[3] = {-1e300, -1e300, -1e300};
    for (int i = start; i < start + count; ++i) {
        int p = (*c.order)[i];
        for (int k = 0; k < 3; ++k) {
            cmin[k] = std::min(cmin[k], c.centers[3 * p + k]);
            cmax[k] = std::max(cmax[k], c.centers[3 * p + k]);
        }
    }
    int axis = 0;
    for (int k = 1; k < 3; ++k)
        if (cmax[k] - cmin[k] > cmax[axis] - cmin[axis]) axis = k;
    int mid = start + count / 2;
    std::nth_element(c.order->begin() + start, c.order->begin() + mid, c.order->begin() + start + count,
                     [&](int a, int b) { return c.centers[3 * a + axis] < c.centers[3 * b + axis]; });
    int l = build_rec(c, start, mid - start);
    int r = build_rec(c, mid, start + count - mid);
    (*c.nodes)[id].left = l;
    (*c.nodes)[id].right = r;
    (*c.nodes)[id].start = 0;
    (*c.nodes)[id].count = 0;
    return id;
}

inline double box_dist2(const BvhNode& n, const V3& p) {
    double d = 0;
    for (int k = 0; k < 3; ++k) {
        double v = p[k];
        if (v < n.bmin[k]) v = n.bmin[k] - v;
        else if (v > n.bmax[k]) v = v - n.bmax[k];
        else v = 0;
        d += v * v;
    }
    return d;
}

}  // namespace

void Bvh::build(const std::vector<double>& boxes, int leaf_size) {
    nodes.clear();
    order.clear();
    int n = (int)(boxes.size() / 6);
    if (n == 0) return;
    order.resize(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    BuildCtx c;
    c.boxes = &boxes;
    c.order = &order;
    c.nodes = &nodes;
    c.leaf = leaf_size;
    c.centers.resize(3 * (size_t)n);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) c.centers[3 * i + k] = 0.5 * (boxes[6 * i + k] + boxes[6 * i + 3 + k]);
    nodes.reserve(2 * n / std::max(1, leaf_size) + 16);
    build_rec(c, 0, n);
}

void TriBvh::build(const std::vector<V3>& P_, const std::vector<int>& F_) {
    P = &P_;
    F = &F_;
    int nf = (int)(F_.size() / 3);
    std::vector<double> boxes(6 * (size_t)nf);
    parallel_each(nf, [&](i64 f) {
        double mn[3] = {1e300, 1e300, 1e300}, mx[3] = {-1e300, -1e300, -1e300};
        for (int c = 0; c < 3; ++c) {
            const V3& p = P_[F_[3 * f + c]];
            for (int k = 0; k < 3; ++k) { mn[k] = std::min(mn[k], p[k]); mx[k] = std::max(mx[k], p[k]); }
        }
        for (int k = 0; k < 3; ++k) { boxes[6 * f + k] = mn[k]; boxes[6 * f + 3 + k] = mx[k]; }
    });
    bvh.build(boxes, 4);
}

int TriBvh::nearest(const V3& p, V3& out, double& d2, double* bary) const {
    if (bvh.empty()) return -1;
    int best = -1;
    double bd = 1e300;
    double bb[3] = {0, 0, 0};
    int stack[96];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        const BvhNode& n = bvh.nodes[ni];
        if (box_dist2(n, p) >= bd) continue;
        if (n.left < 0) {
            for (int i = n.start; i < n.start + n.count; ++i) {
                int f = bvh.order[i];
                double u, v, w;
                V3 q = closest_point_triangle(p, (*P)[(*F)[3 * f]], (*P)[(*F)[3 * f + 1]], (*P)[(*F)[3 * f + 2]], u, v, w);
                double d = norm2(q - p);
                if (d < bd) { bd = d; best = f; out = q; bb[0] = u; bb[1] = v; bb[2] = w; }
            }
        } else {
            double dl = box_dist2(bvh.nodes[n.left], p), dr = box_dist2(bvh.nodes[n.right], p);
            if (dl < dr) { stack[sp++] = n.right; stack[sp++] = n.left; }
            else { stack[sp++] = n.left; stack[sp++] = n.right; }
        }
    }
    d2 = bd;
    if (bary) { bary[0] = bb[0]; bary[1] = bb[1]; bary[2] = bb[2]; }
    return best;
}

int TriBvh::nearest_if(const V3& p, V3& out, double& d2, const std::function<bool(int)>& accept) const {
    if (bvh.empty()) return -1;
    int best = -1;
    double bd = 1e300;
    int stack[96];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        const BvhNode& n = bvh.nodes[ni];
        if (box_dist2(n, p) >= bd) continue;
        if (n.left < 0) {
            for (int i = n.start; i < n.start + n.count; ++i) {
                int f = bvh.order[i];
                if (!accept(f)) continue;
                double u, v, w;
                V3 q = closest_point_triangle(p, (*P)[(*F)[3 * f]], (*P)[(*F)[3 * f + 1]], (*P)[(*F)[3 * f + 2]], u, v, w);
                double d = norm2(q - p);
                if (d < bd) { bd = d; best = f; out = q; }
            }
        } else {
            double dl = box_dist2(bvh.nodes[n.left], p), dr = box_dist2(bvh.nodes[n.right], p);
            if (dl < dr) { stack[sp++] = n.right; stack[sp++] = n.left; }
            else { stack[sp++] = n.left; stack[sp++] = n.right; }
        }
    }
    d2 = bd;
    return best;
}

void SegBvh::build(const std::vector<V3>& a, const std::vector<V3>& b) {
    A = a;
    B = b;
    int n = (int)A.size();
    std::vector<double> boxes(6 * (size_t)n);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            boxes[6 * i + k] = std::min(A[i][k], B[i][k]);
            boxes[6 * i + 3 + k] = std::max(A[i][k], B[i][k]);
        }
    bvh.build(boxes, 4);
}

int SegBvh::nearest(const V3& p, V3& out, double& d2) const {
    if (bvh.empty()) return -1;
    int best = -1;
    double bd = 1e300;
    int stack[96];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        const BvhNode& n = bvh.nodes[ni];
        if (box_dist2(n, p) >= bd) continue;
        if (n.left < 0) {
            for (int i = n.start; i < n.start + n.count; ++i) {
                int s = bvh.order[i];
                double t;
                V3 q = closest_point_segment(p, A[s], B[s], t);
                double d = norm2(q - p);
                if (d < bd) { bd = d; best = s; out = q; }
            }
        } else {
            double dl = box_dist2(bvh.nodes[n.left], p), dr = box_dist2(bvh.nodes[n.right], p);
            if (dl < dr) { stack[sp++] = n.right; stack[sp++] = n.left; }
            else { stack[sp++] = n.left; stack[sp++] = n.right; }
        }
    }
    d2 = bd;
    return best;
}

}  // namespace qf
