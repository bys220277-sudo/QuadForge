// Изотропный ремешинг (Botsch–Kobbelt) с сохранением фичерей и полем размеров.
#pragma once

#include "qf_bvh.h"
#include "qf_mesh.h"

namespace qf {

struct RemeshOptions {
    double target_len = 0.1;                         // целевая длина ребра рабочей сетки
    const std::vector<double>* src_size = nullptr;   // длина на треугольник исходника (иначе target_len)
    int iterations = 8;
    double smooth_lambda = 0.7;
    double corner_angle = 0.7;  // рад: поворот цепочки фичеров больше — «угол» (фиксируется)
    int max_faces = 6000000;
};

struct RemeshStats {
    int splits = 0, collapses = 0, flips = 0;
    int faces_in = 0, faces_out = 0;
};

// src должен быть построен (build). feat — 3*nF флагов фичер-рёбер (включая границы).
// out/out_feat — результат (топология построена внутри out.build()).
bool isotropic_remesh(const TriMesh& src, const std::vector<uint8_t>& feat, const RemeshOptions& opt,
                      TriMesh& out, std::vector<uint8_t>& out_feat, RemeshStats& stats, Ctx& ctx);

}  // namespace qf
