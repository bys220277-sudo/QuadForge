// Поле направлений 4-RoSy (кросс-поле) на треугольной сетке.
#pragma once

#include "qf_mesh.h"

namespace qf {

struct FieldInput {
    // Жёсткие ограничения: |z|=1 — направление зафиксировано (z = e^{i4θ}); 0 — свободно.
    std::vector<C2> fixed;
    // Мягкий приор (кривизна/guides): z_soft и вес (безразмерный, ~ «число рёбер»).
    std::vector<C2> soft_z;
    std::vector<double> soft_w;
    int outer_iters = 6;
    double eps = 1e-3;     // проксимальный член (доля средней диагонали)
    int cg_max_iters = 600;
    double cg_tol = 1e-5;
    unsigned seed = 1;
};

// Угол переноса (из грани g в грань f) через ребро i грани f.
// θ_f ≈ θ_g + transport[3f+i] (mod 2π/4 для 4-RoSy).
struct FieldGeometry {
    std::vector<double> transport;  // 3*nF, NaN на границе
    std::vector<double> weight;     // 3*nF вес ребра
};
void compute_field_geometry(const TriMesh& M, FieldGeometry& G);

// Возвращает z (единичные комплексные) на грань. Ctx — для логов.
void compute_cross_field(const TriMesh& M, const FieldGeometry& G, const FieldInput& in,
                         std::vector<C2>& z, Ctx& ctx);

// Нелинейная релаксация (|z|=1) — расщепляет сингулярности высокого индекса.
void polish_field(const TriMesh& M, const FieldGeometry& G, const FieldInput& in,
                  std::vector<C2>& z, int iterations, Ctx& ctx, double noise = 0.08);

// Истинный индекс сингулярности в четвертях оборота на внутренние вершины
// (0 — регулярная; ±1 — валентность 3/5; ±2, ±4 — нежелательные высокие индексы).
void compute_field_index(const TriMesh& M, const std::vector<C2>& z, std::vector<int>& index);

// Перенос поля с грубой сетки C на мелкую W: направление ближайшей грани C проецируется на грань W.
void prolong_field(const TriMesh& C, const std::vector<C2>& zc, const TriMesh& W, std::vector<C2>& zw);

// Из |z| и сингулярностей: индекс сингулярности (в четвертях) для вершин: k mod 4.
// theta[f] — представитель (из z) в базисе грани.
void field_angles_from_z(const std::vector<C2>& z, std::vector<double>& theta);

}  // namespace qf
