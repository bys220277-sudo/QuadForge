// Целочисленная (seamless) параметризация по полю направлений: MIQ.
#pragma once

#include "qf_field.h"

namespace qf {

struct ParamInput {
    const TriMesh* M = nullptr;
    const FieldGeometry* G = nullptr;
    const std::vector<C2>* z = nullptr;          // кросс-поле (единичные z на грань)
    const std::vector<double>* size = nullptr;   // целевая длина ребра квада на грань
    // Рёбра, которые должны лечь на линии решётки (фичери/границы/guides-hard).
    // 3*nF, симметрично для обеих граней ребра.
    const std::vector<uint8_t>* feature = nullptr;
    int cg_max_iters = 3000;
    bool use_ic = true;
    bool skip_rounding = false;  // отладка: только релаксированное решение
    // Конформная поправка размеров: поле, поворачивающееся на кривой поверхности, требует
    // изменения масштаба (гармонически сопряжённого к углу поля). 0 — выкл., 1 — полностью.
    // Интегрируемость: из поворота поля вычитается его «точная» часть (поворот граней на α,
    // α — решение уравнения Пуассона на дуальном графе). 0 — выкл., 1 — полностью.
    double integrable = 0.0;
    const std::vector<uint8_t>* fixed_faces = nullptr;  // грани с жёстко заданным направлением
    double conformal = 0.0;
    double conformal_clamp = 1.2;  // |log масштаба| не больше этого
    double cg_tol = 1e-6;
    int max_rounds = 40;
    int repair_rounds = 3;
    double round_tol = 0.18;   // сразу фиксируем переменные ближе этого к целому
    double round_frac = 0.3;   // и не менее этой доли оставшихся за раунд
    double cg_tol_round = 1e-4;   // сразу фиксируем переменные ближе этого к целому
    bool verbose = false;
    // Прямой решатель (разреженный Холецкий + дополнение Шура по целым); CG — запасной путь.
    bool direct = true;
    int direct_max_g = 4000;        // предел размера плотного дополнения Шура
    double sliver_clamp = 0.02;     // нижняя граница 2*площади/Σe² для оператора градиента
    int dense_row = 512;            // строки H длиннее — в плотную часть
    i64 direct_max_fill = 80000000; // предел заполнения разреженного фактора
    // Предел памяти локальных систем граней. Тысячи особых точек поля (шумные фичери, мелкие части)
    // раздувают выражения переменных через цепочки разреза — без предела это терабайты и зависание.
    i64 max_local_bytes = 1500000000;
};

struct ParamResult {
    std::vector<double> uv;          // 6*nF: (u,v) трёх углов граней
    std::vector<double> rep;         // представитель угла поля (после комбинга)
    std::vector<uint8_t> kjump;      // 3*nF период. скачок
    std::vector<double> phi;         // конформный лог-масштаб на грань (если считался)
    std::vector<int> foreign;        // отладка
    std::vector<uint8_t> cut;        // 3*nF рёбра разреза
    std::vector<int8_t> vindex;      // индекс сингулярности (0..3) на вершину, -1 для границы
    int nsingular = 0;
    int ncut_edges = 0;
    int n_free = 0, n_int = 0;
    int flipped = 0;
    std::vector<int> bad_faces;      // грани с неположительной площадью в uv
    double min_area_ratio = 0;
    double energy = 0;
    bool ok = false;
    bool too_complex = false;        // система слишком велика (см. max_local_bytes) — нужно поле проще
    std::string error;
};

bool compute_param(const ParamInput& in, ParamResult& out, Ctx& ctx);

}  // namespace qf
