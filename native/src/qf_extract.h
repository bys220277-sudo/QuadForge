// Извлечение квад-сетки из целочисленной параметризации (arrangement iso-линий).
#pragma once

#include "qf_mesh.h"

namespace qf {

struct ExtractResult {
    std::vector<V3> V;                       // вершины решётки
    std::vector<std::array<int, 4>> quads;   // CCW относительно нормалей поверхности
    std::vector<std::vector<int>> others;    // не-квадовые грани (дефекты)
    std::vector<uint8_t> vfeature;           // вершина лежит на фичере/границе (0/1)
    int dead_ends = 0;                       // узлы-тупики (дефекты)
    int boundary_loops = 0;                  // отброшенные внешние циклы
    int nonquad = 0;
    int skipped_triangles = 0;
    bool ok = false;
};

// uv: 6*nF (углы граней в собственной карте). feature: 3*nF флаги рёбер-фичерей (может быть nullptr).
bool extract_quads(const TriMesh& M, const std::vector<double>& uv, const std::vector<uint8_t>* feature,
                   ExtractResult& out, Ctx& ctx);

}  // namespace qf
