// Постобработка квад-сетки: ремонт дефектов (заполнение дыр квадами), релаксация, проекция.
#pragma once

#include "qf_bvh.h"
#include "qf_mesh.h"

namespace qf {

struct PolyMesh {
    std::vector<V3> V;
    std::vector<std::vector<int>> faces;
    std::vector<uint8_t> vfeat;  // вершина на фичере/границе (из извлечения)
};

struct RepairStats {
    int removed_faces = 0;
    int holes_filled = 0;
    int added_quads = 0;
    int added_tris = 0;
    int holes_left = 0;
};

// Удаляет дефектные грани (не-квады, неманифольдные, с противоречивой ориентацией)
// и заполняет образовавшиеся дыры квадами (шаблоны «веер»/рекурсивное деление).
void repair_quad_mesh(PolyMesh& pm, RepairStats& st, Ctx& ctx);
// Удаляет «дублеты»: внутренние вершины валентности 2 между двумя квадами (сливает квады).
int remove_doublets(PolyMesh& pm);

struct RelaxOptions {
    int iterations = 4;
    double lambda = 0.5;
    bool project = true;
};

// Тангенциальная релаксация с проекцией на исходную поверхность и фичер-кривые.
void relax_and_project(PolyMesh& pm, const TriMesh& src, const std::vector<uint8_t>& src_feat, double L,
                       const RelaxOptions& opt, Ctx& ctx);

// Удаляет неиспользуемые вершины.
void compact_mesh(PolyMesh& pm);

struct QuadQuality {
    int faces = 0, quads = 0, tris = 0, ngons = 0;
    int verts = 0, irregular = 0, boundary_edges = 0, nonmanifold_edges = 0;
    double min_angle = 0, max_angle = 0, angle_dev_mean = 0;  // градусы
    double edge_cv = 0;  // коэфф. вариации длин рёбер
};
QuadQuality evaluate_quad_mesh(const PolyMesh& pm);

}  // namespace qf
