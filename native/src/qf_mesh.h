// Треугольная сетка с топологией (смежность граней, веера вершин) и геометрия.
#pragma once

#include "qf_common.h"

namespace qf {

struct TriMesh {
    std::vector<V3> P;
    std::vector<int> F;   // 3*nF индексов вершин, ориентация согласована
    // Смежность: ребро i грани f — (F[3f+i], F[3f+(i+1)%3]).
    std::vector<int> FF;  // сосед через ребро (или -1)
    std::vector<int> FE;  // индекс совпадающего ребра в соседе
    std::vector<V3> FN;   // нормали
    std::vector<double> FA;  // площади
    // Веера вершин в порядке обхода: из грани f пересекаем ребро i (от v) в соседа.
    std::vector<int> vfan_off;        // nV+1
    std::vector<int> vfan_face;       // грани в порядке обхода
    std::vector<uint8_t> vfan_corner; // индекс угла вершины в этой грани
    std::vector<uint8_t> vbound;      // 1 если вершина на границе (веер разомкнут)

    int nV() const { return (int)P.size(); }
    int nF() const { return (int)(F.size() / 3); }
    int vert(int f, int i) const { return F[3 * f + i]; }

    // Строит смежность, нормали, веера. Возвращает false и сообщение при
    // неманифолдности. Согласует ориентацию граней (BFS), если нужно.
    bool build(std::string* err = nullptr);
    void compute_face_geometry();
    bool build_adjacency(std::string* err);
    bool build_fans(std::string* err);
    bool orient_consistently(std::string* err);

    double total_area() const;
    // Цикл граней вокруг вершины: callback(face, cornerIndex).
    template <class Fn>
    void for_fan(int v, Fn&& fn) const {
        for (int k = vfan_off[v]; k < vfan_off[v + 1]; ++k) fn(vfan_face[k], (int)vfan_corner[k]);
    }
    // Базис грани (e1 вдоль первого ребра, e2 = n × e1).
    void face_frame(int f, V3& e1, V3& e2) const {
        e1 = normalized(P[F[3 * f + 1]] - P[F[3 * f]]);
        e2 = cross(FN[f], e1);
    }
};

// Очистка «супа»: удаляет вырожденные/дублирующие грани, неиспользуемые вершины,
// чинит неманифольдные рёбра/вершины. Возвращает отображение старых граней.
struct CleanResult {
    int removed_degenerate = 0;
    int removed_duplicate = 0;
    int removed_nonmanifold = 0;
    int split_vertices = 0;
    int welded_vertices = 0;   // слито совпадающих вершин
    int tjunction_splits = 0;  // разбито граней на T-стыках
    std::vector<int> vert_origin;  // новая вершина -> индекс в исходных
    std::vector<int> weld_rep;     // исходная вершина -> исходный индекс её представителя после сварки
};
CleanResult clean_triangle_soup(std::vector<V3>& P, std::vector<int>& F);

// Кривизна: на грань — главные кривизны и направление главной оси (k1 по d1).
struct CurvatureData {
    std::vector<double> k1, k2;   // k1 >= k2 по модулю? нет: k1 >= k2 со знаком
    std::vector<double> theta;    // угол направления k1 в базисе грани (face_frame)
};
void compute_face_curvature(const TriMesh& M, CurvatureData& out, int smooth_iters = 3);

// Угол вектора w в базисе грани.
inline double angle_in_face(const TriMesh& M, int f, const V3& w) {
    V3 e1, e2;
    M.face_frame(f, e1, e2);
    return std::atan2(dot(w, e2), dot(w, e1));
}

// Нормаль вершин (по углам).
void compute_vertex_normals(const TriMesh& M, std::vector<V3>& VN);

}  // namespace qf
