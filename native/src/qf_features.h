// Фичери: острые рёбра, границы, пользовательские пометки.
#pragma once

#include "qf_mesh.h"

#include <unordered_set>

namespace qf {

struct FeatureOptions {
    double sharp_angle = 0.6;   // рад; 0 = не искать по углу
    bool detect_sharp = true;
    bool boundary = true;
    const int* user_edges = nullptr;  // пары вершин (a,b) ИСХОДНОЙ нумерации, n_user пар
    int n_user = 0;
    const std::vector<int>* origin = nullptr;  // вершина сетки -> исходный индекс
};

// Возвращает 3*nF флагов (симметрично для обеих граней ребра).
void detect_features(const TriMesh& M, const FeatureOptions& opt, std::vector<uint8_t>& feat);

// Пучки близких параллельных фичер-линий (фаски, скругления, тонкие стенки уже fac·L) дают
// несовместные целочисленные ограничения. Из каждого пучка оставляем одну (самую длинную) цепочку,
// остальные снимаем. Граничные рёбра не трогаем. Возвращает число снятых рёбер.
int drop_close_parallel_features(const TriMesh& M, std::vector<uint8_t>& feat, double L, double fac = 0.5);

// Короткие оторванные цепочки фичеров (оба конца свободны, длина < min_len·L) — шум скульпта/скана:
// каждая такая цепочка рождает особые точки поля. Граничные рёбра не трогаем. Возвращает число снятых рёбер.
int drop_short_feature_chains(const TriMesh& M, std::vector<uint8_t>& feat, double L, double min_len = 1.0);

// Жёсткие ограничения поля по фичерям: z=e^{i4α} для граней, смежных с фичер-ребром.
void field_constraints_from_features(const TriMesh& M, const std::vector<uint8_t>& feat,
                                     std::vector<C2>& fixed);

}  // namespace qf
