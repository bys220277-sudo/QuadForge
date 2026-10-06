// C ABI ядра QuadForge (ctypes). Все структуры — POD, раскладка зафиксирована.
#pragma once

#include <stdint.h>

#if defined(_WIN32)
#define QF_API extern "C" __declspec(dllexport)
#else
#define QF_API extern "C" __attribute__((visibility("default")))
#endif

#define QF_API_VERSION 4

#pragma pack(push, 8)
typedef struct QFParams {
    int32_t struct_size;
    int32_t target_faces;        // желаемое число квадов (0 — взять из target_edge_length)
    double target_edge_length;   // 0 — взять из target_faces
    double sharp_angle;          // рад
    int32_t detect_sharp;
    int32_t preserve_boundary;
    double adaptive;             // 0..1 адаптивная плотность по кривизне
    double curvature_align;      // 0..1 выравнивание рёбер по главным кривизнам
    double guide_strength;       // 0..1 влияние направляющих рёбер
    double density_strength;     // 0..1 влияние карты плотности
    int32_t seed;
    int32_t threads;             // 0 = авто
    int32_t use_gpu;
    int32_t quality;             // 0 быстро .. 3 макс
    int32_t relax_iterations;    // релаксация результата (-1 = нет, 0 = авто)
    int32_t debug;
    int32_t verbose;
    int32_t pad0;
    double working_scale;        // k: длина ребра рабочей сетки = L/k (0 = авто)
    double field_scale;          // шаг сетки поля = L*field_scale (0 = авто)
    int32_t i_reserved[16];
    double d_reserved[16];       // [11] — предел градиента размера квадов (0 = 0.3)
} QFParams;

typedef struct QFInput {
    int32_t struct_size;
    int32_t nv;
    const float* V;              // float32[nv*3]
    int32_t nf;
    int32_t n_sharp;
    const int32_t* F;            // int32[nf*3] треугольники
    const int32_t* sharp_edges;  // пары вершин (n_sharp пар) — жёсткие рёбра
    int32_t n_guide;
    int32_t pad0;
    const int32_t* guide_edges;  // пары вершин (n_guide пар) — направляющие
    const float* vdensity;       // плотность на вершину (1 = обычная), может быть NULL
    int32_t reserved[8];
} QFInput;

typedef struct QFStatus {
    volatile int32_t cancel;
    volatile float progress;
    volatile int32_t stage;
    volatile int32_t done;
    char message[192];
} QFStatus;
#pragma pack(pop)

QF_API int qf_api_version(void);
QF_API int qf_default_params(QFParams* p);
QF_API void qf_set_threads(int n);
QF_API int qf_get_threads(void);

// Запуск конвейера. Возвращает дескриптор результата (всегда не NULL); проверять qf_result_ok.
// Функция потокобезопасна относительно других вызовов (без глобального состояния, кроме пула).
QF_API void* qf_run(const QFParams* p, const QFInput* in, QFStatus* st);

QF_API int qf_result_ok(void* r);
QF_API const char* qf_result_error(void* r);
QF_API const char* qf_result_log(void* r);
// Итоговая полигональная сетка: float32 вершины, face_sizes[nfaces], индексы подряд.
QF_API int qf_result_mesh(void* r, const float** verts, int* nverts, const int** face_sizes,
                          int* nfaces, const int** indices, int* nindices);
// Именованные массивы для диагностики: dtype 0 = float64, 1 = int32.
QF_API int qf_result_array(void* r, const char* name, const void** data, int64_t* n, int* dtype);
QF_API void qf_result_free(void* r);
