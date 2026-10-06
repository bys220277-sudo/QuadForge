// QuadForge native core — базовые типы, математика, параллелизм, прогресс.
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

namespace qf {

using i64 = int64_t;
using u64 = uint64_t;
constexpr double kPi = 3.14159265358979323846;
constexpr double kHalfPi = kPi * 0.5;

// ---------------------------------------------------------------------------
// Векторы
// ---------------------------------------------------------------------------
struct V3 {
    double x = 0, y = 0, z = 0;
    V3() = default;
    V3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    double& operator[](int i) { return (&x)[i]; }
    double operator[](int i) const { return (&x)[i]; }
    V3 operator+(const V3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    V3 operator-(const V3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    V3 operator*(double s) const { return {x * s, y * s, z * s}; }
    V3 operator/(double s) const { return {x / s, y / s, z / s}; }
    V3 operator-() const { return {-x, -y, -z}; }
    V3& operator+=(const V3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    V3& operator-=(const V3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    V3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};
inline V3 operator*(double s, const V3& v) { return v * s; }
inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(const V3& a, const V3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm2(const V3& a) { return dot(a, a); }
inline double norm(const V3& a) { return std::sqrt(dot(a, a)); }
inline V3 normalized(const V3& a, const V3& fallback = V3(0, 0, 1)) {
    double n = norm(a);
    return n > 1e-300 ? a / n : fallback;
}

struct V2 {
    double x = 0, y = 0;
    V2() = default;
    V2(double x_, double y_) : x(x_), y(y_) {}
    V2 operator+(const V2& o) const { return {x + o.x, y + o.y}; }
    V2 operator-(const V2& o) const { return {x - o.x, y - o.y}; }
    V2 operator*(double s) const { return {x * s, y * s}; }
};
inline double dot(const V2& a, const V2& b) { return a.x * b.x + a.y * b.y; }

// Комплексное число (без <complex>, чтобы не тянуть лишнего и держать предсказуемую
// арифметику).
struct C2 {
    double re = 0, im = 0;
    C2() = default;
    C2(double r, double i) : re(r), im(i) {}
    C2 operator+(const C2& o) const { return {re + o.re, im + o.im}; }
    C2 operator-(const C2& o) const { return {re - o.re, im - o.im}; }
    C2 operator*(const C2& o) const { return {re * o.re - im * o.im, re * o.im + im * o.re}; }
    C2 operator*(double s) const { return {re * s, im * s}; }
    C2& operator+=(const C2& o) { re += o.re; im += o.im; return *this; }
    C2& operator-=(const C2& o) { re -= o.re; im -= o.im; return *this; }
    double abs() const { return std::sqrt(re * re + im * im); }
    double abs2() const { return re * re + im * im; }
    double arg() const { return std::atan2(im, re); }
};
inline C2 conj(const C2& a) { return {a.re, -a.im}; }
inline C2 polar(double angle) { return {std::cos(angle), std::sin(angle)}; }

// ---------------------------------------------------------------------------
// Параллелизм
// ---------------------------------------------------------------------------
int num_threads();
void set_num_threads(int n);

// body(begin, end) вызывается для непересекающихся диапазонов [begin,end).
// Разбиение зависит только от n и grain (не от числа потоков) — результат
// редукций детерминирован.
void parallel_for(i64 n, const std::function<void(i64, i64)>& body, i64 grain = 2048);

template <class F>
inline void parallel_each(i64 n, F&& f, i64 grain = 2048) {
    parallel_for(n, [&](i64 b, i64 e) { for (i64 i = b; i < e; ++i) f(i); }, grain);
}

// Детерминированная сумма: порядок сложения не зависит от числа потоков.
double parallel_sum(i64 n, const std::function<double(i64, i64)>& partial, i64 grain = 4096);

// ---------------------------------------------------------------------------
// Прогресс / отмена / лог. Структура зеркалируется в Python через ctypes.
// ---------------------------------------------------------------------------
struct Status {
    volatile int32_t cancel;      // 1 — запрос отмены (пишет Python)
    volatile float progress;      // 0..1
    volatile int32_t stage;       // номер этапа
    volatile int32_t done;        // 1 — завершено
    char message[192];            // текущий этап
};

struct Cancelled {};  // исключение для досрочного выхода

struct Ctx {
    Status* st = nullptr;
    std::string log;
    bool verbose = false;
    void progress(float f, int stage, const char* msg) {
        if (st) {
            st->progress = f;
            st->stage = stage;
            std::snprintf(st->message, sizeof(st->message), "%s", msg);
        }
        check_cancel();
    }
    void check_cancel() const {
        if (st && st->cancel) throw Cancelled();
    }
    void logf(const char* fmt, ...) __attribute__((format(printf, 2, 3)));
};

double now_seconds();

}  // namespace qf
