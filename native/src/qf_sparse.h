// Разреженный Холецкий (up-looking, CSparse-подобный) + упорядочение вложенным рассечением.
#pragma once

#include "qf_common.h"

namespace qf {

// Симметричная матрица в CSR (полная: обе половины), строки отсортированы.
struct SymCsr {
    int n = 0;
    std::vector<i64> off;
    std::vector<int> col;
    std::vector<double> val;
};

// Упорядочение рекурсивной координатной бисекцией с разделителями (nested dissection).
// coords: 3*n координат вершины для каждой переменной. Возвращает perm[new] = old.
std::vector<int> nested_dissection_order(const SymCsr& A, const std::vector<double>& coords, int leaf = 64);

// Супернодальный мультифронтальный Холецкий: упорядочение (perm) + постпорядок дерева
// исключения, ослабленная амальгамация супернодов, плотные блочные ядра во фронтах.
class SparseCholesky {
public:
    // A — полная симметричная матрица, perm[new]=old (может быть пустым).
    // max_nnz > 0: отказ (false), если заполнение L превышает лимит.
    bool factor(const SymCsr& A, const std::vector<int>& perm, double reg = 0.0, i64 max_nnz = 0);
    // Решение A x = b (в исходной нумерации).
    void solve(const std::vector<double>& b, std::vector<double>& x) const;
    void solve_inplace(double* bx, double* work) const;  // b -> x, work размера n
    i64 nnz() const { return nnz_; }
    int size() const { return n; }
    int supernodes() const { return ns; }
    bool ok() const { return ok_; }

private:
    int n = 0, ns = 0;
    bool ok_ = false;
    i64 nnz_ = 0;
    std::vector<int> perm_, iperm_;
    std::vector<int> sn_c0;      // ns+1: первый столбец супернода
    std::vector<i64> sn_roff;    // ns+1: смещение структуры строк
    std::vector<int> sn_rows;    // строки супернода (сначала собственные столбцы), по возрастанию
    std::vector<i64> sn_xoff;    // ns+1: смещение блока m x k (по столбцам) в Lx
    std::vector<double> Lx;
};

// Плотный Холецкий (блочный, параллельный): A — n x n по строкам, на выходе нижний треугольник L.
bool dense_cholesky(double* A, int n);
// Решение L L^T x = b (b -> x на месте).
void dense_cholesky_solve(const double* L, int n, double* b);

}  // namespace qf
