// Tests for SymmBlockTriMat: construction, block accessors, to_dense and the dense LU solve.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <iostream>

#include "linalg.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

Mat2 mat2(double a, double b, double c, double d)
{
    Mat2 m;
    m << a, b,
         c, d;
    return m;
}

// Three blocks where every stored entry is distinct, so a block in the wrong place
// or a missing transpose changes the dense matrix.
// The diagonal blocks are symmetric and dominate, so the matrix is invertible.
SymmBlockTriMat three_block_matrix()
{
    SymmBlockTriMat A(3);
    A.diag(0) = mat2(10, 1, 1, 20);
    A.diag(1) = mat2(11, 2, 2, 21);
    A.diag(2) = mat2(12, 3, 3, 22);
    A.upper(0) = mat2(1, 2, 3, 4);   // couples blocks 0 and 1
    A.upper(1) = mat2(5, 6, 7, 8);   // couples blocks 1 and 2
    return A;
}

void test_construction()
{
    const SymmBlockTriMat A(4);
    CHECK(A.num_blocks() == 4);
    CHECK(A.n() == 8);

    const Mat M = A.to_dense();
    CHECK(M.rows() == 8 && M.cols() == 8);
    CHECK(M.isZero());
}

void test_accessors()
{
    SymmBlockTriMat A(2);
    A.diag(1) = mat2(1, 2, 2, 3);
    A.upper(0) += mat2(4, 5, 6, 7);
    A.upper(0) += mat2(1, 1, 1, 1);   // writes accumulate

    // The const overloads read the same stored blocks
    const SymmBlockTriMat& cA = A;
    CHECK(cA.diag(0).isZero());
    CHECK(cA.diag(1) == mat2(1, 2, 2, 3));
    CHECK(cA.upper(0) == mat2(5, 6, 7, 8));
}

void test_to_dense()
{
    const SymmBlockTriMat A = three_block_matrix();

    Mat expected(6, 6);
    expected << 10,  1,  1,  2,  0,  0,
                 1, 20,  3,  4,  0,  0,
                 1,  3, 11,  2,  5,  6,
                 2,  4,  2, 21,  7,  8,
                 0,  0,  5,  7, 12,  3,
                 0,  0,  6,  8,  3, 22;

    const Mat M = A.to_dense();
    CHECK(M == expected);
    CHECK(M == M.transpose());

    // A single block has no off-diagonal blocks
    SymmBlockTriMat B(1);
    B.diag(0) = mat2(1, 2, 2, 3);
    CHECK(B.to_dense() == mat2(1, 2, 2, 3));
}

void test_solve()
{
    const SymmBlockTriMat A = three_block_matrix();
    const Mat M = A.to_dense();

    // Manufactured solution: b = M x_true, then solve recovers x_true
    Vec x_true(6);
    x_true << 1.0, -2.0, 0.5, 3.0, -1.5, 2.0;
    const Vec b = M * x_true;

    Vec x(6);
    A.solve(x, b);
    CHECK((x - x_true).norm() <= 1e-12 * x_true.norm());

    // The right-hand side may be an expression, and x may be a view into a larger vector
    Vec y = Vec::Zero(8);
    A.solve(y.segment(1, 6), -b);
    CHECK((y.segment(1, 6) + x_true).norm() <= 1e-12 * x_true.norm());
    CHECK(y(0) == 0.0 && y(7) == 0.0);
}

} // namespace

int main()
{
    test_construction();
    test_accessors();
    test_to_dense();
    test_solve();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}
