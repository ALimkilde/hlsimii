// Tests for BlockTriMat: construction, block accessors, to_dense, setZero, the dense LU solve, row scaling
// and the matrix-vector product. Also tests the block Thomas solver BlockTriLU against the dense solve.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <cstdlib>
#include <iostream>
#include <limits>

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

// Three blocks where every stored entry is distinct, so a block in the wrong place,
// or a lower block mixed up with a transposed upper block, changes the dense matrix.
// The diagonal blocks dominate their rows, so the matrix is invertible.
BlockTriMat three_block_matrix()
{
    BlockTriMat A(3);
    A.diag(0) = mat2(100, 1, 1, 200);
    A.diag(1) = mat2(110, 2, 2, 210);
    A.diag(2) = mat2(120, 3, 3, 220);
    A.upper(0) = mat2(1, 2, 3, 4);       // block (0, 1)
    A.upper(1) = mat2(5, 6, 7, 8);       // block (1, 2)
    A.lower(0) = mat2(31, 32, 33, 34);   // block (1, 0)
    A.lower(1) = mat2(35, 36, 37, 38);   // block (2, 1)
    return A;
}

void test_construction()
{
    const BlockTriMat A(4);
    CHECK(A.num_blocks() == 4);
    CHECK(A.n() == 8);

    const Mat M = A.to_dense();
    CHECK(M.rows() == 8 && M.cols() == 8);
    CHECK(M.isZero());
}

void test_accessors()
{
    BlockTriMat A(2);
    A.diag(1) = mat2(1, 2, 2, 3);
    A.upper(0) += mat2(4, 5, 6, 7);
    A.upper(0) += mat2(1, 1, 1, 1);   // writes accumulate
    A.lower(0) = mat2(9, 8, 7, 6);

    // The const overloads read the same stored blocks
    const BlockTriMat& cA = A;
    CHECK(cA.diag(0).isZero());
    CHECK(cA.diag(1) == mat2(1, 2, 2, 3));
    CHECK(cA.upper(0) == mat2(5, 6, 7, 8));
    CHECK(cA.lower(0) == mat2(9, 8, 7, 6));
}

void test_to_dense()
{
    const BlockTriMat A = three_block_matrix();

    Mat expected(6, 6);
    expected << 100,   1,   1,   2,   0,   0,
                  1, 200,   3,   4,   0,   0,
                 31,  32, 110,   2,   5,   6,
                 33,  34,   2, 210,   7,   8,
                  0,   0,  35,  36, 120,   3,
                  0,   0,  37,  38,   3, 220;

    const Mat M = A.to_dense();
    CHECK(M == expected);

    // A single block has no off-diagonal blocks
    BlockTriMat B(1);
    B.diag(0) = mat2(1, 2, 2, 3);
    CHECK(B.to_dense() == mat2(1, 2, 2, 3));
}

// setZero clears every block, including the lower ones
void test_set_zero()
{
    BlockTriMat A = three_block_matrix();
    A.setZero();
    CHECK(A.to_dense().isZero());
}

void test_solve()
{
    const BlockTriMat A = three_block_matrix();
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

// multiply_rowwise(f) takes one factor per block and scales both rows of block row b
// by f(b), i.e. A <- diag(f0, f0, f1, f1, ...) A.
// Every block row gets a different factor, so scaling columns instead of rows, using
// the factor of the wrong block row, or skipping a block all change the result.
// The factors are small integers, so the products are exact and == can be used.
void test_multiply_rowwise()
{
    BlockTriMat A = three_block_matrix();
    const Mat M = A.to_dense();

    Vec f(3);
    f << 2, 5, 11;
    Vec f_rows(6);                   // f with each factor repeated for both rows of its block
    f_rows << 2, 2, 5, 5, 11, 11;

    A.multiply_rowwise(f);
    CHECK(A.to_dense() == f_rows.asDiagonal() * M);

    // Spot checks against hand-computed entries, one per kind of block
    const Mat S = A.to_dense();
    CHECK(S(1, 0) == 2 * 1.0);      // diag(0), block row 0
    CHECK(S(1, 3) == 2 * 4.0);      // upper(0), block row 0
    CHECK(S(2, 0) == 5 * 31.0);     // lower(0), block row 1
    CHECK(S(5, 3) == 11 * 38.0);    // lower(1), block row 2

    // A single block has only a diagonal block
    BlockTriMat B(1);
    B.diag(0) = mat2(1, 2, 3, 4);
    Vec g(1);
    g << 3;
    B.multiply_rowwise(g);
    CHECK(B.to_dense() == mat2(3, 6, 9, 12));

    // The factors may be a view into a larger vector
    BlockTriMat C = three_block_matrix();
    Vec h = Vec::Zero(5);
    h.segment(1, 3) = f;
    C.multiply_rowwise(h.segment(1, 3));
    CHECK(C.to_dense() == f_rows.asDiagonal() * M);
}

// A * x must equal the dense product: block row b is
// diag(b) x_b + upper(b) x_{b+1} + lower(b-1) x_{b-1}.
// The entries are small integers, so the products are exact and == can be used.
void test_multiply_vector()
{
    const BlockTriMat A = three_block_matrix();
    const Mat M = A.to_dense();

    Vec x(6);
    x << 1, -2, 3, -4, 5, -6;
    const Vec y = A * x;
    CHECK(y.size() == 6);
    CHECK(y == M * x);

    // Hand-computed middle row, which uses all three blocks:
    // lower(0) row 0 . x_0 + diag(1) row 0 . x_1 + upper(1) row 0 . x_2
    // = (31 - 64) + (330 - 8) + (25 - 36) = 278
    CHECK(y(2) == 278.0);

    // A single block has only a diagonal block, so no off-diagonal term may be read
    BlockTriMat B(1);
    B.diag(0) = mat2(1, 2, 3, 4);
    Vec x1(2);
    x1 << 5, -1;
    CHECK(B * x1 == mat2(1, 2, 3, 4) * x1);

    // Two blocks: each block row has exactly one off-diagonal neighbour
    BlockTriMat C(2);
    C.diag(0) = mat2(1, 2, 3, 4);
    C.diag(1) = mat2(5, 6, 7, 8);
    C.upper(0) = mat2(9, 10, 11, 12);
    C.lower(0) = mat2(13, 14, 15, 16);
    Vec x2(4);
    x2 << 1, 2, -3, 4;
    CHECK(C * x2 == C.to_dense() * x2);

    // x may be a view into a larger vector
    Vec z = Vec::Zero(8);
    z.segment(1, 6) = x;
    CHECK(A * z.segment(1, 6) == M * x);
}

// Random block tridiagonal matrix whose diagonal blocks dominate their rows,
// so it is invertible and block Thomas (which only pivots inside a block) is stable.
BlockTriMat random_matrix(int num_blocks)
{
    BlockTriMat A(num_blocks);
    for (int b = 0; b < num_blocks; b++) {
        A.diag(b) = Mat2::Random() + 10.0 * Mat2::Identity();
        if (b < num_blocks - 1) {
            A.upper(b) = Mat2::Random();
            A.lower(b) = Mat2::Random();
        }
    }
    return A;
}

// Factorizes A with BlockTriLU and checks that solve recovers a manufactured solution.
// x starts out as NaN, so any block the solver forgets to write makes the check fail.
void check_lu_solve(const BlockTriMat& A)
{
    const int n = A.n();
    const Vec x_true = Vec::LinSpaced(n, -1.0, 2.0);
    const Vec b = A.to_dense() * x_true;

    BlockTriLU lu(A.num_blocks());
    lu.compute(A);

    Vec x = Vec::Constant(n, std::numeric_limits<double>::quiet_NaN());
    lu.solve(x, b);
    CHECK((x - x_true).norm() <= 1e-12 * x_true.norm());
    CHECK((A * x - b).norm() <= 1e-12 * b.norm());

    // Same answer as the dense reference solve
    Vec x_dense(n);
    A.solve(x_dense, b);
    CHECK((x - x_dense).norm() <= 1e-12 * x_true.norm());
}

void test_lu_solve()
{
    // Hand-built matrix with distinct, non-symmetric blocks
    check_lu_solve(three_block_matrix());

    // Edge cases: one block (no off-diagonals) and two blocks (one of each)
    std::srand(1);
    check_lu_solve(random_matrix(1));
    check_lu_solve(random_matrix(2));

    // A longer chain, closer to a real line mesh
    check_lu_solve(random_matrix(50));
}

void test_lu_solve_views_and_reuse()
{
    std::srand(2);
    const BlockTriMat A = random_matrix(5);
    const Mat M = A.to_dense();
    const Vec x_true = Vec::LinSpaced(10, 3.0, -1.0);
    const Vec b = M * x_true;

    BlockTriLU lu(5);
    lu.compute(A);

    // The right-hand side may be an expression, and x may be a view into a larger vector
    Vec y = Vec::Zero(12);
    lu.solve(y.segment(1, 10), -b);
    CHECK((y.segment(1, 10) + x_true).norm() <= 1e-12 * x_true.norm());
    CHECK(y(0) == 0.0 && y(11) == 0.0);

    // In-place solve: x and b are the same vector
    Vec v = b;
    lu.solve(v, v);
    CHECK((v - x_true).norm() <= 1e-12 * x_true.norm());

    // Repeated solves with the same factorization give the same answer
    Vec x1(10), x2(10);
    lu.solve(x1, b);
    lu.solve(x2, b);
    CHECK(x1 == x2);

    // Refactorizing the same object with a different matrix uses only the new matrix
    const BlockTriMat B = random_matrix(5);
    const Vec c = B.to_dense() * x_true;
    lu.compute(B);
    Vec x3(10);
    lu.solve(x3, c);
    CHECK((x3 - x_true).norm() <= 1e-12 * x_true.norm());
}

} // namespace

int main()
{
    test_construction();
    test_accessors();
    test_to_dense();
    test_set_zero();
    test_solve();
    test_multiply_rowwise();
    test_multiply_vector();
    test_lu_solve();
    test_lu_solve_views_and_reuse();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}
