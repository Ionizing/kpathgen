#define __DEBUG_ON
#include <matops.h>

void test_mat_mat33_trans(const double mat[3][3]) {
  debug_puts("Before transpose:");
  debug_print_matd33(mat);
  double res[3][3] = {}; 
  memcpy(res, mat, 9 * sizeof(double));
  mat_mat33_trans(res);
  debug_puts("After transpose:");
  debug_print_matd33(res);
  puts("");
}

void test_mat_bra3_dot_mat33(const double bra[3],
                             const double mat[3][3]) {
  debug_print_vecd3(bra);
  debug_print_matd33(mat);
  double res[3] = {0};
  mat_bra3_dot_mat33(bra, mat, res);
  debug_print_vecd3(res);
  puts("");
}

void test_mat_matx3_dot_mat33(const double matA[][3],
                              const int    rowA,
                              const double matB[3][3]) {
  debug_print_matdx3(matA, rowA);
  debug_print_matd33(matB);
  double (*res)[3];
  res = (double (*)[3]) calloc(rowA, sizeof(double[3]));
  mat_matx3_dot_mat33(matA, rowA, matB, res);
  debug_print_matdx3(res, rowA);
  puts("");
}

void test_mat_mat33_det(const double mat[3][3]) {
  debug_print_matd33(mat);
  debug_printf("Determinant of this matrix is: %lf\n", mat_mat33_det(mat));
  puts("");
}

void test_mat_mat33_inverse(const double amat[3][3]) {
  debug_print_matd33(amat);
  double bmat[3][3] = {};
  mat_mat33_inverse(amat, bmat);
  debug_print_matd33(bmat);
  puts("");
}


int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);

  warning_puts("test_mat_mat33_trans.");
  const double mat33[3][3] = {{1, 2, 0}, {0, 1, 1}, {2, 0, 1}};
  test_mat_mat33_trans(mat33);
  puts("");

  warning_puts("test_mat_bra3_dot_mat33.");
  test_mat_bra3_dot_mat33(mat33[0], mat33);
  test_mat_bra3_dot_mat33(mat33[1], mat33);
  test_mat_bra3_dot_mat33(mat33[2], mat33);
  puts("");

  warning_puts("test_mat_matx3_dot_mat33.");
  test_mat_matx3_dot_mat33(mat33, 3, mat33);
  puts("");

  warning_puts("test_mat_mat33_det.");
  test_mat_mat33_det(mat33);
  puts("");

  warning_puts("test_mat_mat33_inverse");
  test_mat_mat33_inverse(mat33);
  puts("");
  return 0;
}
