#include "../include/matops.h"

ERROR_CODE mat_bra3_dot_mat33(const double bra[3],
                              const double mat[3][3],
                                    double res[3]) {
  memset(res, 0, sizeof(double) * 3);
  /*  // Slow
   * for (int j=0; j!=3; ++j) {
   *   for (int k=0; k!=3; ++k) {
   *     temp[k] += a[j] * b[j][k];
   *   }
   * }
   */ 
  res[0] = bra[0] * mat[0][0]  +  bra[1] * mat[1][0]  +  bra[2] * mat[2][0];
  res[1] = bra[0] * mat[0][1]  +  bra[1] * mat[1][1]  +  bra[2] * mat[2][1];
  res[2] = bra[0] * mat[0][2]  +  bra[1] * mat[1][2]  +  bra[2] * mat[2][2];
  return SUCCESS;
}

ERROR_CODE mat_mat33_inverse(const double amat[3][3],
                             double bmat[3][3]) {

/*
 *  bcell = adjoint(acell) / det(A)
 */
  double det = mat_mat33_det(amat);
  if (fabs(det) < 1E-5) {
    error_puts("Cannot evaluate the inverse of current matrix: Singular Matrix.");
    return MATRIX_SINGULAR_MATRIX;
  }
  double invdet = 1 / mat_mat33_det(amat);
  bmat[0][0] = (amat[1][1] * amat[2][2]  -  amat[2][1] * amat[1][2]) * invdet;
  bmat[0][1] = (amat[0][2] * amat[2][1]  -  amat[0][1] * amat[2][2]) * invdet;
  bmat[0][2] = (amat[0][1] * amat[1][2]  -  amat[0][2] * amat[1][1]) * invdet;
  bmat[1][0] = (amat[1][2] * amat[2][0]  -  amat[1][0] * amat[2][2]) * invdet;
  bmat[1][1] = (amat[0][0] * amat[2][2]  -  amat[0][2] * amat[2][0]) * invdet;
  bmat[1][2] = (amat[1][0] * amat[0][2]  -  amat[0][0] * amat[1][2]) * invdet;
  bmat[2][0] = (amat[1][0] * amat[2][1]  -  amat[2][0] * amat[1][1]) * invdet;
  bmat[2][1] = (amat[2][0] * amat[0][1]  -  amat[0][0] * amat[2][1]) * invdet;
  bmat[2][2] = (amat[0][0] * amat[1][1]  -  amat[1][0] * amat[0][1]) * invdet;
  return SUCCESS;
}

inline double mat_mat33_det(const double mat[3][3]) {
  return mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
         mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
         mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

#define __SWAP_DOUBLE(A, B) do {\
  double tmp = (A);\
  (A) = (B);\
  (B) = tmp;\
} while (0)

ERROR_CODE mat_mat33_trans(double mat[3][3]) {
  __SWAP_DOUBLE(mat[0][1], mat[1][0]);
  __SWAP_DOUBLE(mat[0][2], mat[2][0]);
  __SWAP_DOUBLE(mat[1][2], mat[2][1]);
  return SUCCESS;
}
#undef __SWAP_DOUBLE

ERROR_CODE mat_matx3_dot_mat33(const double matA[][3],
                               const int    rowA,
                               const double matB[3][3],
                                     double res[3][3]) {
  if (rowA <= 0) {
    error_puts("Invalid matrix A shape: number of rows.");
    return MATRIX_INVALID_NROW;
  } else {  }

  for (int k=0; k!=3; ++k) {
    for (int i=0; i!=rowA; ++i) {
      /*
       * for (int j=0; j!=3; ++j) {
       *   res[i][j] += matA[i][k] * matB[k][j];
       * }
       */
      res[i][0] += matA[i][k] * matB[k][0];
      res[i][1] += matA[i][k] * matB[k][1];
      res[i][2] += matA[i][k] * matB[k][2];
    }
  }

  return SUCCESS;
}


inline void mat_mat33_multiply(double mat[3][3],
                               double a) {
  mat[0][0] *= a;
  mat[0][1] *= a;
  mat[0][2] *= a;
  mat[1][0] *= a;
  mat[1][1] *= a;
  mat[1][2] *= a;
  mat[2][0] *= a;
  mat[2][1] *= a;
  mat[2][2] *= a;
}


