#pragma once

#include "base.h"
#include "error.h"
#include "debug.h"

/*
 * Matrix Operations, specialized for 3x3 matrix
 */


ERROR_CODE mat_bra3_dot_mat33(const double bra[3], 
                              const double mat[3][3], 
                                    double res[3]);

ERROR_CODE mat_mat33_inverse(const double amat[3][3],
                                   double bmat[3][3]);

double mat_mat33_det(const double mat[3][3]);

ERROR_CODE mat_mat33_trans(double mat[3][3]);

ERROR_CODE mat_matx3_dot_mat33(const double matA[][3],
                               const int    rowA,
                               const double matB[3][3],
                                     double res[][3]);

void mat_mat33_multiply(double mat[3][3],
                              double a);
