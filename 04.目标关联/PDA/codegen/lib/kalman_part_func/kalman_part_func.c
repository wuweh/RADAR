/*
 * File: kalman_part_func.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "kalman_part_func.h"
#include "inv.h"
#include "kalman_part_func_data.h"

/* Function Definitions */

/*
 * Arguments    : const double x_in[4]
 *                const double p_in[16]
 *                double P_predic[16]
 *                double Z_predic[4]
 *                double S[16]
 *                double K[16]
 * Return Type  : void
 */
void kalman_part_func(const double x_in[4], const double p_in[16], double
                      P_predic[16], double Z_predic[4], double S[16], double K
                      [16])
{
  int i0;
  int i1;
  double dv4[16];
  int i2;
  double d0;
  double b_P_predic[16];
  for (i0 = 0; i0 < 4; i0++) {
    for (i1 = 0; i1 < 4; i1++) {
      dv4[i0 + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        dv4[i0 + (i1 << 2)] += A[i0 + (i2 << 2)] * p_in[i2 + (i1 << 2)];
      }
    }

    Z_predic[i0] = 0.0;
    for (i1 = 0; i1 < 4; i1++) {
      d0 = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        d0 += dv4[i0 + (i2 << 2)] * A[i1 + (i2 << 2)];
      }

      P_predic[i0 + (i1 << 2)] = d0 + Q[i0 + (i1 << 2)];
      Z_predic[i0] += C[i0 + (i1 << 2)] * x_in[i1];
    }
  }

  for (i0 = 0; i0 < 4; i0++) {
    for (i1 = 0; i1 < 4; i1++) {
      dv4[i0 + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        dv4[i0 + (i1 << 2)] += C[i0 + (i2 << 2)] * P_predic[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      d0 = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        d0 += dv4[i0 + (i2 << 2)] * C[i1 + (i2 << 2)];
      }

      S[i0 + (i1 << 2)] = d0 + R[i0 + (i1 << 2)];
    }
  }

  invNxN(S, dv4);
  for (i0 = 0; i0 < 4; i0++) {
    for (i1 = 0; i1 < 4; i1++) {
      b_P_predic[i0 + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        b_P_predic[i0 + (i1 << 2)] += P_predic[i0 + (i2 << 2)] * C[i1 + (i2 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      K[i0 + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        K[i0 + (i1 << 2)] += b_P_predic[i0 + (i2 << 2)] * dv4[i2 + (i1 << 2)];
      }
    }
  }
}

/*
 * File trailer for kalman_part_func.c
 *
 * [EOF]
 */
