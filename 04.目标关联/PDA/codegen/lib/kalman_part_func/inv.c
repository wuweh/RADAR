/*
 * File: inv.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "kalman_part_func.h"
#include "inv.h"

/* Function Definitions */

/*
 * Arguments    : const double x[16]
 *                double y[16]
 * Return Type  : void
 */
void invNxN(const double x[16], double y[16])
{
  int i3;
  double b_x[16];
  int j;
  signed char ipiv[4];
  int c;
  int iy;
  int k;
  signed char p[4];
  int ix;
  double smax;
  double s;
  int jy;
  int i;
  for (i3 = 0; i3 < 16; i3++) {
    y[i3] = 0.0;
    b_x[i3] = x[i3];
  }

  for (i3 = 0; i3 < 4; i3++) {
    ipiv[i3] = (signed char)(1 + i3);
  }

  for (j = 0; j < 3; j++) {
    c = j * 5;
    iy = 0;
    ix = c;
    smax = fabs(b_x[c]);
    for (k = 2; k <= 4 - j; k++) {
      ix++;
      s = fabs(b_x[ix]);
      if (s > smax) {
        iy = k - 1;
        smax = s;
      }
    }

    if (b_x[c + iy] != 0.0) {
      if (iy != 0) {
        ipiv[j] = (signed char)((j + iy) + 1);
        ix = j;
        iy += j;
        for (k = 0; k < 4; k++) {
          smax = b_x[ix];
          b_x[ix] = b_x[iy];
          b_x[iy] = smax;
          ix += 4;
          iy += 4;
        }
      }

      i3 = (c - j) + 4;
      for (i = c + 1; i < i3; i++) {
        b_x[i] /= b_x[c];
      }
    }

    iy = c;
    jy = c + 4;
    for (i = 1; i <= 3 - j; i++) {
      smax = b_x[jy];
      if (b_x[jy] != 0.0) {
        ix = c + 1;
        i3 = (iy - j) + 8;
        for (k = 5 + iy; k < i3; k++) {
          b_x[k] += b_x[ix] * -smax;
          ix++;
        }
      }

      jy += 4;
      iy += 4;
    }
  }

  for (i3 = 0; i3 < 4; i3++) {
    p[i3] = (signed char)(1 + i3);
  }

  for (k = 0; k < 3; k++) {
    if (ipiv[k] > 1 + k) {
      iy = p[ipiv[k] - 1];
      p[ipiv[k] - 1] = p[k];
      p[k] = (signed char)iy;
    }
  }

  for (k = 0; k < 4; k++) {
    c = p[k] - 1;
    y[k + ((p[k] - 1) << 2)] = 1.0;
    for (j = k; j + 1 < 5; j++) {
      if (y[j + (c << 2)] != 0.0) {
        for (i = j + 1; i + 1 < 5; i++) {
          y[i + (c << 2)] -= y[j + (c << 2)] * b_x[i + (j << 2)];
        }
      }
    }
  }

  for (j = 0; j < 4; j++) {
    iy = j << 2;
    for (k = 3; k >= 0; k--) {
      jy = k << 2;
      if (y[k + iy] != 0.0) {
        y[k + iy] /= b_x[k + jy];
        for (i = 0; i < k; i++) {
          y[i + iy] -= y[k + iy] * b_x[i + jy];
        }
      }
    }
  }
}

/*
 * File trailer for inv.c
 *
 * [EOF]
 */
