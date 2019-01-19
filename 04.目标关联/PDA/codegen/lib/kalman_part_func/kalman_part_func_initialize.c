/*
 * File: kalman_part_func_initialize.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

/* Include Files */
#include <string.h>
#include "rt_nonfinite.h"
#include "kalman_part_func.h"
#include "kalman_part_func_initialize.h"
#include "kalman_part_func_data.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void kalman_part_func_initialize(void)
{
  static const double dv0[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.010000000000000002,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.010000000000000002 };

  static const double dv1[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const double dv2[16] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0,
    0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.01 };

  static const double dv3[16] = { 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0 };

  rt_InitInfAndNaN(8U);
  memcpy(&R[0], &dv0[0], sizeof(double) << 4);
  memcpy(&C[0], &dv1[0], sizeof(double) << 4);
  memcpy(&Q[0], &dv2[0], sizeof(double) << 4);
  memcpy(&A[0], &dv3[0], sizeof(double) << 4);
}

/*
 * File trailer for kalman_part_func_initialize.c
 *
 * [EOF]
 */
