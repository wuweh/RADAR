/*
 * File: kalman_part_func.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

#ifndef KALMAN_PART_FUNC_H
#define KALMAN_PART_FUNC_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "kalman_part_func_types.h"

/* Function Declarations */
extern void kalman_part_func(const double x_in[4], const double p_in[16], double
  P_predic[16], double Z_predic[4], double S[16], double K[16]);

#endif

/*
 * File trailer for kalman_part_func.h
 *
 * [EOF]
 */
