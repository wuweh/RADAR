/*
 * File: _coder_kalman_part_func_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

#ifndef _CODER_KALMAN_PART_FUNC_API_H
#define _CODER_KALMAN_PART_FUNC_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_kalman_part_func_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void kalman_part_func(real_T x_in[4], real_T p_in[16], real_T P_predic[16],
  real_T Z_predic[4], real_T S[16], real_T K[16]);
extern void kalman_part_func_api(const mxArray * const prhs[2], int32_T nlhs,
  const mxArray *plhs[4]);
extern void kalman_part_func_atexit(void);
extern void kalman_part_func_initialize(void);
extern void kalman_part_func_terminate(void);
extern void kalman_part_func_xil_terminate(void);

#endif

/*
 * File trailer for _coder_kalman_part_func_api.h
 *
 * [EOF]
 */
