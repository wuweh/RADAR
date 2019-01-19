/*
 * File: _coder_kalman_part_func_mex.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

/* Include Files */
#include "_coder_kalman_part_func_api.h"
#include "_coder_kalman_part_func_mex.h"

/* Function Declarations */
static void kalman_part_func_mexFunction(int32_T nlhs, mxArray *plhs[4], int32_T
  nrhs, const mxArray *prhs[2]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[4]
 *                int32_T nrhs
 *                const mxArray *prhs[2]
 * Return Type  : void
 */
static void kalman_part_func_mexFunction(int32_T nlhs, mxArray *plhs[4], int32_T
  nrhs, const mxArray *prhs[2])
{
  const mxArray *outputs[4];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        16, "kalman_part_func");
  }

  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 16,
                        "kalman_part_func");
  }

  /* Call the function. */
  kalman_part_func_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  kalman_part_func_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(kalman_part_func_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  kalman_part_func_initialize();

  /* Dispatch the entry-point. */
  kalman_part_func_mexFunction(nlhs, plhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_kalman_part_func_mex.c
 *
 * [EOF]
 */
