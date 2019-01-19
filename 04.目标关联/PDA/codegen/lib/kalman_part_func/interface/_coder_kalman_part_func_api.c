/*
 * File: _coder_kalman_part_func_api.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_kalman_part_func_api.h"
#include "_coder_kalman_part_func_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131466U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "kalman_part_func",                  /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4];
static const mxArray *b_emlrt_marshallOut(const real_T u[16]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *p_in,
  const char_T *identifier))[16];
static const mxArray *c_emlrt_marshallOut(const real_T u[4]);
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[16];
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *A, const
  char_T *identifier, real_T y[16]);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *x_in,
  const char_T *identifier))[4];
static const mxArray *emlrt_marshallOut(const real_T u[16]);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[16]);
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[16];
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[16]);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[4]
 */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4]
{
  real_T (*y)[4];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const real_T u[16]
 * Return Type  : const mxArray *
 */
  static const mxArray *b_emlrt_marshallOut(const real_T u[16])
{
  const mxArray *y;
  const mxArray *m1;
  static const int32_T iv1[2] = { 0, 0 };

  static const int32_T iv2[2] = { 4, 4 };

  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv1, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m1, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m1, iv2, 2);
  emlrtAssign(&y, m1);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *p_in
 *                const char_T *identifier
 * Return Type  : real_T (*)[16]
 */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *p_in,
  const char_T *identifier))[16]
{
  real_T (*y)[16];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(p_in), &thisId);
  emlrtDestroyArray(&p_in);
  return y;
}
/*
 * Arguments    : const real_T u[4]
 * Return Type  : const mxArray *
 */
  static const mxArray *c_emlrt_marshallOut(const real_T u[4])
{
  const mxArray *y;
  const mxArray *m2;
  static const int32_T iv3[1] = { 0 };

  static const int32_T iv4[1] = { 4 };

  y = NULL;
  m2 = emlrtCreateNumericArray(1, iv3, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m2, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m2, iv4, 1);
  emlrtAssign(&y, m2);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[16]
 */
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[16]
{
  real_T (*y)[16];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *A
 *                const char_T *identifier
 *                real_T y[16]
 * Return Type  : void
 */
  static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *A, const
  char_T *identifier, real_T y[16])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(A), &thisId, y);
  emlrtDestroyArray(&A);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *x_in
 *                const char_T *identifier
 * Return Type  : real_T (*)[4]
 */
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *x_in,
  const char_T *identifier))[4]
{
  real_T (*y)[4];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(x_in), &thisId);
  emlrtDestroyArray(&x_in);
  return y;
}
/*
 * Arguments    : const real_T u[16]
 * Return Type  : const mxArray *
 */
  static const mxArray *emlrt_marshallOut(const real_T u[16])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 4, 4 };

  real_T *pData;
  int32_T i0;
  int32_T i;
  int32_T b_i;
  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m0);
  i0 = 0;
  for (i = 0; i < 4; i++) {
    for (b_i = 0; b_i < 4; b_i++) {
      pData[i0] = u[b_i + (i << 2)];
      i0++;
    }
  }

  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[16]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[16])
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[4]
 */
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4]
{
  real_T (*ret)[4];
  static const int32_T dims[1] = { 4 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[16]
 */
  static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[16]
{
  real_T (*ret)[16];
  static const int32_T dims[2] = { 4, 4 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[16])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[16]
 * Return Type  : void
 */
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[16])
{
  static const int32_T dims[2] = { 4, 4 };

  int32_T i1;
  int32_T i2;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  for (i1 = 0; i1 < 4; i1++) {
    for (i2 = 0; i2 < 4; i2++) {
      ret[i2 + (i1 << 2)] = (*(real_T (*)[16])emlrtMxGetData(src))[i2 + (i1 << 2)];
    }
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray * const prhs[2]
 *                int32_T nlhs
 *                const mxArray *plhs[4]
 * Return Type  : void
 */
void kalman_part_func_api(const mxArray * const prhs[2], int32_T nlhs, const
  mxArray *plhs[4])
{
  real_T (*P_predic)[16];
  real_T (*Z_predic)[4];
  real_T (*S)[16];
  real_T (*K)[16];
  real_T (*x_in)[4];
  real_T (*p_in)[16];
  const mxArray *tmp;
  real_T A[16];
  real_T Q[16];
  real_T C[16];
  real_T R[16];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  P_predic = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  Z_predic = (real_T (*)[4])mxMalloc(sizeof(real_T [4]));
  S = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  K = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));

  /* Marshall function inputs */
  x_in = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "x_in");
  p_in = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "p_in");

  /* Marshall in global variables */
  tmp = emlrtGetGlobalVariable("A");
  if (tmp != NULL) {
    e_emlrt_marshallIn(&st, tmp, "A", A);
  }

  tmp = emlrtGetGlobalVariable("Q");
  if (tmp != NULL) {
    e_emlrt_marshallIn(&st, tmp, "Q", Q);
  }

  tmp = emlrtGetGlobalVariable("C");
  if (tmp != NULL) {
    e_emlrt_marshallIn(&st, tmp, "C", C);
  }

  tmp = emlrtGetGlobalVariable("R");
  if (tmp != NULL) {
    e_emlrt_marshallIn(&st, tmp, "R", R);
  }

  /* Invoke the target function */
  kalman_part_func(*x_in, *p_in, *P_predic, *Z_predic, *S, *K);

  /* Marshall out global variables */
  emlrtPutGlobalVariable("A", emlrt_marshallOut(A));
  emlrtPutGlobalVariable("Q", emlrt_marshallOut(Q));
  emlrtPutGlobalVariable("C", emlrt_marshallOut(C));
  emlrtPutGlobalVariable("R", emlrt_marshallOut(R));

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(*P_predic);
  if (nlhs > 1) {
    plhs[1] = c_emlrt_marshallOut(*Z_predic);
  }

  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(*S);
  }

  if (nlhs > 3) {
    plhs[3] = b_emlrt_marshallOut(*K);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void kalman_part_func_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  kalman_part_func_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void kalman_part_func_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void kalman_part_func_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_kalman_part_func_api.c
 *
 * [EOF]
 */
