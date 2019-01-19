/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 19-Jan-2019 14:01:44
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "kalman_part_func.h"
#include "main.h"
#include "kalman_part_func_terminate.h"
#include "kalman_part_func_initialize.h"

/* Function Declarations */
static void argInit_4x1_real_T(double result[4]);
static void argInit_4x4_real_T(double result[16]);
static double argInit_real_T(void);
static void main_kalman_part_func(void);

/* Function Definitions */

/*
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_4x1_real_T(double result[4])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[16]
 * Return Type  : void
 */
static void argInit_4x4_real_T(double result[16])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    for (idx1 = 0; idx1 < 4; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + (idx1 << 2)] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_kalman_part_func(void)
{
  double dv5[4];
  double dv6[16];
  double P_predic[16];
  double Z_predic[4];
  double S[16];
  double K[16];

  /* Initialize function 'kalman_part_func' input arguments. */
  /* Initialize function input argument 'x_in'. */
  /* Initialize function input argument 'p_in'. */
  /* Call the entry-point 'kalman_part_func'. */
  argInit_4x1_real_T(dv5);
  argInit_4x4_real_T(dv6);
  kalman_part_func(dv5, dv6, P_predic, Z_predic, S, K);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  kalman_part_func_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_kalman_part_func();

  /* Terminate the application.
     You do not need to do this more than one time. */
  kalman_part_func_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
