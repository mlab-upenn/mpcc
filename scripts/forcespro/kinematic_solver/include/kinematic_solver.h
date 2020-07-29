/*
kinematic_solver : A fast customized optimization solver.

Copyright (C) 2013-2020 EMBOTECH AG [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

/* Generated by FORCES PRO v3.1.0 on Wednesday, July 29, 2020 at 10:39:38 AM */

#ifndef SOLVER_STDIO_H
#define SOLVER_STDIO_H
#include <stdio.h>
#endif

#ifndef kinematic_solver_H
#define kinematic_solver_H

/* DATA TYPE ------------------------------------------------------------*/
typedef double kinematic_solver_float;

typedef double kinematic_solverinterface_float;

#ifndef SOLVER_STANDARD_TYPES
#define SOLVER_STANDARD_TYPES

typedef signed char solver_int8_signed;
typedef unsigned char solver_int8_unsigned;
typedef char solver_int8_default;
typedef signed short int solver_int16_signed;
typedef unsigned short int solver_int16_unsigned;
typedef short int solver_int16_default;
typedef signed int solver_int32_signed;
typedef unsigned int solver_int32_unsigned;
typedef int solver_int32_default;
typedef signed long long int solver_int64_signed;
typedef unsigned long long int solver_int64_unsigned;
typedef long long int solver_int64_default;

#endif


/* SOLVER SETTINGS ------------------------------------------------------*/

/* MISRA-C compliance */
#ifndef MISRA_C_kinematic_solver
#define MISRA_C_kinematic_solver (0)
#endif

/* restrict code */
#ifndef RESTRICT_CODE_kinematic_solver
#define RESTRICT_CODE_kinematic_solver (0)
#endif

/* print level */
#ifndef SET_PRINTLEVEL_kinematic_solver
#define SET_PRINTLEVEL_kinematic_solver    (2)
#endif

/* timing */
#ifndef SET_TIMING_kinematic_solver
#define SET_TIMING_kinematic_solver    (1)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define SET_MAXIT_kinematic_solver			(20)	

/* scaling factor of line search (FTB rule) */
#define SET_FLS_SCALE_kinematic_solver		(kinematic_solver_float)(0.99)      

/* maximum number of supported elements in the filter */
#define MAX_FILTER_SIZE_kinematic_solver	(20) 

/* maximum number of supported elements in the filter */
#define MAX_SOC_IT_kinematic_solver			(4) 

/* desired relative duality gap */
#define SET_ACC_RDGAP_kinematic_solver		(kinematic_solver_float)(0.0001)

/* desired maximum residual on equality constraints */
#define SET_ACC_RESEQ_kinematic_solver		(kinematic_solver_float)(1E-06)

/* desired maximum residual on inequality constraints */
#define SET_ACC_RESINEQ_kinematic_solver	(kinematic_solver_float)(1E-06)

/* desired maximum violation of complementarity */
#define SET_ACC_KKTCOMPL_kinematic_solver	(kinematic_solver_float)(1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define OPTIMAL_kinematic_solver      (1)

/* maximum number of iterations has been reached */
#define MAXITREACHED_kinematic_solver (0)

/* wrong number of inequalities error */
#define INVALID_NUM_INEQ_ERROR_kinematic_solver  (-4)

/* factorization error */
#define FACTORIZATION_ERROR_kinematic_solver   (-5)

/* NaN encountered in function evaluations */
#define BADFUNCEVAL_kinematic_solver  (-6)

/* no progress in method possible */
#define NOPROGRESS_kinematic_solver   (-7)

/* invalid values in parameters */
#define PARAM_VALUE_ERROR_kinematic_solver   (-11)

/* licensing error - solver not valid on this machine */
#define LICENSE_ERROR_kinematic_solver  (-100)

/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct
{
    /* vector of size 180 */
    kinematic_solver_float x0[180];

    /* vector of size 6 */
    kinematic_solver_float xinit[6];

    /* vector of size 240 */
    kinematic_solver_float all_parameters[240];


} kinematic_solver_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct
{
    /* vector of size 9 */
    kinematic_solver_float x01[9];

    /* vector of size 9 */
    kinematic_solver_float x02[9];

    /* vector of size 9 */
    kinematic_solver_float x03[9];

    /* vector of size 9 */
    kinematic_solver_float x04[9];

    /* vector of size 9 */
    kinematic_solver_float x05[9];

    /* vector of size 9 */
    kinematic_solver_float x06[9];

    /* vector of size 9 */
    kinematic_solver_float x07[9];

    /* vector of size 9 */
    kinematic_solver_float x08[9];

    /* vector of size 9 */
    kinematic_solver_float x09[9];

    /* vector of size 9 */
    kinematic_solver_float x10[9];

    /* vector of size 9 */
    kinematic_solver_float x11[9];

    /* vector of size 9 */
    kinematic_solver_float x12[9];

    /* vector of size 9 */
    kinematic_solver_float x13[9];

    /* vector of size 9 */
    kinematic_solver_float x14[9];

    /* vector of size 9 */
    kinematic_solver_float x15[9];

    /* vector of size 9 */
    kinematic_solver_float x16[9];

    /* vector of size 9 */
    kinematic_solver_float x17[9];

    /* vector of size 9 */
    kinematic_solver_float x18[9];

    /* vector of size 9 */
    kinematic_solver_float x19[9];

    /* vector of size 9 */
    kinematic_solver_float x20[9];


} kinematic_solver_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct
{
    /* iteration number */
    solver_int32_default it;

	/* number of iterations needed to optimality (branch-and-bound) */
	solver_int32_default it2opt;
	
    /* inf-norm of equality constraint residuals */
    kinematic_solver_float res_eq;
	
    /* inf-norm of inequality constraint residuals */
    kinematic_solver_float res_ineq;

	/* norm of stationarity condition */
    kinematic_solver_float rsnorm;

	/* max of all complementarity violations */
    kinematic_solver_float rcompnorm;

    /* primal objective */
    kinematic_solver_float pobj;	
	
    /* dual objective */
    kinematic_solver_float dobj;	

    /* duality gap := pobj - dobj */
    kinematic_solver_float dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    kinematic_solver_float rdgap;		

    /* duality measure */
    kinematic_solver_float mu;

	/* duality measure (after affine step) */
    kinematic_solver_float mu_aff;
	
    /* centering parameter */
    kinematic_solver_float sigma;
	
    /* number of backtracking line search steps (affine direction) */
    solver_int32_default lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    solver_int32_default lsit_cc;
    
    /* step size (affine direction) */
    kinematic_solver_float step_aff;
    
    /* step size (combined direction) */
    kinematic_solver_float step_cc;    

	/* solvertime */
	kinematic_solver_float solvetime;   

	/* time spent in function evaluations */
	kinematic_solver_float fevalstime;  

} kinematic_solver_info;







/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* Time of Solver Generation: (UTC) Wednesday, July 29, 2020 10:39:40 AM */
/* User License expires on: (UTC) Thursday, January 21, 2021 10:00:00 PM (approx.) (at the time of code generation) */
/* Solver Static License expires on: (UTC) Thursday, January 21, 2021 10:00:00 PM (approx.) */
/* Solver Generation Request Id: eea18b0f-2b3c-4558-ba29-ed24a55507fb */
/* examine exitflag before using the result! */
#ifdef __cplusplus
extern "C" {
#endif		

typedef void (*kinematic_solver_extfunc)(kinematic_solver_float* x, kinematic_solver_float* y, kinematic_solver_float* lambda, kinematic_solver_float* params, kinematic_solver_float* pobj, kinematic_solver_float* g, kinematic_solver_float* c, kinematic_solver_float* Jeq, kinematic_solver_float* h, kinematic_solver_float* Jineq, kinematic_solver_float* H, solver_int32_default stage, solver_int32_default iterations);

extern solver_int32_default kinematic_solver_solve(kinematic_solver_params *params, kinematic_solver_output *output, kinematic_solver_info *info, FILE *fs, kinematic_solver_extfunc evalextfunctions_kinematic_solver);	





#ifdef __cplusplus
}
#endif

#endif
