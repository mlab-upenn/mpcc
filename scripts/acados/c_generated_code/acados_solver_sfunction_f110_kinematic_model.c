/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#define S_FUNCTION_NAME   acados_solver_sfunction_f110_kinematic_model
#define S_FUNCTION_LEVEL  2

#define MDL_START

// acados
#include "acados/utils/print.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
#include "f110_kinematic_model_model/f110_kinematic_model_model.h"
#include "acados_solver_f110_kinematic_model.h"

#include "simstruc.h"

//#define SAMPLINGTIME -1
#define SAMPLINGTIME 0.05

static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

              

    // specify the number of input ports
    if ( !ssSetNumInputPorts(S, 8) )
        return;

    // specify the number of output ports
    if ( !ssSetNumOutputPorts(S, 6) )
        return;

    // specify dimension information for the input ports
    // x0
    ssSetInputPortVectorDimension(S, 0, 7);
    // parameters
    ssSetInputPortVectorDimension(S, 1, (20+1) * 12);
    // lbx
    ssSetInputPortVectorDimension(S, 2, 4);
    // ubx
    ssSetInputPortVectorDimension(S, 3, 4);
    // lbu
    ssSetInputPortVectorDimension(S, 4, 3);
    // ubu
    ssSetInputPortVectorDimension(S, 5, 3);
    // lh
    ssSetInputPortVectorDimension(S, 6, 1);
    // uh
    ssSetInputPortVectorDimension(S, 7, 1);

    // specify dimension information for the output ports
    ssSetOutputPortVectorDimension(S, 0, 3 ); // optimal input
    ssSetOutputPortVectorDimension(S, 1, 1 ); // solver status
    ssSetOutputPortVectorDimension(S, 2, 1 ); // KKT residuals
    ssSetOutputPortVectorDimension(S, 3, 7 ); // first state
    ssSetOutputPortVectorDimension(S, 4, 1); // computation times
    ssSetOutputPortVectorDimension(S, 5, 1 ); // sqp iter

    // specify the direct feedthrough status
    // should be set to 1 for all inputs used in mdlOutputs
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    ssSetInputPortDirectFeedThrough(S, 4, 1);
    ssSetInputPortDirectFeedThrough(S, 5, 1);
    ssSetInputPortDirectFeedThrough(S, 6, 1);
    ssSetInputPortDirectFeedThrough(S, 7, 1);

    // one sample time
    ssSetNumSampleTimes(S, 1);
}


#if defined(MATLAB_MEX_FILE)

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

#endif /* MATLAB_MEX_FILE */


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLINGTIME);
    ssSetOffsetTime(S, 0, 0.0);
}


static void mdlStart(SimStruct *S)
{
    acados_create();
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    InputRealPtrsType in_sign;
    

    // local buffer
    real_t buffer[12];


    /* go through inputs */
    // initial condition
    in_sign = ssGetInputPortRealSignalPtrs(S, 0);
    for (int i = 0; i < 7; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", buffer);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", buffer);







    // parameters - stage-variant !!!
    in_sign = ssGetInputPortRealSignalPtrs(S, 1);

    // update value of parameters
    for (int ii = 0; ii <= 20; ii++) 
    {
        for (int jj = 0; jj < 12; jj++)
            buffer[jj] = (double)(*in_sign[ii*12+jj]);
        acados_update_params(ii, buffer, 12);
    }



    // lbx
    in_sign = ssGetInputPortRealSignalPtrs(S, 2);

    for (int i = 0; i < 4; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int ii = 1; ii < 20; ii++)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, ii, "lbx", buffer);

    // ubx
    in_sign = ssGetInputPortRealSignalPtrs(S, 3);

    for (int i = 0; i < 4; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int ii = 1; ii < 20; ii++)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, ii, "ubx", buffer);


    // lbu
    in_sign = ssGetInputPortRealSignalPtrs(S, 4);

    for (int i = 0; i < 3; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int ii = 0; ii < 20; ii++)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, ii, "lbu", buffer);

    // ubu
    in_sign = ssGetInputPortRealSignalPtrs(S, 5);

    for (int i = 0; i < 3; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int ii = 0; ii < 20; ii++)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, ii, "ubu", buffer);




    // lh
    in_sign = ssGetInputPortRealSignalPtrs(S, 6);

    for (int i = 0; i < 1; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int ii = 0; ii < 20; ii++)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, ii, "lh", buffer);

    // uh
    in_sign = ssGetInputPortRealSignalPtrs(S, 7);

    for (int i = 0; i < 1; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int ii = 0; ii < 20; ii++)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, ii, "uh", buffer);

    /* call solver */
    int rti_phase = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_phase", &rti_phase);
    int acados_status = acados_solve();


    /* set outputs */
    // assign pointers to output signals
    real_t *out_u0, *out_status, *out_sqp_iter, *out_KKT_res, *out_x1, *out_cpu_time;
    int tmp_int;

    out_u0          = ssGetOutputPortRealSignal(S, 0);
    out_status      = ssGetOutputPortRealSignal(S, 1);
    out_KKT_res     = ssGetOutputPortRealSignal(S, 2);
    out_x1          = ssGetOutputPortRealSignal(S, 3);
    out_cpu_time    = ssGetOutputPortRealSignal(S, 4);
    out_sqp_iter    = ssGetOutputPortRealSignal(S, 5);

    // extract solver info
    *out_status = (real_t) acados_status;
    *out_KKT_res = (real_t) nlp_out->inf_norm_res;
//    *out_cpu_time = (real_t) nlp_out->total_time;
    
    // get solution time
    ocp_nlp_get(nlp_config, nlp_solver, "time_tot", (void *) out_cpu_time);

    // get sqp iter
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", (void *) &tmp_int);
    *out_sqp_iter = (real_t) tmp_int;
//    *out_sqp_iter = (real_t) nlp_out->sqp_iter;

    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", (void *) out_u0);

    // get next state
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 1, "x", (void *) out_x1);

}

static void mdlTerminate(SimStruct *S)
{
    acados_free();
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
