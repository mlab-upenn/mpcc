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

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
#include "f110_dynamic_model_model/f110_dynamic_model_model.h"



#include "f110_dynamic_model_constraints/f110_dynamic_model_h_constraint.h"


#include "f110_dynamic_model_cost/f110_dynamic_model_external_cost.h"


#include "acados_solver_f110_dynamic_model.h"

#define NX     9
#define NZ     0
#define NU     3
#define NP     13
#define NBX    3
#define NBX0   9
#define NBU    3
#define NSBX   0
#define NSBU   0
#define NSH    1
#define NSG    0
#define NSPHI  0
#define NSHN   0
#define NSGN   0
#define NSPHIN 0
#define NSBXN  0
#define NS     1
#define NSN    0
#define NG     0
#define NBXN   0
#define NGN    0
#define NY     0
#define NYN    0
#define N      20
#define NH     1
#define NPHI   0
#define NHN    0
#define NPHIN  0
#define NR     0


// ** global data **
ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_plan * nlp_solver_plan;
ocp_nlp_config * nlp_config;
ocp_nlp_dims * nlp_dims;
external_function_param_casadi * forw_vde_casadi;
external_function_param_casadi * expl_ode_fun;



external_function_param_casadi * nl_constr_h_fun;
external_function_param_casadi * nl_constr_h_fun_jac;


external_function_param_casadi nl_constr_h_e_fun_jac;
external_function_param_casadi nl_constr_h_e_fun;
external_function_param_casadi * ext_cost_fun;
external_function_param_casadi * ext_cost_fun_jac;
external_function_param_casadi * ext_cost_fun_jac_hess;


int acados_create()
{
    int status = 0;

    /************************************************
    *  plan & config
    ************************************************/
    nlp_solver_plan = ocp_nlp_plan_create(N);
    nlp_solver_plan->nlp_solver = SQP;
    

    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
    for (int i = 0; i < N; i++)
        nlp_solver_plan->nlp_cost[i] = EXTERNAL;

    nlp_solver_plan->nlp_cost[N] = LINEAR_LS;

    for (int i = 0; i < N; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = ERK;
    }

    for (int i = 0; i < N; i++)
    {
        nlp_solver_plan->nlp_constraints[i] = BGH;
    }
    nlp_solver_plan->nlp_constraints[N] = BGH;
    nlp_config = ocp_nlp_config_create(*nlp_solver_plan);


    /************************************************
    *  dimensions
    ************************************************/
    int nx[N+1];
    int nu[N+1];
    int nbx[N+1];
    int nbu[N+1];
    int nsbx[N+1];
    int nsbu[N+1];
    int nsg[N+1];
    int nsh[N+1];
    int nsphi[N+1];
    int ns[N+1];
    int ng[N+1];
    int nh[N+1];
    int nphi[N+1];
    int nz[N+1];
    int ny[N+1];
    int nr[N+1];
    int nbxe[N+1];

    for (int i = 0; i < N+1; i++)
    {
        // common
        nx[i]     = NX;
        nu[i]     = NU;
        nz[i]     = NZ;
        ns[i]     = NS;
        // cost
        ny[i]     = NY;
        // constraints
        nbx[i]    = NBX;
        nbu[i]    = NBU;
        nsbx[i]   = NSBX;
        nsbu[i]   = NSBU;
        nsg[i] = NSG;
        nsh[i]    = NSH;
        nsphi[i]  = NSPHI;
        ng[i]     = NG;
        nh[i]     = NH;
        nphi[i]   = NPHI;
        nr[i]     = NR;
        nbxe[i]   = 0;
    }

    // for initial state
    nbx[0]  = NBX0;
    nsbx[0] = 0;
    ns[0] = NS - NSBX;
    nbxe[0] = 9;

    // terminal - common
    nu[N]   = 0;
    nz[N]   = 0;
    ns[N]   = NSN;
    // cost
    ny[N]   = NYN;
    // constraint
    nbx[N]   = NBXN;
    nbu[N]   = 0;
    ng[N]    = NGN;
    nh[N]    = NHN;
    nphi[N]  = NPHIN;
    nr[N]    = 0;

    nsbx[N]  = NSBXN;
    nsbu[N]  = 0;
    nsg[N]   = NSGN;
    nsh[N]   = NSHN;
    nsphi[N] = NSPHIN;

    /* create and set ocp_nlp_dims */
    nlp_dims = ocp_nlp_dims_create(nlp_config);

    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "ns", ns);

    for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbx", &nsbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbu", &nsbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsg", &nsg[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbxe", &nbxe[i]);
    }

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nh", &nh[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nsh", &nsh[N]);
    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, N, "ny", &ny[N]);



    /************************************************
    *  external functions
    ************************************************/
    nl_constr_h_fun_jac = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        nl_constr_h_fun_jac[i].casadi_fun = &f110_dynamic_model_constr_h_fun_jac_uxt_zt;
        nl_constr_h_fun_jac[i].casadi_n_in = &f110_dynamic_model_constr_h_fun_jac_uxt_zt_n_in;
        nl_constr_h_fun_jac[i].casadi_n_out = &f110_dynamic_model_constr_h_fun_jac_uxt_zt_n_out;
        nl_constr_h_fun_jac[i].casadi_sparsity_in = &f110_dynamic_model_constr_h_fun_jac_uxt_zt_sparsity_in;
        nl_constr_h_fun_jac[i].casadi_sparsity_out = &f110_dynamic_model_constr_h_fun_jac_uxt_zt_sparsity_out;
        nl_constr_h_fun_jac[i].casadi_work = &f110_dynamic_model_constr_h_fun_jac_uxt_zt_work;
        external_function_param_casadi_create(&nl_constr_h_fun_jac[i], 13);
    }
    nl_constr_h_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        nl_constr_h_fun[i].casadi_fun = &f110_dynamic_model_constr_h_fun;
        nl_constr_h_fun[i].casadi_n_in = &f110_dynamic_model_constr_h_fun_n_in;
        nl_constr_h_fun[i].casadi_n_out = &f110_dynamic_model_constr_h_fun_n_out;
        nl_constr_h_fun[i].casadi_sparsity_in = &f110_dynamic_model_constr_h_fun_sparsity_in;
        nl_constr_h_fun[i].casadi_sparsity_out = &f110_dynamic_model_constr_h_fun_sparsity_out;
        nl_constr_h_fun[i].casadi_work = &f110_dynamic_model_constr_h_fun_work;
        external_function_param_casadi_create(&nl_constr_h_fun[i], 13);
    }
    
    


    // explicit ode
    forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        forw_vde_casadi[i].casadi_fun = &f110_dynamic_model_expl_vde_forw;
        forw_vde_casadi[i].casadi_n_in = &f110_dynamic_model_expl_vde_forw_n_in;
        forw_vde_casadi[i].casadi_n_out = &f110_dynamic_model_expl_vde_forw_n_out;
        forw_vde_casadi[i].casadi_sparsity_in = &f110_dynamic_model_expl_vde_forw_sparsity_in;
        forw_vde_casadi[i].casadi_sparsity_out = &f110_dynamic_model_expl_vde_forw_sparsity_out;
        forw_vde_casadi[i].casadi_work = &f110_dynamic_model_expl_vde_forw_work;
        external_function_param_casadi_create(&forw_vde_casadi[i], 13);
    }

    expl_ode_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        expl_ode_fun[i].casadi_fun = &f110_dynamic_model_expl_ode_fun;
        expl_ode_fun[i].casadi_n_in = &f110_dynamic_model_expl_ode_fun_n_in;
        expl_ode_fun[i].casadi_n_out = &f110_dynamic_model_expl_ode_fun_n_out;
        expl_ode_fun[i].casadi_sparsity_in = &f110_dynamic_model_expl_ode_fun_sparsity_in;
        expl_ode_fun[i].casadi_sparsity_out = &f110_dynamic_model_expl_ode_fun_sparsity_out;
        expl_ode_fun[i].casadi_work = &f110_dynamic_model_expl_ode_fun_work;
        external_function_param_casadi_create(&expl_ode_fun[i], 13);
    }


    // external cost
    ext_cost_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        ext_cost_fun[i].casadi_fun = &f110_dynamic_model_ext_cost_fun;
        ext_cost_fun[i].casadi_n_in = &f110_dynamic_model_ext_cost_fun_n_in;
        ext_cost_fun[i].casadi_n_out = &f110_dynamic_model_ext_cost_fun_n_out;
        ext_cost_fun[i].casadi_sparsity_in = &f110_dynamic_model_ext_cost_fun_sparsity_in;
        ext_cost_fun[i].casadi_sparsity_out = &f110_dynamic_model_ext_cost_fun_sparsity_out;
        ext_cost_fun[i].casadi_work = &f110_dynamic_model_ext_cost_fun_work;

        external_function_param_casadi_create(&ext_cost_fun[i], 13);
    }

    ext_cost_fun_jac = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        // residual function
        ext_cost_fun_jac[i].casadi_fun = &f110_dynamic_model_ext_cost_fun_jac;
        ext_cost_fun_jac[i].casadi_n_in = &f110_dynamic_model_ext_cost_fun_jac_n_in;
        ext_cost_fun_jac[i].casadi_n_out = &f110_dynamic_model_ext_cost_fun_jac_n_out;
        ext_cost_fun_jac[i].casadi_sparsity_in = &f110_dynamic_model_ext_cost_fun_jac_sparsity_in;
        ext_cost_fun_jac[i].casadi_sparsity_out = &f110_dynamic_model_ext_cost_fun_jac_sparsity_out;
        ext_cost_fun_jac[i].casadi_work = &f110_dynamic_model_ext_cost_fun_jac_work;

        external_function_param_casadi_create(&ext_cost_fun_jac[i], 13);
    }

    ext_cost_fun_jac_hess = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        // residual function
        ext_cost_fun_jac_hess[i].casadi_fun = &f110_dynamic_model_ext_cost_fun_jac_hess;
        ext_cost_fun_jac_hess[i].casadi_n_in = &f110_dynamic_model_ext_cost_fun_jac_hess_n_in;
        ext_cost_fun_jac_hess[i].casadi_n_out = &f110_dynamic_model_ext_cost_fun_jac_hess_n_out;
        ext_cost_fun_jac_hess[i].casadi_sparsity_in = &f110_dynamic_model_ext_cost_fun_jac_hess_sparsity_in;
        ext_cost_fun_jac_hess[i].casadi_sparsity_out = &f110_dynamic_model_ext_cost_fun_jac_hess_sparsity_out;
        ext_cost_fun_jac_hess[i].casadi_work = &f110_dynamic_model_ext_cost_fun_jac_hess_work;

        external_function_param_casadi_create(&ext_cost_fun_jac_hess[i], 13);
    }

    /************************************************
    *  nlp_in
    ************************************************/
    nlp_in = ocp_nlp_in_create(nlp_config, nlp_dims);

    double time_steps[N];
    time_steps[0] = 0.05;
    time_steps[1] = 0.05;
    time_steps[2] = 0.05;
    time_steps[3] = 0.05;
    time_steps[4] = 0.05;
    time_steps[5] = 0.05;
    time_steps[6] = 0.05;
    time_steps[7] = 0.05;
    time_steps[8] = 0.05;
    time_steps[9] = 0.05;
    time_steps[10] = 0.05;
    time_steps[11] = 0.05;
    time_steps[12] = 0.05;
    time_steps[13] = 0.05;
    time_steps[14] = 0.05;
    time_steps[15] = 0.05;
    time_steps[16] = 0.05;
    time_steps[17] = 0.05;
    time_steps[18] = 0.05;
    time_steps[19] = 0.05;

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_in_set(nlp_config, nlp_dims, nlp_in, i, "Ts", &time_steps[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "scaling", &time_steps[i]);
    }

    /**** Dynamics ****/
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_vde_forw", &forw_vde_casadi[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_ode_fun", &expl_ode_fun[i]);
    
    }


    /**** Cost ****/
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun", &ext_cost_fun[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun_jac", &ext_cost_fun_jac[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun_jac_hess", &ext_cost_fun_jac_hess[i]);
    }



    double Zl[NS];
    double Zu[NS];
    double zl[NS];
    double zu[NS];
    
    Zl[0] = 1000;

    
    Zu[0] = 1000;

    
    zl[0] = 1000;

    
    zu[0] = 1000;

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zl", Zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zu", Zu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zl", zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zu", zu);
    }


    // terminal cost





    /**** Constraints ****/

    // bounds for initial stage

    // x0
    int idxbx0[9];
    
    idxbx0[0] = 0;
    idxbx0[1] = 1;
    idxbx0[2] = 2;
    idxbx0[3] = 3;
    idxbx0[4] = 4;
    idxbx0[5] = 5;
    idxbx0[6] = 6;
    idxbx0[7] = 7;
    idxbx0[8] = 8;

    double lbx0[9];
    double ubx0[9];
    
    lbx0[0] = 0;
    ubx0[0] = 0;
    lbx0[1] = 0;
    ubx0[1] = 0;
    lbx0[2] = 0;
    ubx0[2] = 0;
    lbx0[3] = 1;
    ubx0[3] = 1;
    lbx0[4] = 0.01;
    ubx0[4] = 0.01;
    lbx0[5] = 0;
    ubx0[5] = 0;
    lbx0[6] = 0;
    ubx0[6] = 0;
    lbx0[7] = 0;
    ubx0[7] = 0;
    lbx0[8] = 0;
    ubx0[8] = 0;

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);


    // idxbxe_0
    int idxbxe_0[9];
    
    idxbxe_0[0] = 0;
    idxbxe_0[1] = 1;
    idxbxe_0[2] = 2;
    idxbxe_0[3] = 3;
    idxbxe_0[4] = 4;
    idxbxe_0[5] = 5;
    idxbxe_0[6] = 6;
    idxbxe_0[7] = 7;
    idxbxe_0[8] = 8;
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbxe", idxbxe_0);


    /* constraints that are the same for initial and intermediate */



    // u
    int idxbu[NBU];
    
    idxbu[0] = 0;
    idxbu[1] = 1;
    idxbu[2] = 2;
    double lbu[NBU];
    double ubu[NBU];
    
    lbu[0] = -10;
    ubu[0] = 10;
    lbu[1] = -2;
    ubu[1] = 2;
    lbu[2] = -0.1;
    ubu[2] = 5;

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxbu", idxbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lbu", lbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ubu", ubu);
    }







    // set up soft bounds for nonlinear constraints
    int idxsh[NSH];
    
    idxsh[0] = 0;
    double lsh[NSH];
    double ush[NSH];
    
    lsh[0] = 0.1;
    ush[0] = -0.1;

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsh", idxsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsh", lsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ush", ush);
    }





    // x
    int idxbx[NBX];
    
    idxbx[0] = 6;
    idxbx[1] = 7;
    idxbx[2] = 8;
    double lbx[NBX];
    double ubx[NBX];
    
    lbx[0] = 0;
    ubx[0] = 100;
    lbx[1] = -1;
    ubx[1] = 5;
    lbx[2] = -0.4;
    ubx[2] = 0.4;

    for (int i = 1; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxbx", idxbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lbx", lbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ubx", ubx);
    }





    // set up nonlinear constraints for stage 0 to N-1 
    double lh[NH];
    double uh[NH];

    
    lh[0] = -10;

    
    uh[0] = 0;
    
    for (int i = 0; i < N; i++)
    {
        // nonlinear constraints for stages 0 to N-1
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac",
                                     &nl_constr_h_fun_jac[i]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun",
                                    &nl_constr_h_fun[i]);
        
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }




    /* terminal constraints */

















    /************************************************
    *  opts
    ************************************************/

    nlp_opts = ocp_nlp_solver_opts_create(nlp_config, nlp_dims);



    int num_steps_val = 1;
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &num_steps_val);

    int ns_val = 4;
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &ns_val);

    int newton_iter_val = 3;
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_iter", &newton_iter_val);

    bool tmp_bool = false;
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &tmp_bool);

    double nlp_solver_step_length = 0.05;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "step_length", &nlp_solver_step_length);

    double levenberg_marquardt = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "levenberg_marquardt", &levenberg_marquardt);

    /* options QP solver */
    int qp_solver_cond_N;
    // NOTE: there is no condensing happening here!
    qp_solver_cond_N = N;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_cond_N", &qp_solver_cond_N);


    int qp_solver_iter_max = 50;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_iter_max", &qp_solver_iter_max);
    // set SQP specific options
    double nlp_solver_tol_stat = 0.0001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_stat", &nlp_solver_tol_stat);

    double nlp_solver_tol_eq = 0.0001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_eq", &nlp_solver_tol_eq);

    double nlp_solver_tol_ineq = 0.0001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_ineq", &nlp_solver_tol_ineq);

    double nlp_solver_tol_comp = 0.0001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_comp", &nlp_solver_tol_comp);

    int nlp_solver_max_iter = 100;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "max_iter", &nlp_solver_max_iter);

    int initialize_t_slacks = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "initialize_t_slacks", &initialize_t_slacks);

    int print_level = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "print_level", &print_level);


    int ext_cost_num_hess = 0;
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "cost_numerical_hessian", &ext_cost_num_hess);
    }


    /* out */
    nlp_out = ocp_nlp_out_create(nlp_config, nlp_dims);

    // initialize primal solution
    double x0[9];

    // initialize with x0
    
    x0[0] = 0;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = 1;
    x0[4] = 0.01;
    x0[5] = 0;
    x0[6] = 0;
    x0[7] = 0;
    x0[8] = 0;


    double u0[NU];
    
    u0[0] = 0.0;
    u0[1] = 0.0;
    u0[2] = 0.0;

    for (int i = 0; i < N; i++)
    {
        // x0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x0);
        // u0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
    }
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, N, "x", x0);
    
    nlp_solver = ocp_nlp_solver_create(nlp_config, nlp_dims, nlp_opts);



    // initialize parameters to nominal value
    double p[13];
    
    p[0] = 0;
    p[1] = 0;
    p[2] = 0;
    p[3] = 0;
    p[4] = 0;
    p[5] = 0;
    p[6] = 0;
    p[7] = 0;
    p[8] = 0;
    p[9] = 0;
    p[10] = 0;
    p[11] = 0;
    p[12] = 0;


    for (int ii = 0; ii < N; ii++)
    {
        forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
        expl_ode_fun[ii].set_param(expl_ode_fun+ii, p);
    }


    // cost
    for (int ii = 0; ii < N; ii++)
    {
        ext_cost_fun[ii].set_param(ext_cost_fun+ii, p);
        ext_cost_fun_jac[ii].set_param(ext_cost_fun_jac+ii, p);
        ext_cost_fun_jac_hess[ii].set_param(ext_cost_fun_jac_hess+ii, p);
    }

    // constraints

    for (int ii = 0; ii < N; ii++)
    {
        nl_constr_h_fun_jac[ii].set_param(nl_constr_h_fun_jac+ii, p);
        nl_constr_h_fun[ii].set_param(nl_constr_h_fun+ii, p);
    }

    status = ocp_nlp_precompute(nlp_solver, nlp_in, nlp_out);

    if (status != ACADOS_SUCCESS)
    {
        printf("\nocp_precompute failed!\n\n");
        exit(1);
    }

    return status;
}


int acados_update_params(int stage, double *p, int np)
{
    int solver_status = 0;

    int casadi_np = 13;
    if (casadi_np != np) {
        printf("acados_update_params: trying to set %i parameters for external functions."
            " External function has %i parameters. Exiting.\n", np, casadi_np);
        exit(1);
    }
    if (stage < 20)
    {
        forw_vde_casadi[stage].set_param(forw_vde_casadi+stage, p);
        expl_ode_fun[stage].set_param(expl_ode_fun+stage, p);
    

        // constraints
    
        nl_constr_h_fun_jac[stage].set_param(nl_constr_h_fun_jac+stage, p);
        nl_constr_h_fun[stage].set_param(nl_constr_h_fun+stage, p);

        // cost
        ext_cost_fun[stage].set_param(ext_cost_fun+stage, p);
        ext_cost_fun_jac[stage].set_param(ext_cost_fun_jac+stage, p);
        ext_cost_fun_jac_hess[stage].set_param(ext_cost_fun_jac_hess+stage, p);

    }
    else // stage == N
    {
        // terminal shooting node has no dynamics
        // cost
        // constraints
    
    }


    return solver_status;
}



int acados_solve()
{
    // solve NLP 
    int solver_status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);

    return solver_status;
}


int acados_free()
{
    // free memory
    ocp_nlp_solver_opts_destroy(nlp_opts);
    ocp_nlp_in_destroy(nlp_in);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_solver_destroy(nlp_solver);
    ocp_nlp_dims_destroy(nlp_dims);
    ocp_nlp_config_destroy(nlp_config);
    ocp_nlp_plan_destroy(nlp_solver_plan);

    /* free external function */
    // dynamics
    for (int i = 0; i < 20; i++)
    {
        external_function_param_casadi_free(&forw_vde_casadi[i]);
        external_function_param_casadi_free(&expl_ode_fun[i]);
    }
    free(forw_vde_casadi);
    free(expl_ode_fun);

    // cost
    for (int i = 0; i < 20; i++)
    {
        external_function_param_casadi_free(&ext_cost_fun[i]);
        external_function_param_casadi_free(&ext_cost_fun_jac[i]);
        external_function_param_casadi_free(&ext_cost_fun_jac_hess[i]);
    }
    free(ext_cost_fun);
    free(ext_cost_fun_jac_hess);

    // constraints
    for (int i = 0; i < 20; i++)
    {
        external_function_param_casadi_free(&nl_constr_h_fun_jac[i]);
        external_function_param_casadi_free(&nl_constr_h_fun[i]);
    }
    free(nl_constr_h_fun_jac);
    free(nl_constr_h_fun);

    return 0;
}

ocp_nlp_in * acados_get_nlp_in() { return  nlp_in; }
ocp_nlp_out * acados_get_nlp_out() { return  nlp_out; }
ocp_nlp_solver * acados_get_nlp_solver() { return  nlp_solver; }
ocp_nlp_config * acados_get_nlp_config() { return  nlp_config; }
void * acados_get_nlp_opts() { return  nlp_opts; }
ocp_nlp_dims * acados_get_nlp_dims() { return  nlp_dims; }
ocp_nlp_plan * acados_get_nlp_plan() { return  nlp_solver_plan; }


void acados_print_stats()
{
    int sqp_iter, stat_m, stat_n, tmp_int;
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);
    ocp_nlp_get(nlp_config, nlp_solver, "stat_n", &stat_n);
    ocp_nlp_get(nlp_config, nlp_solver, "stat_m", &stat_m);

    
    double stat[1000];
    ocp_nlp_get(nlp_config, nlp_solver, "statistics", stat);

    int nrow = sqp_iter+1 < stat_m ? sqp_iter+1 : stat_m;

    printf("iter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\n");
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < stat_n + 1; j++)
        {
            if (j == 0 || j > 4)
            {
                tmp_int = (int) stat[i + j * nrow];
                printf("%d\t", tmp_int);
            }
            else
            {
                printf("%e\t", stat[i + j * nrow]);
            }
        }
        printf("\n");
    }
}
