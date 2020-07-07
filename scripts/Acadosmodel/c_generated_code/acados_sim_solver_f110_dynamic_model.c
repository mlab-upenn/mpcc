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
#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/print.h"


// example specific
#include "f110_dynamic_model_model/f110_dynamic_model_model.h"
#include "acados_sim_solver_f110_dynamic_model.h"


// ** global data **
sim_config  * f110_dynamic_model_sim_config;
sim_in      * f110_dynamic_model_sim_in;
sim_out     * f110_dynamic_model_sim_out;
void        * f110_dynamic_model_sim_dims;
sim_opts    * f110_dynamic_model_sim_opts;
sim_solver  * f110_dynamic_model_sim_solver;


external_function_param_casadi * sim_forw_vde_casadi;
external_function_param_casadi * sim_expl_ode_fun_casadi;



int f110_dynamic_model_acados_sim_create()
{
    // initialize
    int nx = 9;
    int nu = 3;
    int nz = 0;

    
    double Tsim = 0.02;

    
    // explicit ode
    sim_forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_expl_ode_fun_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));

    sim_forw_vde_casadi->casadi_fun = &f110_dynamic_model_expl_vde_forw;
    sim_forw_vde_casadi->casadi_n_in = &f110_dynamic_model_expl_vde_forw_n_in;
    sim_forw_vde_casadi->casadi_n_out = &f110_dynamic_model_expl_vde_forw_n_out;
    sim_forw_vde_casadi->casadi_sparsity_in = &f110_dynamic_model_expl_vde_forw_sparsity_in;
    sim_forw_vde_casadi->casadi_sparsity_out = &f110_dynamic_model_expl_vde_forw_sparsity_out;
    sim_forw_vde_casadi->casadi_work = &f110_dynamic_model_expl_vde_forw_work;
    external_function_param_casadi_create(sim_forw_vde_casadi, 13);

    sim_expl_ode_fun_casadi->casadi_fun = &f110_dynamic_model_expl_ode_fun;
    sim_expl_ode_fun_casadi->casadi_n_in = &f110_dynamic_model_expl_ode_fun_n_in;
    sim_expl_ode_fun_casadi->casadi_n_out = &f110_dynamic_model_expl_ode_fun_n_out;
    sim_expl_ode_fun_casadi->casadi_sparsity_in = &f110_dynamic_model_expl_ode_fun_sparsity_in;
    sim_expl_ode_fun_casadi->casadi_sparsity_out = &f110_dynamic_model_expl_ode_fun_sparsity_out;
    sim_expl_ode_fun_casadi->casadi_work = &f110_dynamic_model_expl_ode_fun_work;
    external_function_param_casadi_create(sim_expl_ode_fun_casadi, 13);

    

    // sim plan & config
    sim_solver_plan plan;
    plan.sim_solver = ERK;

    // create correct config based on plan
    f110_dynamic_model_sim_config = sim_config_create(plan);

    // sim dims
    f110_dynamic_model_sim_dims = sim_dims_create(f110_dynamic_model_sim_config);
    sim_dims_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims, "nx", &nx);
    sim_dims_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims, "nu", &nu);
    sim_dims_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims, "nz", &nz);


    // sim opts
    f110_dynamic_model_sim_opts = sim_opts_create(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims);
    int tmp_int = 4;
    sim_opts_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_opts, "num_stages", &tmp_int);
    tmp_int = 1;
    sim_opts_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_opts, "num_steps", &tmp_int);
    tmp_int = 3;
    sim_opts_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_opts, "newton_iter", &tmp_int);
    bool tmp_bool = false;
    sim_opts_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_opts, "jac_reuse", &tmp_bool);



    // sim in / out
    f110_dynamic_model_sim_in  = sim_in_create(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims);
    f110_dynamic_model_sim_out = sim_out_create(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims);
    sim_in_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims,
               f110_dynamic_model_sim_in, "T", &Tsim);

    // model functions
    f110_dynamic_model_sim_config->model_set(f110_dynamic_model_sim_in->model,
                 "expl_vde_for", sim_forw_vde_casadi);
    f110_dynamic_model_sim_config->model_set(f110_dynamic_model_sim_in->model,
                 "expl_ode_fun", sim_expl_ode_fun_casadi);

    // sim solver
    f110_dynamic_model_sim_solver = sim_solver_create(f110_dynamic_model_sim_config,
                                               f110_dynamic_model_sim_dims, f110_dynamic_model_sim_opts);

    /* initialize parameter values */
    
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
    sim_forw_vde_casadi[0].set_param(sim_forw_vde_casadi, p);
    sim_expl_ode_fun_casadi[0].set_param(sim_expl_ode_fun_casadi, p);
    

    /* initialize input */
    // x
    double x0[9];
    for (int ii = 0; ii < 9; ii++)
        x0[ii] = 0.0;

    sim_in_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims,
               f110_dynamic_model_sim_in, "x", x0);


    // u
    double u0[3];
    for (int ii = 0; ii < 3; ii++)
        u0[ii] = 0.0;

    sim_in_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims,
               f110_dynamic_model_sim_in, "u", u0);

    // S_forw
    double S_forw[108];
    for (int ii = 0; ii < 108; ii++)
        S_forw[ii] = 0.0;
    for (int ii = 0; ii < 9; ii++)
        S_forw[ii + ii * 9 ] = 1.0;


    sim_in_set(f110_dynamic_model_sim_config, f110_dynamic_model_sim_dims,
               f110_dynamic_model_sim_in, "S_forw", S_forw);

    int status = sim_precompute(f110_dynamic_model_sim_solver, f110_dynamic_model_sim_in, f110_dynamic_model_sim_out);

    return status;
}


int f110_dynamic_model_acados_sim_solve()
{
    // integrate dynamics using acados sim_solver
    int status = sim_solve(f110_dynamic_model_sim_solver,
                           f110_dynamic_model_sim_in, f110_dynamic_model_sim_out);
    if (status != 0)
        printf("error in f110_dynamic_model_acados_sim_solve()! Exiting.\n");

    return status;
}


int f110_dynamic_model_acados_sim_free()
{
    // free memory
    sim_solver_destroy(f110_dynamic_model_sim_solver);
    sim_in_destroy(f110_dynamic_model_sim_in);
    sim_out_destroy(f110_dynamic_model_sim_out);
    sim_opts_destroy(f110_dynamic_model_sim_opts);
    sim_dims_destroy(f110_dynamic_model_sim_dims);
    sim_config_destroy(f110_dynamic_model_sim_config);

    // free external function
    external_function_param_casadi_free(sim_forw_vde_casadi);
    external_function_param_casadi_free(sim_expl_ode_fun_casadi);

    return 0;
}


int f110_dynamic_model_acados_sim_update_params(double *p, int np)
{
    int status = 0;
    int casadi_np = 13;
    sim_forw_vde_casadi[0].set_param(sim_forw_vde_casadi, p);
    sim_expl_ode_fun_casadi[0].set_param(sim_expl_ode_fun_casadi, p);

    return status;
}

/* getters pointers to C objects*/
sim_config * f110_dynamic_model_acados_get_sim_config()
{
    return f110_dynamic_model_sim_config;
};

sim_in * f110_dynamic_model_acados_get_sim_in()
{
    return f110_dynamic_model_sim_in;
};

sim_out * f110_dynamic_model_acados_get_sim_out()
{
    return f110_dynamic_model_sim_out;
};

void * f110_dynamic_model_acados_get_sim_dims()
{
    return f110_dynamic_model_sim_dims;
};

sim_opts * f110_dynamic_model_acados_get_sim_opts()
{
    return f110_dynamic_model_sim_opts;
};

sim_solver  * f110_dynamic_model_acados_get_sim_solver()
{
    return f110_dynamic_model_sim_solver;
};

