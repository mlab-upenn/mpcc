/*
 * CasADi to FORCES Template - missing information to be filled in by createCasadi.m 
 * (C) embotech AG, Zurich, Switzerland, 2013-2020. All rights reserved.
 *
 * This file is part of the FORCES client, and carries the same license.
 */ 

#ifdef __cplusplus
extern "C" {
#endif
    
#include "include/kinematic_solver.h"

#define casadi_real kinematic_solver_float

#include "kinematic_solver_model.h" 
#include "kinematic_solver_model_helpers.h"

/* prototyes for models */
    

/* copies data from sparse matrix into a dense one */
static void sparse2fullcopy(solver_int32_default nrow, solver_int32_default ncol, const solver_int32_default *colidx, const solver_int32_default *row, kinematic_solver_float *data, kinematic_solver_float *out)
{
    solver_int32_default i, j;
    
    /* copy data into dense matrix */
    for(i=0; i<ncol; i++)
    {
        for( j=colidx[i]; j < colidx[i+1]; j++ )
        {
            out[i*nrow + row[j]] = data[j];
        }
    }
}

/* CasADi - FORCES interface */
extern void kinematic_solver_casadi2forces(kinematic_solver_float *x,        /* primal vars                                         */
                                 kinematic_solver_float *y,        /* eq. constraint multiplers                           */
                                 kinematic_solver_float *l,        /* ineq. constraint multipliers                        */
                                 kinematic_solver_float *p,        /* parameters                                          */
                                 kinematic_solver_float *f,        /* objective function (scalar)                         */
                                 kinematic_solver_float *nabla_f,  /* gradient of objective function                      */
                                 kinematic_solver_float *c,        /* dynamics                                            */
                                 kinematic_solver_float *nabla_c,  /* Jacobian of the dynamics (column major)             */
                                 kinematic_solver_float *h,        /* inequality constraints                              */
                                 kinematic_solver_float *nabla_h,  /* Jacobian of inequality constraints (column major)   */
                                 kinematic_solver_float *hess,     /* Hessian (column major)                              */
                                 solver_int32_default stage,     /* stage number (0 indexed)                            */
								 solver_int32_default iteration /* iteration number of solver                          */)
{
    /* CasADi input and output arrays */
    const kinematic_solver_float *in[4];
    kinematic_solver_float *out[7];
    
    /* temporary storage for casadi sparse output */
    kinematic_solver_float this_f;
    kinematic_solver_float nabla_f_sparse[6];
    kinematic_solver_float h_sparse[1];
    kinematic_solver_float nabla_h_sparse[3];
    kinematic_solver_float c_sparse[6];
    kinematic_solver_float nabla_c_sparse[23];
            
    
    /* pointers to row and column info for 
     * column compressed format used by CasADi */
    solver_int32_default nrow, ncol;
    const solver_int32_default *colind, *row;
    
    /* set inputs for CasADi */
    in[0] = x;
    in[1] = p; /* maybe should be made conditional */
    in[2] = l; /* maybe should be made conditional */     
    in[3] = y; /* maybe should be made conditional */

    	/* set outputs for CasADi */
	out[0] = &this_f;
	out[1] = nabla_f_sparse;
	if ((0 <= stage && stage <= 18))
	{
		/* set inputs */
		out[2] = h_sparse;
		out[3] = nabla_h_sparse;
		out[4] = c_sparse;
		out[5] = nabla_c_sparse;
		
		/* call CasADi */
		kinematic_solver_model_0(in, out);
		
		/* copy to dense */
		if( nabla_f )
		{
			kinematic_solver_model_0_sparsity(3, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, nabla_f_sparse, nabla_f);
		}
		if( c )
		{
			kinematic_solver_model_0_sparsity(6, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, c_sparse, c);
		}
		if( nabla_c )
		{
			kinematic_solver_model_0_sparsity(7, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, nabla_c_sparse, nabla_c);
		}
		if( h )
		{
			kinematic_solver_model_0_sparsity(4, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, h_sparse, h);
		}
		if( nabla_h )
		{
			kinematic_solver_model_0_sparsity(5, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, nabla_h_sparse, nabla_h);
		}
		
	}
	/* set outputs for CasADi */
	out[0] = &this_f;
	out[1] = nabla_f_sparse;
	if ((19 == stage))
	{
		/* set inputs */
		out[2] = h_sparse;
		out[3] = nabla_h_sparse;
		/* call CasADi */
		kinematic_solver_model_1(in, out);
		/* copy to dense */
		if( nabla_f )
		{
			kinematic_solver_model_1_sparsity(3, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, nabla_f_sparse, nabla_f);
		}
		if( h )
		{
			kinematic_solver_model_1_sparsity(4, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, h_sparse, h);
		}
		if( nabla_h )
		{
			kinematic_solver_model_1_sparsity(5, &nrow, &ncol, &colind, &row);
			sparse2fullcopy(nrow, ncol, colind, row, nabla_h_sparse, nabla_h);
		}
		
	}         
    
    /* add to objective */
    if( f )
    {
        *f += this_f;
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif