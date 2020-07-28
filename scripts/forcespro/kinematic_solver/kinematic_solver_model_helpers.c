#ifdef __cplusplus
extern "C" {
#endif

#include "include/kinematic_solver.h"
#include "kinematic_solver_model.h"

void kinematic_solver_model_0_sparsity(solver_int32_default i, solver_int32_default *nrow, solver_int32_default *ncol,const solver_int32_default **colind, const solver_int32_default **row)
{
const solver_int32_default *s;
if (i < 2) 
{
s = kinematic_solver_model_0_inner_sparsity_in(i);
} else {
s = kinematic_solver_model_0_inner_sparsity_out(i-2);
}
*nrow = s[0];
*ncol = s[1];
*colind = s + 2;
*row = s + 2 + (*ncol + 1);
}
void kinematic_solver_model_1_sparsity(solver_int32_default i, solver_int32_default *nrow, solver_int32_default *ncol,const solver_int32_default **colind, const solver_int32_default **row)
{
const solver_int32_default *s;
if (i < 2) 
{
s = kinematic_solver_model_1_inner_sparsity_in(i);
} else {
s = kinematic_solver_model_1_inner_sparsity_out(i-2);
}
*nrow = s[0];
*ncol = s[1];
*colind = s + 2;
*row = s + 2 + (*ncol + 1);
}

int kinematic_solver_model_0(const casadi_real **arg, casadi_real **res)
{
	return kinematic_solver_model_0_inner(arg, res, NULL, NULL, 0);
}
int kinematic_solver_model_1(const casadi_real **arg, casadi_real **res)
{
	return kinematic_solver_model_1_inner(arg, res, NULL, NULL, 0);
}

#ifdef __cplusplus
}
#endif