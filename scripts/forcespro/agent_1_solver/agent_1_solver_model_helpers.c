#ifdef __cplusplus
extern "C" {
#endif

#include "include/agent_1_solver.h"
#include "agent_1_solver_model.h"

void agent_1_solver_model_0_sparsity(solver_int32_default i, solver_int32_default *nrow, solver_int32_default *ncol,const solver_int32_default **colind, const solver_int32_default **row)
{
const solver_int32_default *s;
if (i < 2) 
{
s = agent_1_solver_model_0_inner_sparsity_in(i);
} else {
s = agent_1_solver_model_0_inner_sparsity_out(i-2);
}
*nrow = s[0];
*ncol = s[1];
*colind = s + 2;
*row = s + 2 + (*ncol + 1);
}
void agent_1_solver_model_1_sparsity(solver_int32_default i, solver_int32_default *nrow, solver_int32_default *ncol,const solver_int32_default **colind, const solver_int32_default **row)
{
const solver_int32_default *s;
if (i < 2) 
{
s = agent_1_solver_model_1_inner_sparsity_in(i);
} else {
s = agent_1_solver_model_1_inner_sparsity_out(i-2);
}
*nrow = s[0];
*ncol = s[1];
*colind = s + 2;
*row = s + 2 + (*ncol + 1);
}

int agent_1_solver_model_0(const casadi_real **arg, casadi_real **res)
{
	return agent_1_solver_model_0_inner(arg, res, NULL, NULL, 0);
}
int agent_1_solver_model_1(const casadi_real **arg, casadi_real **res)
{
	return agent_1_solver_model_1_inner(arg, res, NULL, NULL, 0);
}

#ifdef __cplusplus
}
#endif