void FORCES_NLP_solver_model_0_sparsity(solver_int32_default i, solver_int32_default *nrow, solver_int32_default *ncol, const solver_int32_default **colind, const solver_int32_default **row);
void FORCES_NLP_solver_model_1_sparsity(solver_int32_default i, solver_int32_default *nrow, solver_int32_default *ncol, const solver_int32_default **colind, const solver_int32_default **row);
int FORCES_NLP_solver_model_0(const casadi_real **arg, casadi_real **res);
int FORCES_NLP_solver_model_1(const casadi_real **arg, casadi_real **res);