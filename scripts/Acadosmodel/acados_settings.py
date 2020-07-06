#!/usr/bin/python3
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from dynamic_model import dynamic_model
import scipy.linalg
import numpy as np


def acados_settings(Tf, N):
    # create render arguments
    ocp = AcadosOcp()

    # load model, (add in param file later)
    model, constraint = dynamic_model(0)

    # define acados ODE
    model_ac = AcadosModel()
    model_ac.f_impl_expr = model.f_impl_expr
    model_ac.f_expl_expr = model.f_expl_expr
    model_ac.x = model.x
    model_ac.xdot = model.xdot
    #inputvector
    model_ac.u = model.u
    #parameter vector
    model_ac.p = model.p
    model_ac.name = model.name
    ocp.model = model_ac



    # set dimensions
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    ocp.dims.nx = nx
    ocp.dims.np = 0
    ocp.dims.ny = ny
    ocp.dims.ny_e = ny_e
    ocp.dims.nbx = 1
    ocp.dims.nsbx = 0
    ocp.dims.nbu = nu
    ocp.dims.nu = nu
    ocp.dims.N = N
    ocp.dims.ns = 2

    # set cost
    CHANGE THIS
    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    unscale = N / Tf

    ocp.cost.W = unscale * scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Qe / unscale

    Vx = np.zeros((ny, nx))
    Vx[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx = Vx

    Vu = np.zeros((ny, nu))
    Vu[6, 0] = 1.0
    Vu[7, 1] = 1.0
    ocp.cost.Vu = Vu

    Vx_e = np.zeros((ny_e, nx))
    Vx_e[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx_e = Vx_e

    ocp.cost.zl = 100 * np.ones((ocp.dims.ns,))
    ocp.cost.zu = 100 * np.ones((ocp.dims.ns,))
    ocp.cost.Zl = 0 * np.ones((ocp.dims.ns,))
    ocp.cost.Zu = 0 * np.ones((ocp.dims.ns,))

    # set intial references
    ocp.cost.yref = np.array([1, 0, 0, 0, 0, 0, 0, 0])
    ocp.cost.yref_e = np.array([0, 0, 0, 0, 0, 0])

    # setting constraints
    ocp.constraints.lbx = np.array([-12])
    ocp.constraints.ubx = np.array([12])
    ocp.constraints.idxbx = np.array([1])
    ocp.constraints.lbu = np.array([model.dthrottle_min, model.ddelta_min])
    ocp.constraints.ubu = np.array([model.dthrottle_max, model.ddelta_max])
    ocp.constraints.idxbu = np.array([0, 1])
    # ocp.constraints.lsbx=np.zero s([1])
    # ocp.constraints.usbx=np.zeros([1])
    # ocp.constraints.idxsbx=np.array([1])
    ocp.constraints.lh = np.array(
        [
            constraint.along_min,
            constraint.alat_min,
            model.n_min,
            model.throttle_min,
            model.delta_min,
        ]
    )
    ocp.constraints.uh = np.array(
        [
            constraint.along_max,
            constraint.alat_max,
            model.n_max,
            model.throttle_max,
            model.delta_max,
        ]
    )
    ocp.constraints.lsh = np.zeros(ocp.dims.nsh)
    ocp.constraints.ush = np.zeros(ocp.dims.nsh)
    ocp.constraints.idxsh = np.array([0, 2])

    # set intial condition
    ocp.constraints.x0 = model.x0

    # set QP solver and integration
    ocp.solver_options.tf = Tf
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.nlp_solver_type = "SQP_RTI"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "ERK"
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.sim_method_num_steps = 3
    # ocp.solver_options.nlp_solver_step_length = 0.05
    # ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.tol = 1e-4
    # ocp.solver_options.nlp_solver_tol_comp = 1e-1

    # create solver
    acados_solver = AcadosOcpSolver(ocp, json_file="acados_ocp.json")

    print("solver created returning to main")
    return constraint, model, acados_solver

def main():
    """ Main program """
    Tf = 1
    N = 50
    acados_settings(Tf, N)
    return 0

if __name__ == "__main__":
    main()
