#!/usr/bin/python3
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from dynamic_model import dynamic_model
import scipy.linalg
import numpy as npy


def acados_settings(Tf, N, track_lu_table, modelparams):
    # create render arguments
    ocp = AcadosOcp()

    #load model
    model, constraints = dynamic_model(modelparams)

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
    #parameter vector
    model_ac.z = model.z
    #external cost function
    model_ac.cost_expr_ext_cost = model.stage_cost
    #model_ac.cost_expr_ext_cost_e = 0
    #model_ac.con_h_expr = model.con_h_expr
    model_ac.name = model.name
    ocp.model = model_ac



    # set dimensions

    nx  = model.x.size()[0]
    nu  = model.u.size()[0]
    nz  = model.z.size()[0]
    np  = model.p.size()[0]
    #ny  = nu + nx
    #ny_e = nx


    ocp.dims.nx   = nx
    ocp.dims.nz   = nz
    #ocp.dims.ny   = ny
    #ocp.dims.ny_e  = ny_e
    ocp.dims.nu   = nu
    ocp.dims.np   = np
    #ocp.dims.nh = 2
    #ocp.dims.ny = ny
    #ocp.dims.ny_e = ny_e
    #ocp.dims.nbx = nx
    #ocp.dims.nsbx = 0
    ocp.dims.nbu = nu

    #number of soft on h constraints
    #ocp.dims.nsh = 2

    ocp.dims.N = N


    # set cost to casadi expression defined above
    ocp.cost.cost_type = "EXTERNAL"
    #ocp.cost.cost_type_e = "EXTERNAL"

    #not sure if needed
    unscale = N / Tf

    #constraints
    #stagewise halfspace constraints for tracks with slack
    #ocp.constraints.uh = npy.zeros(2)
    #ocp.constraints.lh = -1e9*npy.ones(2)
    #ocp.constraints.Jsh = 1

    # boxconstraints
    # change these later
    #ocp.constraints.lbx = -30 * npy.ones(nx)
    #ocp.constraints.ubx = 30 * npy.ones(nx)
    #ocp.constraints.idxbx = npy.arange(nx)
    #print("ocp.constraints.idxbx: ",ocp.constraints.idxbx)

    ocp.constraints.lbu = npy.array([model.ddot_min, model.deltadot_min, model.thetadot_min])
    ocp.constraints.ubu = npy.array([model.ddot_max, model.deltadot_max, model.thetadot_max])
    ocp.constraints.idxbu = npy.array([0, 1, 2])
    # ocp.constraints.lsbx=npy.zero s([1])
    # ocp.constraints.usbx=npy.zeros([1])
    # ocp.constraints.idxsbx=npy.array([1])

    #ocp.constraints.lsh = npy.zeros(ocp.dims.nsh)
    #ocp.constraints.ush = npy.zeros(ocp.dims.nsh)
    #ocp.constraints.idxsh = npy.array([0, 2])

    # set intial condition
    ocp.constraints.x0 = model.x0

    # set QP solver and integration
    ocp.solver_options.tf = Tf
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.nlp_solver_type ="SQP"#"SQP_RTI" #
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "ERK"
    ocp.parameter_values = npy.zeros(np)
    #ocp.solver_options.sim_method_num_stages = 4
    #ocp.solver_options.sim_method_num_steps = 3
    #ocp.solver_options.nlp_solver_step_length = 0.05
    ocp.solver_options.nlp_solver_max_iter = 3
    ocp.solver_options.tol = 1e-4
    #ocp.solver_options.print_level = 1
    # ocp.solver_options.nlp_solver_tol_comp = 1e-1


    # create solver
    acados_solver = AcadosOcpSolver(ocp, json_file="acados_ocp2.json")

    print("solver created returning to main")
    return constraints, model, acados_solver, ocp
