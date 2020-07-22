#!/usr/bin/python3
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from models import dynamic_model, kinematic_model
import scipy.linalg
import numpy as npy


def acados_settings_dyn(Tf, N, modelparams):
    #create render arguments
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
    model_ac.con_h_expr = model.con_h_expr
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
    ocp.dims.nh = 1
    #ocp.dims.ny = ny
    #ocp.dims.ny_e = ny_e
    ocp.dims.nbx = 3
    ocp.dims.nsbx = 0
    ocp.dims.nbu = nu

    #number of soft on h constraints
    ocp.dims.nsh = 1
    ocp.dims.ns = 1

    ocp.dims.N = N


    # set cost to casadi expression defined above
    ocp.cost.cost_type = "EXTERNAL"
    #ocp.cost.cost_type_e = "EXTERNAL"
    ocp.cost.zu = 1000 * npy.ones((ocp.dims.ns,))
    ocp.cost.Zu = 1000 * npy.ones((ocp.dims.ns,))
    ocp.cost.zl = 1000 * npy.ones((ocp.dims.ns,))
    ocp.cost.Zl = 1000 * npy.ones((ocp.dims.ns,))
    #not sure if needed
    #unscale = N / Tf

    #constraints
    #stagewise  constraints for tracks with slack
    ocp.constraints.uh = npy.array([0.00])
    ocp.constraints.lh = npy.array([-10])
    ocp.constraints.lsh = 0.1*npy.ones(ocp.dims.nsh)
    ocp.constraints.ush = -0.1*npy.ones(ocp.dims.nsh)
    ocp.constraints.idxsh = npy.array([0])
    #ocp.constraints.Jsh = 1

    # boxconstraints
    ocp.constraints.lbx = npy.array([model.theta_min, model.d_min, model.delta_min])
    ocp.constraints.ubx = npy.array([model.theta_max, model.d_max, model.delta_max])
    ocp.constraints.idxbx = npy.array([6,7,8])
    #ocp.constraints.lsbx= -0.1 * npy.ones(ocp.dims.nbx)
    #ocp.constraints.usbx= 0.1 * npy.ones(ocp.dims.nbx)
    #ocp.constraints.idxsbx= npy.array([6,7,8])
    #print("ocp.constraints.idxbx: ",ocp.constraints.idxbx)

    ocp.constraints.lbu = npy.array([model.ddot_min, model.deltadot_min, model.thetadot_min])
    ocp.constraints.ubu = npy.array([model.ddot_max, model.deltadot_max, model.thetadot_max])
    ocp.constraints.idxbu = npy.array([0, 1, 2])


    # set intial condition
    ocp.constraints.x0 = model.x0

    # set QP solver and integration
    ocp.solver_options.tf = Tf
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.nlp_solver_type = "SQP"#"SQP_RTI" #
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "ERK"
    ocp.parameter_values = npy.zeros(np)
    #ocp.solver_options.sim_method_num_stages = 4
    #ocp.solver_options.sim_method_num_steps = 3
    ocp.solver_options.nlp_solver_step_length = 0.05
    ocp.solver_options.nlp_solver_max_iter = 100
    ocp.solver_options.tol = 1e-4
    #ocp.solver_options.print_level = 1
    # ocp.solver_options.nlp_solver_tol_comp = 1e-1

    # create solver
    acados_solver = AcadosOcpSolver(ocp, json_file="acados_ocp2.json")

    print("solver created returning to main")
    return constraints, model, acados_solver, ocp


def acados_settings_kin(Tf, N, modelparams):
    # create render arguments
    ocp = AcadosOcp()

    #load model
    model, constraints = kinematic_model(modelparams)

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
    model_ac.con_h_expr = model.con_h_expr
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
    ocp.dims.nh = 1
    #ocp.dims.ny = ny
    #ocp.dims.ny_e = ny_e
    ocp.dims.nbx = 4
    #ocp.dims.nsbx = 4
    ocp.dims.nbu = nu

    #number of soft on h constraints
    #ocp.dims.nsh = 1
    #ocp.dims.ns = 1

    ocp.dims.N = N


    # set cost to casadi expression defined above
    ocp.cost.cost_type = "EXTERNAL"
    #ocp.cost.cost_type_e = "EXTERNAL"
    #ocp.cost.zu = 1000 * npy.ones((ocp.dims.ns,))
    #ocp.cost.zl = 1000 * npy.ones((ocp.dims.ns,))
    #ocp.cost.Zu = 1000 * npy.ones((ocp.dims.ns,))
    #ocp.cost.Zl = 1000 * npy.ones((ocp.dims.ns,))
    #not sure if needed
    #unscale = N / Tf

    #constraints
    #stagewise  constraints for tracks with slack
    ocp.constraints.uh = npy.array([0.00])
    ocp.constraints.lh = npy.array([-10])
    #ocp.constraints.lsh = 0.1*npy.ones(ocp.dims.nsh)
    #ocp.constraints.ush = 0.001*npy.ones(ocp.dims.nsh)
    #ocp.constraints.idxsh = npy.array([0])
    #ocp.constraints.Jsh = 1

    # boxconstraints
    # change these later
    ocp.constraints.lbx = npy.array([model.vx_min, model.theta_min, model.d_min, model.delta_min])
    ocp.constraints.ubx = npy.array([model.vx_max, model.theta_max, model.d_max, model.delta_max])
    ocp.constraints.idxbx = npy.array([3,4,5,6])
    #ocp.constraints.lsbx= npy.zeros(ocp.dims.nbx)
    #ocp.constraints.usbx= npy.zeros(ocp.dims.nbx)
    #ocp.constraints.idxsbx= npy.array([3,4,5,6])
    #print("ocp.constraints.idxbx: ",ocp.constraints.idxbx)

    ocp.constraints.lbu = npy.array([model.ddot_min, model.deltadot_min, model.thetadot_min])
    ocp.constraints.ubu = npy.array([model.ddot_max, model.deltadot_max, model.thetadot_max])
    ocp.constraints.idxbu = npy.array([0, 1, 2])


    # set intial condition
    ocp.constraints.x0 = model.x0

    # set QP solver and integration
    ocp.solver_options.tf = Tf
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.nlp_solver_type = "SQP"#"SQP_RTI" #
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "ERK"
    ocp.parameter_values = npy.zeros(np)
    #ocp.solver_options.sim_method_num_stages = 4
    #ocp.solver_options.sim_method_num_steps = 3
    ocp.solver_options.nlp_solver_step_length = 0.05
    ocp.solver_options.nlp_solver_max_iter = 100
    ocp.solver_options.tol = 1e-4
    #ocp.solver_options.print_level = 1
    # ocp.solver_options.nlp_solver_tol_comp = 1e-1


    # create solver
    acados_solver = AcadosOcpSolver(ocp, json_file="acados_ocp2.json")

    print("solver created returning to main")
    return constraints, model, acados_solver, ocp
