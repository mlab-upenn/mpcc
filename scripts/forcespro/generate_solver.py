import casadi
import yaml
import numpy as np
import forcespro.nlp


def get_forces_solver(N, Tf, modelparams = "modelparams.yaml"):

    #load constant model parameters
    with open(modelparams) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lwb = lf + lr

    #initial example with kinematic_model
    model = forcespro.nlp.SymbolicModel()

    Ts = Tf/N

    #set dimensions
    model.N = N
    model.nvar = 9 #stage variables z = [u, x]'
    model.neq = 6 #number of equality constraints
    model.nh = 1 #number of inequality constraints
    model.npar = 12 #

    #let z = [u, x] = [vxdot, deltadot, thetadot, posx, posy, phi, vx, delta, theta]

    #define objective
    def stage_cost(z, p):
        #extract parameters
        xt = p[0]
        yt = p[1]
        phit = p[2]
        sin_phit = p[3]
        cos_phit = p[4]
        theta_hat = p[5]
        Qc = p[6]
        Ql = p[7]
        Q_theta = p[8]
        R_d = p[9]
        R_delta = p[10]
        #r = p[11] not needed here

        #extract states
        posx = z[3]
        posy = z[4]
        theta = z[8]

        #extract inputs
        vxdot = z[0]
        deltadot = z[1]
        thetadot = z[2]

        #compute approximate linearized contouring and lag error
        xt_hat = xt + cos_phit * ( theta - theta_hat)
        yt_hat = yt + sin_phit * ( theta - theta_hat)

        e_cont = sin_phit * (xt_hat - posx) - cos_phit *(yt_hat - posy)
        e_lag = cos_phit * (xt_hat - posx) + sin_phit *(yt_hat - posy)

        cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + vxdot * R_d * vxdot + deltadot * R_delta * deltadot

        return cost

    model.objective = lambda z, p: stage_cost(z,p)

    def continuous_dynamics(x, u, p):

        phi = x[2]
        vx = x[3]
        delta = x[4]
        theta = x[5]


        vxdot = u[0]
        deltadot = u[1]
        thetadot = u[2]

        statedot = np.array([
            vx * casadi.cos(phi),       #posxdot
            vx * casadi.sin(phi),       #posydot
            vx/lwb * casadi.tan(delta), #phidot
            vxdot,                          #vxdot
            deltadot,
            thetadot
            ])
        return statedot

    #set model to continuous dynamics mode
    model.continuous_dynamics = continuous_dynamics

    #dynamics only in state Variables
    model.E = np.concatenate([np.zeros((6,3)),np.eye(6)], axis=1)

    #nonlinear constraints
    def nonlinear_ineq(z, p):
        #extract parameters
        #extract parameters
        xt = p[0]
        yt = p[1]
        phit = p[2]
        sin_phit = p[3]
        cos_phit = p[4]
        theta_hat = p[5]
        r = p[11]

        #extract relevant states
        posx = z[3]
        posy = z[4]
        theta = z[8]

        #compute approximate linearized contouring and lag error
        xt_hat = xt + cos_phit * ( theta - theta_hat)
        yt_hat = yt + sin_phit * ( theta - theta_hat)

        hval = (xt_hat-posx)**2 + (yt_hat-posy)**2 - r**2
        return hval

    model.ineq = lambda z, p: nonlinear_ineq(z, p)
    model.hu = np.array([0.01])
    model.hl = np.array([-10])


    #boxconstraints
    vxdot_min = -10.0
    vxdot_max = 10.0

    delta_min = -0.40  # minimum steering angle [rad]
    delta_max = 0.40  # maximum steering angle [rad]

    deltadot_min = -2  # minimum steering angle cahgne[rad/s]
    deltadot_max = 2 # maximum steering angle cahgne[rad/s]

    thetadot_min = 0.1  # minimum adv param speed [m/s]
    thetadot_max = 5 # maximum adv param speed [m/s]

    theta_min = 0.00  # minimum adv param [m]
    theta_max = 100 # maximum adv param  [m]

    vx_max = 2 # max long vel [m/s]
    vx_min = -1 # min long vel [m/s]

    #vars=[vxdot, deltadot, thetadot, posx, posy, phi, vx, delta, theta]

    model.ub = np.array([vxdot_max, deltadot_max, thetadot_max, 100, 100, 1000, vx_max, delta_max, theta_max])
    model.lb = np.array([vxdot_min, deltadot_min, thetadot_min , -100, -100, -1000, vx_min, delta_min, theta_min])

    #put initial condition on all state variables x
    model.xinitidx = 3 + np.arange(model.nvar -3)

    # Set solver options
    codeoptions = forcespro.CodeOptions()
    codeoptions.nlp.integrator.type = 'ERK4'
    codeoptions.nlp.integrator.Ts = Ts
    codeoptions.nlp.integrator.nodes = 3 #intermediate integration nodes

    codeoptions.maxit = 40  # Maximum number of iterations
    codeoptions.printlevel = 2  # Use printlevel = 2 to print progress (but not for timings)
    codeoptions.optlevel = 1  # 0 no optimization, 1 optimize for size, 2 optimize for speed, 3 optimize for size & speed
    codeoptions.nlp.stack_parambounds = 1
    #codeoptions.noVariableElimination = True
    # Creates code for symbolic model formulation given above, then contacts server to generate new solver
    solver = model.generate_solver(codeoptions)
    return solver

if __name__ == "__main__":
    N = 20
    Tf = 1

    solver = get_forces_solver(N, Tf)
