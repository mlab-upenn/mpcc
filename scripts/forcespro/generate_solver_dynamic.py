import casadi
import yaml
import numpy as np
import forcespro.nlp


def get_forces_solver_dynamic(N, Tf, modelparams = "modelparams.yaml"):
    #load global constant model parameters
    with open(modelparams) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)

    m = params['m'] #[kg]
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    Iz = params['Iz'] #[kg*m^3]

    #pajecka and motor coefficients
    Bf = params['Bf']
    Br = params['Br']
    Cf = params['Cf']
    Cr = params['Cr']
    Cm1 = params['Cm1']
    Cm2 = params['Cm2']
    Croll = params['Croll']
    Cd = params['Cd']
    Df = params['Df']
    Dr = params['Dr']

    #forces model
    model = forcespro.nlp.SymbolicModel()

    #compute sampling time for integration of continuous dynamics
    Ts = Tf/N

    #set dimensions
    model.N = N
    model.nvar = 12 #stage variables z = [u, x]'
    model.neq = 9 #number of equality constraints
    model.nh = 1 #number of inequality constraints
    model.npar = 12 #
    ninputs = 3

    #let z = [u, x] = [ddot, deltadot, thetadot, posx, posy, phi, vx, vy, omega, d, delta, theta]
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r']

    #define objective
    def stage_cost(z, p):
        #extract parameters
        xt = p[pvars.index('xt')]
        yt = p[pvars.index('yt')]
        phit = p[pvars.index('phit')]
        sin_phit = p[pvars.index('sin_phit')]
        cos_phit = p[pvars.index('cos_phit')]
        theta_hat = p[pvars.index('theta_hat')]
        Qc = p[pvars.index('Qc')]
        Ql = p[pvars.index('Ql')]
        Q_theta = p[pvars.index('Q_theta')]
        R_d = p[pvars.index('R_d')]
        R_delta = p[pvars.index('R_delta')]

        #extract states
        posx = z[zvars.index('posx')]
        posy = z[zvars.index('posy')]
        theta = z[zvars.index('theta')]

        #extract inputs
        ddot = z[zvars.index('ddot')]
        deltadot = z[zvars.index('deltadot')]
        thetadot = z[zvars.index('thetadot')]

        #compute approximate linearized contouring and lag error
        xt_hat = xt + cos_phit * ( theta - theta_hat)
        yt_hat = yt + sin_phit * ( theta - theta_hat)

        e_cont = sin_phit * (xt_hat - posx) - cos_phit *(yt_hat - posy)
        e_lag = cos_phit * (xt_hat - posx) + sin_phit *(yt_hat - posy)

        cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + ddot * R_d * ddot + deltadot * R_delta * deltadot

        return cost

    model.objective = lambda z, p: stage_cost(z,p)

    def continuous_dynamics(x, u, p):
        #extract states and inputs
        posx = x[zvars.index('posx')-ninputs]
        posy = x[zvars.index('posy')-ninputs]
        phi = x[zvars.index('phi')-ninputs]
        vx = x[zvars.index('vx')-ninputs]
        vy = x[zvars.index('vy')-ninputs]
        omega = x[zvars.index('omega')-ninputs]
        d = x[zvars.index('d')-ninputs]
        delta = x[zvars.index('delta')-ninputs]
        theta = x[zvars.index('theta')-ninputs]

        ddot = u[zvars.index('ddot')]
        deltadot = u[zvars.index('deltadot')]
        thetadot = u[zvars.index('thetadot')]

        #build CasADi expressions for dynamic model
        #front lateral tireforce
        alphaf = -casadi.atan2((omega*lf + vy), vx) + delta
        Ffy = Df*casadi.sin(Cf*casadi.atan(Bf*alphaf))

        #rear lateral tireforce
        alphar = casadi.atan2((omega*lr - vy),vx)
        Fry = Dr*casadi.sin(Cr*casadi.atan(Br*alphar))

        #rear longitudinal forces
        Frx = (Cm1-Cm2*vx) * d - Croll -Cd*vx*vx

        #let z = [u, x] = [ddot, deltadot, thetadot, posx, posy, phi, vx, vy, omega, d, delta, theta]

        statedot = np.array([
                vx * casadi.cos(phi) - vy * casadi.sin(phi),        #posxdot
                vx * casadi.sin(phi) + vy * casadi.cos(phi),        #posydot
                omega,                                              #phidot
                1/m * (Frx - Ffy*casadi.sin(delta) + m*vy*omega),   #vxdot
                1/m * (Fry + Ffy*casadi.cos(delta) - m*vx*omega),   #vydot
                1/Iz * (Ffy*lf*casadi.cos(delta) - Fry*lr),         #omegadot
                ddot,
                deltadot,
                thetadot
                ])
        return statedot

    #set model to continuous dynamics mode
    model.continuous_dynamics = continuous_dynamics

    #dynamics only in state Variables
    model.E = np.concatenate([np.zeros((9,3)),np.eye(9)], axis=1)

    #nonlinear constraints
    def nonlinear_ineq(z, p):
        #extract parameters
        xt = p[pvars.index('xt')]
        yt = p[pvars.index('yt')]
        phit = p[pvars.index('phit')]
        sin_phit = p[pvars.index('sin_phit')]
        cos_phit = p[pvars.index('cos_phit')]
        theta_hat = p[pvars.index('theta_hat')]
        r = p[pvars.index('r')]

        #extract relevant states
        posx = z[zvars.index('posx')]
        posy = z[zvars.index('posy')]
        theta = z[zvars.index('theta')]

        #compute approximate linearized contouring and lag error
        xt_hat = xt + cos_phit * ( theta - theta_hat)
        yt_hat = yt + sin_phit * ( theta - theta_hat)

        hval = (xt_hat-posx)**2 + (yt_hat-posy)**2 - r**2
        return hval

    model.ineq = lambda z, p: nonlinear_ineq(z, p)
    model.hu = np.array([0.0000])
    model.hl = np.array([-10])

    #boxconstraints
    ddot_min = -10.0 #min change in d [-]
    ddot_max = 10.0  #max change in d [-]

    d_min = -0.1 #min d [-]
    d_max = 1 #max d [-]

    delta_min = -0.40  # minimum steering angle [rad]
    delta_max = 0.40  # maximum steering angle [rad]

    deltadot_min = -2  # minimum steering angle cahgne[rad/s]
    deltadot_max = 2 # maximum steering angle cahgne[rad/s]

    omega_min = -100 # minimum yawrate [rad/sec]
    omega_max = 100 # maximum yawrate [rad/sec]

    thetadot_min = 0.05  # minimum adv param speed [m/s]
    thetadot_max = 5 # maximum adv param speed [m/s]

    theta_min = 0.00  # minimum adv param [m]
    theta_max = 1000 # maximum adv param  [m]

    vx_max = 3.5 # max long vel [m/s]
    vx_min = -0.5 #0.05 # min long vel [m/s]

    vy_max = 3 # max lat vel [m/s]
    vy_min = -3 # min lat vel [m/s]

    #Note: z = [u, x] = [vxdot, deltadot, thetadot, posx, posy, phi, vx, vy, omega, d, delta, theta]
    model.ub = np.array([ddot_max, deltadot_max, thetadot_max, 10, 10, 100, vx_max, vy_max, omega_max, d_max, delta_max, theta_max])
    model.lb = np.array([ddot_min, deltadot_min, thetadot_min , -10, -10, -100, vx_min, vy_min, omega_min, d_min, delta_min, theta_min])

    #put initial condition on all state variables x
    model.xinitidx = 3 + np.arange(model.nvar -3)
    # Set solver options
    codeoptions = forcespro.CodeOptions("dynamic_solver")
    codeoptions.nlp.integrator.type = 'ERK4'
    codeoptions.nlp.integrator.Ts = Ts
    codeoptions.nlp.integrator.nodes = 3 #intermediate integration nodes

    codeoptions.maxit = 15  # Maximum number of iterations
    codeoptions.printlevel = 2  # Use printlevel = 2 to print progress (but not for timings)
    codeoptions.optlevel = 0  # 0 no optimization, 1 optimize for size, 2 optimize for speed, 3 optimize for size & speed
    codeoptions.nlp.stack_parambounds = 2
    #codeoptions.noVariableElimination = True
    # Creates code for symbolic model formulation given above, then contacts server to generate new solver
    solver = model.generate_solver(codeoptions)
    return solver

if __name__ == "__main__":
    N = 20
    Tf = 1

    solver = get_forces_solver_dynamic(N, Tf)
