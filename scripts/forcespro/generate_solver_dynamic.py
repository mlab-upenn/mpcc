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
    Cr = params['Cr']
    Cd = params['Cd']
    Df = params['Df']
    Dr = params['Dr']

    #forces model
    model = forcespro.nlp.SymbolicModel()

    #set dimensions
    model.N = N
    model.nvar = 12 #stage variables z = [u, x]'
    model.neq = 9 #number of equality constraints
    model.nh = 1 #number of inequality constraints
    model.npar = 12 #

    #let z = [u, x] = [vxdot, deltadot, thetadot, posx, posy, phi, vx, delta, theta]
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

        cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + vxdot * R_d * vxdot + deltadot * R_delta * deltadot

        return cost

    def continuous_dynamics(x, u, p):
        #extract
        posx = z[zvars.index('posx')]
        posy = z[zvars.index('posy')]
        phi = z[zvars.index('phi')]
        vx = z[zvars.index('vx')]
        vy = z[zvars.index('vy')]
        omega = z[zvars.index('omega')]
        d = z[zvars.index('d')]
        delta = z[zvars.index('delta')]
        theta = z[zvars.index('theta')]

        #front lateral tireforce
        alphaf = -casadi.atan((omega*lf + vy)/ vx) + delta
        Ffy = Df*casadi.sin(Cf*casadi.atan(Bf*alphaf))

        #rear lateral tireforce
        alphar = casadi.atan((omega*lr - vy)/vx)
        Fry = Dr*casadi.sin(Cf*casadi.atan(Br*alphar))

        #rear longitudinal forces
        Frx = (Cm1-Cm2*vx) * d - Cr -Cd*vx**2

        statedot = np.array([
                vx*cos(phi) - vy * sin(phi), #posxdot
                vx*sin(phi) + vy * cos(phi), #posxdot
                omega,                       #phidot
                1/m * (Frx - Ffy*sin(delta) + m*vy*omega), #vxdot
                1/m * (Fry + Ffy*cos(delta) - m*vx*omega), #vydot
                1/Iz * (Ffy*lf*cos(delta) - Fry*lr),       #omegadot
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

    model.ub = np.array([ddot_max, deltadot_max, thetadot_max, 100, 100, 1000, vx_max, delta_max, theta_max])
    model.lb = np.array([ddot_min, deltadot_min, thetadot_min , -100, -100, -1000, vx_min, delta_min, theta_min])

    #put initial condition on all state variables x
    model.xinitidx = 3 + np.arange(model.nvar -3)
