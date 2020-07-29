import casadi
import yaml
import numpy as np

def dynamic_model(modelparams):

    # define casadi struct
    model = casadi.types.SimpleNamespace()
    constraints = casadi.types.SimpleNamespace()

    model_name = "f110_dynamic_model"
    model.name = model_name


    #loadparameters
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

    print("CasADi model created with the following parameters: filename: ",modelparams,"\n values:", params)

    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    uvars = ['ddot', 'deltadot', 'thetadot']
    pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r']
    #parameter vector, contains linearization poitns
    xt =  casadi.SX.sym("xt")
    yt =  casadi.SX.sym("yt")
    phit = casadi.SX.sym("phit")
    sin_phit = casadi.SX.sym("sin_phit")
    cos_phit = casadi.SX.sym("cos_phit")
    #stores linearization point
    theta_hat = casadi.SX.sym("theta_hat")
    Qc = casadi.SX.sym("Qc")
    Ql = casadi.SX.sym("Ql")
    Q_theta = casadi.SX.sym("Q_theta")

    #cost on smoothnes of motorinput
    R_d = casadi.SX.sym("R_d")
    #cost on smoothness of steering_angle
    R_delta = casadi.SX.sym("R_delta")
    #trackwidth
    r = casadi.SX.sym("r")

    #pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r']
    p = casadi.vertcat(xt, yt, phit, sin_phit, cos_phit, theta_hat, Qc, Ql, Q_theta, R_d, R_delta, r)


    #single track model with pajecka tireforces as in  Optimization-Based Autonomous Racing of 1:43 Scale RC Cars Alexander Liniger, Alexander Domahidi and Manfred Morari
    #pose
    posx = casadi.SX.sym("posx")
    posy = casadi.SX.sym("posy")


    #vel (long and lateral)
    vx = casadi.SX.sym("vx")
    vy = casadi.SX.sym("vy")

    #body angular Rate
    omega = casadi.SX.sym("omega")
    #heading
    phi = casadi.SX.sym("phi")
    #steering_angle
    delta = casadi.SX.sym("delta")
    #motorinput
    d = casadi.SX.sym("d")

    #dynamic forces
    Frx = casadi.SX.sym("Frx")
    Fry = casadi.SX.sym("Fry")
    Ffx = casadi.SX.sym("Ffx")
    Ffy = casadi.SX.sym("Ffy")

    #arclength progress
    theta = casadi.SX.sym("theta")

    #temporal derivatives
    posxdot = casadi.SX.sym("xdot")
    posydot = casadi.SX.sym("ydot")
    vxdot = casadi.SX.sym("vxdot")
    vydot = casadi.SX.sym("vydot")
    phidot = casadi.SX.sym("phidot")
    omegadot = casadi.SX.sym("omegadot")
    deltadot = casadi.SX.sym("deltadot")
    thetadot = casadi.SX.sym("thetadot")
    ddot = casadi.SX.sym("ddot")

    #inputvector
    #uvars = ['ddot', 'deltadot', 'thetadot']
    u = casadi.vertcat(ddot, deltadot, thetadot)

    #car state Dynamics
    #xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    x = casadi.vertcat(
        posx,
        posy,
        phi,
        vx,
        vy,
        omega,
        d,
        delta,
        theta
        )

    xdot = casadi.vertcat(
        posxdot,
        posydot,
        phidot,
        vxdot,
        vydot,
        omegadot,
        ddot,
        deltadot,
        thetadot
        )

    #build CasADi expressions for dynamic model
    #front lateral tireforce
    alphaf = -casadi.atan2((omega*lf + vy), vx) + delta
    Ffy = Df*casadi.sin(Cf*casadi.atan(Bf*alphaf))

    #rear lateral tireforce
    alphar = casadi.atan2((omega*lr - vy),vx)

    Fry = Dr*casadi.sin(Cr*casadi.atan(Br*alphar))

    #rear longitudinal forces

    Frx = (Cm1-Cm2*vx) * d - Croll -Cd*vx*vx

    f_expl = casadi.vertcat(
        vx * casadi.cos(phi) - vy * casadi.sin(phi),        #posxdot
        vx * casadi.sin(phi) + vy * casadi.cos(phi),        #posydot
        omega,                                              #phidot
        1/m * (Frx - Ffy*casadi.sin(delta) + m*vy*omega),   #vxdot
        1/m * (Fry + Ffy*casadi.cos(delta) - m*vx*omega),   #vydot
        1/Iz * (Ffy*lf*casadi.cos(delta) - Fry*lr),         #omegadot
        ddot,
        deltadot,
        thetadot
        )

    # algebraic variables
    z = casadi.vertcat([])

    model.f_expl_expr = f_expl
    model.f_impl_expr = xdot - f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p
    model.z = z

    #boxconstraints
    model.ddot_min = -10.0 #min change in d [-]
    model.ddot_max = 10.0  #max change in d [-]

    model.d_min = -0.1 #min d [-]
    model.d_max = 1 #max d [-]

    model.delta_min = -0.40  # minimum steering angle [rad]
    model.delta_max = 0.40  # maximum steering angle [rad]

    model.deltadot_min = -2  # minimum steering angle cahgne[rad/s]
    model.deltadot_max = 2 # maximum steering angle cahgne[rad/s]

    model.omega_min = -100 # minimum yawrate [rad/sec]
    model.omega_max = 100 # maximum yawrate [rad/sec]

    model.thetadot_min = 0.05  # minimum adv param speed [m/s]
    model.thetadot_max = 5 # maximum adv param speed [m/s]

    model.theta_min = 0.00  # minimum adv param [m]
    model.theta_max = 1000 # maximum adv param  [m]

    model.vx_max = 3.5 # max long vel [m/s]
    model.vx_min = -0.5 #0.05 # min long vel [m/s]

    model.vy_max = 3 # max lat vel [m/s]
    model.vy_min = -3 # min lat vel [m/s]


    model.x0 = np.array([0, 0, 0, 1, 0.01, 0, 0, 0, 0])


    #compute approximate linearized contouring and lag error
    xt_hat = xt + cos_phit * ( theta - theta_hat)
    yt_hat = yt + sin_phit * ( theta - theta_hat)

    e_cont = sin_phit * (xt_hat - posx) - cos_phit *(yt_hat - posy)
    e_lag = cos_phit * (xt_hat - posx) + sin_phit *(yt_hat - posy)

    cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + ddot * R_d * ddot + deltadot * R_delta * deltadot

    #error = casadi.vertcat(e_cont, e_lag)
    #set up stage cost
    #Q = diag(vertcat(Qc, Ql))
    model.con_h_expr = (xt_hat-posx)**2 + (yt_hat-posy)**2 - r**2
    model.stage_cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + ddot * R_d * ddot + deltadot * R_delta * deltadot

    return model, constraints

def kinematic_model(modelparams):
    # define casadi struct
    # define casadi struct
    model = casadi.types.SimpleNamespace()
    constraints = casadi.types.SimpleNamespace()

    model_name = "f110_kinematic_model"
    model.name = model_name


    #loadparameters
    with open(modelparams) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)


    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lwb = lf+lr

    print("CasADi model created with the following parameters: filename: ",modelparams,"\n values:", params)

    #parameter vector, contains linearization poitns
    xt =  casadi.SX.sym("xt")
    yt =  casadi.SX.sym("yt")
    phit = casadi.SX.sym("phit")
    sin_phit = casadi.SX.sym("sin_phit")
    cos_phit = casadi.SX.sym("cos_phit")
    gt_upper = casadi.SX.sym("gt_upper")
    gt_lower = casadi.SX.sym("gt_lower")
    #stores linearization point
    theta_hat = casadi.SX.sym("theta_hat")
    Qc = casadi.SX.sym("Qc")
    Ql = casadi.SX.sym("Ql")
    Q_theta = casadi.SX.sym("Q_theta")

    #cost on smoothnes of motorinput
    R_d = casadi.SX.sym("R_d")

    #cost on smoothness of steering_angle
    R_delta = casadi.SX.sym("R_delta")

    r = casadi.SX.sym("r")
    p = casadi.vertcat(xt, yt, phit, sin_phit, cos_phit, theta_hat, Qc, Ql, Q_theta, R_d, R_delta, r)


    #single track model
    #pose
    posx = casadi.SX.sym("posx")
    posy = casadi.SX.sym("posy")

    #longitudinal velocity
    vx = casadi.SX.sym("vx")
    #heading
    phi = casadi.SX.sym("phi")
    #steering_angle
    delta = casadi.SX.sym("delta")
    #motorinput
    d = casadi.SX.sym("d")
    #arclength progress
    theta = casadi.SX.sym("theta")
    #temporal derivatives
    posxdot = casadi.SX.sym("xdot")
    posydot = casadi.SX.sym("ydot")
    vxdot = casadi.SX.sym("vxdot")
    phidot = casadi.SX.sym("phidot")
    deltadot = casadi.SX.sym("deltadot")
    thetadot = casadi.SX.sym("thetadot")
    ddot = casadi.SX.sym("ddot")

    #inputvector
    u = casadi.vertcat(ddot, deltadot, thetadot)

    #car state Dynamics
    x = casadi.vertcat(
        posx,
        posy,
        phi,
        vx,
        theta,
        d,
        delta
        )

    xdot = casadi.vertcat(
        posxdot,
        posydot,
        phidot,
        vxdot,
        thetadot,
        ddot,
        deltadot
        )

    f_expl = casadi.vertcat(
        vx*casadi.cos(phi),
        vx*casadi.sin(phi),
        vx/lwb * casadi.tan(delta),
        d,
        thetadot,
        ddot,
        deltadot
        )

    # algebraic variables
    z = casadi.vertcat([])

    model.f_expl_expr = f_expl
    model.f_impl_expr = xdot - f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p
    model.z = z
    #boxconstraints
    model.d_min = -3.0
    model.d_max = 5.0

    model.ddot_min = -10.0
    model.ddot_max = 10.0

    model.delta_min = -0.40  # minimum steering angle [rad]
    model.delta_max = 0.40  # maximum steering angle [rad]

    model.deltadot_min = -2  # minimum steering angle cahgne[rad/s]
    model.deltadot_max = 2 # maximum steering angle cahgne[rad/s]

    model.thetadot_min = -0.1  # minimum adv param speed [m/s]
    model.thetadot_max = 5 # maximum adv param speed [m/s]

    model.theta_min = 0.00  # minimum adv param [m]
    model.theta_max = 100 # maximum adv param  [m]

    model.vx_max = 2 # max long vel [m/s]
    model.vx_min = -1 # min long vel [m/s]

    model.x0 = np.array([0, 0, 0, 1, 0, 0, 0])

    #halfspace constraints on x capturing the track at each stage
    #n = casadi.vertcat(-sin_phit, cos_phit)
    #constraints.h_upper = casadi.vertcat(0,0)
    #g_upper = casadi.vertcat(gt_upper, -gt_lower)
    #halfspace constriants for track boundaries, con_expr <= 0
    #model.con_h_expr = casadi.vertcat(n[0]*x[0]+n[1]*x[1]-g_upper[0], -n[0]*x[0]-n[1]*x[1]-g_upper[1])

    #compute approximate linearized contouring and lag error
    xt_hat = xt + cos_phit * ( theta - theta_hat)
    yt_hat = yt + sin_phit * ( theta - theta_hat)

    e_cont = sin_phit * (xt_hat - posx) - cos_phit *(yt_hat - posy)
    e_lag = cos_phit * (xt_hat - posx) + sin_phit *(yt_hat - posy)

    model.con_h_expr = e_cont**2+e_lag**2-(r)**2
    #virtual turn radius
    #error = casadi.vertcat(e_cont, e_lag)
    #set up stage cost
    #Q = diag(vertcat(Qc, Ql))
    model.stage_cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + ddot * R_d * ddot + deltadot * R_delta * deltadot
    return model, constraints
