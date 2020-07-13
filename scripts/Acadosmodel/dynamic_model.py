from casadi import *
import yaml


def dynamic_model(modelparams):

    # define casadi struct
    model = types.SimpleNamespace()
    constraints = types.SimpleNamespace()

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
    Cr = params['Cr']
    Cd = params['Cd']
    Df = params['Df']
    Dr = params['Dr']

    print("CasADi model created with the following parameters: filename: ",modelparams,"\n values:", params)

    #parameter vector
    xt =  SX.sym("xt")
    yt =  SX.sym("yt")
    phit = SX.sym("phit")
    sin_phit = SX.sym("sin_phit")
    cos_phit = SX.sym("cos_phit")
    gt_upper = SX.sym("gt_upper")
    gt_lower = SX.sym("gt_lower")
    #stores linearization point
    theta_hat = SX.sym("theta_hat")
    Qc = SX.sym("Qc")
    Ql = SX.sym("Ql")
    Q_theta = SX.sym("Q_theta")

    #cost on smoothnes of motorinput
    R_d = SX.sym("R_d")

    #cost on smoothness of steering_angle
    R_delta = SX.sym("R_delta")

    p = vertcat(xt, yt, phit, sin_phit, cos_phit, gt_upper, gt_lower, theta_hat, Qc, Ql, Q_theta, R_d, R_delta)


    #single track model with pajecka tireforces as in  Optimization-Based Autonomous Racing of 1:43 Scale RC Cars Alexander Liniger, Alexander Domahidi and Manfred Morari
    #pose
    posx = SX.sym("posx")
    posy = SX.sym("posy")


    #vel (long and lateral)
    vx = SX.sym("vx")
    vy = SX.sym("vy")


    #body angular Rate
    omega = SX.sym("omega")
    #heading
    phi = SX.sym("phi")
    #steering_angle
    delta = SX.sym("delta")
    #motorinput
    d = SX.sym("d")




    #dynamic forces
    Frx = SX.sym("Frx")
    Fry = SX.sym("Fry")
    Ffx = SX.sym("Ffx")
    Ffy = SX.sym("Ffy")

    #arclength progress
    theta = SX.sym("theta")

    #temporal derivatives
    posxdot = SX.sym("xdot")
    posydot = SX.sym("ydot")
    vxdot = SX.sym("vxdot")
    vydot = SX.sym("vydot")
    phidot = SX.sym("phidot")
    omegadot = SX.sym("omegadot")
    deltadot = SX.sym("deltadot")
    thetadot = SX.sym("thetadot")
    ddot = SX.sym("ddot")
    #inputvector
    u = vertcat(ddot, deltadot, thetadot)

    #car state Dynamics
    x = vertcat(
        posx,
        posy,
        phi,
        vx,
        vy,
        omega,
        theta,
        d,
        delta
        )

    xdot = vertcat(
        posxdot,
        posydot,
        phidot,
        vxdot,
        vydot,
        omegadot,
        thetadot,
        ddot,
        deltadot
        )

    #front lateral tireforce
    alphaf = -atan((omega*lf + vy)/vx) + delta
    Ffy = Df*sin(Cf*atan(Bf*alphaf))

    #rear lateral tireforce
    alphar = atan((omega*lr - vy)/vx)
    Fry = Dr*sin(Cf*atan(Br*alphar))

    #rear longitudinal forces
    Frx = (Cm1-Cm2*vx) * d - Cr -Cd*vx**2

    f_expl = vertcat(
        vx*cos(phi) - vy * sin(phi),
        vx*sin(phi) + vy * cos(phi),
        omega,
        1/m * (Frx - Ffy*sin(delta) + m*vy*omega),
        1/m * (Fry + Ffy*cos(delta) - m*vx*omega),
        1/Iz * (Ffy*lf*cos(delta) - Fry*lr),
        thetadot,
        ddot,
        deltadot
        )

    # algebraic variables
    z = vertcat([])

    model.f_expl_expr = f_expl
    model.f_impl_expr = xdot - f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p
    model.z = z
    #boxconstraints
    model.d_min = -3.0
    model.d_max = 3.0

    model.ddot_min = -10.0
    model.ddot_max = 10.0

    model.delta_min = -0.40  # minimum steering angle [rad]
    model.delta_max = 0.40  # maximum steering angle [rad]

    model.deltadot_min = -10  # minimum steering angle cahgne[rad/s]
    model.deltadot_max = 10 # maximum steering angle cahgne[rad/s]

    model.thetadot_min = 0  # minimum adv param speed [m/s]
    model.thetadot_max = 10 # maximum adv param speed [m/s]
    '''
    vars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    xt0 = track_lu_table[0,vars.index('xtrack')]
    yt0 = track_lu_table[0,vars.index('ytrack')]
    phit0 = track_lu_table[0,vars.index('phitrack')]
    theta_hat0 = track_lu_table[0,vars.index('sval')]
    '''
    model.x0 = np.array([0, 0, 0, 1, 0.01, 0, 0, 0, 0])

    #halfspace constraints on x capturing the track at each stage
    n = vertcat(-sin_phit, cos_phit)
    constraints.h_upper = vertcat(0,0)
    g_upper = vertcat(gt_upper, -gt_lower)
    #con_expr <= 0
    model.con_h_expr = vertcat(n[0]*x[0]+n[1]*x[1]-g_upper[0], -n[0]*x[0]-n[1]*x[1]-g_upper[1])

    #compute approximate linearized contouring and lag error
    xt_hat = xt + cos_phit * ( theta - theta_hat)
    yt_hat = yt + sin_phit * ( theta - theta_hat)

    e_cont = sin_phit * (posx - xt_hat) - cos_phit * (posy - yt_hat)
    e_lag = -cos_phit * (posx - xt_hat) - sin_phit * (posy - yt_hat)

    error = vertcat(e_cont, e_lag)
    #set up stage cost
    Q = diag(vertcat(Qc, Ql))
    model.stage_cost = bilin(Q, error, error) - Q_theta * thetadot + bilin(R_d , ddot, ddot) + bilin(R_delta , deltadot, deltadot)

    return model, constraints
