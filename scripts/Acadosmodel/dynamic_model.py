from casadi import *
import Bezier


def dynamic_model(Tf,N):
    '''
    #load track
    waypoints = Bezier.getwaypoints(track)
    #abez,bbez coeffs
    a_bez, b_bez = Bezier.interpolate(waypoints)
    order_inverse = 4
    norm_inverse = 2
    coeffs = Bezier.fit_st(order_inverse, norm_inverse, waypoints)
    '''
    # define casadi struct
    model = types.SimpleNamespace()
    constraints = types.SimpleNamespace()

    model_name = "f110_dynamic_model"
    #loadparameters
    m = 2 #[kg]
    lf = 0.1 #[m]
    lr = 0.1 #[m]
    Iz = 1 #[kg*m^3]

    #pajecka and motor coefficients
    Bf = 1
    Br = 1
    Cf = 1
    Cr = 1
    Cm1 = 1
    Cm2 = 2
    Cr = 1
    Cd = 1
    Df = 1
    Dr = 1


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
    R = SX.sym("R")

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

    #virtual input projected track velocity
    vt = SX.sym("vt")

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
    thetadot = SX.sym("thetadot")
    deltadot = SX.sym("deltadot")

    #car state Dynamics
    x = vertcat(
        posx,
        posy,
        phi,
        vx,
        vy,
        omega,
        theta
        )

    xdot = vertcat(
        posxdot,
        posydot,
        phidot,
        vxdot,
        vydot,
        omegadot,
        thetadot
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
        1/m * (Frx - Fry*sin(delta) + m*vy*omega),
        1/m * (Fry + Fry*cos(delta) - m*vx*omega),
        1/Iz * (Ffy*lf*cos(delta) - Fry*lr),
        vt
        )


    model.f_expl_expr = f_expl
    model.f_impl_expr = xdot - f_expl
    model.x = x

    #halfspace constraints on x capturing the track
    n = vertcat(-sin_phit, cos_phit)
    constraints.C = np.matrix([[n[0], n[1], 0, 0, 0, 0, 0],  [-n[0], -n[1], 0, 0, 0, 0, 0]])
    constraints.ge_upper = vertcat(gt_upper, -gt_lower)


    #compute approximate linearized contouring and lag error
    xt_hat = xt + cos_phit * ( theta - theta_hat)
    yt_hat = yt + sin_phit * ( theta - theta_hat)

    e_cont = sin_phit * (posx - xt_hat) - cos_phit * (posy - yt_hat)
    e_lag = -cos_phit * (posx - xt_hat) - sin_phit * (posy - yt_hat)

    error = vertcat(e_cont, e_lag)
    #set up stage cost
    Q = diag(vertcat(Qc, Ql))
    model.stage_cost = bilin(Q, error, error)

    return model, constraints


if __name__=="__main__":
    model = dynamic_model(0,0,0)
