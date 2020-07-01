from casadi import *
import


def dynamic_model(Tf,N,track):

    # define casadi struct
    model = types.SimpleNamespace()

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

    #laod track to use for constraints


    #single track model with pajecka tireforces as in  Optimization-Based Autonomous Racing of 1:43 Scale RC Cars Alexander Liniger, Alexander Domahidi and Manfred Morari
    #pose
    posx = MX.sym("posx")
    posy = MX.sym("posy")

    #vel (long and lateral)
    vx = MX.sym("vx")
    vy = MX.sym("vy")

    #body angular Rate
    omega = MX.sym("omega")

    #heading
    phi = MX.sym("phi")

    #steering_angle
    delta = MX.sym("delta")

    #motorinput
    d = MX.sym("d")

    #dynamic forces
    Frx = MX.sym("Frx")
    Fry = MX.sym("Fry")
    Ffx = MX.sym("Ffx")
    Ffy = MX.sym("Ffy")

    #arclength progress
    theta = MX.sym("theta")

    #temporal derivatives
    posxdot = MX.sym("xdot")
    posydot = MX.sym("ydot")
    vxdot = MX.sym("vxdot")
    vydot = MX.sym("vydot")
    phidot = MX.sym("phidot")
    omegadot = MX.sym("omegadot")
    thetadot = MX.sym("thetadot")
    deltadot = MX.sym("deltadot")

    #car state Dynamics
    x = vertcat(
        posxdot,
        posydot,
        phidot,
        vxdot,
        vydot,
        omegadot,
        )

    xdot = vertcat(
        posxdot,
        posydot,
        phidot,
        vxdot,
        vydot,
        omegadot,
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
        )

    model.f_expl_expr = f_expl
    model.f_impl_expr = xdot - f_expl
    model.x = x
    return model
