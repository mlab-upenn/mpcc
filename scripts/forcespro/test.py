import casadi
import yaml
import numpy as np
import forcespro.nlp

#initial example with kinematic_model
model = forcespro.nlp.SymbolicModel()

N = 20
Tf = 1
Ts = Tf/N

#set dimensions
model.N = N
model.nvar = 10 #stage variables z = [u, x]'
model.neq = 7 #number of equality constraints
model.nh = 1 #number of inequality constraints
model.npar = 13 #

#let z = [x,u] = [ddot,deltadot,thetadot,posx, posy,phi,vx,theta,d,delta]

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
    #lwb = p[12] not needed here

    #extract states
    posx = z[3]
    posy = z[4]
    theta = z[7]
    ddot = z[0]
    deltadot = z[1]
    thetadot = z[2]

    #compute approximate linearized contouring and lag error
    xt_hat = xt + cos_phit * ( theta - theta_hat)
    yt_hat = yt + sin_phit * ( theta - theta_hat)

    e_cont = sin_phit * (xt_hat - posx) - cos_phit *(yt_hat - posy)
    e_lag = cos_phit * (xt_hat - posx) + sin_phit *(yt_hat - posy)

    cost = e_cont * Qc * e_cont + e_lag * Qc * e_lag - Q_theta * thetadot + ddot * R_d * ddot + deltadot * R_delta * deltadot

    return cost

model.objective = lambda z, p: stage_cost(z,p)

def continuous_dynamics(x, u, p):
    '''
    return np.array([
            x[3] * casadi.cos(x[2]),
            x[3] * casadi.sin(x[2]),
            x[3]/p[12] * casadi.tan(x[6]),
            x[5],
            u[2],
            u[0],
            u[2]
            ])
    '''
    phi = x[2]
    vx = x[3]
    theta = x[4]
    d = x[5]
    delta = x[6]

    ddot = u[0]
    deltadot = u[1]
    thetadot = u[2]

    lwb = p[12]

    statedot = np.array([
        vx * casadi.cos(phi),
        vx * casadi.sin(phi),
        vx/lwb * casadi.tan(delta),
        d,
        thetadot,
        ddot,
        deltadot
        ])

    return statedot

#set model to continuous dynamics mode
model.continuous_dynamics = continuous_dynamics

#dynamics only in state Variables
model.E = np.concatenate([np.zeros((7,3)),np.eye(7)], axis=1)

#nonlinear constraints
def nonlinear_ineq(z, p):
    #extract parameters
    xt = p[0]
    yt = p[1]
    r = p[11]

    #extract relevant states
    posx = z[3]
    posy = z[4]

    hval = (xt-posx)**2 + (yt-posy)**2 - r**2

    return hval

model.ineq = lambda z, p: nonlinear_ineq(z, p)
model.hu = np.array([0])
model.hl = np.array([-10000])


#boxconstraints
d_min = -3.0
d_max = 5.0

ddot_min = -10.0
ddot_max = 10.0

delta_min = -0.40  # minimum steering angle [rad]
delta_max = 0.40  # maximum steering angle [rad]

deltadot_min = -2  # minimum steering angle cahgne[rad/s]
deltadot_max = 2 # maximum steering angle cahgne[rad/s]

thetadot_min = -0.1  # minimum adv param speed [m/s]
thetadot_max = 5 # maximum adv param speed [m/s]

theta_min = 0.00  # minimum adv param [m]
theta_max = 100 # maximum adv param  [m]

vx_max = 2 # max long vel [m/s]
vx_min = -1 # min long vel [m/s]

model.ub = np.array([100, 100, 1000, vx_max, theta_max, d_max, delta_max, ddot_max, deltadot_max, thetadot_max ])
model.lb = np.array([-100, -100, -1000, vx_min, theta_min, d_min, delta_min, ddot_min, deltadot_min, thetadot_min ])

#put initial condition on all z variables
model.xinitidx = range(model.nvar)

# Set solver options
codeoptions = forcespro.CodeOptions()
codeoptions.nlp.integrator.type = 'ERK4'
codeoptions.nlp.integrator.Ts = Ts
codeoptions.nlp.integrator.nodes = 5 #intermediate integration nodes

codeoptions.maxit = 200  # Maximum number of iterations
codeoptions.printlevel = 2  # Use printlevel = 2 to print progress (but not for timings)
codeoptions.optlevel = 0  # 0 no optimization, 1 optimize for size, 2 optimize for speed, 3 optimize for size & speed
#codeoptions.nlp.stack_parambounds = 1

# Creates code for symbolic model formulation given above, then contacts server to generate new solver
solver = model.generate_solver(codeoptions)
