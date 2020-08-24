import casadi
import numpy as np
import yaml

class dynamics_simulator():
    def __init__(self, modelparams, Ts, x0, nodes):

        with open(modelparams) as file:
            params = yaml.load(file, Loader= yaml.FullLoader)

        self.xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
        self.uvars = ['ddot', 'deltadot', 'thetadot']

        m = params['m'] #[kg]
        lf = params['lf'] #[m]
        lr = params['lr'] #[m]
        Iz = params['Iz'] #[kg*m^3]
        lencar = lf+lr
        widthcar = lencar/2

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


        x = casadi.SX.sym("x", len(self.xvars))
        u = casadi.SX.sym("u", len(self.uvars))

        #extract states and inputs
        posx = x[self.xvars.index('posx')]
        posy = x[self.xvars.index('posy')]
        phi = x[self.xvars.index('phi')]
        vx = x[self.xvars.index('vx')]
        vy = x[self.xvars.index('vy')]
        omega = x[self.xvars.index('omega')]
        d = x[self.xvars.index('d')]
        delta = x[self.xvars.index('delta')]
        theta = x[self.xvars.index('theta')]

        ddot = u[self.uvars.index('ddot')]
        deltadot = u[self.uvars.index('deltadot')]
        thetadot = u[self.uvars.index('thetadot')]


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

        #xdot = f(x,u)
        self.f = casadi.Function('f',[x,u], [statedot])

        #state of system
        self.x = x0
        #sampling timestep
        self.Ts = Ts
        #integration nodes
        self.nodes = nodes

    def tick(self, u):
        T_int = self.Ts/self.nodes
        xtemp = self.x
        for idx in range(self.nodes):
            xtemp = self._integrate(T_int, u)
        return self.x

    def set_theta(self, theta):
        self.x[self.xvars.index('theta')] = theta

    def wrap_phi(self):
        if self.x[self.xvars.index('phi')] > 2 * 3.14159:
            self.x[self.xvars.index('phi')] -= 2 * 3.14159
            wrapdir = 1
        elif self.x[self.xvars.index('phi')] < -2 * 3.14159:
            self.x[self.xvars.index('phi')] += 2 * 3.14159
            wrapdir = -1
        else:
             wrapdir = 0
        return wrapdir
    #RK4 integration
    def _integrate(self, Ts, u):

        k1 = self.f(self.x, u).__array__().reshape(-1,)
        k2 = self.f(self.x + Ts/2 * k1, u).__array__().reshape(-1,)
        k3 = self.f(self.x + Ts/2 * k2, u).__array__().reshape(-1,)
        k4 = self.f(self.x + Ts * k3, u).__array__().reshape(-1,)
        self.x = self.x + 1/6 * Ts*(k1 + 2* k2 + 2*k3 + k4).reshape(-1,)
        return self.x
