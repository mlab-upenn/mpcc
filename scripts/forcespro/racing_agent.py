'''
// MIT License

// Copyright (c) 2020 Peter Werner

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
'''

import numpy as np
from generate_col_avoid_solver import get_col_avoid_solver
import forcespro.nlp
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import Bezier
import yaml
import sys


class racer():
    def __init__(self, name = "default_racer", modelparams = "modelparams.yaml", solverparams = "solverparams", track, Tsim):

        self.name = name
        #load global constant model parameters
        with open(modelparams) as file:
            params = yaml.load(file, Loader= yaml.FullLoader)
        self.lf = params['lf'] #[m]
        self.lr = params['lr'] #[m]
        #pajecka and motor coefficients
        self.Bf = params['Bf']
        self.Br = params['Br']
        self.Cf = params['Cf']
        self.Cr = params['Cr']
        self.Cm1 = params['Cm1']
        self.Cm2 = params['Cm2']
        self.Croll = params['Croll']
        self.Cd = params['Cd']
        self.Df = params['Df']
        self.Dr = params['Dr']
        self.lencar = (lf+lr)
        self.widthcar = lencar/2

        #solver params
        with open(solverparams) as file:
            params = yaml.load(file, Loader= yaml.FullLoader)

        self.Tf = params['Tf']
        self.N = params['N']
        self.Nsim = np.int(np.floor(self.N/self.Tf*Tsim))

        self.Qc = params['Qc']
        self.Ql = params['Ql']
        self.Q_theta = params['Q_theta']
        self.R_d = params['R_d']
        self.R_delta = params['R_delta']

        self.r = track['r']
        self.smax = track['smax']
        self.track_lu_table = blah

        self.solver = get_col_avoid_solver(self.N, self.Tf, modelparams)

        self.trackvars =['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
        self.xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
        self.uvars = ['ddot', 'deltadot', 'thetadot']
        self.pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r', 'x_ob', 'y_ob', 'phi_ob', 'l_ob', 'w_ob']
        self.zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']

        self.z_current = np.zeros((N, len(self.zvars) ))
        self.theta_current = np.zeros((N,))

        #list to store all visited states
        self.zinit_vals = np.zeros(Nsim, len(self.zvars))
        #list containing also the prediction horizons
        self.z_data = np.zeros(Nsim, N, len(self.zvars))
        #sim step trackign
        self.simidx = 0
        self.laps = 0

    def initialize_trajectory(self, xinit, enemyinfo, startidx):
        x_ob = enemyinfo['x_ob']
        y_ob = enemyinfo['y_ob']
        phi_ob = enemyinfo['phi_ob']
        l_ob = enemyinfo['l_ob']
        w_ob = enemyinfo['w_ob']

        #initialization for theta values
        iter = 15

        self.zinit = np.concatenate([np.array([0,0,0]), xinit])
        self.z_current = np.tile(zinit,(N,1))

        #arbitrarily set theta  values and
        theta_old = zinit[self.zvars.index('theta')]*np.ones((N,)) + 0.01*np.arange(N)
        self.z_current[:,self.zvars.index('theta')] = theta_old
        index_lin_points = 100 * theta_old
        index_lin_points = index_lin_points.astype(np.int32)
        track_lin_points = track_lu_table[index_lin_points,:]

        #initialize x values on track
        self.z_current[:,3] = track_lin_points[:,self.trackvars.index('xtrack')]
        self.z_current[:,4] = track_lin_points[:,self.trackvars.index('ytrack')]
        self.z_current[:,5] = track_lin_points[:,self.trackvars.index('phitrack')]

        for idx in range(iter):
            print("theta values", theta_old)
            all_parameters = []
            #get track linearization
            index_lin_points = 100 * theta_old
            index_lin_points = index_lin_points.astype(np.int32)
            track_lin_points = track_lu_table[index_lin_points,:]

            for stageidx in range(self.N):
                p_val = np.array([track_lin_points[stageidx,self.trackvars.index('xtrack')],
                                    track_lin_points[stageidx,self.trackvars.index('ytrack')],
                                    track_lin_points[stageidx,self.trackvars.index('phitrack')],
                                    track_lin_points[stageidx,self.trackvars.index('sin(phi)')],
                                    track_lin_points[stageidx,self.trackvars.index('cos(phi)')],
                                    track_lin_points[stageidx,self.trackvars.index('sval')],  #aka theta_hat
                                    Qc,
                                    Ql,
                                    Q_theta,
                                    R_d,
                                    R_delta,
                                    r,
                                    x_ob,
                                    y_ob,
                                    phi_ob,
                                    l_ob,
                                    w_ob
                                    ])
                all_parameters.append(p_val)

            all_parameters = np.array(all_parameters)

            #problem dictionary, arrays have to be flattened
            problem = {"x0": self.z_current.reshape(-1,),
                       "xinit": xinit,
                       "all_parameters": all_parameters.reshape(-1,)}
            #solve problem
            output, exitflag, info = self.solver.solve(problem)
            #print(info)
            #print(exitflag)
            #input("hit [enter] to continue.")

            #extract theta values
            idx_sol = 0
            for key in output:
                #print(key)
                zsol = output[key]
                usol = zsol[0:3]
                self.z_current[idx_sol,:] = zsol
                idx_sol = idx_sol+1

            self.theta_current = z_current[:, zvars.index('theta')]

            #compute difference
            theta_diff = np.sum(np.abs(theta_current-theta_old))
            print(self.name +": theta init difference: ", theta_diff)
            print("theta values", self.theta_current)
            theta_old = self.theta_current

        return self.z_current

    def update(self, enemyinfo):

        x_ob = enemyinfo['x_ob']
        y_ob = enemyinfo['y_ob']
        phi_ob = enemyinfo['phi_ob']
        l_ob = enemyinfo['l_ob']
        w_ob = enemyinfo['w_ob']

        all_parameters = []

        theta_old = self.theta_current
        #get track linearization
        index_lin_points = 100 * theta_old
        index_lin_points = index_lin_points.astype(np.int32)
        print("track linearized around entries:", index_lin_points )
        track_lin_points = self.track_lu_table[index_lin_points,:]

        #######################################################################
        #set params and warmstart
        for stageidx in range(self.N-1):
            p_val = np.array([track_lin_points[stageidx,self.trackvars.index('xtrack')],
                                track_lin_points[stageidx,self.trackvars.index('ytrack')],
                                track_lin_points[stageidx,self.trackvars.index('phitrack')],
                                track_lin_points[stageidx,self.trackvars.index('sin(phi)')],
                                track_lin_points[stageidx,self.trackvars.index('cos(phi)')],
                                track_lin_points[stageidx,self.trackvars.index('sval')],  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r,
                                x_ob,
                                y_ob,
                                phi_ob,
                                l_ob,
                                w_ob
                                ])
            #create parameter matrix
            all_parameters.append(p_val)
            #stack state initializations, last state is copied
            self.z_current[stageidx,:] = self.z_current[stageidx+1,:]

        #last stage copy old solution for init
        stageidx = self.N-1
        p_val = np.array([track_lin_points[stageidx,self.trackvars.index('xtrack')],
                            track_lin_points[stageidx,self.trackvars.index('ytrack')],
                            track_lin_points[stageidx,self.trackvars.index('phitrack')],
                            track_lin_points[stageidx,self.trackvars.index('sin(phi)')],
                            track_lin_points[stageidx,self.trackvars.index('cos(phi)')],
                            track_lin_points[stageidx,self.trackvars.index('sval')],  #aka theta_hat
                            Qc,
                            Ql,
                            Q_theta,
                            R_d,
                            R_delta,
                            r,
                            x_ob,
                            y_ob,
                            phi_ob,
                            l_ob,
                            w_ob
                            ])
        all_parameters.append(p_val)
        all_parameters = np.array(all_parameters)
        #last state of z_current is already copied.

        #######################################################################
        #problem dictionary, arrays have to be flattened
        problem = {"x0": self.z_current.reshape(-1,),
                   "xinit": self.xinit,
                   "all_parameters": all_parameters.reshape(-1,)}
        #solve problem
        output, exitflag, info = self.solver.solve(problem)

        #extract solution
        idx_sol = 0
        for key in output:
            zsol = output[key]
            self.z_current[idx_sol, :] = zsol
            idx_sol = idx_sol+1

        #log solution
        self.z_data[simidx,:,:] = self.z_current
        '''
        print("theta: ", zinit[self.zvars.index('theta')])
        print("vx: ", zinit[self.zvars.index('vx')])
        print("vy: ", zinit[self.zvars.index('vy')])
        print("omega: ", zinit[self.zvars.index('omega')])
        print("phi: ", zinit[self.zvars.index('phi')]*180/3.1415)
        print("d: ", zinit[self.zvars.index('d')])
        print("delta: ", zinit[self.zvars.index('delta')])

        vx = zinit[self.zvars.index('vx')]
        vy = zinit[self.zvars.index('vy')]
        delta = zinit[self.zvars.index('delta')]
        omega = zinit[self.zvars.index('omega')]
        d = zinit[self.zvars.index('d')]
        alphaf = -np.arctan2((omega*lf + vy), vx) + delta
        Ffy = Df*np.sin(Cf*np.arctan(Bf*alphaf))
        alphar = np.arctan2((omega*lr - vy),vx)
        Fry = Dr*np.sin(Cr*np.arctan(Br*alphar))
        Frx = (Cm1-Cm2*vx) * d - Cr -Cd*vx*vx

        print("Ffy: ", Ffy)
        print("Fry: ", Fry)
        print("Frx: ", Frx)
        '''
        self.theta_current = self.z_current[:,self.zvars.index('theta')]
        self.zinit = self.z_current[0,:]
        self.xinit = self.zinit[3:]
        self.zinit_vals[simidx,:] = self.zinit

        '''
        #plotting result
        trk_plt.plot_horizon(theta_vals, step_sol_z_arr[:, 3:6])
        trk_plt.plot_input_state_traj(step_sol_z_arr, zvars)
        #plt.show()

        plt.pause(0.1)
        input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()
        '''
        #preparation for next timestep
        self.theta_current = np.hstack((self.z_current[1:, self.zvars.index('theta')], self.z_current[-1, self.zvars.index('theta')]+0.1))

        if self.theta_current[0] > self.smax :
            print("#################################RESET###############################")
            laps = laps + 1
            print("lap:", laps)
            theta_vals = theta_vals - smax
            self.z_current[:,self.zvars.index('theta')] = theta_vals

        self.simidx = self.simidx + 1
        return self.z_current

    def return_sim_data():
        return self.zinit_vals, self.z_data
