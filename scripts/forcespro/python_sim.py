import numpy as np
from generate_solver_kinematic import get_forces_solver_kinematic
import forcespro.nlp
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import Bezier
import sys


def main_kin():
    np.set_printoptions(precision=4)
    np.set_printoptions(threshold=sys.maxsize)
    # model parameters
    # model parameters
    paramfile = "modelparams.yaml"
    #load global constant model parameters
    with open(paramfile) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lencar = lf+lr
    #sim parameters
    Tsim = 4
    Tf = 1
    N = 20
    Qc = 0.1
    Ql = 1000
    Q_theta = 10
    R_d = 0.1
    R_delta = 0.1
    Nsim = np.int(np.floor(N/Tf*Tsim))
    r = 0.3 #trackwidth

    solver = get_forces_solver_kinematic(N, Tf, paramfile)

    track_lu_table, smax = Bezier.generatelookuptable("tracks/simpleoval")
    trk_plt = plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()

    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx = 10

    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    uvars = ['ddot', 'deltadot', 'thetadot']
    pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r']

    car_soln = []
    xt0 = track_lu_table[startidx,vars.index('xtrack')]
    yt0 = track_lu_table[startidx,vars.index('ytrack')]
    phit0 = track_lu_table[startidx,vars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,vars.index('sval')]

    #initial condition
    xinit = np.array([xt0, yt0, phit0, 1, 0, theta_hat0])
    zinit = np.concatenate([np.array([0,0,0]),xinit])
    ############################################################################
    #initialization for theta values
    iter = 100
    z_current = np.tile(zinit,(N,1))
    #arbitrarily set theta  values and
    theta_old = theta_hat0*np.ones((N,)) + 0.1*np.arange(N)
    z_current[:,8] = theta_old
    index_lin_points = 100 * theta_old
    index_lin_points = index_lin_points.astype(np.int32)
    track_lin_points = track_lu_table[index_lin_points,:]

    #initialize x values on track
    z_current[:,3] = track_lin_points[:,vars.index('xtrack')]
    z_current[:,4] = track_lin_points[:,vars.index('ytrack')]
    z_current[:,5] = track_lin_points[:,vars.index('phitrack')]

    for idx in range(iter):

        all_parameters = []
        #get track linearization
        index_lin_points = 100 * theta_old
        index_lin_points = index_lin_points.astype(np.int32)
        track_lin_points = track_lu_table[index_lin_points,:]

        for stageidx in range(N):
            p_val = np.array([track_lin_points[stageidx,vars.index('xtrack')],
                                track_lin_points[stageidx,vars.index('ytrack')],
                                track_lin_points[stageidx,vars.index('phitrack')],
                                track_lin_points[stageidx,vars.index('sin(phi)')],
                                track_lin_points[stageidx,vars.index('cos(phi)')],
                                track_lin_points[stageidx,vars.index('sval')],  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-0.05
                                ])
            all_parameters.append(p_val)

        all_parameters = np.array(all_parameters)

        #problem dictionary, arrays have to be flattened
        problem = {"x0": z_current.reshape(-1,),
                   "xinit": xinit,
                   "all_parameters": all_parameters.reshape(-1,)}
        #solve problem
        output, exitflag, info = solver.solve(problem)
        #print(info)
        #print(output)
        #input("hit [enter] to continue.")

        #extract theta values
        idx_sol = 0
        for key in output:
            #print(key)
            zsol = output[key]
            usol = zsol[0:3]
            z_current[idx_sol,:] = zsol
            idx_sol = idx_sol+1

        theta_current = z_current[:,8]

        #compute difference
        theta_diff = np.sum(np.abs(theta_current-theta_old))
        print("theta init difference: ", theta_diff)
        print("theta values", theta_current)
        theta_old = theta_current
    ############################################################################
    #setup using estimated initial trajectory
    theta_vals = theta_current

    step_sol_z_arr = z_current

    print("plotting initialization")
    trk_plt.plot_horizon(theta_vals, step_sol_z_arr[:, 3:6])
    trk_plt.plot_input_state_traj(step_sol_z_arr, zvars)
    plt.pause(0.1)
    input("hit [enter] to continue.")
    plt.pause(0.1)
    trk_plt.clear_horizion()
    trk_plt.clear_input_state_traj()

    #list storing visited states
    zinit_vals = []
    laps = 0
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):
        step_sol_u = []
        all_parameters = []
        theta_old = theta_vals
        #get track linearization
        index_lin_points = 100 * theta_old - 100*laps*smax
        index_lin_points = index_lin_points.astype(np.int32)
        print("track linearized around entries:", index_lin_points )
        track_lin_points = track_lu_table[index_lin_points,:]

        #######################################################################
        #set params and warmstart
        for stageidx in range(N-1):
            p_val = np.array([track_lin_points[stageidx,vars.index('xtrack')],
                                track_lin_points[stageidx,vars.index('ytrack')],
                                track_lin_points[stageidx,vars.index('phitrack')],
                                track_lin_points[stageidx,vars.index('sin(phi)')],
                                track_lin_points[stageidx,vars.index('cos(phi)')],
                                track_lin_points[stageidx,vars.index('sval')] + laps*smax,  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-0.05
                                ])
            #create parameter matrix
            all_parameters.append(p_val)
            #stack state initializations
            z_current[stageidx,:] = step_sol_z_arr[stageidx+1,:]

        #last stage copy old solution for init
        stageidx = N-1
        p_val = np.array([track_lin_points[stageidx,vars.index('xtrack')],
                            track_lin_points[stageidx,vars.index('ytrack')],
                            track_lin_points[stageidx,vars.index('phitrack')],
                            track_lin_points[stageidx,vars.index('sin(phi)')],
                            track_lin_points[stageidx,vars.index('cos(phi)')],
                            track_lin_points[stageidx,vars.index('sval')] + laps*smax,  #aka theta_hat
                            Qc,
                            Ql,
                            Q_theta,
                            R_d,
                            R_delta,
                            r
                            ])
        all_parameters.append(p_val)
        all_parameters = np.array(all_parameters)
        z_current[N-1,:] = step_sol_z_arr[N-1,:]

        #######################################################################
        #problem dictionary, arrays have to be flattened
        problem = {"x0": z_current.reshape(-1,),
                   "xinit": xinit,
                   "all_parameters": all_parameters.reshape(-1,)}
        #solve problem
        output, exitflag, info = solver.solve(problem)

        #extract solution
        idx_sol = 0
        for key in output:
            zsol = output[key]
            usol = zsol[0:3]
            step_sol_u.append(usol)
            step_sol_z_arr[idx_sol,:] = zsol
            idx_sol = idx_sol+1

        theta_vals = step_sol_z_arr[:,8]
        zinit = step_sol_z_arr[1,:]
        xinit = zinit[3:]
        zinit_vals.append(zinit)

        step_sol_u_arr = np.array(step_sol_u)

        objective = compute_objective(Tf/float(N),
                    Qc,
                    Ql,
                    Q_theta,
                    R_d,
                    R_delta,
                    theta_vals,
                    theta_old,
                    step_sol_z_arr[:, 3:5],
                    step_sol_u_arr,
                    track_lin_points[:,vars.index('xtrack'):vars.index('ytrack')+1],
                    track_lin_points[:,vars.index('phitrack')]
                    )


        #plotting result

        print("objective value", objective)
        trk_plt.plot_horizon(theta_vals, step_sol_z_arr[:, 3:6])
        trk_plt.plot_input_state_traj(step_sol_z_arr, zvars)
        #plt.show()

        plt.pause(0.1)
        input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()

        #preparation for next timestep
        theta_vals = np.hstack((step_sol_z_arr[1:, 8], step_sol_z_arr[-1, 8]+0.1))

        if theta_vals[0] > (laps+1)*smax :
            print("#################################RESET###############################")
            laps = laps + 1


    ###############################/SIMULATION##################################
    trk_plt.plot_traj(np.array(x0vals))
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

def main_dyn():
    print("inactive")
    return 0

if __name__ == "__main__":
    arg = sys.argv[1]

    if arg == "kin":
        print("executing with kinematic_model")
        main_kin()
    elif arg == "dyn":
        print("executing with dynamic_model")
        main_dyn()
    else:
        print("Please run model with valid modeltype: [kin, dyn]")
        print("e.g. python3 python_sim.py kin")
