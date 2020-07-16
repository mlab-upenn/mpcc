import numpy as np
from acados_settings import *
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import Bezier



def main():
    # model parameters
    paramfile = "modelparams.yaml"
    #sim parameters
    Tsim = 2
    Tf = 1
    N = 20
    Qc = 100
    Ql = 10000
    Q_theta = 10
    R_d = 0.1
    R_delta = 0.1

    Nsim = np.int(np.floor(N/Tf*Tsim))

    track_lu_table, smax = Bezier.generatelookuptable("tracks/simpleoval")
    trk_plt = plotter(track_lu_table, smax)
    trk_plt.plot_track()

    constraints, model, acados_solver, ocp = acados_settings(Tf, N, track_lu_table, paramfile)
    #plot_pajecka(paramfile)

    plt.show(block=False)
    plt.pause(0.001)
    input("hit [enter] to continue.")

    startidx = 0
    vars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    car_soln = []
    xt0 = track_lu_table[startidx,vars.index('xtrack')]
    yt0 = track_lu_table[startidx,vars.index('ytrack')]
    phit0 = track_lu_table[startidx,vars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,vars.index('sval')]
    x0 = np.array([xt0, yt0, phit0, 1, 0.01, 0, theta_hat0, 0, 0])

    #initialization for theta values
    initial_theta_spacing = 0.01
    theta_vals = np.arange(0,N*initial_theta_spacing,initial_theta_spacing)

    x0vals = []
    x0vals.append(x0)
    full_sol_x = []
    full_sol_u = []

    for simidx in range(Nsim):
        step_sol_x = []
        step_sol_u = []

        #get track linearization
        index_lin_points = 100 * theta_vals
        index_lin_points = index_lin_points.astype(np.int32)
        print("track linearized around entries:", index_lin_points )
        track_lin_points = track_lu_table[index_lin_points,:]

        #set params
        for stageidx in range(N):
            p_val = np.array([track_lin_points[stageidx,vars.index('xtrack')],
                                track_lin_points[stageidx,vars.index('ytrack')],
                                track_lin_points[stageidx,vars.index('phitrack')],
                                track_lin_points[stageidx,vars.index('sin(phi)')],
                                track_lin_points[stageidx,vars.index('cos(phi)')],
                                track_lin_points[stageidx,vars.index('g_upper')],
                                track_lin_points[stageidx,vars.index('g_lower')],
                                track_lin_points[stageidx,vars.index('sval')],  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta
                                ])
            #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
            acados_solver.set(stageidx,"p", p_val)

        #initialize
        acados_solver.set(0, "lbx", x0)
        acados_solver.set(0, "ubx", x0)

        #warmstart stagewise (implement later)
        #for stageidx in range(N-1):
            #acados_solver
        status = acados_solver.solve()
        #acados_solver.print_statistics()
        #print(status)
        print("usol : ", acados_solver.get(0,"u"))
        print("xsol : ", acados_solver.get(0,"x"))

        x0 = acados_solver.get(1,"x")
        x0vals.append(x0)

        for idx_sol in range(N):
            step_sol_x.append(acados_solver.get(idx_sol,"x"))
            step_sol_u.append(acados_solver.get(idx_sol,"u"))

        step_sol_x_arr = np.array(step_sol_x)
        step_sol_u_arr = np.array(step_sol_u)

        theta_old = theta_vals
        theta_vals_sol = step_sol_x_arr[:, 6]

        objective = compute_objective(Tf/float(N),
                    Qc,
                    Ql,
                    Q_theta,
                    R_d,
                    R_delta,
                    theta_vals_sol,
                    theta_old,
                    step_sol_x_arr[:, :2],
                    step_sol_u_arr,
                    track_lin_points[:,vars.index('xtrack'):vars.index('ytrack')+1],
                    track_lin_points[:,vars.index('phitrack')]
                    )

        print("objective value", objective)
        trk_plt.plot_horizon(theta_vals_sol, step_sol_x_arr[:, :2])
        trk_plt.plot_input_state_traj(step_sol_x_arr, step_sol_u_arr)
        plt.pause(0.1)
        input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()


        #preparation for next timestep
        theta_vals = np.hstack((step_sol_x_arr[1:, 6], step_sol_x_arr[-1, 6]+0.1))
        #print("theta vals", theta_vals)
        if theta_vals[0]>smax:
            theta_vals = theta_vals-smax
            print("#################################RESET###############################")



    trk_plt.plot_traj(np.array(x0vals))
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

if __name__ == "__main__":
    main()
