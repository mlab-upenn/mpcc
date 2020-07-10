import numpy as np
from acados_settings import *
import Bezier



def main():
    """ Main program """
    Tsim = 0.1
    Tf = 1
    N = 50
    Qc = 0.1
    Ql = 10
    Q_theta = 10
    R_d = 0.01
    R_delta = 0.01

    Nsim = np.int(np.floor(N/Tf*Tsim))

    track_lu_table = Bezier.generatelookuptable("my_fav_track")
    constraints, model, acados_solver, ocp = acados_settings(Tf, N, track_lu_table)

    vars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    smax = track_lu_table[-1,vars.index('sval')]
    car_soln = []
    xt0 = track_lu_table[0,vars.index('xtrack')]
    yt0 = track_lu_table[0,vars.index('ytrack')]
    phit0 = track_lu_table[0,vars.index('phitrack')]
    theta_hat0 = track_lu_table[0,vars.index('sval')]
    x0 = np.array([xt0, yt0, phit0, 0, 0, 0, theta_hat0, 0, 0])

    initial_theta_spacing = 0.05
    index_lin_points = 100 * np.arange(0,N*initial_theta_spacing,initial_theta_spacing)
    index_lin_points = index_lin_points.astype(np.int32)
    track_lin_points = track_lu_table[index_lin_points,:]

    for simidx in range(Nsim):

        #set params
        for stageidx in range(N):
            p_val = np.array([track_lin_points[stageidx,vars.index('xtrack')],
                                track_lin_points[stageidx,vars.index('ytrack')],
                                track_lin_points[stageidx,vars.index('phitrack')],
                                track_lin_points[stageidx,vars.index('sin(phi)')],
                                track_lin_points[stageidx,vars.index('cos(phi)')],
                                track_lin_points[stageidx,vars.index('g_upper')],
                                track_lin_points[stageidx,vars.index('g_lower')],
                                track_lin_points[stageidx,vars.index('sval')],
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
        acados_solver.print_statistics()
        print(status)
        print("xsol : ", acados_solver.get(0,"x"))
        print("usol : ", acados_solver.get(0,"u"))

        #print(simidx)
    return 0

if __name__ == "__main__":
    main()
