import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

def interpolate(waypoints):
    #interpolates with cubic bezier curves with cyclic boundary condition
    n = len(waypoints)
    M = np.zeros([n,n])

    #build M
    tridiagel = np.matrix([[1, 4, 1]])
    for idx in range(n-2):
        M[idx+1:idx+2, idx:idx+3] = tridiagel

    M[0,0:2]= tridiagel[:,1:3]
    M[-1,-2:]= tridiagel[:,0:2]
    M[0:2,-1] = tridiagel[:,0].reshape(1,-1)
    M[-1,0] = tridiagel[:,0].reshape(1,-1)


    #build sol vector
    s =np.zeros([n,2])
    for idx in range(n-1):
        s[idx,:] = 2*(2*waypoints[idx,:] + waypoints[idx+1,:])
    s[-1:] = 2*(2*waypoints[-1,:] + waypoints[0,:])

    #solve for a & b
    Ax = np.linalg.solve(M,s[:,0])
    Ay = np.linalg.solve(M,s[:,1])

    a = np.vstack([Ax,Ay])
    b = np.zeros([2,n])

    b[:,:-1] = 2*waypoints.T[:,1:] - a[:,1:]
    b[:,-1] = 2*waypoints.T[:,0] - a[:,0]

    return a, b

def compute_t(coef, order, s):
    res = 0
    for idx in range(order):
        res += coef[idx]*np.power(s,order-idx)
    return res

def eval_raw(waypoints, a, b, t):
    segment = np.floor(t)
    segment = np.int(segment)
    n = len(waypoints)
    if segment>=n:
        t =n-0.0001
        segment = n-1
    elif t<0:
        t = 0
    t_val = t-segment
    coords = np.power(1 - t_val, 3) * waypoints.T[:,segment] + 3 * np.power(1 - t_val, 2) * t_val * a[:,segment]\
    + 3 * (1 - t_val) * np.power(t_val, 2) * b[:,segment] + np.power(t_val, 3) * waypoints.T[:,np.int(np.mod(segment+1,n))]

    return coords

def getangle_raw(waypoints, a, b, t):
    der = eval_raw(waypoints, a, b, t+0.1) - eval_raw(waypoints, a, b, t)
    phi = np.arctan2(der[1],der[0])
    return phi

def fit_st(order, norm, waypoints, a, b):

    #fit  the s-t rel.
    nwp = len(waypoints)
    npoints = 500

    #compute approx distance to arc param
    tvals = np.linspace(0, nwp, npoints)

    coords =[]
    for t in tvals:
        coords.append(eval_raw(waypoints, a, b, t))
    coords = np.array(coords)

    dists = []
    for idx in range(npoints):
        dists.append(np.sqrt(np.sum(np.square(coords[idx,:]-coords[np.mod(idx+1,npoints-1),:]))))

    dists = np.cumsum(np.array(dists))
    coef = cp.Variable((order,1))

    #create regressor
    A = dists.reshape(-1,1)
    for idx in range(order-1):
        #print('ord conc %f', idx+2)
        A = np.concatenate((np.power(dists,idx+2).reshape(-1,1),A), axis =1)
    y = tvals.reshape(-1,1)

    objective = cp.Minimize(cp.norm(A@coef-y,norm))
    constraint =[(A@coef)[-1] == y[-1]]
    prob = cp.Problem(objective, constraint)
    result = prob.solve(solver = 'ECOS', verbose=True)

    coeffs = coef.value
    smax = dists[-1]
    svals = np.linspace(0, smax, npoints)
    t_corr = compute_t(coeffs,order,svals)

    plt.figure()
    plt.plot(tvals, dists)
    plt.plot(t_corr, svals)
    plt.xlabel("t (Bezier param) [-]")
    plt.ylabel("s (approx. distance traveled) [m] ")

    return coeffs, smax

def getwaypoints(track):
    #placeholder for now
    trackx = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, \
                  0.8 ,0.8 ,0.8 ,0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0, 0.0, 0.0, 0.0 ])
    tracky = np.array([0.05, 0.3, 0.4, 0.2, 0.2, 0.0, 0.0, 0.0, 0.05, \
                    0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1 ])
    waypoints = np.vstack([trackx,tracky]).T
    return waypoints


def generatelookuptable(track):
    #load track
    waypoints = getwaypoints(track)
    #trackwidth
    r = 0.01
    #abez,bbez coeffs
    a, b = interpolate(waypoints)
    order_inverse = 8
    norm_inverse = 2
    coeffs, smax = fit_st(order_inverse, norm_inverse, waypoints, a, b)

    lutable_density = 100 #[p/m]

    npoints = np.int(np.floor(smax * lutable_density))


    svals = np.linspace(0, smax, npoints)
    tvals = compute_t(coeffs, order_inverse, svals)

    #  entries :
    names_table = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    table = []
    for idx in range(npoints):
        track_point = eval_raw(waypoints, a, b, tvals[idx])
        phi = getangle_raw(waypoints, a, b, tvals[idx])
        n = [-np.sin(phi), np.cos(phi)]
        g_upper = r + track_point[0]*n[0] + track_point[1]*n[1]
        g_lower = -r + track_point[0]*n[0] + track_point[1]*n[1]
        table.append([svals[idx], tvals[idx], track_point[0], track_point[1], phi, np.cos(phi), np.sin(phi), g_upper, g_lower])

    table = np.array(table)
    plot_track(table)
    print("stored as names_table = ", names_table)
    return table

def plot_track (table):
    #downsample
    downsampling = 8
    coords = table[::downsampling, 2:4]
    phis = table[::downsampling, 4]
    svals = table[::downsampling, 0]
    tvals = table[::downsampling, 1]
    cos_phi = table[::downsampling, 5]
    sin_phi = table[::downsampling, 6]
    gvals = table[::downsampling, 7]


    dists = []
    npoints = len(coords)
    for idx in range(npoints):
        dists.append(np.sqrt(np.sum(np.square(coords[idx,:]-coords[np.mod(idx+1,npoints-1),:]))))
    dists = np.cumsum(np.array(dists))

    plt.figure()
    plt.plot(svals, dists)
    plt.plot([0, 4], [0, 4])
    plt.xlabel("t (Bezier param corrected) [m]")
    plt.ylabel("s (approx. distance traveled) [m] ")
    plt.legend(["arclength vs t_corr","x=y"])

    plt.figure()
    len_indicator = 0.05
    plt.plot(table[:,2],table[:,3])
    for idx in range(npoints):
        n = [-sin_phi[idx], cos_phi[idx]]
        g = gvals[idx]
        lm = 0.001
        baseupper = [(coords[idx,0]-lm ), (g-n[0]*(coords[idx,0]-lm )) / n[1]]
        endupper = [(coords[idx,0]+lm ), (g-n[0]*(coords[idx,0]+lm )) / n[1]]
        base = coords[idx,:]
        end = len_indicator * np.array([cos_phi[idx], sin_phi[idx]]) + base
        plt.plot([base[0], end[0]],[base[1], end[1]], color = 'r')
        #plt.plot([baseupper[0], endupper[0]],[baseupper[1], endupper[1]], color = 'g')
