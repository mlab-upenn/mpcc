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
#import cvxpy as cp
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import yaml

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
    n = len(waypoints)
    t = np.mod(t, n)
    segment = np.floor(t)
    segment = np.int(segment)

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

def fit_st(waypoints, a, b):
    #using two revolutions to account for horizon overshooting end of lap

    #fit  the s-t rel.
    nwp = len(waypoints)
    npoints = 20 * nwp
    #compute approx max distance
    tvals = np.linspace(0, nwp, npoints+1)
    coords =[]
    for t in tvals:
        coords.append(eval_raw(waypoints, a, b, t))
    coords = np.array(coords)
    dists = []
    dists.append(0)
    for idx in range(npoints):
        dists.append(np.sqrt(np.sum(np.square(coords[idx,:]-coords[np.mod(idx+1,npoints-1),:]))))
    dists = np.cumsum(np.array(dists))
    smax = dists[-1]

    #--------fit  the s-t rel. to two track revolutions------
    npoints = 2 * 20 * nwp

    #compute approx distance to arc param
    tvals = np.linspace(0, 2*nwp, npoints+1)

    coords =[]
    for t in tvals:
        coords.append(eval_raw(waypoints, a, b, np.mod(t, nwp)))
    coords = np.array(coords)

    distsr = []
    distsr.append(0)
    for idx in range(npoints):
        distsr.append(np.sqrt(np.sum(np.square(coords[idx,:]-coords[np.mod(idx+1,npoints-1),:]))))
    dists = np.cumsum(np.array(distsr))

    ts_inverse = CubicSpline(dists, tvals)
    svals = np.linspace(0, 2*smax, npoints)
    t_corr = ts_inverse(svals)
    #t_corr = compute_t(coeffs,order,svals)

    #plt.figure()
    #plt.plot(tvals, dists, linestyle = '--')
    #plt.plot(t_corr, svals)
    #plt.xlabel("t (Bezier param) [-]")
    #plt.ylabel("s (approx. distance traveled) [m] ")

    return ts_inverse, smax

def getwaypoints(track):
    #placeholder for now
    scaler = 10
    trackx = scaler*np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, \
                  0.8 ,0.8 ,0.8 ,0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0, 0.0, 0.0, 0.0 ])
    tracky = scaler*np.array([0.05, 0.3, 0.4, 0.2, 0.2, 0.0, 0.0, 0.0, 0.05, \
                    0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1 ])
    waypoints = np.vstack([trackx,tracky]).T

    waypoints = np.genfromtxt(track + '.csv', delimiter=',')
    return waypoints


def generatelookuptable(track):
    #load track
    waypoints = getwaypoints(track)
    #plt.scatter(waypoints[:,0], waypoints[:,1])
    #trackwidth
    r = 0.1
    #abez,bbez coeffs
    a, b = interpolate(waypoints)
    ts_inverse, smax = fit_st(waypoints, a, b)

    lutable_density = 100 #[p/m]

    npoints = np.int(np.floor(2 * smax * lutable_density))
    print("table generated with npoints = ", npoints)
    svals = np.linspace(0, 2*smax, npoints)
    tvals = ts_inverse(svals)

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
    #plot_track(table)
    print("Variables stored in following order = ", names_table)
    np.savetxt(str(track) + '_lutab.csv', table, delimiter = ', ')

    dict = {'smax': float(smax), 'ppm' : lutable_density}
    with open(r''+track+'_params.yaml', 'w') as file:
        documents = yaml.dump(dict, file)
    return table, smax


def plot_track (table):
    #downsample
    downsampling = 20
    coords = table[:, 2:4]
    phis = table[::downsampling, 4]
    svals = table[::downsampling, 0]
    tvals = table[::downsampling, 1]
    cos_phi = table[::downsampling, 5]
    sin_phi = table[::downsampling, 6]
    gvals = table[::downsampling, 7]


    dists = []
    dists.append(0)
    npoints = len(coords)
    for idx in range(npoints-1):
        dists.append(np.sqrt(np.sum(np.square(coords[idx,:]-coords[np.mod(idx+1,npoints-1),:]))))
    dists = np.cumsum(np.array(dists))
    dists = dists[::downsampling]
    coords = coords[::downsampling]
    npoints = len(coords)

    plt.figure()
    plt.plot(svals, dists)
    plt.plot([0, svals[-1]], [0, svals[-1]])
    plt.xlabel("t (Bezier param corrected) [m]")
    plt.ylabel("s (approx. distance traveled) [m] ")
    plt.legend(["arclength vs t_corr","x=y"])

    plt.figure()
    len_indicator = 0.05
    downsampling = 10
    plt.plot(table[:,2],table[:,3])
    plt.scatter(table[::downsampling,2],table[::downsampling,3], marker = 'o' )
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
    plt.show()
