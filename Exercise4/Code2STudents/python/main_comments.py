# -*- coding: utf-8 -*-
"""
POSITIONING & LOCATION BASED SERVICES
AA 2023/2024

EX 4: Point Positioning

@author: Marianna Alghisi
"""

import numpy as np
#import functions4 as fn
from variables import *
from functions2students import CartesianToGeodetic, topocent, tropo_correction, iono_correction, R_0

# Initialization of the parameters
max_iter = 10
convergence = 0.1
X_start = np.array([[0], [0], [0], [0]])
e = 0.0818191908426215
a = 6378137

# Least Squares iteration --> tip: use a loop!
for i in range(max_iter):
    # Initialize the structures for the Least Squares (A, b)
    A = np.zeros((11, 4))
    b = np.zeros((11, 1))

    # For each satellite compute:
    for j in range(len(sat_ids)):
        # - topocentric positions of satellite (approx distance, elevation and azimuth)
        xyz_sat_j = np.array([xyz_sat[j,:]]).T
        ro_approx, az, el = topocent(X_start[0:3], xyz_sat_j)
        # - tropospheric and ionospheric corrections
        phi_r, lam_r, h_r = CartesianToGeodetic(X_start[0][0], X_start[1][0], X_start[2][0], a, e)
        tropo_corr = tropo_correction(h_r, el)
        iono_corr = iono_correction(phi_r, lam_r, az, el, time_rx, alpha + beta)
        # - b = LS known term vector
        b[j] = ro_approx - c*dtS[j] + tropo_corr + iono_corr
        # - A = LS design matrix
        e_j = (X_start[0:3] - xyz_sat_j)/ro_approx
        A[j] = np.array([e_j[0][0], e_j[1][0], e_j[2][0], 1])

    # Implement LS solution and estimate the corrections to receiver pos
    N = A.T @ A
    dP_r_s = pr_C1 - b
    est_corr = np.linalg.inv(N) @ A.T @ dP_r_s
    # Estimated coord = approximate + estimated correction
    x_r_hat = X_start[0:3] + est_corr[0:3]

    #Estimated corrections, delta x, delta y, delta z

    # Update the coordinates and the clock offset of the receiver
    X_start = np.array([x_r_hat[0], x_r_hat[1], x_r_hat[2],est_corr[3]/c])
    
    # Check convergence of the resut, in case break the loop
    if max(abs(est_corr[0:3])) <= convergence:
        print('Break at cycle: ', i+1)
        break
    # Check at the end that convergence didn't fail
    if i == (max_iter-1):
        print('Convergence failed')

X = X_start
xr, yr, zr = X[0][0], X[1][0], X[2][0]
dtr = X[3][0]

# Final estimated unknowns
# LS residuals and sigma2
y = A @ X + b
v = pr_C1 - y
sigma2 = (v.T @ v)/(len(sat_ids)-4)
# Covariance matrix of the estimated coordinates
Cxx = sigma2 * np.linalg.inv(N)
# PDOP computeation
Q_geom = np.linalg.inv(N)[0:3, 0:3]
phi_r, lam_r, h_r = CartesianToGeodetic(xr, yr, zr, a, e)
Q = R_0(np.deg2rad(phi_r), np.deg2rad(lam_r)) @ Q_geom @ R_0(np.deg2rad(phi_r), np.deg2rad(lam_r)).T
PDOP = np.sqrt(Q[0,0] + Q[1,1] + Q[2,2])

print('Coordinates of the receiver:')
print('X: ', xr, '\nY: ', yr, '\nZ: ', zr)
print('Clock offset of the receiver: \n', dtr)
print('PDOP: \n', PDOP)


print('\n\n----------------------------------------')
'''Repeat with cut-off angle of 5°'''
# Initialization of the parameters
X_start = X # result previously obtained
cutoff = 5

# Loop over all the available satellites:
new_sat = []
for i in range(len(sat_ids)):
    # - Compute the elevations
    xyz_sat_i = np.array([xyz_sat[i,:]]).T
    ro_approx, az, el = topocent(X_start[0:3], xyz_sat_i)
    # - Discard all the satellites with elevation < 5°
    if el >= cutoff:
        new_sat.append(i)


# Implement LS solution for the new staellite configuration
for i in range(max_iter):
    # Initialize the structures for the Least Squares (A, b)
    A = []
    b = []

    # For each satellite compute:
    for j in new_sat:
        # - topocentric positions of satellite (approx distance, elevation and azimuth)
        xyz_sat_j = np.array([xyz_sat[j,:]]).T
        ro_approx, az, el = topocent(X_start[0:3], xyz_sat_j)
        # - tropospheric and ionospheric corrections
        phi_r, lam_r, h_r = CartesianToGeodetic(X_start[0][0], X_start[1][0], X_start[2][0], a, e)
        tropo_corr = tropo_correction(h_r, el)
        iono_corr = iono_correction(phi_r, lam_r, az, el, time_rx, alpha + beta)
        # - b = LS known term vector
        b.append(ro_approx - c*dtS[j] + tropo_corr + iono_corr)
        # - A = LS design matrix
        e_j = (X_start[0:3] - xyz_sat_j)/ro_approx
        A.append([e_j[0][0], e_j[1][0], e_j[2][0], 1])

    A = np.array(A)
    b = np.array(b)
    # Implement LS solution and estimate the corrections to receiver pos
    N = A.T @ A
    dP_r_s = pr_C1[new_sat] - b
    est_corr = np.linalg.inv(N) @ A.T @ dP_r_s
    # Estimated coord = approximate + estimated correction
    x_r_hat = X_start[0:3] + est_corr[0:3]

    #Estimated corrections, delta x, delta y, delta z

    # Update the coordinates and the clock offset of the receiver
    X_start = np.array([x_r_hat[0], x_r_hat[1], x_r_hat[2],est_corr[3]/c])
    
    # Check convergence of the resut, in case break the loop
    if max(abs(est_corr[0:3])) <= convergence:
        print('Break at cycle: ', i+1)
        break
    # Check at the end that convergence didn't fail
    if i == (max_iter-1):
        print('Convergence failed')

X = X_start
xr, yr, zr = X[0][0], X[1][0], X[2][0]
dtr = X[3][0]

# Final estimated unknowns
# LS residuals and sigma2
y = A @ X + b # predicted values
v = pr_C1[new_sat] - y # residuals
sigma2 = (v.T @ v)/(len(new_sat)-4)
# Covariance matrix of the estimated coordinates
Cxx = sigma2 * np.linalg.inv(N)
# PDOP computeation
Q_geom = np.linalg.inv(N)[0:3, 0:3]
phi_r, lam_r, h_r = CartesianToGeodetic(X[0][0], X[1][0], X[2][0], a, e)
Q = R_0(np.deg2rad(phi_r), np.deg2rad(lam_r)) @ Q_geom @ R_0(np.deg2rad(phi_r), np.deg2rad(lam_r)).T
PDOP = np.sqrt(Q[0,0] + Q[1,1] + Q[2,2])


print('In view satellites: \n', len(new_sat))
print('Coordinates of the receiver:')
print('X: ', xr, '\nY: ', yr, '\nZ: ', zr)
print('Clock offset of the receiver: \n', dtr)
print('PDOP: \n', PDOP)