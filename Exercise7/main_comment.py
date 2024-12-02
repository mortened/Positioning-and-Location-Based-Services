# -*- coding: utf-8 -*-
"""
POSITIONING & LOCATION BASED SERVICES
AA 2023/2024

EX 7: Inertial Navigation
@author: Marianna Alghisi

-------------------------------------------------------------------------------
Guidelines: a vehicle that is moving in the 2D X-Y inertial planar system with 
two accelerometers and a gyroscope mounted on it. According to the provided data
in Body System compute the position of the vehicle in Iertial System

Input data: both .dat files contain:
    - col1: Epoch [s]
    - col2: Acceleration in X [m/s]
    - col3: Acceleration in Y [m/s]
    - col4: Angular velocity [rad/s]
'Inertial_data.dat' --> observations simulated without errors
'Inertial_data_ni.dat' --> observations simulated with errors

Workflow:
    1) Import the data from both datasets (with and without errors)
       HINT: numpy.loadtxt() function can be used to load .dat files as arrays
    2) Use calcTrajectory function to compute the position of the vehicle
       in each epoch
    3) Plot the results 
-------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt

#1) Import the data
# Load the data from the .dat files using numpy.loadtxt() function

def load_data(filename):
    data = np.loadtxt(filename)
    return data[:,0], data[:,1], data[:,2], data[:,3]

epoch, a_x, a_y, omega = load_data('Inertial_data.dat')
epoch_error, a_x_error, a_y_error, omega_error = load_data('Inertial_data_ni.dat')

#Rotation matrix
R = lambda alpha: np.array([
    [np.cos(alpha), np.sin(alpha)],
    [-np.sin(alpha), np.cos(alpha)]
])

#2) Use the function to compute the trajectory in Inertial RS for both datasets
def calcTrajectory(a_x, a_y, omegaz, epoch):
    # Initialize and compute velocities and delta positions in body frame
    alpha = np.array([0])
    vx = np.array([0])            # Initial velocity in body frame vx
    vy = np.array([0])            # Initial velocity in body frame vy
    XY = np.array([100, 100])           # Initial position in body frame X, Y 
    delta_t = np.diff(epoch)[0]   # Delta t, given that all time steps are equal              
    delta_X = np.array([vx[0]*delta_t + 0.5*a_x[0]*(delta_t**2)])
    # Implement a loop to compute X velocities and delta positions in body frame
    for i in range(1, len(epoch)):   
        vxi = vx[i-1] + a_x[i]*delta_t
        vx = np.append(vx, vxi)
        delta_Xi = vxi * delta_t + 0.5 * a_x[i] * (delta_t**2)
        delta_X = np.append(delta_X, delta_Xi)
    
    # y is constrained to zero because the cart is on a rail: only centrifugal, no skidding
    # clean apparent centrifugal from Y acceleration
    a_y_centr = omegaz * vx
    a_y_clean = a_y  - a_y_centr
    delta_Y = np.array([vy[0] *delta_t + 0.5 * a_y_clean[0] * (delta_t**2)])

    # compute, in case, skidding velocity and displacement
    for i in range(1, len(epoch)):
        vyi = vy[i-1] + a_y_clean[i] * delta_t
        vy = np.append(vy, vyi)
        delta_Yi = vyi * delta_t + 0.5 * a_y_clean[i] * (delta_t**2)
        delta_Y = np.append(delta_Y, delta_Yi)
    
    # Implement a loop to compute for each epoch alpha, R(alpha), rotate Dx from body to
    for i in range(1, len(epoch)):
        alpha = np.append(alpha, alpha[i-1] + omega[i] * delta_t)

    # inertial and update intertial coordinates
        delta_XY = R(alpha[i]) @ np.array([delta_X[i], delta_Y[i]]).flatten()
        XY = np.vstack((XY, XY[i-1] + delta_XY))
    # Return the structure containing the
    return XY

trajectory_errorless = calcTrajectory(a_x, a_y, omega, epoch)
trajectory_error = calcTrajectory(a_x_error, a_y_error, omega_error, epoch_error)

#3) Plot the results
plt.figure(figsize=(10,6))
plt.plot(trajectory_errorless[:,0], trajectory_errorless[:,1], label='No noise', color='Red')
plt.plot(trajectory_error[:,0], trajectory_error[:,1], label='With noise', color='Blue', linestyle='--')
plt.title('Trajectory 2D')
plt.ylabel('y[m]')
plt.xlabel('x[m]')
plt.legend()
plt.show()
