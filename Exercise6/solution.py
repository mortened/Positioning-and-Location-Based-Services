import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Read csv data
data = pd.read_csv('data.csv', header=None).values
x = data[:,0]
y = data[:,1]

#Standard deviations
sigma_xy = 1 #m
sigma_xy_dot = 1 #m/s
sigma_obs = 25 #m

#Covariance matrices
C_model = np.diag([sigma_xy**2, sigma_xy**2, sigma_xy_dot**2, sigma_xy_dot**2])
C_obs = np.diag([sigma_obs**2, sigma_obs**2])

#Initial covariance error
C_error = np.diag([10**2, 10**2, 1**2, 1**2])

#Initial state
xi = np.array([x[0], y[0], 0, 0])

#State transition matrix
T = np.array([[1, 0, 1, 0],
             [0, 1, 0, 1],
             [0, 0, 1, 0],
             [0, 0, 0, 1]])
             

#Observation matrix
A = np.array([[1, 0, 1, 0],
              [0, 1, 0, 1]])

#Kalman filter outputs
x_kalman, y_kalman = [], []

for i in range(1, len(x)):
    #State prediction
    Xi = T @ xi
    print(Xi)

    #Update prediction
    K = C_model + T @ C_error @ T.T
    G = K @ A.T @ np.linalg.inv(A @ K @ A.T + C_obs)
    I = np.eye(4)
    yi_plus_1 = np.array([x[i], y[i]])
    xi = G @ yi_plus_1 + (I - G @ A) @ Xi

    #Update error covariance
    C_error = (I - G @ A) @ K

    #Store the results
    x_kalman.append(xi[0])
    y_kalman.append(xi[1])

#Plotting data against Kalman filter output
plt.figure(figsize=(8,6))
plt.plot(x, y, 'o', label='Data points', color='blue', markersize=3)
plt.plot(x_kalman, y_kalman, '-', label='Kalman filter', color='red', linewidth=4)
plt.xlabel('X')
plt.ylabel('Y')
plt.xlim(100, 900)
plt.ylim(450, 850)
plt.legend()
plt.show()