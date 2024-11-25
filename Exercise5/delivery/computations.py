import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Read csv data
data = pd.read_csv('coordinates.csv')

#Fetch the coordinates
x = np.array(data['X'].values)
y = np.array(data['Y'].values)
z = np.array(data['Z'].values)

#Compute mean coordinates for the week of accusation
x_mean = np.mean(x)
y_mean = np.mean(y)
z_mean = np.mean(z)

#Compute N, E, U coordinates from the mean coordinates
e = 0.0818191908426215
a = 6378137
Rn = lambda a, e, phi: a/np.sqrt(1-e**2*np.sin(phi)**2) # Curvature radius
def CartesianToGeodetic(X, Y, Z, a, e):
    r = np.sqrt(X**2 + Y**2)
    b = a * np.sqrt(1-e**2)
    e_b_2 = (a**2 - b**2)/b**2
    psi = np.arctan2(Z,(r*np.sqrt(1-e**2)))
    lam = np.arctan2(Y,X) # Longitude
    phi = np.arctan2(Z + e_b_2* b * np.sin(psi)**3, r - e**2*a*np.cos(psi)**3) # Latitude
    Rn_ = Rn(a, e, phi)
    h = r/(np.cos(phi)) - Rn_
    # Convert lat, lon from radians to degrees
    return np.rad2deg(phi), np.rad2deg(lam), h

phi, lam, h = CartesianToGeodetic(x_mean, y_mean, z_mean, a, e)

#Geocentrix baselines from the mean coordinates
delta_x = x - x_mean
delta_y = y - y_mean
delta_z = z - z_mean
delta_xyz = np.vstack((delta_x, delta_y, delta_z))

#Rotation matrix from local cartesian to global cartesian
def R_0(phi, lam):
    return np.array([[-np.sin(lam), np.cos(lam), 0],
                     [-np.sin(phi)*np.cos(lam), -np.sin(phi)*np.sin(lam), np.cos(phi)],
                     [np.cos(phi)*np.cos(lam), np.cos(phi)*np.sin(lam), np.sin(phi)]])

R = R_0(np.deg2rad(phi), np.deg2rad(lam))
ENU = R @ delta_xyz

E = ENU[0]
N = ENU[1]
U = ENU[2]
Dates = np.array(data['Date'].values)

print(f'Mean coordinates: \nX: {x_mean} m\nY: {y_mean} m\nZ: {z_mean} m')


plt.figure(figsize=(15,3))
plt.plot(Dates, E, label='East', color='blue')
plt.xlabel('Date')
plt.ylabel('East [m]')
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(15,3))
plt.plot(Dates, N, label='North', color='green')
plt.xlabel('Date')
plt.ylabel('North [m]')
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(15,3))
plt.plot(Dates, U, label='Up', color='red')
plt.xlabel('Date')
plt.ylabel('Up [m]')
plt.legend()
plt.grid()
plt.show()