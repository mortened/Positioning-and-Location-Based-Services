'''
Positioning and Location Based Services
A.A. 2024/2025
Exercise 2:  GPS orbits computation

Marianna Alghisi

References:
    - http://www.navipedia.net/index.php/Clock_Modelling
    - http://www.ajgeomatics.com/SatelliteOrbitsRev0828-2010.pdf
    - http://www.navipedia.net/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation
'''

# Import all the required libraries and external script containing the functions

import math
import numpy as np
import matplotlib.pyplot as plt
import functions as fn

# Load Almanac of satellite SVN 63, PRN 01 (Block IIR) for 2016
# The almanac contains the following information:
dt0 = -7.661711424589e-05
dt1 = -3.183231456205e-12
dt2 = 0.000000000000e+00
square_root_a = 5.153650835037e+03 #square root of meters
e = 3.841053112410e-03
M0 = 1.295004883409e+00
Omega0 = -2.241692424630e-01 #radians
Omegadot = -8.386063598924e-09 #radians/sec
i0 = 9.634782624741e-01 #Radians
idot = -7.286017777600e-11 #Radians/sec
w0 = 9.419793734505e-01 #Radians
wdot = 0.000000000000e+00 #Radians/sec

GMe = 3.986005e14 # m^3/s^2
OmegaEdot = 7.2921151467e-5 #radians

# Create a list containing all the epochs
t_start = 1 
t_end = 60*60*24 # 1 day in seconds
t_step = 1 # 1 second
t = list(range(t_start, t_end, t_step))

'''
1) Compute clock offsets and plot it
'''
# Tip: use a loop on t and append the elements to a list
clock_offsets = []
for dt in t:
    clock_offset = dt0 + dt1*dt + dt2*dt**2
    clock_offsets.append(clock_offset)

# plot
fig, ax = plt.subplots(figsize=(10,6))
ax.set(xlabel='Seconds in a day', ylabel='Clock-offset', title = 'Clock-offset variations in one day')
ax.plot(t, clock_offsets, '-', color='blue')
#plt.show()

'''
2) Compute the position of the satellite in Cartesian ITRF [X, Y, Z]
'''
coord_ORS = []
coord_ITRF = []

n = np.sqrt((GMe)/(square_root_a**6)) # Mean motion

#Mean anomaly M(t) = M0 + n(t - t0)
M = [M0 + n * (dt - t_start) for dt in t]  

# Eccentric anomaly
eta = []

# Computing eccentric anomaly
def eccentric_anomaly(M, e, tol):
    #Initial guess
    eta = M
    #Iterate until convergence
    while True:
        eta_next = M + e*np.sin(eta)
        if abs(eta_next - eta) < tol:
            break
        eta = eta_next
    return eta

for i in range(len(t)):
    eta_val = eccentric_anomaly(M[i], e, 1e-10)
    eta.append(eta_val)

# Rotation matrices
def Rz(gamma):
    return np.array([[np.cos(gamma), np.sin(gamma), 0], [-np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])
def Rx(alpha):
    return np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])

#Storing ITRF coordinates
coord_ITRF = []

for i, dt in enumerate(t):
	# Compute psi
    psi = np.atan2(
        (np.sqrt(1-e**2)*np.sin(eta[i])),(np.cos(eta[i])-e)
    )
	# Compute radius r
    r = (square_root_a**2 *(1-e**2))/(1+e*np.cos(psi))
	# Compute the coordinates of the satellite in ORS and store it in coord_ORS
    xORS = r * np.cos(psi)
    yORS = r * np.sin(psi)
    zORS = 0
    # Compute rotation angles omega, i, OMEGA
    OMEGA = Omega0 + ((Omegadot - OmegaEdot) * (dt - t_start))
    i = i0 + ((idot)*(dt - t_start))
    omega = w0 + (wdot *(dt - t_start))
    # Compute the rotation matrices required to transform from ORS to ITRF
    # R3(-OMEGA(t))
    ROMEGA = Rz(-OMEGA)
    # R1(-i(t))
    Ri = Rx(-i)
    # R3(-omega(t))
    Romega = Rz(-omega)
    # Final rotation matrix R
    R = ROMEGA @ Ri @ Romega
    # Compute the coordinates of the satellites in ITRF and store it in coord_ORS
    xITRF, yITRF, zITRF = R @ np.array([xORS, yORS, zORS])
    coord_ITRF.append([xITRF, yITRF, zITRF])


'''
3) Convert satellite's coordinates from global cartesian [X, Y, Z] to godetic [latitude, longitude, height]
'''
# Tip: You can use the cartesian to geodetic function that you implemented for the first exercise

#Function from exercise 1
Rn = lambda a, e, phi: a/np.sqrt(1-e**2*np.sin(phi)**2) # Curvature radius
def CartesianToGeodetic(X, Y, Z, a, e):
    r = np.sqrt(X**2 + Y**2)
    b = a * np.sqrt(1-e**2)
    e_b_2 = (a**2 - b**2)/b**2
    psi = np.arctan2(Z,(r*np.sqrt(1-e**2)))
    lam = np.arctan2(Y,X) # Longitude
    phi = np.arctan2(Z + e_b_2*b*np.sin(psi)**3, r - e**2*a*np.cos(psi)**3) # Latitude
    Rn_ = Rn(a, e, phi)
    h = r/(np.cos(phi)) - Rn_

    # Convert lat, lon from radians to degrees
    phi, lam = phi*180/np.pi, lam*180/np.pi
    return phi, lam, h

# Create 3 lists containing the values of Latitude, Longitude and Height
lat = []
lon = []
h = []
for i in range(len(t)):
    xITRF, yITRF, zITRF = coord_ITRF[i]
    lat_, lon_, h_ = CartesianToGeodetic(xITRF, yITRF, zITRF, square_root_a**2, e)
    lat.append(lat_)
    lon.append(lon_)
    h.append(h_)
# Compute the average radius
r_av = np.mean(h)

'''
4) Plot satellite's daily trajectory with basemap
'''
# EXAMPLE OF PLOTTING THE GROUNDTRACK
# REQUIRED LIBRARIES FOR THIS SOLUTION: pandas and geopandas
import pandas as pd
import geopandas as gpd

# Realize a dataframe containing satellite coordinates
df = pd.DataFrame()  # Create an empty DataFrame
df['time'] = t  # Assign the time list to the 'time' column
df['Latitude'] = lat  # Assign latitude values
df['Longitude'] = lon  # Assign longitude values


# Transform the DataFrame in a GeoDataFrame
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['Longitude'], df['Latitude']), crs = 3857)

# Load the basemap
world = gpd.read_file('world/world.shp')

# Plot the trajectory with world basemap
fig, ax = plt.subplots (figsize = (8,6))
world.plot(ax=ax)
ax.set(xlabel='Longitude', ylabel='Latitude', title='Satellite daily trajectory')
gdf.plot(ax = ax, marker='o', color='red')
plt.show()
'''
# 5) Plot height of the satellite 
# '''
# mean_h = #average value of the sat heights in KM
# fig, ax = plt.subplots(figsize=(10,6))
# ax.set(xlabel='seconds in one day (00:00 - 23:59 = 86400 sec)', ylabel='[km]', title='ellipsoidic height variations [km] around mean height = '+str(mean_h)+' [km]')
# ax.plot(t, your_variables_heights_in_km, '-', color='blue')
# '''
# 6) Print on a text file the results
# '''




