# -*- coding: utf-8 -*-
"""
Positioning and Location Based Services
A.A. 2023/2024
3rd Exercise: Ionospheric delay computation

@author: Marianna Alghisi
"""

'''
Goals: 
   1) 4 zenithal maps of ionospheric error corrections:
   - Elevation: 90°
   - Latitude: [-80°, 80°, step = 0.5°]
   - Longitude: [-180°, 180°, step = 0.5°]
   - time: [00:00, 06:00, 12:00, 18:00]
   
   2) 2 polar maps of ionospheric error corrections:
    - Observer located in Milan
    - Elevation: [0, 90°, step = 0.5°]
    - Azimuth: [-180°, 180°, step = 0.5°]
    time: [00:00, 12:00]
'''

# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import ionoCorrection as ic
import pandas as pd
import geopandas as gpd
from ionoCorrection import iono_correction

# Ionospheric correction parameters:
alpha = [7.4506*10**(-9), 1.4901*10**(-8), -5.9605*10**(-8), -1.1921*10**(-7)]
beta = [9.2160*10**(4), 1.3107*10**(5), -6.5536*10**(4), -5.2429*10**(5)]
ionoparams = alpha + beta
'''
1) ZENITHAL MAPS
'''
# Initialization of the parameters: define inputs for the Zenithal maps (elevation, azimuth, time, lat, lon)

# Loop on time, latitude and longitude --> compute for each point the Ionospheric delay
# TIP: store latitude, longitude and iono_delay in list objects!

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '06:00', '12:00', '18:00']
latitudes = np.arange(-80, 80, 0.5)
longitudes = np.arange(-180, 180, 0.5)
elevation = 90
iono_delay, latitude, longitude = np.array([]), np.array([]), np.array([])
azi = 0

for i in range(4):
    time_in_sec = int(timePrint[i][0:2])*3600
    for lat in latitudes:
        for lon in longitudes:
            IONO = iono_correction(lat, lon, azi, elevation, time_in_sec, ionoparams)
            iono_delay = np.append(iono_delay, IONO)
            latitude = np.append(latitude, lat)
            longitude = np.append(longitude, lon)
    #IONO =   #list of iono delays for considered epoch
    results = pd.DataFrame()
    results['latitude'] = latitude #list of latitude
    results['longitude'] = longitude #list of longitude
    results['iono_delay'] = iono_delay #IONO

    #clear for next iteration
    iono_delay, latitude, longitude = np.array([]), np.array([]), np.array([])

    gdf = gpd.GeoDataFrame(results, geometry=gpd.points_from_xy(results.longitude, results.latitude), crs = 3857)
    world = gpd.read_file('world/world.shp')
    fig, ax = plt.subplots (figsize = (8,6))
    world.boundary.plot(ax=ax, color='black')
    ax.set(xlabel='Longitude', ylabel='Latitude', title='Zenithal map of ionospheric delay at '+ str(timePrint[i]))
    gdf.plot(column='iono_delay', ax = ax, marker='o', legend=True, cmap='jet')
    plt.show()

'''
2) POLAR MAPS
'''
# Definition of Milano's coordinates
lat_mi = 45 + 28/60 + 38.28/60**2
lon_mi = 9 + 10/60 + 53.4/60**2

# Inizialization of the parameters for the loop: time, elevation, azimuth

# Loop on the parameters ---> compute for each point the Ionospheric delay

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '12:00']
elevations = np.arange(0, 90, 0.5)
azimuths = np.arange(-180, 180, 0.5)
iono_delay, azimuth, elevation = np.array([]), np.array([]), np.array([])

for i in [0,1]:
    t = int(timePrint[i][0:2])*3600
    for el in elevations:
        for az in azimuths:
            IONO = iono_correction(lat_mi, lon_mi, az, el, t, ionoparams)
            iono_delay = np.append(iono_delay, IONO)
            elevation = np.append(elevation, 90-el)
            azimuth = np.append(azimuth, np.deg2rad(az))

    fig = plt.figure(figsize = (8,6))
    ax = fig.add_subplot(projection='polar')
    plt.scatter(azimuth, elevation, c=iono_delay, cmap='jet',  alpha=0.75, label=iono_delay)
    ax.set_title('Ionospheric Error Polar Map for Milan Observer time = '+str(t))
    plt.colorbar(label='Ionospheric Delay')
    #rotate the plot to have 0° at the bottom
    ax.set_theta_zero_location('S')
    plt.show()

    #clear for next iteration
    iono_delay, azimuth, elevation = np.array([]), np.array([]), np.array([])