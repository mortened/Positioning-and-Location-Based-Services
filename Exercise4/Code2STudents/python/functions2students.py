"""
Positioning & Location Based Services
AA: 2023/2024

@author: Marianna Alghisi
"""

import numpy as np
import math

e = 0.0818191908426215
a = 6378137

#Rotation matrix from local cartesian to global cartesian
def R_0(phi, lam):
    return np.array([[-np.sin(lam), np.cos(lam), 0],
                     [-np.sin(phi)*np.cos(lam), -np.sin(phi)*np.sin(lam), np.cos(phi)],
                     [np.cos(phi)*np.cos(lam), np.cos(phi)*np.sin(lam), np.sin(phi)]])

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

def topocent(X, sat):
    '''
    Parameters 
    ----------
    X : np.array of the receiver position in global cartesian.
    sat : np.array of the satellite position in global cartesian.

    Returns
    -------
    ro_approx : distance between the satellite and the receiver [m].
    el : elevation of the satellite with respect the receiver [rad].
    az : azimuth of the satellite with respect the receiver [rad].

    '''
    
    delta = sat - X
    ro_approx = np.sqrt(delta[0][0]**2 + delta[1][0]**2 + delta[2][0]**2)

    X_g = CartesianToGeodetic(X[0][0], X[1][0], X[2][0], a, e)
    lat0 = np.deg2rad(X_g[0])
    lon0 = np.deg2rad(X_g[1])
    
    R = np.array([ [-np.sin(lon0), np.cos(lon0), 0],
                    [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
                    [np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
                    ])
    
    ENU = np.dot(R, delta)
    E, N, U = ENU[0][0], ENU[1][0], ENU[2][0]
    
    d_h = np.sqrt(E**2 + N**2)
    
    if d_h < 0.1:
        az = 0
        el = 90
    else:
        az = np.arctan2(E, N)
        el = np.arctan2(U, d_h)
    
        #CONVERT AZIMUTH AND ELEVATION TO DEGREES
        az = np.rad2deg(az)
        el = np.rad2deg(el)

    return [ro_approx, az, el]    


def tropo_correction (h, el):
    '''
    Computation of the pseudorange correction due to tropospheric refraction.
    Saastamoinen model.
    
    Parameters
    ----------
    h : height f the receiver [m].
    el : elevation of the satellite with respect the receiver [degrees].

    Returns
    -------
    tropoDelay : tropospheric effect [m].

    '''
    if (h > -500) and (h < 5000):
        eta = np.deg2rad(el)
        #eta satellite elevation in radians
        Po = 1013.25 #mBar
        To = 291.15 #degree Kelvin
        Ho = 50/100
        ho = 0
        height = h - ho # h is the ellipsoidal height of the receiver
        Pr = Po * (1-0.0000226 * height)**5.225 #pressure
        Tr = To - 0.0065 * height #temperature
        Hr = Ho * math.exp(-0.0006396 * height)
        er = Hr * math.exp(-37.2465 + (0.213166 * Tr) - 0.000256908 * (Tr)**2) #humidity
        tropoDelay = 0.002277 / np.sin(eta) * (Pr + er * (1255/Tr + 0.05) - (np.tan(eta)) ** -2)
        return tropoDelay
    else:
        tropoDelay = 0
        return tropoDelay

def iono_correction(lat, lon, az, el, time, ionoparams):
    # Initialization
    c = 299792458

    # Ionospheric parameters
    a0, a1, a2, a3, b0, b1, b2, b3 = ionoparams

    # Elevation from 0 to 90 degrees
    el = np.abs(el)

    # Conversion to semicircles
    lat = lat / 180
    lon = lon / 180
    az = az / 180
    el = el / 180

    psi = (0.0137 / (el + 0.11)) - 0.022

    phi = lat + psi * np.cos(az * np.pi)

    if phi > 0.416:
        phi = 0.416
    elif phi < -0.416:
        phi = -0.416

    # Geodetic longitude of the earth projection of the ionospheric intersection point
    lambda_ = lon + (psi * np.sin(az * np.pi)) / np.cos(phi * np.pi)

    # Geomagnetic latitude of the earth projection of the ionospheric intersection point
    ro = phi + 0.064 * np.cos((lambda_ - 1.617) * np.pi)

    # Local time in seconds
    t = lambda_ * 43200 + time

    if t >= 86400:
        t -= 86400
    elif t < 0:
        t += 86400

    # Obliquity factor
    f = 1 + 16 * (0.53 - el) ** 3

    a = a0 + a1 * ro + a2 * ro ** 2 + a3 * ro ** 3
    a = np.where(a < 0, 0, a)

    p = b0 + b1 * ro + b2 * ro ** 2 + b3 * ro ** 3
    # p = np.where(p < 72000, 72000, p)

    x = (2 * np.pi * (t - 50400)) / p

    # Ionospheric delay
    if abs(x) < 1.57:
        delay = c * f * (5e-9 + a * (1 - (x ** 2) / 2 + (x ** 4) / 24))

    else: 
        delay = c * f * 5e-9

    return delay