# -*- coding: utf-8 -*-
"""
Positioning and Location Based Services
A.A. 2023/2024
3rd Exercise: Ionospheric delay computation

@author: Marianna Alghisi
"""

import numpy as np 

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
