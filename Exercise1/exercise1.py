import numpy as np
pi = np.pi

#Data
#Origin in ITRF Geodetic coordinates
X_0_ITRF = np.array([(44+23/60+24/3600)*pi/180, (8+56/60+20/3600)*pi/180, 70])
sigma_0 = 0.1 #m

#Converting Origin to Global Cartesian ITRF (X,Y,Z)

Rn = lambda a, e, phi: a/np.sqrt(1-e**2*np.sin(phi)**2) # Curvature radius

def GeodeticToGlobalCartesian(phi, lam, h, a, e):
    Rn_ = Rn(a, e, phi)
    X = (Rn_ + h)*np.cos(phi)*np.cos(lam)
    Y = (Rn_ + h)*np.cos(phi)*np.sin(lam)
    Z = (Rn_*(1-e**2)+h)*np.sin(phi)
    return np.array([X, Y, Z])

a = 6378137 #semimajor axis in meters
e = 0.0818191908426215 #eccentricity
#Origin in Global Cartesian ITRF
X_0_ITRF_Cartesian = GeodeticToGlobalCartesian(X_0_ITRF[0], X_0_ITRF[1], X_0_ITRF[2], a, e)

#A,B,C in Boody Frame with respect to origin O:
X_A_BF = np.array([0, 30, 0])
sigma_a = 0.02 #m
X_B_BF = np.array([0, -30, 0])
sigma_b = 0.02 #m
X_C_BF = np.array([200, 0, 0])
sigma_c = 0.1 #m

#Rotation angles body frame/local cartesian frame
alpha = (10.23/3600)*pi/180 #rad
beta = (9.5/3600)*pi/180 #rad
gamma = ((30)+(27/60)+(18/3600))*pi/180 #rad

#From body frame to local cartesian frame

#Rotation matrix
def Rx(alpha):
    return np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
def Ry(beta):
    return np.array([[np.cos(beta), 0, -np.sin(beta)], [0, 1, 0], [np.sin(beta), 0, np.cos(beta)]])
def Rz(gamma):
    return np.array([[np.cos(gamma), np.sin(gamma), 0], [-np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])

def LocalCartesianToBodyFrame(alpha, beta, gamma):
    return Rz(gamma) @ Ry(beta) @ Rx(alpha)

def BodyFrameToLocalCartesian(alpha, beta, gamma):
    return LocalCartesianToBodyFrame(alpha, beta, gamma).T

#In local cartesian frame
X_A_LC = BodyFrameToLocalCartesian(alpha, beta, gamma) @ X_A_BF
X_B_LC = BodyFrameToLocalCartesian(alpha, beta, gamma) @ X_B_BF
X_C_LC = BodyFrameToLocalCartesian(alpha, beta, gamma) @ X_C_BF
#in meters
print("Local cartesian coordinates of points A, B, C")
print(f'Point A:\n N = {X_A_LC[0]:.3f}m \n E = {X_A_LC[1]:.3f}m \n U = {X_A_LC[2]:.3f}m')
print(f'Point B:\n N = {X_B_LC[0]:.3f}m \n E = {X_B_LC[1]:.3f}m \n U = {X_B_LC[2]:.3f}m')
print(f'Point C:\n N = {X_C_LC[0]:.3f}m \n E = {X_C_LC[1]:.3f}m \n U = {X_C_LC[2]:.3f}m')

#Local Cartesian to Global Cartesian

def LocalCartesianToGlobalCartesian(X_0_GC, R_0, X_LC):
    X_GC = X_0_GC + R_0.T @ X_LC
    return X_GC

#Rotation matrix from local cartesian to global cartesian
def R_0(phi, lam):
    return np.array([[-np.sin(lam), np.cos(lam), 0],
                     [-np.sin(phi)*np.cos(lam), -np.sin(phi)*np.sin(lam), np.cos(phi)],
                     [np.cos(phi)*np.cos(lam), np.cos(phi)*np.sin(lam), np.sin(phi)]])

#In Global Cartesian
X_A_GC = LocalCartesianToGlobalCartesian(X_0_ITRF_Cartesian, R_0(X_0_ITRF[0], X_0_ITRF[1]), X_A_LC)
X_B_GC = LocalCartesianToGlobalCartesian(X_0_ITRF_Cartesian, R_0(X_0_ITRF[0], X_0_ITRF[1]), X_B_LC)
X_C_GC = LocalCartesianToGlobalCartesian(X_0_ITRF_Cartesian, R_0(X_0_ITRF[0], X_0_ITRF[1]), X_C_LC)

print("ITRF global cartesian coordinates of points A,B,C ")
print(f'Point A:\n X = {X_A_GC[0]}m \n Y = {X_A_GC[1]}m \n Z = {X_A_GC[2]}m')
print(f'Point B:\n X = {X_B_GC[0]}m \n Y = {X_B_GC[1]}m \n Z = {X_B_GC[2]}m')
print(f'Point C:\n X = {X_C_GC[0]}m \n Y = {X_C_GC[1]}m \n Z = {X_C_GC[2]}m')

#Converted coordinates to ETRF2014 on on September 1st 2022 using ETRF/ITRF Coordinate Transformation Tool
#doi:10.24414/ROB-EUREF-ECTT. 

#Point A
X_A_ETRF = np.array([4509854.8132, 709344.7332, 4439228.7612])

#Points you get when only using three decimals in latitude and longitude from above (i.e. the ones in results.txt)
# X_A_ETRF = np.array([4509854.8129, 709344.7326, 4439228.7612])

#Point B
X_B_ETRF = np.array([4509885.8305, 709380.3976, 4439191.8020])
#Point C
X_C_ETRF = np.array([4509773.4716, 709521.8569, 4439282.7126])

#Converting ETRF Global Cartesian to Geodetic

def CartesianToGeodetic(X, Y, Z, a, e):
    r = np.sqrt(X**2 + Y**2)
    b = a*np.sqrt(1-e**2)
    e_b_2 = (a**2 - b**2)/b**2
    psi = np.arctan2(Z,(r*np.sqrt(1-e**2)))
    lam = np.arctan2(Y,X)
    phi = np.arctan2(Z + e_b_2*b*np.sin(psi)**3, r - e**2*a*np.cos(psi)**3)
    Rn_ = Rn(a, e, phi)
    h = r/(np.cos(phi)) - Rn_
    return np.array([phi, lam, h])

print("Geodetic coordinates of points A, B, C in ETRF")
def DegreestoDMS(deg):
    d = int(deg)
    m = (deg-d)*60
    s = (m-int(m))*60
    return d, int(m), s

phi_A, lam_A, h_A = CartesianToGeodetic(X_A_ETRF[0], X_A_ETRF[1], X_A_ETRF[2], a, e)
phi_B, lam_B, h_B = CartesianToGeodetic(X_B_ETRF[0], X_B_ETRF[1], X_B_ETRF[2], a, e)
phi_C, lam_C, h_C = CartesianToGeodetic(X_C_ETRF[0], X_C_ETRF[1], X_C_ETRF[2], a, e)

#Converting to degrees, minutes and seconds
phi_A = DegreestoDMS(np.degrees(phi_A))
lam_A = DegreestoDMS(np.degrees(lam_A))
phi_B = DegreestoDMS(np.degrees(phi_B))
lam_B = DegreestoDMS(np.degrees(lam_B))
phi_C = DegreestoDMS(np.degrees(phi_C))
lam_C = DegreestoDMS(np.degrees(lam_C))


print(f'Point A:\n lat = {phi_A}')
print(f' long = {lam_A}')
print(f' h = {h_A}m')
print(f'Point B:\n lat = {phi_B}')
print(f' long = {lam_B}')
print(f' h = {h_B}m')
print(f'Point C:\n lat = {phi_C}')
print(f' long = {lam_C}')
print(f' h = {h_C}m')


#Covariance Propagation
#Covariance matrix of each point in body frame
sigma_a, sigma_b, sigma_c, sigma_0 = sigma_a*100, sigma_b*100, sigma_c*100, sigma_0*100 #m to cm
C_BF_A = np.diag([sigma_a**2, sigma_a**2, sigma_a**2])
C_BF_B = np.diag([sigma_b**2, sigma_b**2, sigma_b**2])
C_BF_C = np.diag([sigma_c**2, sigma_c**2, sigma_c**2])

#Covariance matrix of each point in local cartesian frame
C_LC_A = BodyFrameToLocalCartesian(alpha, beta, gamma) @ C_BF_A @ BodyFrameToLocalCartesian(alpha, beta, gamma).T
C_LC_B = BodyFrameToLocalCartesian(alpha, beta, gamma) @ C_BF_B @ BodyFrameToLocalCartesian(alpha, beta, gamma).T
C_LC_C = BodyFrameToLocalCartesian(alpha, beta, gamma) @ C_BF_C @ BodyFrameToLocalCartesian(alpha, beta, gamma).T

#Covariance matrix of P0 in global cartesian frame
C_GC_P0 = np.diag([sigma_0**2, sigma_0**2, sigma_0**2])

#Covariance matrix of each point in global cartesian frame
C_GC_Delta_X_A = R_0(X_0_ITRF[0], X_0_ITRF[1]).T @ C_LC_A @ R_0(X_0_ITRF[0], X_0_ITRF[1])
C_GC_A = C_GC_P0 + C_GC_Delta_X_A
C_GC_Delta_X_B = R_0(X_0_ITRF[0], X_0_ITRF[1]).T @ C_LC_B @ R_0(X_0_ITRF[0], X_0_ITRF[1])
C_GC_B = C_GC_P0 + C_GC_Delta_X_B
C_GC_Delta_X_C = R_0(X_0_ITRF[0], X_0_ITRF[1]).T @ C_LC_C @ R_0(X_0_ITRF[0], X_0_ITRF[1])
C_GC_C = C_GC_P0 + C_GC_Delta_X_C

STD_A = np.sqrt(np.diag(C_GC_A))
STD_B = np.sqrt(np.diag(C_GC_B))
STD_C = np.sqrt(np.diag(C_GC_C))

print("Standard deviations of points A, B, C in cm")
print(f'Point A:\n N = {STD_A[0]:.3f}cm \n E = {STD_A[1]:.3f}cm \n U = {STD_A[2]:.3f}cm')
print(f'Point B:\n N = {STD_B[0]:.3f}cm \n E = {STD_B[1]:.3f}cm \n U = {STD_B[2]:.3f}cm')
print(f'Point C:\n N = {STD_C[0]:.3f}cm \n E = {STD_C[1]:.3f}cm \n U = {STD_C[2]:.3f}cm')