import numpy as np
import math
import time

# Constants
ALIVE = 1
DEAD = 0
NN = 100
THRESHOLD = 0.01
CHANCE = 0.1
PI = np.pi

# Functions for Monte Carlo Simulation

def random_gen():
    return np.random.rand()

def rotSphi(S, phi):
    """ Rotate Stokes vector by azimuthal angle phi """
    cos2phi = np.cos(2 * phi)
    sin2phi = np.sin(2 * phi)
    S2 = np.zeros(4)
    S2[0] = S[0]
    S2[1] = S[1] * cos2phi + S[2] * sin2phi
    S2[2] = -S[1] * sin2phi + S[2] * cos2phi
    S2[3] = S[3]
    return S2

def rotateYYZZ(YY, ZZ, phi, theta):
    """ Rotate coordinates in local photon frame """
    YY2 = np.zeros(3)
    ZZ2 = np.zeros(3)
    
    ct = np.cos(phi)
    st = np.sin(phi)
    vt = 1 - ct
    
    YY2[0] = (ZZ[0] * ZZ[0] * vt + ct) * YY[0] + (ZZ[0] * ZZ[1] * vt - ZZ[2] * st) * YY[1] + (ZZ[0] * ZZ[2] * vt + ZZ[1] * st) * YY[2]
    YY2[1] = (ZZ[0] * ZZ[1] * vt + ZZ[2] * st) * YY[0] + (ZZ[1] * ZZ[1] * vt + ct) * YY[1] + (ZZ[1] * ZZ[2] * vt - ZZ[0] * st) * YY[2]
    YY2[2] = (ZZ[0] * ZZ[2] * vt - ZZ[1] * st) * YY[0] + (ZZ[1] * ZZ[2] * vt + ZZ[0] * st) * YY[1] + (ZZ[2] * ZZ[2] * vt + ct) * YY[2]
    
    temp = 1 / np.sqrt(YY2[0]**2 + YY2[1]**2 + YY2[2]**2)
    YY2 *= temp

    ct = np.cos(theta)
    st = np.sin(theta)
    vt = 1 - ct

    ZZ2[0] = (YY2[0] * YY2[0] * vt + ct) * ZZ[0] + (YY2[0] * YY2[1] * vt - YY2[2] * st) * ZZ[1] + (YY2[0] * YY2[2] * vt + YY2[1] * st) * ZZ[2]
    ZZ2[1] = (YY2[0] * YY2[1] * vt + YY2[2] * st) * ZZ[0] + (YY2[1] * YY2[1] * vt + ct) * ZZ[1] + (YY2[1] * YY2[2] * vt - YY2[0] * st) * ZZ[2]
    ZZ2[2] = (YY2[0] * YY2[2] * vt - YY2[1] * st) * ZZ[0] + (YY2[1] * YY2[2] * vt + YY2[0] * st) * ZZ[1] + (YY2[2] * YY2[2] * vt + ct) * ZZ[2]
    
    temp = 1 / np.abs(np.sqrt(ZZ2[0]**2 + ZZ2[1]**2 + ZZ2[2]**2))
    ZZ2 *= temp

    return YY2, ZZ2

# Main Monte Carlo simulation function
def monte_carlo_simulation():
    # Parameters
    radius = 2.0 / 2  # microns
    lambda_ = 0.6328  # microns
    rho = 1.152e-4  # Dilution
    Nphotons = int(1e5)
    mua = 0.0  # absorption coefficient
    nre_p = 1.59
    nre_med = 1.33
    nim_med = 0.0
    nangles = 1000

    # Initialize matrices for Stokes parameters (intensity and polarization)
    MM = NN - 1
    IR = np.zeros((NN, NN))
    QR = np.zeros((NN, NN))
    UR = np.zeros((NN, NN))
    VR = np.zeros((NN, NN))

    # Allocate Stokes vectors and direction vectors
    S0 = np.zeros(4)
    S = np.zeros(4)
    S2 = np.zeros(4)
    XX = np.array([1.0, 0.0, 0.0])
    YY = np.array([0.0, 1.0, 0.0])
    ZZ = np.array([0.0, 0.0, 1.0])

    start_time = time.time()

    # Launch different polarization states (H, V, P, R)
    for jjj in range(1, 5):
        if jjj == 1:
            S0 = [1, 1, 0, 0]  # Horizontal polarization
            print("launch H")
        elif jjj == 2:
            S0 = [1, -1, 0, 0]  # Vertical polarization
            print("launch V")
        elif jjj == 3:
            S0 = [1, 0, 1, 0]  # Plus 45° polarization
            print("launch P")
        elif jjj == 4:
            S0 = [1, 0, 0, 1]  # Right circular polarization
            print("launch R")

        # Launch Nphotons
        for i_photon in range(Nphotons):
            x = y = z = 0.0  # Initial photon position
            W = 1  # Photon weight
            photon_status = ALIVE

            # Photon Stokes vector
            S = np.copy(S0)

            while photon_status == ALIVE:
                # Step size calculation
                rnd = random_gen()
                while rnd == 0:
                    rnd = random_gen()
                s = -np.log(random_gen()) / (rho)

                x += ZZ[0] * s
                y += ZZ[1] * s
                z += ZZ[2] * s

                # Absorption
                absorb = W * (1 - rho)
                W -= absorb

                if z <= 0:
                    IR[int(NN/2)][int(NN/2)] += S[0]  # Store intensity in central pixel
                    photon_status = DEAD

        # Write results to file
        with open(f'outH{jjj}.dat', 'w') as f:
            for row in IR:
                f.write('\t'.join(f"{val:.5f}" for val in row) + '\n')

    finish_time = time.time()
    print(f"Elapsed Time = {finish_time - start_time} seconds")

if __name__ == "__main__":
    monte_carlo_simulation()
