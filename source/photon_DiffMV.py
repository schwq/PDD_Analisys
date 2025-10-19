""" Script to generate various photon beam of different MV values (MegaVolts)
The code uses a custom exponential function with some build-up, dose builds up as the electron that have been liberated 
deposit their dose, reaching a dmax value. The dmax value, for small MV values, can be estimated by dmax = MV / 4 (cm)"""

import numpy as np
import matplotlib.pyplot as plt
from enum import Enum

# All suported energies level for the photon beam 
class PhotonEnergy(Enum):
    MV2 = 2
    MV4 = 4
    MV6 = 6
    MV8 = 8
    MV10 = 10
    MV12 = 12
    MV15 = 15
    # If used, this should go alone 
    MVALL = 0

# Beam parameters (MV : [a, b])
PhotonsEnergiesParameters = {
    PhotonEnergy.MV2:  [3.40, 0.08],
    PhotonEnergy.MV4:  [2.80, 0.07],
    PhotonEnergy.MV6:  [2.20, 0.06],
    PhotonEnergy.MV8:  [1.80, 0.055],
    PhotonEnergy.MV10: [1.60, 0.045],
    PhotonEnergy.MV12: [1.40, 0.042],
    PhotonEnergy.MV15: [1.20, 0.040],
}

def _Photon_ExpFunc(a, b, z):
    D = (1 - np.exp(-a*z)) * np.exp(-b*z)
    D_norm = 100 * D / D.max()  # normalize to 100%
    dMax = np.argmax(D)
    return D, D_norm, dMax

# Plots a single photon beam of a certain energy
def Plot_SinglePhotonBeam(z, energy):
    a, b = PhotonsEnergiesParameters[energy]
    D, D_norm, dMax = _Photon_ExpFunc(a, b, z)
    plt.plot(z, D_norm, label=f"{energy.value} MV (a={a:.2f}, b={b:.3f}) dmax={z[dMax]:.2f}cm")
    return D_norm

# Plots all photon beams of the requested energies
def Plot_PhotonBeam(energies : tuple[PhotonEnergy]):

    z = np.linspace(0, 25, 500)  # depth in cm
    plt.figure(figsize=(9,5))
    
    

    en = []

    if energies[0] == PhotonEnergy.MVALL:
        for MV, _ in PhotonsEnergiesParameters.items():
            en.append(Plot_SinglePhotonBeam(z, MV))
    else:
        for energy in energies:
            if energy not in PhotonsEnergiesParameters:
                print(f"Photon Beam energy of {energy.value}MV not supported!")
                continue 
            en.append(Plot_SinglePhotonBeam(z, energy))

    plt.xlabel("Depth in water (cm)")
    plt.ylabel("Relative dose (%)")
    plt.title("Approximate PDD Curves for Different Photon Beam Energies")
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 110)
    plt.xlim(0, 25)
    plt.legend()
    return z, en
