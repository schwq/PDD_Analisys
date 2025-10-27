"""Script to generate various photon beam of different MV values (MegaVolts)
The code uses a custom exponential function with some build-up, dose builds up as the electron that have been liberated
deposit their dose, reaching a dmax value. The dmax value, for small MV values, can be estimated by dmax = MV / 4 (cm)

The article for the dose-depth aproximation for photon energies, is:
Li X-J, Ye Y-C, Zhang Y-S, Wu J-M (2022) Empirical modeling of the percent depth dose for megavoltage photon beams. PLoS ONE
17(1): e0261042. https://doi.org/10.1371/journal.pone.0261042
"""

import numpy as np
import matplotlib.pyplot as plt
from enum import Enum


# All suported energies level for the photon beam
class PhotonEnergy(Enum):
    MV4 = 4
    MV6 = 6
    MV10 = 10
    MV18 = 18
    # If used, this should go alone
    MVALL = 0


# Beam parameters (MV : [n, u])
PhotonsEnergiesParameters = {
    PhotonEnergy.MV4: [0.17, 0.0605],
    PhotonEnergy.MV6: [0.208, 0.0515],
    PhotonEnergy.MV10: [0.495, 0.0458],
    PhotonEnergy.MV18: [1.2, 0.0422],
}

# Plots a single photon beam of a certain energy
def Plot_SinglePhotonBeam(z, energy):
    n, u = PhotonsEnergiesParameters[energy]
    D, D_norm, dMax = EMPDDMV_PhotonBeam(z, n, u)
    plt.plot(
        z,
        D_norm,
        label=f"{energy.value} MV (n={n:.2f}, u={u:.3f}) dmax={z[dMax]:.2f}cm",
    )
    return D_norm

#Empirical modeling of the percent depth dose for megavoltage photon beams
def EMPDDMV_PhotonBeam(depth, n, u):
    d = depth / np.sqrt((depth ** 2) + n) * np.exp(-u * depth)
    dnorm = d = 100 * d / np.max(d)
    dmax = np.argmax(d)
    return d, dnorm, dmax

# Plots all photon beams of the requested energies
def Plot_PhotonBeam(energies: tuple[PhotonEnergy]):

    z = np.linspace(0, 25, 500)  # depth in cm
    plt.figure(figsize=(9, 5))

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
