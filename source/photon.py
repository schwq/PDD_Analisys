"""
The article for the dose-depth aproximation for photon energies, is:
Li X-J, Ye Y-C, Zhang Y-S, Wu J-M (2022) Empirical modeling of the percent depth dose for megavoltage photon beams. PLoS ONE
17(1): e0261042. https://doi.org/10.1371/journal.pone.0261042
"""

import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
import json

# Plus-Minus unicode 
pm = '\u00B1'

class PhotonEnergy(Enum):
    MV4 = 4
    MV6 = 6
    MV10 = 10
    MV18 = 18
    MVALL = 0


# Beam parameters (MV : [n, u])
PhotonsEnergiesParameters = {
    PhotonEnergy.MV4: [0.17, 0.0605],
    PhotonEnergy.MV6: [0.208, 0.0515],
    PhotonEnergy.MV10: [0.495, 0.0458],
    PhotonEnergy.MV18: [1.2, 0.0422],
}

# Get all stored data from real measures in photon therapy
def read_measured_data(show_energies):

    measured_data = []
    
    with open("./data/measured_pdd.json", 'r') as file: data = json.load(file)

    for key in PhotonsEnergiesParameters.keys():

        if show_energies[0] is not PhotonEnergy.MVALL: 
            if key not in show_energies: continue 

        beam = data[f"{key.value}mv"]
        depths = []
        doses = []
        
        for val in beam:
            depth = val['d']
            dose = val['v']
            depths.append(depth)
            doses.append(dose)
        
        measured_data.append((key, depths, doses))
        
        plt.plot(depths, doses, '+', label=f'Measured data {key.value}MV')

    return measured_data

# Plots a graph about a photon with given energy 
def plt_single_photon(depth, energy):
    n, u = PhotonsEnergiesParameters[energy]
    D, D_norm, dMax = empddmv_photon(depth, n, u)
    plt.plot(
        depth,
        D_norm,
        label=f"{energy.value} MV (n={n:.2f}, u={u:.3f}) dmax={depth[dMax]:.2f}cm",
    )
    return (energy, D, D_norm)

# Empirical modeling of the percent depth dose for megavoltage photon beams
def empddmv_photon(depth, n, u):
    d = depth / np.sqrt((depth ** 2) + n) * np.exp(-u * depth)
    dnorm = d = 100 * d / np.max(d)
    dmax = np.argmax(d)
    return d, dnorm, dmax

# Given the measured data and the theotical data, compare the error between them 
def compare_measured_vs_model(measured_data, theoretical_data, depth):
    for energy, depths, doses_measured in measured_data:
        
        _, _, dnorm = next((d for d in theoretical_data if d[0] == energy), None)
        doses_model_interp = np.interp(depths, depth, dnorm)
        doses_measured = np.array(doses_measured)
        diff = doses_measured - doses_model_interp
        mae = np.mean(np.abs(diff))

        print(f"Photon Aproximation: Mean absolute error for {energy.value}MV: {pm}{mae:.2f}")

#Plots all photon PDD graphs of the given energies, use MVALL alone for all energies
def plt_photons(energies: tuple[PhotonEnergy]):

    depths = np.linspace(0, 25, 500) 
    plt.figure(figsize=(9, 5))

    theoretical_data = []

    measured_data = read_measured_data(energies)

    if energies[0] == PhotonEnergy.MVALL:
        for MV, _ in PhotonsEnergiesParameters.items():
            theoretical_data.append(plt_single_photon(depths, MV))
    else:
        for energy in energies:
            if energy not in PhotonsEnergiesParameters:
                print(f"Photon Beam energy of {energy.value}MV not supported!")
                continue
            theoretical_data.append(plt_single_photon(depths, energy))

    compare_measured_vs_model(measured_data, theoretical_data, depths)

    plt.xlabel("Depth in water (cm)")
    plt.ylabel("Relative dose (%)")
    plt.title("Approximate PDD Curves for Different Photon Beam Energies")
    
    plt.grid(True, alpha=0.3)
    
    plt.ylim(0, 110)
    plt.xlim(0, 25)
    plt.legend()

    return depths, theoretical_data, measured_data
