import numpy as np
from scipy import signal, integrate
import matplotlib.pyplot as plt

p = 1.5
alpha = 1.9e-3 

ap = alpha ** (1/p)
dp = 1 - (1/p)

def DBP(d):
    return np.piecewise(d, [d<=0, d>0], [0, lambda d: 1 / (p*ap*d**dp)])

#For protons with energies between 10 and 200 MeV one finds p ≈ 1.8 (Evans 1982, Raju 1980). 
def R(E0): 
    if E0 < 10 or E0 > 200: print("Not supported energy ", E0)
    return alpha * (E0 ** p)

def E(d, R):
    return ((R-d)/alpha) ** (1/p)

#W (R) for the Bragg peaks such that the superposition results in a flat SOBP of height D0 within an interval [da , db] 
def W(R, d0, da, db):
    t1 = p * np.sin(np.pi / p) * ap 
    t2 = np.pi * (db - R)**(1/p)
    return np.piecewise(R, [R < da or R > db, da <= R < db], [0, lambda R: d0 * (t1/t2)])

def D_SOBD(d, d0, da, db):
    r = (da - d)/(db - da)
    rhat = r**(1/3)
    bp = lambda d: d0*(0.75+np.sqrt(3)/(4*np.pi)*np.log((1+rhat)**2/(1-rhat+rhat**2))-3.0/(2.0*np.pi)*np.arctan((2*rhat-1)/np.sqrt(3))) 
    return np.piecewise(d, [(0<=d)and(d<da), (da<=d)and(d<=db)], [bp, d0, 0])

da = 10
db = 15 
d0 = 1 


r = R(100) # Mev 
z = np.linspace(0, 25, 500)
plt.figure(figsize=(9,5))

plt.xlabel("Depth in water (cm)")
plt.ylabel("Relative dose (%)")
plt.title("Approximate PDD Curves for Different Photon Beam Energies")
plt.grid(True, alpha=0.3)
plt.ylim(0, 110)
plt.xlim(0, 25)
plt.legend()