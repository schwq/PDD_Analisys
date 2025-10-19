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
    return np.piecewise(R,[(da<=R) & (R<db)],[lambda R: d0*p*np.sin(np.pi/p)*ap/(np.pi*(db-R)**(1.0/p)), 0])

def D_SOBD(d, d0, da, db):
    return np.piecewise(d, [(0<=d)&(d<da), (da<=d)&(d<=db)], [lambda d: _D_SOBD_Fun(d, d0, da, db), d0, 0])

def _D_SOBD_Fun(d, d0, da, db):
    r = (da - d)/(db - da)
    rhat = r**(1/3)
    return d0*(0.75+np.sqrt(3)/(4*np.pi)*np.log((1+rhat)**2/(1-rhat+rhat**2))-3.0/(2.0*np.pi)*np.arctan((2*rhat-1)/np.sqrt(3))) 
    

def add_legend():
    plt.grid(True, alpha=0.3)
    plt.grid('on')
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, loc = 'best')

def set_lim():
    min, max = plt.gca().get_ylim()
    plt.gca().set_ylim(min, max)

def Plot_Proton_SOBP(da, db, d0):

    de = np.linspace(0, 25, 500)
    plt.figure(figsize=(9,5))

    bx = DBP(db - de)
    bx /= bx[0] 
    bx *= 10
    plt.plot(de, bx, 'r', label='Exact')
    plt.title("Bragg peak")
    plt.xlabel('$d$ (cm)')
    plt.ylabel("Relative Dose (%)")
    set_lim()
    add_legend()

    plt.figure(figsize=(9,5))
    wx = W(de, d0, da, db) 
    plt.plot(de, wx, 'b', label='Weighting')
    plt.xlabel('$R$ (cm)')
    plt.ylabel(r'$W\left(R\right)$')
    plt.title('Weighting function')
    set_lim()
    add_legend()

    plt.figure(figsize=(9,5))
    sx = D_SOBD(de, d0, da, db)
    sx *= 100
    plt.plot(de, sx, 'g', label='Exact')
    plt.xlabel('$d$ (cm)')
    plt.ylabel('Relative Dose (%)')
    plt.title('(SOBP) Spread-out Bragg peak')
    set_lim()
    add_legend()
    
    return de, bx, wx, sx
