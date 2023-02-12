#!/usr/bin/env python3

import numpy as np
from scipy.constants import hbar, Boltzmann as bmc
from scipy.special import zeta
import matplotlib.pyplot as plt
import pyvista as pv


def dist(x):
    temp = 1
    return x**2 / (np.exp(hbar * x / temp / bmc) - 1) / (2*zeta(3)*(bmc*temp/hbar)**3)

if __name__ == "__main__":
    a = np.genfromtxt("build/test_dsmc.out.csv", delimiter=",")
    pd = pv.PolyData(a[:,:3])
    pd.plot()

    plt.hist(a[:,3], bins=60, density=True)
    xs = np.linspace(0.01, 300e11, 10000)
    plt.plot(xs, dist(xs))
    plt.show()
    
    plt.hist2d(a[:,0], a[:,1])
    plt.gca().set_aspect("equal")
    plt.colorbar()
    plt.show()
    
    plt.hist2d(a[:,0], a[:,2])
    plt.gca().set_aspect("equal")
    plt.colorbar()
    plt.show()
   
    plt.hist2d(a[:,1], a[:,2])
    plt.gca().set_aspect("equal")
    plt.colorbar()
    plt.show()
