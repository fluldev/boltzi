#!/usr/bin/env python3

import numpy as np
from scipy.constants import hbar, Boltzmann as bmc
from scipy.special import zeta
import matplotlib.pyplot as plt
import pyvista as pv


def dist(x):
    temp = 8
    return x**2 / (np.exp(hbar * x / temp / bmc) - 1) / (2*zeta(3)*(bmc*temp/hbar)**3)

if __name__ == "__main__":
    a = np.genfromtxt("build/test_lin_material.out.csv", delimiter=",")
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


    a = np.genfromtxt("build/test_lin_material_surf.out.csv", delimiter=",")
    a = a[::12]
    a[:, 3:6]=a[:,:3] + a[:, 3:6] / np.linalg.norm(a[:, 3:6],axis=-1)[:, np.newaxis]*.03
    a = a[:,:6].reshape(-1, 3)
    pd = pv.PolyData(a, lines = [[2, i*2, i*2+1] for i in range(len(a)//2)])
    plttr = pv.Plotter()
    plttr.add_mesh(pd)
    plttr.add_mesh(pv.Cube(center=(.5,.5,.5)))
    plttr.add_axes()
    plttr.show()

