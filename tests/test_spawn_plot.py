#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    a = np.genfromtxt("build/test_spawn.out.csv", delimiter=",")
    plt.hist2d(a[:,0], a[:,1])
    plt.gca().set_aspect("equal")
    plt.colorbar()
    plt.show()

