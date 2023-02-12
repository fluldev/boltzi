#!/usr/bin/env python3

import numpy as np
from scipy.special import lambertw, zeta, factorial
from scipy.integrate import quad
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

n_max = 15
samples = 10000
relative_end_intercept = 8+8/3
relative_end_slope = -4/3
relative_end_min = 3

if __name__ == "__main__":
    def to_int(x, n):
        return x**n / (np.exp(x)-1)

    def inf_int(n):
        return factorial(n)*zeta(n+1)

    def term_max(n):
        return (lambertw(-n*np.exp(-n))+n).real

    def fitting_fn(x_d, a, b, c):
        return a*np.tanh(b+c*x_d)
    
    fit_tbl = []

    params = [1, -.2, .3]
    for n in range(2, n_max+1):
        x_ds = np.linspace(term_max(n)*.8, max(relative_end_intercept+relative_end_slope*n, relative_end_min) * term_max(n), samples)
        
        results = []

        for x_d in x_ds:
            v, *_ = quad(to_int, 0, x_d, args=(n,))
            results.append(v / inf_int(n))

        params, *_ = curve_fit(fitting_fn, x_ds, results, p0=params)

        fig, (ax0, ax1) = plt.subplots(1,2)

        xs = np.linspace(0, x_ds[-1]*1.1, 1000)
        ax0.plot(xs, to_int(xs, n)) 
        ax0.vlines([x_ds[0], x_ds[-1]], 0, to_int(term_max(n), n))
        ax1.plot(x_ds, results, ls="-", c="black")
        ax1.plot(x_ds, fitting_fn(x_ds, *params), ls="--", c="red")
        plt.show()

        fit_tbl.append(params)

    print("const double tanh_fitparams[][3] = {") 
    for line in fit_tbl:
        print(f"{{{line[0]}, {line[1]}, {line[2]}}},")
    print("};")






