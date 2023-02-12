#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt 
from scipy.constants import Boltzmann as bmc, hbar
from scipy.special import zeta
import argparse

fit_params = [
    [1, -0.23959032269460867, 0.33623052919856133],
    [1, -0.4903470666751711, 0.2972615358259651],
    [1, -0.6928498360286857, 0.27110901462808296],
    [1, -0.8646777924910959, 0.2511663890388728],
    [1, -1.0160019730937484, 0.23510156523730397],
    [1, -1.1545260510396398, 0.2219574233047828],
    [1, -1.2828630854900795, 0.21085455001207437],
    [1, -1.4034307792700464, 0.20132983821950692],
    [1, -1.5177563940014838, 0.193047103667941],
    [1, -1.626870723141303, 0.18575909008382566],
    [1, -1.7315076947085744, 0.17928051914075824],
    [1, -1.8322190469193815, 0.1734703599426384],
    [1, -1.929437705283456, 0.16821950444688086],
    [1, -2.0235147574865535, 0.16344222128710925]
]


def kappa(temp, length, vel, omega_debye):
    a,b,c = fit_params[1]
    return 3*zeta(4)*bmc**4/hbar**3/(np.pi*vel)**2 * temp**3 * length\
            * a*np.tanh(b+c*hbar*omega_debye/(bmc*temp))/100

def kappa_corr(temp, temp2, length, vel, omega_debye):
    a,b,c = fit_params[1]
    return 3/4*zeta(4)*bmc**4/hbar**3/(np.pi*vel)**2 * (temp**4-temp2**4)/abs(temp-temp2) * length\
            * a*np.tanh(b+c*hbar*omega_debye/(bmc*(temp+temp2)/2))/100


def total_emission(temp, vel, omega_debye):
    a,b,c = fit_params[0]
    return 3/(np.pi**2 * vel**2)*zeta(4) * hbar*(bmc*temp/hbar)**4 / 4\
            * a*np.tanh(b+c*hbar*omega_debye/(bmc*temp))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("temperature", type=float)
    parser.add_argument("length", type=float)
    parser.add_argument("debye_temperature", type=float)
    parser.add_argument("velocity", type=float, nargs="+")
    parser.add_argument("--dev", type=float, default=0)

    args = parser.parse_args()

    debye_freq = bmc/hbar*args.debye_temperature

    print("kappa_max = ", sum([
        kappa(args.temperature, args.length, vel, debye_freq)
        for vel in args.velocity 
    ]).__format__(".3e"), "W*cm^-1*K^-1")
    if args.dev>0:
        print("kappa_corr = ", sum([
            kappa_corr(args.temperature, args.dev, args.length, vel, debye_freq)
            for vel in args.velocity 
        ]).__format__(".3e"), "W*cm^-1*K^-1")
    print("total emission = ", sum([
        total_emission(args.temperature, vel, debye_freq)
        - total_emission(args.dev, vel, debye_freq) 
        for vel in args.velocity
    ]).__format__(".3e"), "m^-2*s^-1")
