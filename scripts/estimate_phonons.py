#!/usr/bin/env python3

import argparse

import numpy as np
from scipy.constants import hbar, Boltzmann as bmc
from scipy.special import zeta, factorial

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

def integral_dist(n, vel, debye_temperature, temperature):
    x= debye_temperature / temperature
    a,b,c = fit_params[n-2]
    return factorial(n)*zeta(n+1)*a*np.tanh(b+c*x)\
           / (2*np.pi**2 * vel**3)\
           * (bmc*temperature / hbar)**(n+1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("debye_temperature", type=str)
    parser.add_argument("velocity", type=str)
    parser.add_argument("temperature", type=str)
    parser.add_argument("--rel", type=str, default=-1)
    parser.add_argument("--mul", type=str, default=1)
    parser.add_argument("--dev", type=str, default=0)

    args = parser.parse_args()

    res = eval(args.mul)*integral_dist(2, eval(args.velocity), eval(args.debye_temperature), eval(args.temperature))
    if args.dev:
        res -= eval(args.mul)*integral_dist(2, eval(args.velocity), eval(args.debye_temperature), eval(args.temperature)+eval(args.dev))
        res = np.abs(res)
    if args.rel:
        d = np.abs(integral_dist(2, eval(args.velocity), eval(args.debye_temperature), eval(args.rel))\
                   - integral_dist(2, eval(args.velocity), eval(args.debye_temperature), eval(args.rel) + eval(args.dev)))
        res /= d
    print(res)

