#!/usr/bin/env python3

import yaml
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar, Boltzmann as bmc

import argparse

def distr(k, temp, v):
    return (v*k)**2 / (2*np.pi**2*v**3) / (np.exp(hbar*(k*v)/(bmc*temp))-1)

def parametric(k, temp, fit_const, k_pow, temp_pow):
    return fit_const * k**k_pow * temp**temp_pow

def parametric_exp(k, temp, fit_const, k_pow, temp_pow, exp_const):
    return fit_const * k**k_pow * temp**temp_pow * np.exp(-exp_const/temp)

if __name__ == "__main__":
    plt.style.use("science")
    parser = argparse.ArgumentParser(prog="boltzi_graphs")
    parser.add_argument("config", type=str)
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config=yaml.safe_load(f)

    dev_temp = eval(f"lambda x,y,z: {str(config['deviational_temperature']).replace('^','**')}")
    temp = dev_temp(0,0,0)
    time_step = float(config["time_step"])
    
    for material in config["materials"]:
        for branch in list(material.values())[0]["branches"]:
            branch = list(branch.values())[0]
            v_g = float(branch["velocity"])
            k_debye = float(branch["debye_frequency"]) / v_g
            ks = np.linspace(k_debye*.001, k_debye, 1000)
            plt.plot(ks, distr(ks, temp, v_g), label="n($|\\vec k|$)", c="black", ls="--")
            plt.ylabel("Phononendichte / $m^{-3}$")
            plt.xlabel("$\\omega$ / $s^{-1}$")
            ax = plt.gca().twinx()
            ax.set_ylabel("$\\Delta t / \\tau(\\omega,T)$")
            for i,prc in enumerate(branch["two_phonon_processes"]):
                name = list(prc.keys())[0]
                prc = list(prc.values())[0]
                match name:
                    case "parametric":
                        ax.plot(ks, time_step*parametric(
                            ks, 
                            temp, 
                            float(prc["fit_constant"]),
                            float(prc["wavevector_power"]),
                            float(prc["temperature_power"]),

                        ), label=f"two-phonon parametric {i+1}")
                    case "parametric_exp":
                        ax.plot(ks, time_step*parametric_exp(
                            ks, 
                            temp, 
                            float(prc["fit_constant"]),
                            float(prc["wavevector_power"]),
                            float(prc["temperature_power"]),
                            float(prc["exp_constant"])
                        ), label=f"two-phonon parametric_exp {i+1}")

            for i,prc in enumerate(branch["three_phonon_processes"]):
                name = list(prc.keys())[0]
                prc = list(prc.values())[0]
                match name:
                    case "parametric":
                        ax.plot(ks, time_step*parametric(
                            ks, 
                            temp, 
                            float(prc["fit_constant"]),
                            float(prc["wavevector_power"]),
                            float(prc["temperature_power"]),

                        ), label=f"three-phonon parametric {i+1}")
                    case "parametric_exp":
                        ax.plot(ks, time_step*parametric_exp(
                            ks, 
                            temp, 
                            float(prc["fit_constant"]),
                            float(prc["wavevector_power"]),
                            float(prc["temperature_power"]),
                            float(prc["exp_constant"])
                        ), label=f"three-phonon parametric_exp {i+1}")
            plt.legend()
            plt.show()


