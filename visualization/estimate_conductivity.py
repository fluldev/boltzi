#!/usr/bin/env python3

import argparse

import numpy as np
from scipy.constants import Boltzmann as bmc, hbar
from scipy.special import factorial, zeta
import matplotlib.pyplot as plt
import yaml

def branch_weight(v):
    return 1/v**3

def general_scatter_rate(omega, temp, v, fit_constant, k_pow, temp_pow, exp_const, debye_frequency):
    return fit_constant * omega**k_pow * temp**temp_pow * np.exp(-exp_const/temp) 

def dist(omega, temp, v):
    return hbar**2/bmc * omega**4  / temp**2 * np.exp(hbar*omega/bmc/temp) / (np.exp(hbar*omega/bmc/temp)-1)**2\
            * v**2/3 / (2*np.pi**2 * v**3)

def dist_n_scatter(omega, temp, v, fit_constant, k_pow, temp_pow, exp_const, debye_frequency):
    return hbar**2/bmc * omega**(4-k_pow)  / temp**2 * np.exp(hbar*omega/bmc/temp) / (np.exp(hbar*omega/bmc/temp)-1)**2\
            * v**2/3 / (2*np.pi**2 * v**3)\
            / (fit_constant * temp**temp_pow) * np.exp(exp_const/temp)

expfitparams = [
    [-0.23959032269460867, 0.33623052919856133],
    [-0.4903470666751711, 0.2972615358259651],
    [-0.6928498360286857, 0.27110901462808296],
    [-0.8646777924910959, 0.2511663890388728],
    [-1.0160019730937484, 0.23510156523730397],
    [-1.1545260510396398, 0.2219574233047828],
    [-1.2828630854900795, 0.21085455001207437],
    [-1.4034307792700464, 0.20132983821950692],
    [-1.5177563940014838, 0.193047103667941],
    [-1.626870723141303, 0.18575909008382566],
    [-1.7315076947085744, 0.17928051914075824],
    [-1.8322190469193815, 0.1734703599426384],
    [-1.929437705283456, 0.16821950444688086],
    [-2.0235147574865535, 0.1634422212871092]
]

#def conductivity_estimator(temp, v, fit_constant, k_pow, temp_pow, exp_const, debye_frequency):
#    x = hbar*debye_frequency/(bmc*temp)
#    b,c = expfitparams[-int(k_pow)+2] 
#    delta = np.tanh(b+c*x)
#        #return 1
#
#    return -1/fit_constant * hbar**2 / (bmc*temp**2)*temp**-temp_pow/(4*np.pi**2*v)*(bmc*temp/hbar)**(5-k_pow)\
#            * (factorial(4-k_pow)*zeta(4-k_pow)*delta-x**(4-k_pow)/(np.exp(x)-1)) * np.exp(exp_const/temp)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="boltzi_graphs")
    parser.add_argument("config", type=str)
    parser.add_argument("-t","--temp_range", nargs=2, type=float, default=(1,300))
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config=yaml.safe_load(f)

    weight_sum = 0
    weights = []
    parameters = []

    temperatures = np.linspace(args.temp_range[0], args.temp_range[1], 1000)

    for material in config["materials"]:
        for branch in list(material.values())[0]["branches"]:
            cur_parameters = []
            branch = list(branch.values())[0]
            v_g = float(branch["velocity"])
            weights.append(branch_weight(v_g))
            weight_sum+=branch_weight(v_g)
            if "debye_frequency" in branch:
                omega_debye = float(branch["debye_frequency"])
            else:
                omega_debye = float(branch["debye_temperature"]) * bmc / hbar
            for i,prc in enumerate(branch["two_phonon_processes"]):
                name = list(prc.keys())[0]
                prc = list(prc.values())[0]
                match name:
                    case "parametric":
                        cur_parameters.append((
                            v_g, 
                            float(prc["fit_constant"]) / v_g**float(prc["wavevector_power"]), 
                            float(prc["wavevector_power"]), 
                            float(prc["temperature_power"]),
                            0,
                            omega_debye
                        ))
                    case "parametric_exp":
                        cur_parameters.append((
                            v_g, 
                            float(prc["fit_constant"]) / v_g**float(prc["wavevector_power"]), 
                            float(prc["wavevector_power"]), 
                            float(prc["temperature_power"]),
                            float(prc["exp_constant"]),
                            omega_debye
                        ))
            
            for i,prc in enumerate(branch["three_phonon_processes"]):
                name = list(prc.keys())[0]
                prc = list(prc.values())[0]
                match name:
                    case "parametric":
                        cur_parameters.append((
                            v_g, 
                            float(prc["fit_constant"]) / v_g**float(prc["wavevector_power"]), 
                            float(prc["wavevector_power"]), 
                            float(prc["temperature_power"]),
                            0,
                            omega_debye
                        ))
                    case "parametric_exp":
                        cur_parameters.append((
                            v_g, 
                            float(prc["fit_constant"]) / v_g**float(prc["wavevector_power"]), 
                            float(prc["wavevector_power"]), 
                            float(prc["temperature_power"]),
                            float(prc["exp_constant"]),
                            omega_debye
                        ))
            parameters.append(cur_parameters)
            
    branch_cond = []
    for param in parameters: 
        #def conductivity_integrand(omega, temp):
        #    rate = np.zeros_like(omega)
        #    for cur_param in param:
        #        rate+=general_scatter_rate(omega, temp, *cur_param)
        #    scatter_times = 1/rate

        #    return dist(omega, temp, param[0][0]) * scatter_times
        
        n_iterations = 100000

        #xs = np.linspace(param[0][-1]*.01,param[0][-1], 1000)
        #plt.plot(xs,dist(xs, 100, param[0][0]))
        #plt.gca().twinx().plot(xs, general_scatter_rate(xs, 100, *param[0]))
        #plt.show()
        
        cur_cond = []
        for temp in temperatures:
            cur_cond.append(np.mean(
                #conductivity_integrand((.01+np.random.random(n_iterations)*.99)*param[0][-1], temp)
                dist_n_scatter((.01+np.random.random(n_iterations)*.99)*param[0][-1], temp, *param[1])
            ) * param[0][-1])
        branch_cond.append(cur_cond)
    branch_cond = np.sum(branch_cond, axis=0)
    plt.plot(temperatures, branch_cond) 
    plt.show()
