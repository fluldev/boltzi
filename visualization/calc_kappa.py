#!/usr/bin/env python3

import argparse
from collections import defaultdict

import numpy as np
import re
import os
import matplotlib.pyplot as plt

if __name__ == "__main__":
    plt.style.use("science")
    parser = argparse.ArgumentParser()
    parser.add_argument("current", type=str)
    parser.add_argument("gradient", type=str)
    parser.add_argument("--range", type=float, nargs=2, default=[])
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--loglog", action="store_true")
    parser.add_argument("--single", action="store_true")
    parser.add_argument("--csv", type=str, default="")
    parser.add_argument("--bugfix", type=float, default=5)

    args = parser.parse_args()

    gradient = eval(args.gradient)

    res = defaultdict(lambda: [])
    for file in filter(lambda x: re.match(args.current, x), os.listdir()):
        m = re.match(args.current, file)
        data = np.genfromtxt(file, delimiter=",")
        mask = np.ones_like(data[:,0], dtype=bool)
        if args.range:
            mask = (args.range[1]>=data[:,0]) & (data[:,0]>=args.range[0])

        if args.single:
            plt.scatter(data[mask, 0], data[mask, 1] / gradient / 100, c="red", marker="x")
            mean = np.mean(data[mask, 1] / gradient / 100)
            std = np.std(data[mask, 1] / gradient / 100)
            mi,ma = np.min(data[mask, 0]), np.max(data[mask, 0]) 
            plt.hlines([mean], ls="-", color="black", xmin=mi, xmax=ma)
            plt.hlines([mean+std,mean-std], ls="--", color="black", xmin=mi, xmax=ma)
            plt.xlabel("$t$ / $s$")
            plt.ylabel("$\\kappa$ / $Wcm^{-1}K^{-1}$")
            plt.show()
            data[:,1]/=gradient * 100
            np.savetxt(args.csv, data[mask], delimiter=",", newline="\n")
            exit(0)
        
        res[m.group(1)].append(np.abs(np.mean(data[mask, 1] / gradient / 100)))  # /100 == cm

    data = []
    errs = []
    temps = []
    
    for temp in sorted(res.keys()):
        mm = np.min(res[temp])
        # tempfix for bug
        while True:
            for ele in res[temp]:
                if ele > mm*args.bugfix:
                    res[temp].remove(ele)
                    break
            else:
                break

        temps.append(float(temp))
        data.append(np.mean(res[temp]))
        errs.append(np.std(res[temp]))

    data = np.array(data)
    errs = np.array(errs)

    if not args.csv:
        if args.log or args.loglog:
            plt.gca().set_yscale("log")
        if args.loglog:
            plt.gca().set_xscale("log")
        plt.errorbar(temps, data, errs, 
            fmt="x", 
            c="red", 
            capsize=5,
            ecolor="black"
        )
        plt.xlabel("$T$ / $K$")
        plt.ylabel("$\\kappa$ / $Wcm^{-1}K^{-1}$")
        plt.show()
    else:
        res = np.vstack([temps,data,errs]).T
        np.savetxt(args.csv, res, delimiter=",", newline="\n")
