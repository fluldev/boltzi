#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import re
import os

from collections import defaultdict


if __name__ == "__main__": 
    plt.style.use("science")
    parser = argparse.ArgumentParser()
    parser.add_argument("file_regex", type=str, help="Regular expression matching files with capture group 1 beeing the temperature.")
    parser.add_argument("--title", type=str, default="")
    parser.add_argument("--ylabel", type=str, default="")
    parser.add_argument("--expr", type=str, default="np.linalg.norm(x)")
    parser.add_argument("--csv", action="store_true")
    parser.add_argument("-v","--vectors", action="store_true", help="Use vector artifacts.")
    parser.add_argument("--log", action="store_true")
    args = parser.parse_args()

    expr = eval(f"lambda x: {args.expr}")

    data_agg = defaultdict(lambda: [])
    files = os.listdir()
    for file in filter(lambda f: re.match(args.file_regex, f), files):
        data = []
        with open(file, "r") as f:
            next(f)
            next(f)
            next(f)
            if args.vectors:
                for k,v1,*v in zip(f,f,f,f):
                    data.append([
                        *list(map(float, k.removesuffix("\n").split(","))), 
                        *list(map(float, v1.removesuffix("\n").replace("\t","").split(","))),
                        *list(map(lambda x: float(x.removesuffix("\n").removeprefix("\t")), v))
                    ])
                    next(f)
            else:
                for k,v1 in zip(f,f):
                    data.append([
                        *list(map(float, k.removesuffix("\n").split(","))), 
                        *list(map(float, v1.removesuffix("\n").split(","))),
                    ])

        entries = len(data)
        temp = float(re.match(args.file_regex, file)[1])
        data = np.array(data).reshape(entries, -1)
        data_sum = np.sum(np.apply_along_axis(expr, 1, data[:, 4:])*data[:,3])
        data_agg[temp].append(data_sum / np.sum(data[:,3]))

    res = []
    errs = []
    temps = sorted(data_agg.keys())
    for key in temps:
        res.append(np.mean(data_agg[key])) 
        errs.append(np.std(data_agg[key]))

    if not args.csv:
        fig, ax = plt.subplots(1,1)
        if args.log:
            ax.set_yscale("log")
        ax.errorbar(temps, res, errs, 
            fmt="x", 
            c="red", 
            capsize=5,
            ecolor="black"
        )
        ax.set_xlabel("$T$ / $K$")
        ax.set_ylabel(args.ylabel)
        ax.set_title(args.title)
        plt.show()
    else:
        for t, v in zip(temps, res):
            print(t,v,sep=",")
