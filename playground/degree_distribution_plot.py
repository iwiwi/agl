#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import json
import sys
import re
from matplotlib.ticker import LogLocator


if __name__ == "__main__":

    log = open(sys.argv[1], 'r')
    json_data = json.load(log)
    log.close()

    plt.rcParams['font.size'] = 20
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams['lines.linewidth'] = 2
    # fig = plt.figure(figsize=(10, 6))
    labelfontsize = 27

    for rad in json_data["distribution"][0]:
        if(int(rad) >= 20):
            continue
        m = {}
        for deg in json_data["distribution"][0][rad][0]:
            tx = int(deg)
            ty = json_data["distribution"][0][rad][0][deg]
            m[tx] = ty
        x = []
        y = []
        for k, v in sorted(m.items()):
            x.append(k + 1)
            y.append(v)

        for i in reversed(xrange(0, len(y) - 1)):
            y[i] += y[i + 1]

        # fig = plt.figure()
        # axes = fig.add_subplot(1, 1, 1)
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel("$p(k)$", labelpad=5, fontsize=labelfontsize)
        plt.xlabel("$k$", fontsize=labelfontsize)

        plt.tight_layout(pad=0.2)
        plt.xlim(1.0, 10**6)
        plt.ylim(1.0, 10**7)
        # axes.legend(loc='upper right',  frameon=False, fontsize=labelfontsize)
        plt.plot(x, y,

                 )
    plt.savefig("rad_" + rad + ".png")
    plt.close()
