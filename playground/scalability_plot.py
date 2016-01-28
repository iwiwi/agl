#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import statsmodels.api as sm
# import seaborn as sns
import json
import sys
import re


def plot(v, t, label):
    plt.plot(v[0:len(t)], t, label=label)

if __name__ == "__main__":
    v = np.array([0, 42, 287, 2002, 14007, 98042, 686287])

    coloring_t = np.array([0, 0.001055, 0.0404019, 3.50213, 862.684])
    sketch_t = np.array([0, 0.000427008,  0.031154182, 1.33151496,
                         29.5330592,  598.807429,  9605.42106])
    burning_t = np.array([0, 0.01176904])
    memb_t = np.array(
        [0, 0.000296832, 0.00783086,  0.42593, 21.814944,   1768.2418])

    plot(v, sketch_t, "Sketch")
    plot(v, memb_t, "MEMB")
    plot(v, burning_t, "Burning")
    plot(v, coloring_t, "Coloring")

    plt.legend(loc='best', fontsize=18)

    # plt.xscale("log")
    # plt.yscale("log")
    plt.xlabel("vertices")
    plt.ylabel("time (sec)")
    plt.savefig("scalability.png")
