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


def plot(v, t, symbol, label):
    axes.plot(v[0:len(t)], t, symbol, label=label)

if __name__ == "__main__":
    v = np.array([4, 12, 44, 172, 684, 2732, 10924,
                  43692, 174764, 699052, 2796204])

    coloring_t = np.array([2.60E-05, 0.000100136, 0.000658989,
                           0.0115712, 0.195711, 4.71543, 171.553, 5726.21])
    sketch_t = np.array([8.30E-05, 0.0001959799, 0.001171112, 0.011091946, 0.16714579,
                         1.77581678, 18.6191851, 119.248128, 950.874183, 7624.6712, 67822.5998])
    burning_t = np.array([4.29E-05, 0.0002357962, 0.006925825, 148.7523158])
    memb_t = np.array(
        [4.51E-05, 8.99E-05, 0.000389099, 0.002979043, 0.03322981, 0.5393008, 8.349336, 145.91113, 4402.5736])

    # fig = plt.figure(figsize=(6, 2))
    fig = plt.figure()
    axes = fig.add_subplot(1, 1, 1)

    plot(v, sketch_t, "ro-", "Sketch")
    plot(v, memb_t, "k^-", "MEMB")
    plot(v, burning_t, "gx-", "Burning")
    plot(v, coloring_t, "bs-", "Coloring")

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.legend(loc='best', fontsize=18)

    axes.set_xlabel("vertices")
    axes.set_ylabel("time (sec)")

    fig.tight_layout()
    fig.savefig("scalability.png")
