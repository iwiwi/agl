#!/usr/bin/env python
# coding: utf-8,
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import statsmodels.api as sm
# import seaborn as sns
import json
import sys
import re
from matplotlib.ticker import LogLocator


def get_v(v, t):
    return v[0:len(t)]


def ba_time():
    coloring = np.array([0.00636888,
                         0.029815,
                         0.0768359,
                         0.386336,
                         1.57891,
                         7.57744,
                         36.6437,
                         154.188,
                         675.788,
                         2508.85,
                         12457.3,
                         59724.3])
    cbb = np.array([0.07580832,
                    0.4595,
                    4.486935,
                    44.851705,
                    403.659349,
                    4397.759203,
                    48481.64863])
    sketch = np.array([0.014270,
                       0.045620,
                       0.185428,
                       0.268524,
                       1.137688,
                       2.692887,
                       6.668913,
                       17.178226,
                       65.146818,
                       161.390212,
                       376.680491,
                       880.483400,
                       2011.067960,
                       6474.374790,
                       15390.866200,
                       36125.046400])
    burning = np.array([0.05513482,
                        54.13326206])
    memb = np.array(
        [0.003478288,
         0.008737798,
         0.03229518,
         0.1059241,
         0.2576659,
         2.7547131,
         17.162784,
         75.558687,
         240.69974,
         773.42738,
         3535.4736])
    return coloring, cbb, sketch, burning, memb


def ba_memory():

    coloring = np.array([1632,
                         1628,
                         1904,
                         2128,
                         2452,
                         3104,
                         4404,
                         6992,
                         12180,
                         22568,
                         43332])
    cbb = np.array([1628,
                    1892,
                    1984,
                    5048,
                    6444,
                    7144,
                    8536])
    sketch = np.array([2584,
                       3636,
                       5408,
                       9088,
                       16708,
                       31616,
                       67724,
                       135660,
                       272128,
                       549548,
                       1101596,
                       2206308,
                       4412232,
                       8802440,
                       17519932,
                       35101632])
    burning = np.array([2028,
                        2480])
    memb = np.array(
        [1896,
         2344,
         3344,
         11592,
         24404,
         113304,
         384312,
         1272624,
         7085720,
         25925588,
         89074480])
    return coloring, cbb, sketch, burning, memb


def flower_time():

    coloring = np.array([0.000787973,
                         0.00949097,
                         0.185456,
                         6.21543,
                         227.902,
                         7318.25])
    cbb = np.array([0.001662016,
                    0.01830321,
                    0.371506,
                    6.1214441,
                    121.821401,
                    1664.88904,
                    31869.9089])
    sketch = np.array([0.00111127,
                       0.013381958,
                       0.17328906,
                       1.78817971,
                       15.3899181,
                       192.399174,
                       982.27542,
                       8628.23664,
                       62138])
    burning = np.array([0.006328338,
                        199.001667])
    memb = np.array(
        [0.000493764,
         0.003545046,
         0.0410938,
         0.7632161,
         10.301574,
         246.40831,
         3959.0801])
    return coloring, cbb, sketch, burning, memb


def flower_memory():

    coloring = np.array([1632,
                         1900,
                         2080,
                         3148,
                         9780,
                         55424])
    cbb = np.array([1636,
                    1900,
                    1948,
                    2756,
                    5444,
                    16196,
                    49824])
    sketch = np.array([1896,
                       2740,
                       6668,
                       24460,
                       87540,
                       374556,
                       1495936,
                       5934680,
                       23660592])
    burning = np.array([1632,
                        2116])
    memb = np.array(
        [1636,
         1896,
         3952,
         32348,
         475696,
         7508300,
         97310308])
    return coloring, cbb, sketch, burning, memb


def flower_v():
    flower_v = np.array([125,
                         250,
                         500,
                         1000,
                         2000,
                         4000,
                         8000,
                         16000,
                         32000,
                         64000,
                         128000,
                         256000,
                         512000,
                         1024000,
                         2048000,
                         4096000,
                         8192000])
    return flower_v


def ba_v():
    ba_v = np.array([125,
                     250,
                     500,
                     1000,
                     2000,
                     4000,
                     8000,
                     16000,
                     32000,
                     64000,
                     128000,
                     256000,
                     512000,
                     1024000,
                     2048000,
                     4096000,
                     8192000])
    return ba_v

if __name__ == "__main__":

    plt.rcParams['font.size'] = 30
    plt.rcParams['lines.markersize'] = 15
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 0
    for x in xrange(0, 4):

        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)

        if x == 0:
            v = flower_v()
            coloring, cbb, sketch, burning, memb = flower_time()
            axes.set_ylabel("time (sec)")
            filename = "flower_time"
        elif x == 1:
            v = flower_v()
            coloring, cbb, sketch, burning, memb = flower_memory()
            axes.set_ylabel("memory (KB)")
            filename = "flower_memory"
        elif x == 2:
            v = ba_v()
            coloring, cbb, sketch, burning, memb = ba_time()
            axes.set_ylabel("time (sec)")
            filename = "ba_time"
        else:
            v = ba_v()
            coloring, cbb, sketch, burning, memb = ba_memory()
            axes.set_ylabel("memory (KB)")
            filename = "ba_memory"

        axes.plot(get_v(v, coloring), coloring, "bs-",
                  label="GC", markerfacecolor='w', markeredgecolor='b', markeredgewidth=2)
        axes.plot(get_v(v, cbb), cbb, "mv-", label="CBB")
        axes.plot(get_v(v, memb), memb, "k^--", label="MEMB")
        axes.plot(get_v(v, sketch), sketch, "ro-", label="Sketch")
        axes.plot(get_v(v, burning), burning, "gx-", label="MVB")

        axes.set_xscale("log")
        axes.set_yscale("log")

        handles, labels = axes.get_legend_handles_labels()
        # reverse the order
        labels = labels[::-1]
        handles = handles[::-1]
        mvb_label = labels[0]
        mvb_handle = handles[0]
        labels = labels[1:]
        handles = handles[1:]
        labels.append(mvb_label)
        handles.append(mvb_handle)


        axes.legend(handles, labels,
                    loc='best', fontsize=24, frameon=False)
        # axes.legend(loc='best', fontsize=24, frameon=False)

        axes.set_xlabel("number of vertices", fontsize=30, labelpad=0.1)
        if x == 0 or x == 2:
            axes.yaxis.set_major_locator(LogLocator(
                base=100.0, subs=[1.0], numdecs=4, numticks=15))
        [i.set_linewidth(2) for i in axes.spines.itervalues()]

        fig.tight_layout(pad=0.1)
        fig.savefig("scalability_" + filename + ".pdf")
