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


def residual_by_k():
    k = np.array([
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        512,
    ])
    residuals = {
        # "flower-223950-2-4": np.array([43870785, 35619544, 31361885, 28768126,
        # 28668242, 28635201, 28609758, 28604738, ]),
        # "(3,3)-flower-223950":
        "(3,3)-flower": np.array([90810, 53379, 14800, 10612, 6823, 1860, 941, 1448, ]),
        # "flower-292970-2-3": np.array([11339405, 13270296, 14289872, 15038191,
        # 15088651, 14992996, 14956696, 14895667, ]),
        # "(4,2)-shm-234376"
        "(4,2)-shm": np.array([48343, 19292, 4576, 3400, 2105, 1906, 965, 1469, ]),
        # "(5,2)-shm-312501"
        "(5,2)-shm": np.array([65755, 21516, 7307, 7199, 3765, 1663, 1394, 2954, ]),

    }
    filename = 'squared_residuals_with_k'

    return residuals, k, filename


def residual_by_p():
    p = np.array([
        0.0625,
        0.125,
        0.25,
        0.5,
        1.0,
        2.0,
        4.0,
    ])
    residuals = {
        # "flower-223950-2-4": np.array([36685465, 32964472, 29519425, 28684532,
        # 28635201, 28604186, 28606774, ]),
        # "(3,3)-flower-223950":
        "(3,3)-flower":   np.array([92300, 35542, 13531, 5829, 1860, 1239, 466, ]),
        # "flower-292970-2-3": np.array([12990132, 14764276, 15122801, 15034396,
        #                                14992996, 14936604, 14889567, ]),
        # "(2,2)-flower-43692": np.array([3774, 2073, 1658, 939, 416, 416, 68, ]),
        # "(4,2)-shm-234376"
        "(4,2)-shm":   np.array([25529, 8870, 7814, 2628, 1906, 1521, 779, ]),
        # "(5,2)-shm-312501"
        "(5,2)-shm":   np.array([34354, 18003, 13753, 3783, 1663, 1409, 657, ]),

    }
    filename = 'squared_residuals_with_p'

    return residuals, p, filename


if __name__ == "__main__":
    markers = {
        "(3,3)-flower":   "ko-",
        "(4,2)-shm":  "r^-",
        "(5,2)-shm":  "bx-",
    }

    # residuals, k, name = residual_by_k()
    residuals, k, name = residual_by_p()

    fig = plt.figure(figsize=(6, 3))
    # fig = plt.figure()
    axes = fig.add_subplot(1, 1, 1)
    axes.set_xscale("log", basex=2)
    axes.set_yscale("log")
    axes.set_ylim(300, 100000)

    for graph_name in residuals.keys():
        axes.plot(k, residuals[graph_name],
                  markers[graph_name], label=graph_name)

    # axes.set_xlabel("k")
    axes.set_xlabel("p")

    axes.set_ylabel("squared residuals")
    axes.legend(loc='best', fontsize=10)

    fig.tight_layout()
    fig.savefig(name + ".pdf")
