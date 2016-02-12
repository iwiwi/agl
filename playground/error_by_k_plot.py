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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.lines as mlines
from matplotlib.markers import MarkerStyle


def residual_by_k():
    k = np.array([16, 32, 64, 128, 256, 512, 1024,
                  ])
    residuals = {
        "(3, 0, 6)-SHM": np.array([1.01303547911, 1.01395349842, 1.01143661609, 1.00769589824, 1.00539704767, 1.00389170897, 1.0, ]),
        "(3, 3, 7)-flower": np.array([1.19903705818, 1.19060581606, 1.23245792876, 1.20965716564, 1.18548440066, 1.17346469622, 1.18185550082, ]),
        "(2, 0, 8)-SHM": np.array([1.03668893673, 1.05144652459, 1.04971979683, 1.05271868918, 1.04560075041, 1.04539523499, 1.03845709571, ]),
    }
    error_bar = {
        "(3, 0, 6)-SHM": np.array([0.0216637613746, 0.0335390685564, 0.020514163035, 0.0138057874779, 0.0122831496339, 0.00900290057152, 0.0, ]),
        "(3, 3, 7)-flower": np.array([0.228632301028, 0.198255216953, 0.257758024212, 0.252781424154, 0.214752879588, 0.188521640451, 0.229089666098, ]),
        "(2, 0, 8)-SHM": np.array([0.0570556230701, 0.0768943364425, 0.0677197822587, 0.0719465259731, 0.0707007343681, 0.061490729673, 0.0710289107696, ]),
    }
    filename = 'average_approximation_ratio_with_k'

    return residuals, error_bar, k, filename


def residual_by_p():
    p = np.array([0.125000, 0.250000, 0.500000, 1.000000, 2.000000, 4.000000, 8.000000,
                  ])
    residuals = {
        "(3, 0, 6)-SHM": np.array([1.0105565343, 1.01015447279, 1.01079628386, 1.00769589824, 1.00769589824, 1.00769589824, 1.00402298851, ]),
        "(3, 3, 7)-flower": np.array([1.22280614666, 1.21983349956, 1.20166263276, 1.20965716564, 1.18612479475, 1.22044334975, 1.21425287356, ]),
        "(2, 0, 8)-SHM": np.array([1.05126241841, 1.05453650188, 1.05550689904, 1.05271868918, 1.04906607577, 1.04272359807, 1.03709453088, ]),
    }
    error_bar = {
        "(3, 0, 6)-SHM": np.array([0.0170628162736, 0.0172689991427, 0.017231431309, 0.0138057874779, 0.0138057874779, 0.0138057874779, 0.0142059851659, ]),
        "(3, 3, 7)-flower": np.array([0.221603380899, 0.249952625884, 0.240899171061, 0.252781424154, 0.222364532626, 0.264948594074, 0.259088458563, ]),
        "(2, 0, 8)-SHM": np.array([0.0670298119502, 0.0692714620566, 0.0740803576669, 0.0719465259731, 0.0703512166304, 0.0686347829218, 0.0701311161932, ]),
    }
    filename = 'average_approximation_ratio_with_a'

    return residuals, error_bar, p, filename


if __name__ == "__main__":

    for x in xrange(0, 2):
        linestyles = {
            "(3, 0, 6)-SHM": "-",
            "(3, 3, 7)-flower": "-",
            "(2, 0, 8)-SHM": "-",
        }
        markers = {
            "(3, 0, 6)-SHM": 'o',
            "(3, 3, 7)-flower": "s",
            "(2, 0, 8)-SHM": "^",
        }
        colors = {
            "(3, 0, 6)-SHM": "b",
            "(3, 3, 7)-flower": "r",
            "(2, 0, 8)-SHM": "g",
        }
        markerfacecolors = {
            "(3, 0, 6)-SHM": "b",
            "(3, 3, 7)-flower": "r",
            "(2, 0, 8)-SHM": "w",
        }
        markeredgecolors = {
            "(3, 0, 6)-SHM": "k",
            "(3, 3, 7)-flower": "k",
            "(2, 0, 8)-SHM": "g",
        }
        displace = {
            "(3, 0, 6)-SHM": 1.0 * (2.0**0.1),
            "(3, 3, 7)-flower": 1.0 / (2.0**0.1),
            "(2, 0, 8)-SHM": 1.0,
        }
        # fig = plt.figure(figsize=(6, 3))
        plt.rcParams['font.size'] = 30
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 8
        plt.rcParams['lines.linewidth'] = 3
        fig = plt.figure(figsize=(8, 6))
        axes = fig.add_subplot(1, 1, 1)
        axes.set_xscale("log", basex=2)
        # axes.set_yscale("log")

        axes.set_ylim(0.9, 1.85)

        if x == 0:
            axes.set_xlabel("$\it{k}$", fontsize=40, labelpad=None)
            residuals, error_bar, k, name = residual_by_k()
            axes.set_xlim(16 / np.sqrt(2), 1024 * np.sqrt(2))
        else:
            axes.set_xlabel("$\it{" + u'\u03b1' + "}$",
                            fontsize=40, labelpad=None)
            residuals, error_bar, k, name = residual_by_p()
            axes.set_xlim(0.125 / np.sqrt(2), 8 * np.sqrt(2))

        keys = ["(3, 0, 6)-SHM", "(3, 3, 7)-flower", "(2, 0, 8)-SHM", ]
        lines = []
        labels = []
        for graph_name in keys:
            (_, caps, _) = axes.errorbar(k * displace[graph_name], residuals[graph_name],
                                         yerr=[
                error_bar[graph_name], error_bar[graph_name],
            ], ls=linestyles[graph_name],
                color=colors[graph_name],
                marker=markers[graph_name], label=graph_name,
                capsize=10, elinewidth=3, markersize=18,
                markerfacecolor=markerfacecolors[graph_name],
                markeredgecolor=markeredgecolors[graph_name],
                markeredgewidth=2)
            for cap in caps:
                cap.set_markeredgewidth(2)

            line = mlines.Line2D([], [], color=colors[graph_name], marker=markers[graph_name],
                                 markersize=18, markeredgewidth=2,
                                 markerfacecolor=markerfacecolors[graph_name],
                                 markeredgecolor=markeredgecolors[graph_name],
                                 label=graph_name)
            lines.append(line)
            labels.append(graph_name)

        axes.set_ylabel("approximation ratio",
                        fontsize=30, labelpad=None)

        axes.legend(lines, labels, loc='upper right',
                    fontsize=23, frameon=False)

        # axes.legend(loc='upper right', fontsize=23, frameon=False)

        axes.yaxis.set_major_locator(MultipleLocator(0.2))
        # axes.yaxis.set_label_coords(1.0, 1.0)
        [i.set_linewidth(2) for i in axes.spines.itervalues()]

        fig.tight_layout(pad=0.2)
        fig.savefig(name + ".pdf")
        fig.savefig(name + ".png")
        plt.close()
