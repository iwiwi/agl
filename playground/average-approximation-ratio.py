#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import statsmodels.api as sm
import matplotlib.pyplot as plt
import json
import sys
import re
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.lines as mlines
from matplotlib.markers import MarkerStyle

analytical = {
    "flower-223950-3-3": {4: 37326,
                          10: 6222,
                          28: 1038,
                          82: 174,
                          244: 30,
                          730: 6,
                          2188: 2},
    "shm-234376-4-2": {2: 46876,
                       8: 9376,
                       26: 1876,
                       80: 376,
                       242: 76,
                       728: 16,
                       2186: 4,
                       6560: 1},
    "shm-312501-5-2": {2: 62501,
                       8: 12501,
                       26: 2501,
                       80: 501,
                       242: 101,
                       728: 21,
                       2186: 5,
                       6560: 1},
    "shm-67229-5-3": {2: 9605,
                      8: 1373,
                      26: 197,
                      80: 29,
                      242: 5,
                      728: 1},
}


def residual_by_k(data):
    a = 1.0
    filename = 'average_approximation_ratio_with_k'

    residuals = {}
    error_bar = {}
    graph_keys = sorted(data)
    for graph in graph_keys:
        if graph not in residuals:
            residuals[graph] = np.array([])
        ks = np.array(sorted(data[graph], key=lambda x: int(x)))
        for k in ks:
            rs = sorted(data[graph][k][a])
            s = 0
            all4std = np.array([])
            for r in rs:
                s += data[graph][k][a][r].mean()
                all4std = np.append(all4std, data[graph][k][a][r])
            residuals[graph] = np.append(residuals[graph], (s / len(rs)))

        if graph not in error_bar:
            error_bar[graph] = np.array([])
        for k in ks:
            rs = sorted(data[graph][k][a])
            all4std = np.array([])
            for r in rs:
                all4std = np.append(all4std, data[graph][k][a][r])
            error_bar[graph] = np.append(error_bar[graph], all4std.std())
    return residuals, error_bar, ks, filename


def residual_by_a(data):
    k = 128
    filename = 'average_approximation_ratio_with_a'

    residuals = {}
    error_bar = {}
    graph_keys = sorted(data)
    for graph in graph_keys:
        if graph not in residuals:
            residuals[graph] = np.array([])
        alphas = np.array(sorted(data[graph][k], key=lambda x: float(x)))
        for a in alphas:
            rs = sorted(data[graph][k][a])
            s = 0
            all4std = np.array([])
            for r in rs:
                s += data[graph][k][a][r].mean()
                all4std = np.append(all4std, data[graph][k][a][r])
            residuals[graph] = np.append(residuals[graph], (s / len(rs)))

        if graph not in error_bar:
            error_bar[graph] = np.array([])
        for a in alphas:
            rs = sorted(data[graph][k][a])
            all4std = np.array([])
            for r in rs:
                all4std = np.append(all4std, data[graph][k][a][r])
            error_bar[graph] = np.append(error_bar[graph], all4std.std())
    return residuals, error_bar, alphas, filename


def plot_approximation_ratio(data):
    for x in xrange(0, 2):
        linestyles = {"(3, 0, 6)-SHM": "-",
                      "(3, 3, 7)-flower": "-", "(2, 0, 8)-SHM": "-", }
        markers = {"(3, 0, 6)-SHM": 'o', "(3, 3, 7)-flower": "s",
                   "(2, 0, 8)-SHM": "^", }
        colors = {"(3, 0, 6)-SHM": "b", "(3, 3, 7)-flower": "r",
                  "(2, 0, 8)-SHM": "g", }
        markerfacecolors = {"(3, 0, 6)-SHM": "b",
                            "(3, 3, 7)-flower": "r", "(2, 0, 8)-SHM": "w", }
        markeredgecolors = {"(3, 0, 6)-SHM": "k",
                            "(3, 3, 7)-flower": "k", "(2, 0, 8)-SHM": "g", }
        displace = {"(3, 0, 6)-SHM": 1.0 * (2.0**0.1),
                    "(3, 3, 7)-flower": 1.0 / (2.0**0.1), "(2, 0, 8)-SHM": 1.0, }
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
            residuals, error_bar, k, name = residual_by_k(data)
            axes.set_xlim(16 / np.sqrt(2), 1024 * np.sqrt(2))
        else:
            axes.set_xlabel("$\it{" + u'\u03b1' + "}$",
                            fontsize=40, labelpad=None)
            residuals, error_bar, k, name = residual_by_a(data)
            axes.set_xlim(0.125 / np.sqrt(2), 8 * np.sqrt(2))

        graph_namemap = {
            "(3, 0, 6)-SHM": "shm-67229-5-3",
            "(3, 3, 7)-flower": "flower-223950-3-3",
            "(2, 0, 8)-SHM": "shm-312501-5-2",
        }
        keys = ["(3, 0, 6)-SHM", "(3, 3, 7)-flower", "(2, 0, 8)-SHM", ]
        lines = []
        labels = []
        for graph_name in keys:
            print graph_name
            print graph_namemap[graph_name]
            (_, caps, _) = axes.errorbar(k * displace[graph_name],
                                         residuals[graph_namemap[graph_name]],
                                         yerr=[
                                             error_bar[
                                                 graph_namemap[graph_name]],
                                             error_bar[
                                                 graph_namemap[graph_name]],
            ],
                ls=linestyles[graph_name],
                color=colors[graph_name],
                marker=markers[graph_name],
                label=graph_name,
                capsize=10,
                elinewidth=3,
                markersize=18,
                markerfacecolor=markerfacecolors[
                                             graph_name],
                markeredgecolor=markeredgecolors[
                                             graph_name],
                markeredgewidth=2)
            for cap in caps:
                cap.set_markeredgewidth(2)

            line = mlines.Line2D([],
                                 [],
                                 color=colors[graph_name],
                                 marker=markers[graph_name],
                                 markersize=18,
                                 markeredgewidth=2,
                                 markerfacecolor=markerfacecolors[graph_name],
                                 markeredgecolor=markeredgecolors[graph_name],
                                 label=graph_name)
            lines.append(line)
            labels.append(graph_name)

        axes.set_ylabel("approximation ratio", fontsize=30, labelpad=None)

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


def xy_from_json(json_data):
    boxSizes = json_data['size']
    if 'radius' in json_data:
        radiuses = json_data['radius']
    if 'diameter' in json_data:
        diameters = json_data['diameter']
    x = []
    y = []
    for i in range(0, len(boxSizes)):
        if radiuses[i] == 0:
            continue
        if 'diameter' in json_data:
            x.append(diameters[i])
        else:
            x.append(radiuses[i] * 2)
        y.append(boxSizes[i])
        if boxSizes[i] == 1:
            break
    x = np.array(x)
    y = np.array(y)
    return x, y


if __name__ == "__main__":

    data = {}
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        px, py = xy_from_json(json_data)

        graph = json_data['graph_info'][0]['graph']
        if graph not in data:
            data[graph] = {}

        k = int(json_data['k'])
        if k not in data[graph]:
            data[graph][k] = {}

        a = float(json_data['upper_param'])
        if a not in data[graph][k]:
            data[graph][k][a] = {}

        for i in range(len(px)):
            if px[i] not in data[graph][k][a]:
                data[graph][k][a][px[i]] = []
            data[graph][k][a][px[i]].append(py[i])

    graph_keys = sorted(data)
    # for graph in graph_keys:
    #     sys.stdout.write(graph + "\n")
    #     ks = sorted(data[graph], key=lambda x: int(x))
    #     for k in ks:
    #         sys.stdout.write(str(k) + ",")
    #     sys.stdout.write("\n")
    #
    #     alphas = sorted(data[graph][k])
    #     for a in alphas:
    #         sys.stdout.write(str(a) + ",")
    #     sys.stdout.write("\n")
    #
    #     for a in alphas:
    #         for k in ks:
    #             rs = sorted(data[graph][k][a])
    #             s = 0
    #             all4std = np.array([])
    #             for r in rs:
    #                 data[graph][k][a][r] = np.array(data[graph][k][a][r])
    #                 # sys.stdout.write(str(data[graph][k][a][r].mean()) + ",")
    #                 data[graph][k][a][r] = (data[graph][k][a][r] / (analytical[graph][r] + 0.0))
    #                 s += data[graph][k][a][r].mean()
    #                 all4std = np.append(all4std, data[graph][k][a][r])
    #             sys.stdout.write(str(s / len(rs)) + ",")
    #     sys.stdout.write("\n")
    #
    #     for a in alphas:
    #         for k in ks:
    #             rs = sorted(data[graph][k][a])
    #             all4std = np.array([])
    #             for r in rs:
    #                 all4std = np.append(all4std, data[graph][k][a][r])
    #             sys.stdout.write(str(all4std.std()) + ",")
    #     sys.stdout.write("\n")
    # sys.stdout.write("\n")

    for graph in data:
        for k in data[graph]:
            for a in data[graph][k]:
                for r in data[graph][k][a]:
                    data[graph][k][a][r] = np.array(data[graph][k][a][r])
                    data[graph][k][a][r] = (
                        data[graph][k][a][r] / (analytical[graph][r] + 0.0))
    plot_approximation_ratio(data)
