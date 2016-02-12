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
from matplotlib.ticker import LogLocator


def xy_from_json(json_data):
    boxSizes = json_data['size']
    if 'radius' in json_data:
        radiuses = json_data['radius']
    if 'diameter' in json_data:
        diameters = json_data['diameter']
    if len(boxSizes) <= 1:
        print "very little"
        return
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
    name = json_data['name']
    return x, y, name


def plotLine(beta):
    qx = np.linspace(0.5, 1000000, 1000000)
    axes.plot(qx, theoreticalValue(beta, qx),
              label="$\propto\  x^{" + "{0:.1f}".format(beta[1]) + "}$",
              color="#FF0000")


def plotExpo(beta):
    qx = np.linspace(0.5, 1000000, 1000000)
    axes.plot(qx, expoValue(beta, qx),
              label="$\propto\  e^{ " + "{0:.1f}".format(beta[1]) + " x}$",
              linestyle="--",  color="#000000")


def linearRegression(x, y):
    X = np.column_stack((np.repeat(1, x.size), np.log(x)))
    model = sm.OLS(np.log(y), X)
    results = model.fit()
    a, b = results.params
    return np.array([a, b])


def linearExpo(x, y):
    X = np.column_stack((np.repeat(1, x.size), x))
    model = sm.OLS(np.log(y), X)
    results = model.fit()
    a, b = results.params
    return np.array([a, b])


def getName(jsonData):
    graph_name = jsonData['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\@]*)\.agl$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    return graph_name


def theoreticalValue(beta, x):
    a = beta[0]
    b = beta[1]
    return (x**b) * (np.exp(1)**a)


def fitFunc(beta, x, y):
    residual = y - theoreticalValue(beta, x)
    return residual


def expoValue(beta, x):
    a = beta[0]
    b = beta[1]
    return (np.exp(1)**(b * x)) * (np.exp(1)**a)


def expoFit(beta, x, y):
    residual = y - expoValue(beta, x)
    return residual


if __name__ == "__main__":
    for ai in range(1, len(sys.argv)):
        print sys.argv[ai]
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        px, py, pname = xy_from_json(json_data)
        # plt.plot(px, py, 'o', label=pname)
        plt.rcParams['font.size'] = 20
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 8
        plt.rcParams['lines.linewidth'] = 4
        # fig = plt.figure(figsize=(10, 6))
        labelfontsize = 27

        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
        axes.set_xlabel("$\it{\ell}$", labelpad=5, fontsize=labelfontsize)
        axes.set_ylabel("$\it{b(\ell)}$", fontsize=labelfontsize)
        [i.set_linewidth(2) for i in axes.spines.itervalues()]

        tmpx = px
        tmpy = py
        axes.set_xlim(1, 10000)

        name = getName(json_data)

        if name == "in-2004":
            px = px[3:]
            py = py[3:]
            axes.set_xlim(1, 100)
        if name == "indochina-2004":
            px = px[4:]
            py = py[4:]
            axes.set_xlim(1, 1000)

        expoResult = scipy.optimize.leastsq(
            expoFit, linearExpo(px, py), args=(px, py))
        fracResult = scipy.optimize.leastsq(
            fitFunc, linearRegression(px, py), args=(px, py))

        plotLine(fracResult[0])
        plotExpo(expoResult[0])
        axes.set_xscale("log")
        axes.set_yscale("log")
        # axes.yaxis.set_major_locator(LogLocator(
        #     base=100.0, subs=[1.0], numdecs=4, numticks=15))
        axes.set_ylim(1, 10**6)
        plt.tight_layout(pad=0.2)
        axes.legend(loc='upper right',  frameon=False, fontsize=labelfontsize)
        axes.plot(tmpx, tmpy, 'o',
                  markersize=9,
                  markerfacecolor='w',
                  markeredgecolor='b',
                  markeredgewidth=2)
        plt.savefig(name + ".pdf")
        plt.savefig(name + ".png")
        plt.close()
