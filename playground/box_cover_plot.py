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
    plt.plot(qx, theoreticalValue(beta, qx), label="$x^{" +
             re.sub(
        r'^(\-*[0-9]*\.[0-9]{0,3})[0-9]*$',
        r"\1",
        str(beta[1])
    )
        + "} e^{" + re.sub(
                 r'^(\-*[0-9]*\.[0-9]{0,3})[0-9]*$',
                 r"\1",
        str(beta[0])
    ) + "}$", color="#FF0000")


def plotExpo(beta):
    qx = np.linspace(0.5, 1000000, 1000000)
    plt.plot(qx, expoValue(beta, qx), label="$e^{ " +
             re.sub(r'^(\-*[0-9]*\.[0-9]{0,3})[0-9]*$',  r"\1", str(beta[1])) +
             " x + " +
             re.sub(r'^(\-*[0-9]*\.[0-9]{0,3})[0-9]*$',
                    r"\1", str(beta[0])) + "}$", linestyle="--",  color="#000000")


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


def saveFig(jsonData):
    plt.xlabel("$\it{l_{B}}$", fontsize=18)
    plt.ylabel("$\it{N_{B}}$", fontsize=18)
    graph_name = jsonData['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.@]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    print graph_name

    vertices = jsonData['graph_info'][0]['vertices']
    name = jsonData['name']

    plt.xlim(xmin=0.5)
    plt.xlim(xmax=100)
    plt.ylim(ymin=1)
    plt.ylim(ymax=1000000)
    plt.tight_layout(pad=0.2)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc='best', fontsize=15)
    plt.savefig(graph_name + "_" + str(vertices) + "_" + name + ".png")
    # plt.savefig(graph_name + "_" + str(vertices) + "_" + name + ".pdf")


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


def plotAnalytical(json_text):
    graph_name = json_text['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    boxSizes = json_text['size']
    if 'radius' in json_text:
        radiuses = json_text['radius']
    if 'diameter' in json_text:
        diameters = json_text['diameter']

    px = np.linspace(0.5, 1000, 10000)

    # plot analytical line
    if 'flower' in graph_name or 'shm' in graph_name:
        m = re.match(
            r"(?P<g>[a-z]+)\-(?P<vs>[0-9]+)\-(?P<u>[0-9]+)\-(?P<v>[0-9]+)", graph_name)
        if 'flower' in graph_name:
            fb = -np.log(int(m.group('u')) + int(m.group('v'))) / \
                np.log(int(m.group('u')))
        else:
            fb = -np.log(2 * int(m.group('v')) + 1) / np.log(3)
        n = int((len(boxSizes) - 1) / 2)
        if 'diameter' in json_text:
            fa = np.log(boxSizes[n] / (diameters[n]**fb))
        else:
            fa = np.log(boxSizes[n] / ((radiuses[n] * 2)**fb))
        plotLine(np.array([fa, fb]))


if __name__ == "__main__":
    for ai in range(1, len(sys.argv)):
        print sys.argv[ai]
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        px, py, pname = xy_from_json(json_data)
        # plt.plot(px, py, 'o', label=pname)
        # plt.plot(px, py, 'o', label=json_data['name'])
        graph_name = json_data['graph_info'][0]['graph'].replace(" ", "_")
        if '/' in graph_name:
            r = re.compile("/([a-zA-Z0-9_\-\.@]*)$")
            m = r.search(graph_name)
            graph_name = m.group(1)

        # px = px[4:]
        # py = py[4:]

        expoResult = scipy.optimize.leastsq(
            expoFit, linearExpo(px, py), args=(px, py))
        fracResult = scipy.optimize.leastsq(
            fitFunc, linearRegression(px, py), args=(px, py))

        plotLine(fracResult[0])
        plotExpo(expoResult[0])
        f = sum(fitFunc(fracResult[0], px, py)**2)
        e = sum(expoFit(expoResult[0], px, py)**2)
        plt.plot(px, py, 'o', label=graph_name +
                 " (" + str(np.log10(e / f)) + ")")
        # plotAnalytical(json_data)
        saveFig(json_data)
        plt.close()
