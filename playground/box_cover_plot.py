#!/usr/bin/env python
# coding: UTF-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import json
import sys
import re


def plot_data(jsonData):
    graph_name = jsonData['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    vertices = jsonData['graph_info'][0]['vertices']
    print graph_name
    boxSizes = jsonData['size']
    radiuses = jsonData['radius']
    if 'diameter' in jsonData:
        diameters = jsonData['diameter']
    name = jsonData['name']

    x = []
    y = []
    for i in range(0, len(boxSizes)):
        if radiuses[i] == 0:
            continue

        if 'diameter' in jsonData:
            x.append(diameters[i])
        else:
            x.append(radiuses[i] * 2)
        y.append(boxSizes[i])
        if boxSizes[i] == 1:
            break

    x = np.array(x)
    y = np.array(y)
    nsample = x.size
    lnx = np.log(x)
    lny = np.log(y)

    X = np.column_stack((np.repeat(1, nsample), lnx))
    model = sm.OLS(lny, X)
    results = model.fit()
    # print summary
    # print results.summary()
    # get estimated params
    a, b = results.params

    X2 = np.column_stack((np.repeat(1, nsample), x))
    model2 = sm.OLS(lny, X2)
    results2 = model2.fit()
    # print summary
    # print results2.summary()
    # get estimated params
    a2, b2 = results2.params
    str1 = results.summary().as_text()
    str2 = results2.summary().as_text()
    outstr = ""

    px = np.linspace(0.5, 1000, 10000)
    fb = -np.log(5) / np.log(3)
    fa = np.log(boxSizes[2] / (diameters[2]**fb))

    # show plot
    plt.plot(x, y, 'o', label=name)

    # plot carved lines
    # if name == "MEMB":
    plt.plot(px, (px**fb) * (np.exp(1)**fa), label="y=x^(log5/log3)" + "*e^" + str(fa))
    plt.plot(px, (px**b) * (np.exp(1)**a), label="y=x^" + str(b) + "*e^" + str(a))
    # plt.plot(px, (np.exp(1)**(b2 * px)) * (np.exp(1)**a2), label="y=e^(" + str(b2) + "x+" + str(a2) + ")")

    outstr += "y=x^" + str(b) + "*e^" + str(a) + "\n" + str1 + "\n\n"
    outstr += "y=e^(" + str(b2) + "x+" + str(a2) + ")" + "\n" + str2 + "\n\n"
    plt.xlim(xmin=0.5)
    plt.xlim(xmax=1000)
    plt.ylim(ymin=1)
    plt.ylim(ymax=10000000)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc='best')
    plt.savefig(graph_name + "_" + name + ".png")
    # plt.close()
    outf = open(graph_name + "_" + name + ".log", "w")
    outf.write(outstr)
    outf.close()

if __name__ == '__main__':
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()
        plot_data(json_data)
