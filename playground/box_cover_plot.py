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
    boxSizes = jsonData['algorithms'][0]['size']
    name = jsonData['algorithms'][0]['name']

    x = []
    y = []
    for i in range(1, len(boxSizes)):
        x.append(2 * i)
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

    px = np.linspace(0.5, 100, 10000)
    fu = 1
    fv = 2
    fb = -np.log(fu + fv) / np.log(fu)
    fa = np.log(vertices)

    # show plot
    plt.plot(x, y, 'o', label=name)

    # plot carved lines
    # if name == "MEMB":
    # plt.plot(px, (px**fb) * (np.exp(1)**fa), label="y=x^" + str(fb) + "*e^" + str(fa))
    # plt.plot(px, (px**b) * (np.exp(1)**a), label="y=x^" + str(b) + "*e^" + str(a))
    # plt.plot(px, (np.exp(1)**(b2 * px)) * (np.exp(1)**a2), label="y=e^(" + str(b2) + "x+" + str(a2) + ")")

    outstr += "y=x^" + str(b) + "*e^" + str(a) + "\n" + str1 + "\n\n"
    outstr += "y=e^(" + str(b2) + "x+" + str(a2) + ")" + "\n" + str2 + "\n\n"
    plt.xlim(xmin=0.5)
    plt.xlim(xmax=100)
    plt.ylim(ymin=1)
    plt.ylim(ymax=10000000)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc='best')
    plt.savefig(graph_name + "_" + name + ".png")
    plt.close()
    outf = open(graph_name + "_" + name + ".log", "w")
    outf.write(outstr)
    outf.close()

if __name__ == '__main__':
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()
        plot_data(json_data)
