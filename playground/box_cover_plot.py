# coding: UTF-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import json
import sys
import re

if __name__ == '__main__':
    log = open(sys.argv[1], 'r')
    jsonData = json.load(log)
    log.close()
    graph_name = jsonData['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    vertices = jsonData['graph_info'][0]['vertices']
    print graph_name
    for j in range(0, len(jsonData['algorithms'])):
        boxSizes = jsonData['algorithms'][j]['size']
        name = jsonData['algorithms'][j]['name']

        # calculate
        # fu = 2
        # fv = 2
        # fw = fu + fv
        # n = 1
        # nv = fw
        # while nv < vertices:
        #     n = n + 1
        #     nv = fw * nv - fw
        # d = fw
        # dx = []
        # dy = []
        # m = 1
        # while m <= n:
        #     d = fu * d + (fv - fu)
        #     box_t = (fw - 2) * (fw**(n - m)) / (fw - 1) - fw / (fw - 1)
        #     if box_t <= 0:
        #         break
        #     dx.append(d)
        #     dy.append(box_t)
        #     m = m + 1
        # dx = np.array(dx)
        # dy = np.array(dy)
        # plt.plot(dx, dy, 'o', label="rironchi")
        # print dx
        # print dy

        x = []
        y = []
        for i in range(0, len(boxSizes)):
            x.append(2 * i + 1)
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
        fu = 2
        fv = 2
        fb = -np.log(fu + fv) / np.log(fu)
        fa = np.log(vertices)
        plt.plot(px, (px**fb) * (np.exp(1)**fa), label="y=x^" + str(fb) + "*e^" + str(fa))

        # show plot
        plt.plot(x, y, 'o', label=name)
        plt.plot(px, (px**b) * (np.exp(1)**a), label="y=x^" + str(b) + "*e^" + str(a))
        outstr += "y=x^" + str(b) + "*e^" + str(a) + "\n" + str1 + "\n\n"
        plt.plot(px, (np.exp(1)**(b2 * px)) * (np.exp(1)**a2), label="y=e^(" + str(b2) + "x+" + str(a2) + ")")
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
