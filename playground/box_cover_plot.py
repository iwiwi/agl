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
    graph_name = jsonData['graph_info'][0]['graph']
    r = re.compile("/([a-zA-Z0-9]*\.[a-zA-Z0-9]*)$")
    m = r.search(graph_name)
    graph_name = m.group(1)
    for j in range(0, len(jsonData['algorithms'])):
        boxSizes = jsonData['algorithms'][j]['size']
        name = jsonData['algorithms'][j]['name'].replace(' ', '_')

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

        # show plot
        plt.plot(x, y, 'o', label=name)
        px = np.linspace(0.5, 100, 10000)
        plt.plot(px, (px**b) * (np.exp(1)**a), label="y=x^" + str(b) + "*e^" + str(a))
        outstr += "y=x^" + str(b) + "*e^" + str(a) + "\n" + str1 + "\n\n"
        plt.plot(px, (np.exp(1)**(b2 * px)) * (np.exp(1)**a2), label="y=e^(" + str(b2) + "x+" + str(a2) + ")")
        outstr += "y=e^(" + str(b2) + "x+" + str(a2) + ")" + "\n" + str2 + "\n\n"
        plt.xlim(xmin=0.5)
        plt.xlim(xmax=25)
        plt.ylim(ymin=1)
        plt.ylim(ymax=100000)
        plt.xscale("log")
        plt.yscale("log")
        plt.legend(loc='best')
        plt.savefig(graph_name + "_" + name + ".png")
        plt.close()
        outf = open(graph_name + "_" + name + ".log", "w")
        outf.write(outstr)
        outf.close()
