# coding: UTF-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import json
import sys

if __name__ == '__main__':
    log = open(sys.argv[1], 'r')
    jsonData = json.load(log)
    log.close()
    graph_name = jsonData['graph_info'][0]['graph'].replace('/', '-')
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
        print results.summary()
        # get estimated params
        a, b = results.params

        X2 = np.column_stack((np.repeat(1, nsample), x))
        model2 = sm.OLS(lny, X2)
        results2 = model2.fit()
        # print summary
        print results2.summary()
        # get estimated params
        a2, b2 = results2.params

        # show plot
        plt.plot(x, y, 'o', label=name)
        px = np.linspace(0.5, 100, 10000)
        plt.plot(px, (px**b) * (np.exp(1)**a), label="y=x^" + str(b) + "*e^" + str(a))
        plt.plot(px, (np.exp(1)**(b2 * px)) * (np.exp(1)**a2), label="y=e^(" + str(b2) + "x+" + str(a2) + ")")
        plt.xlim(xmin=0.5)
        plt.xlim(xmax=25)
        plt.ylim(ymin=1)
        plt.ylim(ymax=100000)
        plt.xscale("log")
        plt.yscale("log")
        plt.legend(loc='best')
        plt.savefig(graph_name + "_" + str(j) + ".png")
        print j
        plt.close()
