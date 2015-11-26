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
    boxSize = jsonData['algorithms'][0]['size']
    x = []
    y = []
    for i in range(0, len(boxSize)):
        x.append(2 * i + 1)
        y.append(boxSize[i])
        if boxSize[i] == 1:
            break

    x = np.array(x)
    y = np.array(y)
    nsample = x.size
    lnx = []
    lny = []
    for i in range(0, x.size):
        lnx.append(np.log10(x[i]))
        lny.append(np.log10(y[i]))

    lnx = np.array(lnx)
    lny = np.array(lny)

    # おまじない (後で解説)
    X = np.column_stack((np.repeat(1, nsample), lnx))
    # 回帰実行
    model = sm.OLS(lny, X)
    results = model.fit()

    # 結果の概要を表示
    print results.summary()

    # パラメータの推定値を取得
    a, b = results.params

    # プロットを表示
    plt.plot(x, y, 'o')
    plt.plot(x, (x**b) * (10**a))

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(xmin=0.5)
    plt.xlim(xmax=25)
    plt.savefig(sys.argv[2])
