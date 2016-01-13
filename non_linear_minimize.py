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


def linearRegression(x, y):
    X = np.column_stack((np.repeat(1, x.size), np.log(x)))
    model = sm.OLS(np.log(y), X)
    results = model.fit()
    a, b = results.params
    return np.array([a, b])


def theoreticalValue(beta, x):
    a = beta[0]
    b = beta[1]
    return (x**b) * (np.exp(1)**a)


def fitFunc(beta, x, y):
    residual = y - theoreticalValue(beta, x)
    return residual


if __name__ == "__main__":
    # サンプルデータの作成
    px = np.array([2,
                   5,
                   16,
                   47,
                   142,
                   425])
    py = np.array([97147,
                   62409,
                   9787,
                   1547,
                   212,
                   19])

    initialValue = linearRegression(px, py)
    result = scipy.optimize.leastsq(fitFunc, initialValue, args=(px, py))

    qx = np.linspace(0.5, 1000, 10000)
    print(result)
    plt.xlim(xmin=0.5)
    plt.xlim(xmax=500)
    plt.ylim(ymin=0.5)
    plt.ylim(ymax=120000)
    plt.plot(px, py, 'o')
    plt.plot(qx, theoreticalValue(result[0], qx))
    plt.savefig("test.png")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(ymax=1000000)
    plt.savefig("test2.png")
    plt.close()
