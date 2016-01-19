#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import statsmodels.api as sm
import json
import sys
import re


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


def linearRegression(x, y):
    X = np.column_stack((np.repeat(1, x.size), x))
    model = sm.OLS(y, X)
    results = model.fit()
    a, b = results.params
    return np.array([a, b])


def fractalValue(beta, x):
    a = beta[0]
    b = beta[1]
    return (x**b) * (np.exp(1)**a)


def fractalFit(beta, x, y):
    residual = y - fractalValue(beta, x)
    return residual


def expoValue(beta, x):
    a = beta[0]
    b = beta[1]
    return (np.exp(1)**(b * x)) * (np.exp(1)**a)


def expoFit(beta, x, y):
    residual = y - expoValue(beta, x)
    return residual


if __name__ == "__main__":
    data = {}
    methods = {}
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        px, py = xy_from_json(json_data)

        frac_init = linearRegression(np.log(px), np.log(py))
        fractal = scipy.optimize.leastsq(fractalFit, frac_init, args=(px, py))
        frac_result = sum(fractalFit(fractal[0], px, py)**2)

        expo_init = linearRegression(px, np.log(py))
        expo = scipy.optimize.leastsq(expoFit, expo_init, args=(px, py))
        expo_result = sum(expoFit(expo[0], px, py)**2)

        # Graph Information
        graph_name = json_data['graph_info'][0]['graph']
        edges = json_data['graph_info'][0]['edges']
        vertices = json_data['graph_info'][0]['vertices']

        # General Results
        time = json_data['run']['time']
        method = json_data['name']

        if method != 'MEMB':
            method = 'Sketch'

        if graph_name not in data:
            data[graph_name] = {
                "edges": edges,
                "vertices": vertices
            }

        # Sketch
        if 'k' in json_data:
            k = str(json_data['k']).zfill(4)
            ub = str(json_data['size_upper']).zfill(8)
            method = "Sketch-k." + k + "-upper_bound." + ub
            data[graph_name][method] = {
                "k": int(k),
                "upper_bound": int(ub),
                "fractal": frac_result,
                "exponential": expo_result,
                "time": time
            }
            methods[method] = {
                "k": int(k),
                "upper_bound": int(ub)
            }
        else:
            data[graph_name][method] = {
                "fractal": frac_result,
                "exponential": expo_result,
                "time": time
            }
            methods[method] = {}

    for key in data.keys():
        for method in methods.keys():
            if method not in data[key]:
                data[key][method] = {}
    json_str = json.dumps(data, sort_keys=True, indent=4)
    print json_str
    # print methods
