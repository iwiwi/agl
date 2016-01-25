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
        sys.stderr.write(sys.argv[ai] + "\n")

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
            ub = str(json_data['upper_param'])
            method = "Sketch (k=" + k + " u=" + ub + ")"
            data[graph_name][method] = {
                "k": int(k),
                "upper_param": float(ub),
                "fractal": frac_result,
                "exponential": expo_result,
                "time": time
            }
            methods[method] = {
                "k": int(k),
                "upper_param": float(ub)
            }
        else:
            data[graph_name][method] = {
                "fractal": frac_result,
                "exponential": expo_result,
                "time": time
            }
            methods[method] = {}

    sorted(methods.items(), key=lambda x: x[0])
    for key in data.keys():
        for method in methods.keys():
            if method not in data[key]:
                data[key][method] = {}
    json_str = json.dumps(data, sort_keys=True, indent=4)
    # print json_str

    sys.stdout.write(",,,")
    for method in methods.keys():
        sys.stdout.write(method + ",,,")
    sys.stdout.write("\n")
    sys.stdout.write("model,vertices,edges,")
    for method in methods.keys():
        sys.stdout.write("fractal,exponential,time,")
    sys.stdout.write("\n")
    for graph_name in data.keys():
        if re.match(r"^flower\-", graph_name) or re.match(r"^shm\-", graph_name):
            u = re.sub(r'^.*\-(\d+)\-\d+$', r"\1", graph_name)
            v = re.sub(r'^.*\-(\d+)$',  r"\1", graph_name)
            model = re.sub(r'^([a-z]*)\-.*$',  r"\1", graph_name)
            model = "(" + u + " " + v + ")-" + model
            sys.stdout.write(model + ",")
        sys.stdout.write(str(data[graph_name]["vertices"]) + ",")
        sys.stdout.write(str(data[graph_name]["edges"]) + ",")
        for method in methods.keys():
            if 'fractal' in data[graph_name][method]:
                sys.stdout.write(
                    str(data[graph_name][method]["fractal"]) + ",")
                sys.stdout.write(
                    str(data[graph_name][method]["exponential"]) + ",")
                sys.stdout.write(str(data[graph_name][method]["time"]) + ",")
            else:
                sys.stdout.write(",,,")
        sys.stdout.write("\n")

    # Spreadsheat
    # sys.stdout.write(",graph,")
    # for graph_name in data.keys():
    #     sys.stdout.write(graph_name + ",")
    # sys.stdout.write("\n")
    # sys.stdout.write(",vertices,")
    # for graph_name in data.keys():
    #     sys.stdout.write(str(data[graph_name]['vertices']) + ",")
    # sys.stdout.write("\n")
    # sys.stdout.write(",edges,")
    # for graph_name in data.keys():
    #     sys.stdout.write(str(data[graph_name]['edges']) + ",")
    # sys.stdout.write("\n")
    # for method in methods.keys():
    #     sys.stdout.write(method + "," + "exp,")
    #     for graph_name in data.keys():
    #         if 'exponential' in data[graph_name][method]:
    #             sys.stdout.write(
    #                 str(data[graph_name][method]['exponential']) + ",")
    #         else:
    #             sys.stdout.write(",")
    #     sys.stdout.write("\n")

    #     sys.stdout.write(",frac,")
    #     for graph_name in data.keys():
    #         if 'fractal' in data[graph_name][method]:
    #             sys.stdout.write(
    #                 str(data[graph_name][method]['fractal']) + ",")
    #         else:
    #             sys.stdout.write(",")
    #     sys.stdout.write("\n")

    # sys.stdout.write(",judge,")
    # for graph_name in data.keys():
    #     if 'exponential' in data[graph_name][method] and 'fractal' in data[graph_name][method]:
    #         e = data[graph_name][method]['exponential']
    #         f = data[graph_name][method]['fractal']
    #         pattern = r"^.*\-1\-[0-9]*$"
    #         if re.search(pattern, graph_name):
    #             if f > e:
    #                 sys.stdout.write(str(f - e) + ",")
    #             else:
    #                 sys.stdout.write(str(f - e) + ",")
    #         elif f < e:
    #             sys.stdout.write(str(e - f) + ",")
    #         else:
    #             sys.stdout.write(str(e - f) + ",")
    #     else:
    #         sys.stdout.write(",")
    # sys.stdout.write("\n")

    # sys.stdout.write(",time,")
    # for graph_name in data.keys():
    #     if 'time' in data[graph_name][method]:
    #         sys.stdout.write(str(data[graph_name][method]['time']) + ",")
    #     else:
    #         sys.stdout.write(",")
    # sys.stdout.write("\n")
