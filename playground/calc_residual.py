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
        # sys.stderr.write(sys.argv[ai] + "\n")

        px, py = xy_from_json(json_data)
        if (len(px) < 3):
            continue

        # Real_data
        # ok = True
        # for i in xrange(1, len(px)):
        #     if px[i] > px[i - 1] + 2:
        #         ok = False
        #     if i is len(px) - 1 and len(px) < 9 and py[i] > 1:
        #         ok = False
        # if not ok:
        #     continue

        frac_init = linearRegression(np.log(px), np.log(py))
        fractal = scipy.optimize.leastsq(fractalFit, frac_init, args=(px, py))
        frac_result = sum(fractalFit(fractal[0], px, py)**2)

        expo_init = linearRegression(px, np.log(py))
        expo = scipy.optimize.leastsq(expoFit, expo_init, args=(px, py))
        expo_result = sum(expoFit(expo[0], px, py)**2)

        # print fractalFit(fractal[0], px, py)
        # print expoFit(expo[0], px, py)
        # print fractalFit(fractal[0], px, py)**2
        # print expoFit(expo[0], px, py)**2
        # print frac_result, expo_result

        # Graph Information
        graph_name = json_data['graph_info'][0]['graph']
        edges = json_data['graph_info'][0]['edges']
        vertices = json_data['graph_info'][0]['vertices']

        # General Results
        time = 0
        for x in xrange(len(json_data['time'])):
            if x > 0 and json_data['size'][x - 1] == 1:
                break
            time += json_data['time'][x]
        method = json_data['name']

        memory = -1
        if "run" in json_data and "memory" in json_data["run"]:
            memory = json_data["run"]["memory"]

        if graph_name not in data:
            data[graph_name] = {"edges": edges, "vertices": vertices}

        # Sketch
        if 'k' in json_data:
            k = str(json_data['k']).zfill(4)

            if 'upper_param' in json_data:
                ub = str(json_data['upper_param'])
            elif 'alpha' in json_data:
                ub = str(json_data['alpha'])

            method = "Sketch (k=" + k + ", u=" + ub + ")"
            data[graph_name][method] = {
                "k": int(k),
                "upper_param": float(ub),
                "fractal": frac_result,
                "exponential": expo_result,
                "time": time,
                "memory": memory
            }
            methods[method] = {"k": int(k), "upper_param": float(ub)}
        else:
            if method in data[graph_name] and data[graph_name][method]["time"] < time:
                continue
            data[graph_name][method] = {
                "fractal": frac_result,
                "exponential": expo_result,
                "time": time,
                "memory": memory
            }
            methods[method] = {}

    sorted(methods.items(), key=lambda x: x[0])
    for key in data.keys():
        for method in methods.keys():
            if method not in data[key]:
                data[key][method] = {}
    json_str = json.dumps(data, sort_keys=True, indent=4)
    # print json_str

    # CSV
    sys.stdout.write("\t\t\t")
    for method in methods.keys():
        sys.stdout.write(method + "\t\t\t")
    sys.stdout.write("\n")
    sys.stdout.write("model\tvertices\tedges\t")
    for method in methods.keys():
        sys.stdout.write("-log10(P/E)\ttime\tmemory\t")
    sys.stdout.write("\n")
    for graph_name in data.keys():
        if re.match(r"^flower\-", graph_name):
            u = re.sub(r'^.*\-(\d+)\-\d+$', r"\1", graph_name)
            v = re.sub(r'^.*\-(\d+)$', r"\1", graph_name)
            model = re.sub(r'^([a-z]*)\-.*$', r"\1", graph_name)
            model = "(" + u + ", " + v + ")-" + model
            sys.stdout.write(model + "\t")
        else:
            sys.stdout.write(graph_name + "\t")
        sys.stdout.write(str(data[graph_name]["vertices"]) + "\t")
        sys.stdout.write(str(data[graph_name]["edges"]) + "\t")
        for method in methods.keys():
            if 'fractal' in data[graph_name][method]:
                sys.stdout.write(str(-np.log10(data[graph_name][method]["fractal"] / data[graph_name][method][
                    "exponential"])) + "\t")
                sys.stdout.write(str(data[graph_name][method]["time"]) + "\t")
                sys.stdout.write(
                    str(data[graph_name][method]["memory"]) + "\t")
            else:
                sys.stdout.write("\t")
                sys.stdout.write("\t")
                sys.stdout.write("\t")
        sys.stdout.write("\n")
