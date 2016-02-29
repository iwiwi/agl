#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import statsmodels.api as sm
import matplotlib.pyplot as plt
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


if __name__ == "__main__":
    analytical = {
        "flower-223950-3-3": {4: 37326,
                              10: 6222,
                              28: 1038,
                              82: 174,
                              244: 30,
                              730: 6,
                              2188: 2},
        "shm-234376-4-2": {2: 46876,
                           8: 9376,
                           26: 1876,
                           80: 376,
                           242: 76,
                           728: 16,
                           2186: 4,
                           6560: 1},
        "shm-312501-5-2": {2: 62501,
                           8: 12501,
                           26: 2501,
                           80: 501,
                           242: 101,
                           728: 21,
                           2186: 5,
                           6560: 1},
        "shm-67229-5-3": {2: 9605,
                          8: 1373,
                          26: 197,
                          80: 29,
                          242: 5,
                          728: 1},
    }
    data = {}
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        px, py = xy_from_json(json_data)

        graph = json_data['graph_info'][0]['graph']
        if graph not in data:
            data[graph] = {}

        k = json_data['k']
        if k not in data[graph]:
            data[graph][k] = {}

        a = json_data['upper_param']
        if a not in data[graph][k]:
            data[graph][k][a] = {}

        for i in range(len(px)):
            if px[i] not in data[graph][k][a]:
                data[graph][k][a][px[i]] = []
            data[graph][k][a][px[i]].append(py[i])

    json_str = json.dumps(data, sort_keys=True, indent=4)
    # print json_str

    # graph_keys = sorted(data)
    # for graph in graph_keys:
    #     sys.stdout.write(graph + "\n")
    #     ks = sorted(data[graph], key=lambda x: int(x))
    #     for k in ks:
    #         alphas = sorted(data[graph][k])
    #         for a in alphas:
    #             sys.stdout.write(str(k) + "\t")
    #             sys.stdout.write(str(a) + "\t")
    #             rs = sorted(data[graph][k][a])
    #             s = 0
    #             all4std = np.array([])
    #             for r in rs:
    #                 data[graph][k][a][r] = np.array(data[graph][k][a][r])
    #                 # sys.stdout.write(str(data[graph][k][a][r].mean()) + "\t")
    #                 data[graph][k][a][r] = (
    #                     data[graph][k][a][r] / (analytical[graph][r] + 0.0))
    #                 s += data[graph][k][a][r].mean()
    #                 all4std = np.append(all4std, data[graph][k][a][r])
    #             sys.stdout.write(str(s / len(rs)) + "\t")
    #             sys.stdout.write(str(all4std.std()) + "\n")
    #     sys.stdout.write("\n")

    graph_keys = sorted(data)
    for graph in graph_keys:
        sys.stdout.write(graph + "\n")
        ks = sorted(data[graph], key=lambda x: int(x))
        for k in ks:
            sys.stdout.write(str(k) + ",")
        sys.stdout.write("\n")

        alphas = sorted(data[graph][k])
        for a in alphas:
            sys.stdout.write(str(a) + ",")
        sys.stdout.write("\n")

        for a in alphas:
            for k in ks:
                rs = sorted(data[graph][k][a])
                s = 0
                all4std = np.array([])
                for r in rs:
                    data[graph][k][a][r] = np.array(data[graph][k][a][r])
                    # sys.stdout.write(str(data[graph][k][a][r].mean()) + ",")
                    data[graph][k][a][r] = (data[graph][k][a][r] / (analytical[graph][r] + 0.0))
                    s += data[graph][k][a][r].mean()
                    all4std = np.append(all4std, data[graph][k][a][r])
                sys.stdout.write(str(s / len(rs)) + ",")
        sys.stdout.write("\n")

        for a in alphas:
            for k in ks:
                rs = sorted(data[graph][k][a])
                all4std = np.array([])
                for r in rs:
                    all4std = np.append(all4std, data[graph][k][a][r])
                sys.stdout.write(str(all4std.std()) + ",")
        sys.stdout.write("\n")
    sys.stdout.write("\n")

    if len(graph_keys) == 1:
        cv = np.array([])
        graph_keys = sorted(data)
        for graph in graph_keys:
            ks = sorted(data[graph], key=lambda x: int(x))
            alphas = sorted(data[graph][k])
            if len(ks) > 1:
                continue
            if len(alphas) > 1:
                continue

            for a in alphas:
                for k in ks:
                    rs = sorted(data[graph][k][a])
                    for r in rs:
                        cv = np.append(cv, data[graph][k][a][r].std() / data[graph][k][a][r].mean(),)
        rads = sorted(analytical[graph])
        print rads, cv
        plt.plot(rads, cv, "o-")
        plt.xscale("log")
        plt.xlabel("$\it{l_{B}}$", fontsize=18)
        plt.ylabel("$CV$", fontsize=18)
        plt.savefig(graph + ".png")
        plt.savefig(graph + ".pdf")
