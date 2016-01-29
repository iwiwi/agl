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

if __name__ == "__main__":
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        px, py = xy_from_json(json_data)
        sys.stdout.write(str(json_data['graph_info'][0]['graph']) + "\t")
        sys.stdout.write(str(json_data['k']) +
                         "\t" + str(json_data['upper_param']) + "\t")
        for i in range(len(px)):
            sys.stdout.write(str(py[i]) + "\t")
        sys.stdout.write("\n")
