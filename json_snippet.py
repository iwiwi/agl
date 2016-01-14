#!/usr/bin/env python
# coding: UTF-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import json
import sys
import re


def xy_from_json(json_data):
    graph_name = json_data['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    vertices = json_data['graph_info'][0]['vertices']
    print graph_name
    boxSizes = json_data['size']
    if 'radius' in json_data:
        radiuses = json_data['radius']
    if 'diameter' in json_data:
        diameters = json_data['diameter']
    name = json_data['name']

    if len(boxSizes) <= 1:
        print "very little"
        return
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

if __name__ == '__main__':
    for ai in range(1, len(sys.argv)):
        print sys.argv[ai]
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()
        px, py = xy_from_json(json_data)
        print px
        print py
