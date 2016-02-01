#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import statsmodels.api as sm
import json
import sys
import re


if __name__ == "__main__":
    combined = {}
    rads = []
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()

        radius = json_data["radius"][0]
        time = json_data['time'][0]
        size = json_data['size'][0]
        combined[radius] = {'size': size, 'time': time, 'radius': radius}
        rads.append(radius)

        graph_info = json_data['graph_info']
        name = json_data['name']
    rads.sort()
    data = {
        "time": [],
        "size": [],
        "radius": [],
        "name": name,
        "graph_info": graph_info,
    }
    for radius in rads:
        data["time"].append(combined[radius]["time"])
        data["radius"].append(combined[radius]["radius"])
        data["size"].append(combined[radius]["size"])
    json_str = json.dumps(data, sort_keys=True, indent=4)
    # print json_str
    graph_name = graph_info[0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.@]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    f = open(name + "-" + graph_name + ".json", 'w')
    f.write(json_str)
    f.close()
