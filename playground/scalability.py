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
    data = {}
    methods = {}
    for ai in range(1, len(sys.argv)):
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()
        sys.stderr.write(sys.argv[ai] + "\n")

        # Graph Information
        graph_name = json_data['graph_info'][0]['graph']
        edges = json_data['graph_info'][0]['edges']
        vertices = json_data['graph_info'][0]['vertices']

        # General Results
        time = 0
        for t in json_data['time']:
            time += t
        method = json_data['name']

        if graph_name not in data:
            data[graph_name] = {
                "edges": edges,
                "vertices": vertices
            }

        data[graph_name][method] = {
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

    # CSV
    sys.stdout.write("\t\t\t")
    for method in methods.keys():
        sys.stdout.write(method + "\t")
    sys.stdout.write("\n")
    sys.stdout.write("model\tvertices\tedges\t")
    for method in methods.keys():
        sys.stdout.write("time\t")
    sys.stdout.write("\n")
    for graph_name in data.keys():
        if re.match(r"^flower\-", graph_name) or re.match(r"^shm\-", graph_name):
            u = re.sub(r'^.*\-(\d+)\-\d+$', r"\1", graph_name)
            v = re.sub(r'^.*\-(\d+)$',  r"\1", graph_name)
            model = re.sub(r'^([a-z]*)\-.*$',  r"\1", graph_name)
            model = "(" + u + ", " + v + ")-" + model
            sys.stdout.write(model + "\t")
        else:
            sys.stdout.write(graph_name + "\t")
        sys.stdout.write(str(data[graph_name]["vertices"]) + "\t")
        sys.stdout.write(str(data[graph_name]["edges"]) + "\t")
        for method in methods.keys():
            if 'time' in data[graph_name][method]:
                sys.stdout.write(str(data[graph_name][method]["time"]) + "\t")
            else:
                sys.stdout.write("\t")
        sys.stdout.write("\n")
