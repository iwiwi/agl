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
from matplotlib.ticker import LogLocator

Vertices = 0
Time = 1
Memory = 2


def get_v(v, t):
    return v[0:len(t)]


def pull_data(data, kind=Vertices, flower=True):
    if flower:
        target_graph = "flower"
    else:
        target_graph = "ba"

    v = []
    vmap = {}
    for graph in data.keys():
        if target_graph not in graph:
            continue
        v.append(data[graph]["vertices"])
        vmap[data[graph]["vertices"]] = graph
    v.sort()

    if kind is Vertices:
        return v

    graph_sorted = [vmap[p] for p in v]
    coloring = []
    cbb = []
    sketch = []
    burning = []
    memb = []
    for graph in graph_sorted:
        if target_graph not in graph:
            continue

        label = "time"
        if kind is not Time:
            label = "memory"

        for method in data[graph].keys():

            if (method is "vertices") or (method is "edges"):
                continue
            if label not in data[graph][method]:
                continue
            if "Burning" in method:
                burning.append(data[graph][method][label])
            elif "MEMB" in method:
                memb.append(data[graph][method][label])
            elif "Coloring" in method:
                coloring.append(data[graph][method][label])
            elif "CBB" in method:
                cbb.append(data[graph][method][label])
            elif "sketch" in method:
                sketch.append(data[graph][method][label])
    return coloring, cbb, sketch, burning, memb


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

        if method not in data[graph_name]:
            data[graph_name][method] = {"time": time}
        data[graph_name][method]["time"] = min(
            data[graph_name][method]["time"], time)
        methods[method] = {}
        if "memory" not in data[graph_name][method]:
            data[graph_name][method]["memory"] = json_data["run"]["memory"]

        data[graph_name][method]["memory"] = min(
            data[graph_name][method]["memory"],
            json_data["run"]["memory"],
        )

    sorted(methods.items(), key=lambda x: x[0])
    for key in data.keys():
        for method in methods.keys():
            if method not in data[key]:
                data[key][method] = {}
    json_str = json.dumps(data, sort_keys=True, indent=4)

    # # CSV
    # sys.stdout.write("\t\t\t")
    # for method in methods.keys():
    #     sys.stdout.write(method + "\t")
    # sys.stdout.write("\n")
    # sys.stdout.write("model\tvertices\tedges\t")
    # for method in methods.keys():
    #     sys.stdout.write("time\t")
    # sys.stdout.write("\n")
    # for graph_name in data.keys():
    #     if re.match(r"^flower\-", graph_name) or re.match(r"^shm\-", graph_name):
    #         u = re.sub(r'^.*\-(\d+)\-\d+$', r"\1", graph_name)
    #         v = re.sub(r'^.*\-(\d+)$',  r"\1", graph_name)
    #         model = re.sub(r'^([a-z]*)\-.*$',  r"\1", graph_name)
    #         model = "(" + u + ", " + v + ")-" + model
    #         sys.stdout.write(model + "\t")
    #     else:
    #         sys.stdout.write(graph_name + "\t")
    #     sys.stdout.write(str(data[graph_name]["vertices"]) + "\t")
    #     sys.stdout.write(str(data[graph_name]["edges"]) + "\t")
    #     for method in methods.keys():
    #         if 'time' in data[graph_name][method]:
    #             sys.stdout.write(str(data[graph_name][method]["time"]) + "\t")
    #         else:
    #             sys.stdout.write("\t")
    #     sys.stdout.write("\n")

    # plot figures
    plt.rcParams['font.size'] = 30
    plt.rcParams['lines.markersize'] = 15
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 0

    for x in xrange(0, 4):
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)

        if x == 0:
            v = pull_data(data=data, kind=Vertices, flower=True)
            coloring, cbb, sketch, burning, memb = pull_data(
                data=data, kind=Time, flower=True)
            axes.set_ylabel("time (sec)")
            filename = "flower_time"
        elif x == 1:
            v = pull_data(data=data, kind=Vertices, flower=True)
            coloring, cbb, sketch, burning, memb = pull_data(
                data=data, kind=Memory, flower=True)
            axes.set_ylabel("memory (KB)")
            filename = "flower_memory"
        elif x == 2:
            v = pull_data(data=data, kind=Vertices, flower=False)
            coloring, cbb, sketch, burning, memb = pull_data(
                data=data, kind=Time, flower=False)
            axes.set_ylabel("time (sec)")
            filename = "ba_time"
        else:
            v = pull_data(data=data, kind=Vertices, flower=False)
            coloring, cbb, sketch, burning, memb = pull_data(
                data=data, kind=Memory, flower=False)
            axes.set_ylabel("memory (KB)")
            filename = "ba_memory"

        axes.plot(get_v(v, coloring), coloring, "bs-",
                  label="GC", markerfacecolor='w', markeredgecolor='b', markeredgewidth=2)
        axes.plot(get_v(v, cbb), cbb, "mv-", label="CBB")
        axes.plot(get_v(v, memb), memb, "k^--", label="MEMB")
        axes.plot(get_v(v, sketch), sketch, "ro-", label="Sketch")
        axes.plot(get_v(v, burning), burning, "gx-", label="MVB")

        axes.set_xscale("log")
        axes.set_yscale("log")

        handles, labels = axes.get_legend_handles_labels()
        # reverse the order
        labels = labels[::-1]
        handles = handles[::-1]
        mvb_label = labels[0]
        mvb_handle = handles[0]
        labels = labels[1:]
        handles = handles[1:]
        labels.append(mvb_label)
        handles.append(mvb_handle)

        axes.legend(handles, labels,
                    loc='best', fontsize=24, frameon=False)
        # axes.legend(loc='best', fontsize=24, frameon=False)

        axes.set_xlabel("number of vertices", fontsize=30, labelpad=0.1)
        if x == 0 or x == 2:
            axes.yaxis.set_major_locator(LogLocator(
                base=100.0, subs=[1.0], numdecs=4, numticks=15))
        [i.set_linewidth(2) for i in axes.spines.itervalues()]

        fig.tight_layout(pad=0.1)
        fig.savefig("scalability_" + filename + ".pdf")
        fig.savefig("scalability_" + filename + ".png")
