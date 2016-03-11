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
    f = open(sys.argv[1], 'r')
    coloring_json = json.load(f)
    f.close()
    f = open(sys.argv[2], 'r')
    template_json = json.load(f)
    f.close()

    assert coloring_json['graph_info'] == template_json['graph_info']

    cleaned_json = {
        "time": coloring_json['time'],
        "size": [],
        "diameter": [],
        "name": coloring_json['name'],
        "graph_info": coloring_json['graph_info'],
    }

    coloring_map = {}
    for i in xrange(0, len(coloring_json["diameter"])):
        coloring_map[coloring_json["diameter"][i]] = coloring_json["size"][i]

    for rad in template_json['radius']:
        d = rad * 2
        if d in coloring_map:
            cleaned_json['diameter'].append(d)
            cleaned_json['size'].append(coloring_map[d])

    f = open(sys.argv[1] + ".cleaned.json", 'w')
    f.write(json.dumps(cleaned_json, sort_keys=True, indent=4))
    f.close()
