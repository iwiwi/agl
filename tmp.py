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


def plot_data(jsonData, filename):
    graph_name = jsonData['graph_info'][0]['graph']
    print graph_name
    prefix = re.compile('^.*\-' + graph_name)
    next_name = graph_name + prefix.sub('', filename).strip()
    m = re.match(r"(?P<f>[a-z]*\-[0-9]*\-[0-9]*\-[0-9]*)", next_name)
    jsonData['graph_info'][0]['graph'] = m.group('f')
    json_str = json.dumps(jsonData, sort_keys=True, indent=4)
    f = open(filename, "w")
    f.write(json_str)
    f.close()

if __name__ == '__main__':
    for ai in range(1, len(sys.argv)):
        print sys.argv[ai]
        log = open(sys.argv[ai], 'r')
        json_data = json.load(log)
        log.close()
        plot_data(json_data, sys.argv[ai])
