# coding: UTF-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import json
import sys
import re

if __name__ == '__main__':
    log = open(sys.argv[1], 'r')
    jsonData = json.load(log)
    log.close()
    graph_name = jsonData['graph_info'][0]['graph'].replace(" ", "_")
    if '/' in graph_name:
        r = re.compile("/([a-zA-Z0-9_\-\.]*)$")
        m = r.search(graph_name)
        graph_name = m.group(1)
    r = re.compile("\.([0-9\-]*)\.[0-9]*$")
    m = r.search(sys.argv[1])
    file_date = m.group(1)
    print graph_name
    for j in range(0, len(jsonData['algorithms'])):
        boxSizes = jsonData['algorithms'][j]['size']
        name = jsonData['algorithms'][j]['name']
        f = open(graph_name + "_" + file_date + ".json", 'w')
        json.dump(jsonData, f, sort_keys=True, indent=4)
        f.close()
