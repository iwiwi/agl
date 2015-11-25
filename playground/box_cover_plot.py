import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import sys

if __name__ == '__main__':
    log = open(sys.argv[1], 'r')
    jsonData = json.load(log)
    log.close()
    boxSize = jsonData['algorithms'][0]['size']
    x = []
    y = []
    for i in range(0, len(boxSize)):
        x.append(2 * i + 1)
        y.append(boxSize[i])
        if boxSize[i] == 1:
            break
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(xmin=0.5)
    plt.xlim(xmax=25)
    plt.scatter(x, y)
    plt.savefig(sys.argv[2])
