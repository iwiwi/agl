import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
import sys

if __name__ == '__main__':
    log = open(sys.argv[1], 'r')
    jsonData = json.load(log)
    log.close()
    boxSize = jsonData['algorithms'][0]['size']
    x = []
    for i in range(0, len(boxSize)):
        x.append(2 * i + 1)
    df = pd.DataFrame({'2*radius+1': x, 'Box Size': boxSize})
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(yscale="log", xscale="log")
    # ax.set(xscale="log", yscale="log")
    sns_plot = sns.pointplot(x="2*radius+1", y="Box Size",
                             data=df, ax=ax, linestyles=['  '])
    fig = sns_plot.get_figure()
    fig.savefig(sys.argv[2])
