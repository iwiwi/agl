import numpy as np
import pandas as pd
import seaborn as sns
import json
import sys

if __name__ == '__main__':
    f = open(sys.argv[1], 'r')
    jsonData = json.load(f)
    f.close()
    boxSize = jsonData['algorithms'][0]['size']
    df = pd.DataFrame({'Rad': range(len(boxSize)), 'Size': boxSize})
    sns_plot = sns.jointplot(x="Rad", y="Size", data=df)
    # sns.plot(boxSize)
    sns_plot.savefig(sys.argv[2])
