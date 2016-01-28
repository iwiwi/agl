#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import statsmodels.api as sm
# import seaborn as sns
import json
import sys
import re


if __name__ == "__main__":
    k = np.array([32, 64, 128, 256, 512])
    residual = np.array([739, 427, 416, 442, 1])

    fig = plt.figure(figsize=(6, 2))
    axes = fig.add_subplot(1, 1, 1)
    axes.plot(k, residual, "ko-")
    axes.set_xscale("log", basex=2)
    axes.set_xlabel("k")
    axes.set_ylabel("squared residuals")
    axes.set_yticks([0, 200, 400, 600, 800])

    fig.tight_layout()
    fig.savefig('squared_residuals_with_k.pdf')
