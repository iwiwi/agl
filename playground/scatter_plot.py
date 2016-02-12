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


def xy_from_json(json_data):
    boxSizes = json_data['size']
    if 'radius' in json_data:
        radiuses = json_data['radius']
    if 'diameter' in json_data:
        diameters = json_data['diameter']
    if len(boxSizes) <= 1:
        print "very little"
        return
    x = []
    y = []
    for i in range(0, len(boxSizes)):
        if radiuses[i] == 0:
            continue
        if 'diameter' in json_data:
            x.append(diameters[i])
        else:
            x.append(radiuses[i] * 2)
        y.append(boxSizes[i])
        if boxSizes[i] == 1:
            break
    x = np.array(x)
    y = np.array(y)
    name = json_data['name']
    return x, y, name


if __name__ == "__main__":
    for x in xrange(0, 2):

        plt.rcParams['font.size'] = 30
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 8
        plt.rcParams['lines.linewidth'] = 3
        fig = plt.figure(figsize=(8, 6))
        # fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
        axes.set_xlabel("$\it{\ell}$", fontsize=40, labelpad=None)
        axes.set_xscale("log")
        axes.set_xlim(1.0, 10000)
        [i.set_linewidth(2) for i in axes.spines.itervalues()]
        if x == 0:
            axes.set_ylabel("$\it{b(\ell)}$", fontsize=40)
            axes.set_yscale("log")

            qx = np.linspace(0.5, 1000000, 1000000)
            axes.plot(qx, qx**(-np.log(6) / np.log(3))
                      * np.exp(1)**12.318, 'r-', label='theory')
            
            for ai in range(1, len(sys.argv)):
                print sys.argv[ai]
                log = open(sys.argv[ai], 'r')
                json_data = json.load(log)
                log.close()

                px, py, pname = xy_from_json(json_data)
                axes.plot(px, py, 'bo',
                          markersize=16,
                          markerfacecolor='w',
                          markeredgewidth=2,
                          )
            plt.ylim(ymin=1)
            plt.ylim(ymax=200000)
            plt.tight_layout(pad=0.2)
            axes.legend(loc='upper right', frameon=False, fontsize=32)
            plt.savefig("scatter.pdf")
            plt.savefig("scatter.png")
            plt.close()
        else:
            axes.set_ylabel("CV", fontsize=35)
            pxs = {}
            for ai in range(1, len(sys.argv)):
                print sys.argv[ai]
                log = open(sys.argv[ai], 'r')
                json_data = json.load(log)
                log.close()

                px, py, pname = xy_from_json(json_data)
                if ai == 1:
                    for r in px:
                        pxs[r] = np.array([])

                for i in range(0, len(py)):
                    pxs[px[i]] = np.append(pxs[px[i]], py[i])

            cv = np.array([])
            for r in px:
                cv = np.append(cv, pxs[r].std() / pxs[r].mean())
            axes.plot(px, cv, 'bo-', markersize=18)
            axes.set_ylim(-0.01, 0.21)
            plt.tight_layout(pad=0.2)
            plt.savefig("cv.pdf")
            plt.savefig("cv.png")
            plt.close()
