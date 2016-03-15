#!/usr/bin/env python
# coding=utf-8

import networkx as nx
import sys
import pandas as pd
import numpy as np
import json
import pylab


def loadTSV(tsv_file):
    es = pd.read_csv(tsv_file, delimiter='\t', header=None)
    print "dataframe loaded"
    G = nx.from_pandas_dataframe(es, 0, 1)
    print "graph generated"
    return G


def getMergedList(G, centers_list, rad):
    N = nx.number_of_nodes(G)
    merged = [-1] * N
    queue = []
    cnt = 0
    for start in centers_list:
        queue.append(start)
        merged[start] = start
        cnt += 1
    print "centers loaded"
    while queue:
        q = queue.pop(0)
        neighbor_list = G.neighbors(q)
        for u in neighbor_list:
            if merged[u] >= 0:
                continue
            merged[u] = merged[q]
            queue.append(u)
            cnt += 1
            if cnt % 5000 == 0:
                print cnt
    return merged


def shrinkGraph(G, centers_list, rad, need_pos=False):
    merged = getMergedList(G, centers_list, rad)
    covered_size = [0] * nx.number_of_nodes(G)
    for m in merged:
        covered_size[m] += 1
    shrinked_edges = []
    for row in G.edges():
        s = row[0]
        t = row[1]
        if s >= len(merged) or t >= len(merged):
            print s, t, len(merged), G.number_of_nodes()
            continue
        if merged[s] == merged[t]:
            continue
        shrinked_edges.append((merged[s], merged[t]))
    shrinked_G = nx.from_pandas_dataframe(pd.DataFrame(shrinked_edges), 0, 1)
    pos = {}
    if need_pos:
        pos = nx.spring_layout(shrinked_G)
        for i in xrange(0, len(merged)):
            if i not in pos:
                pos[i] = pos[merged[i]]
    return (shrinked_G, pos, covered_size)


if __name__ == "__main__":
    G = loadTSV(sys.argv[1])
    f = open(sys.argv[2], 'r')
    centers = json.load(f)["centers"]
    f.close()
    print "graph loaded"

    rad_list = []
    centers_dict = {}
    for pair in centers:
        rad_str = pair.keys()[0]
        rad = int(rad_str)
        rad_list.append(rad)
        centers_dict[rad] = pair[rad_str]
    rad_list.sort()
    print "jlog loaded"

    pos = {}
    for rad in rad_list:
        print rad
        if len(centers_dict[rad]) > 500:
            continue
        if len(pos) == 0:
            shrinked_G, pos, covered_size = shrinkGraph(
                G, centers_dict[rad], rad, True)
        else:
            shrinked_G, _, covered_size = shrinkGraph(
                G, centers_dict[rad], rad)
        nx.draw_networkx_nodes(shrinked_G,
                               pos,
                               nodelist=shrinked_G.nodes(),
                               node_size=[covered_size[v] / 1000 for v in shrinked_G.nodes()])
        nx.draw_networkx_edges(shrinked_G, pos, width=0.1)

        pylab.xlim(-1.1, 1.1)
        pylab.ylim(-1.1, 1.1)
        pylab.savefig(str(rad) + ".png")
        pylab.close()
