#!/usr/bin/env python
# coding: UTF-8
import sys
import re
import json

if __name__ == '__main__':
    log = open(sys.argv[1], 'r')
    lines = log.readlines()
    log.close()
    p = re.compile('^\[[0-9 :]*\] ')
    prefix = re.compile('^.* = ')
    suffix = re.compile(' = .*$')

    data = {}
    data['radius'] = []
    data['diameters'] = []
    data['coverage'] = []
    data['time'] = []
    data['size'] = []
    data['graph_info'] = [{}]
    data_name = ''
    for line in lines:
        if p.search(line):
            line = p.sub('', line)
            key = suffix.sub('', line).strip()
            value = prefix.sub('', line).strip()
            if key in data:
                data[key].append(float(value))
            elif key == 'name':
                data[key] = value
            elif key in {"vertices", "edges"}:
                data['graph_info'][0][key] = int(value)
            elif key == 'graph':
                data['graph_info'][0][key] = value.split(" ")[0]
                data_name = value
    if len(data['radius']) == 0:
        del data['radius']
    if len(data['diameters']) == 0:
        del data['diameters']
    print json.dumps(data, sort_keys=True, indent=4)
