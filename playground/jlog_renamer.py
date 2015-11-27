import sys
import os
import re
import json

if __name__ == '__main__':
    files = os.listdir(sys.argv[1])
    for file in files:
        log = open(file, 'r')
        jsonData = json.load(log)
        log.close()
        boxSize = jsonData['algorithms'][0]['size']
