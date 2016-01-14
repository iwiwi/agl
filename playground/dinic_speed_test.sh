#!/bin/bash

../bin/dinic_speed_test -type=agl -graph=/home/math/data/soc-Slashdot0902.agl -method=oneside --jlog_suppress_log 2> oneside.out
../bin/dinic_speed_test -type=agl -graph=/home/math/data/soc-Slashdot0902.agl -method=twoside --jlog_suppress_log 2> twoside.out
diff oneside.out twoside.out

