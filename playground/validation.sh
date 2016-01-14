#!/bin/bash

cd ..
# 愚直解
# bin/validater -type=agl -graph=/home/math/data/facebook_combined.agl -method=test -validation_data_path=bin/facebook-flow.data 
# 自作解
bin/gomory_hu -type=agl -graph=/home/math/data/facebook_combined.agl -method=all_pair -validation_data_path=bin/facebook_gusfield.data
diff bin/facebook-flow.data bin/facebook_gusfield.data

