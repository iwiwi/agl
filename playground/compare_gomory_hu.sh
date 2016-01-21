#!/bin/bash

bin/gomory_hu -type=agl -graph=/home/math/data/com-dblp.ungraph.agl -method=print_gomory_hu_tree -gomory_fu_builder=Gusfield3 -validation_data_path=bin/com-dblp_gusfield3.data
bin/gomory_hu -type=agl -graph=/home/math/data/com-dblp.ungraph.agl -method=print_gomory_hu_tree -gomory_fu_builder=Gusfield4 -validation_data_path=bin/com-dblp_gusfield4.data

bin/gomory_hu_tree_isomorphic -s1=bin/com-dblp_gusfield3.data -s2=bin/com-dblp_gusfield4.data
