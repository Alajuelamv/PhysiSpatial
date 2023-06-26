#!/bin/sh


model=Fumia2013
#Example  simulations normalized RNA as transitions rates with an amplification factor of 10
sim_case=Visium_META_RNA_norm

python3 Scripts/Simulations/MaBoSS_specific.py $model -sy Linux -p 1 "/results_"$sim_case".txt" -s $sim_case -rb "Results/Profiles/Visium_META_RNA_norm.csv" -rf 10
