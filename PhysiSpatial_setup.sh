#!/bin/sh

python3 clustering_ESTIMATE.py

Rscript Ensemble2Entrez.R

Rscript ESTIMATE_nkmeans.R

python3 input_preprocess_PROFILE.py

Rscript Ensemble2Entrez_2.R

python3 set_up_physi.py

Rscript Fumia_Visium_Profile_onlyRNA.R

cp fromVisium2Simulation.sh ./PROFILE-master

cd PROFILE-master

./fromVisium2Simulation.sh
