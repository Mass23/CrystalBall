# CrystalBall

## 1. Introduction
Code created for the "Crystal ball" study that models glacier-fed stream environmental parameters and bacterial strain abundances, and projects them onto future scenarios of climate change.

## 2. Installation
- 2.1 git clone the repository
- 2.2 install the conda envs in the /env directory

## 3. Run the study
- 3.1 The first script (1_CreateData.py) compiles the data: python3 1_CreateData.py
- 3.2 The second script (2_Analyse.R) creates models, runs all analyses, and generates the plots and stats: Rscript 2_Analyse.R
- 3.3 All results can be accessed in the stats/ and plots/ folders
