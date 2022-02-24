# Stochastic-Notch-Patterns
Here you'll find the code for generating and analyzing Notch-Delta mediated pattern formation. Our results can be found in "Stochastic fluctuations promote ordered pattern formation of cells in the Notch-Delta signaling pathway"

Please cite the preprint at https://arxiv.org/abs/2202.00763 if you use this code. 


Files to generate the data include:

+seedsForSims.txt
+generate_data.py


Files to analyze the data:

+ analyze_data.py
+ aux_functions.py


Files to generate figures in the manuscript:

+ aux_functions.py
+ plot_script_functions.py
+ produce_data_for_figures.py
+ FiguresForMain.ipynb


Files to generate figures in the SI:

+ pseudocount.py
+ get_data_for_si.py
+ get_psuedo_change.py
+ figures_for_SI.ipynb

Files to generate movies for the SI:

+SI_pattSim_rand_nucl.py


---
Will need to create the directories 

data/pattern
data/lattices
data/random
data/dev
data/mistakes


analysis/pattern
analysis/lattices
analysis/random
analysis/dev
analysis/mistakes


movie_si_patt
