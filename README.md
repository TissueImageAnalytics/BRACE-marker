# BRACE marker
H&E image features based breast cancer marker for ER+/HER2-

# Survival analysis

Description of the files:   
/survival_analysis/KM_curves.py: is used to generate survival results (p-value, C-Index, HR) and KM curves.   
It uses the csv file containing the features and the time/events located under /features/Cohort-B.csv   
KM curves (pdf files) and the summary of the survival results (.csv file) are saved under /results/   
This script also generate csv file (/features/Cohort-B_for_forest_plots.csv) which can then be used by the R script under /Forest_plots/forest_plot_R_script.R) to generate forest plots.   

/models/ contain the trained models.


