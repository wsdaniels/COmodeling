# COmodeling

Repo to hold data and code for manuscript "Predicting Fire Season Intensity in Maritime Southeast Asia with Interpretable Models". 
Manuscript can be found at: https://doi.org/10.31223/X59320

**Response variable data**
Contains week averaged CO observations from MOPITT:
MSEA_V8JMOPITT_weeklyanomalies_WEDCEN_nofill.csv


**Predictor variable data**
Contains week averaged climate mode indices (and OLR as a proxy for MJO):
nino34_weekly_avg.csv
aao_weekly_avg.csv
dmi_weekly_avg.csv
tsa_weekly_avg.csv
msea_olr.csv


**R scripts**
Two R scripts are provided to reproduce analysis from manuscript:
"RUN_FIRST-create_data_matrix.R" takes the response and predictor variable data and prepares them for modeling
"RUN_SECOND-find_best_model.R" takes the prepared data and finds a sequence of optimally performing models for the given response region
