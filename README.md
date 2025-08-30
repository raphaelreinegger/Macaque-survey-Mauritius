This repository contains all code and datasheets necessary to (largely) replicate the analyses from the article "Mauritius' macaques under the thermal camera lens: a 40-year overdue update on a controversial non-native species".

Only the macaque capture figures have been redacted (in "Grid_submission.csv" and "siteCovs_submission.csv"), as these data are not owned by the authors. Analyses can be run without this parameter. Note that the main script still calls the capture variable; it should be removed from the relevant sections before running the code.

The main script is Code_macaque_survey.R, which uses the data from:

- "Grid_submission.csv"
- "obsCovs_submission.csv"
- "siteCovs_submission.csv"

Within the main script, the following functions are sourced to conduct power and sample size sensitivity analyses:

- ModifiedHT_power_analysis.R
- ModifiedHT_SS_analysis.R
- NMix_power_analysis.R
- NMix_SS_analysis.R
