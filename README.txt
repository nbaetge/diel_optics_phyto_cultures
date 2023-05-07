This repository was created to share the data and analyses used for :

Baetge, N., Halesy, K., Graff, J.R., Ver Wey, B., Westberry, T., Appel, A.E., Bourdin, G., Begouen Demeaux, C., Boss, E., Behrenfeld, M.J. (in review). Bio-optical properties of cultured phytoplankton over complete diel cycles. 

All code and data provided in this repository are freely available for use. Please do not hesitate to contact me if you have questions or comments or if you find errors/inaccuracies.


R scripts used in analyses are found in the "analyses" folder and can be run using R Studio. The analyses can be replicated as follows:

-Clone repository

-Open R project

-Run scripts in this order:

1. Growth_Rates.Rmd
2. FRR.Rmd
3. POC.Rmd
4. POC_filtrates_FCM.Rmd
5. Culture_and_Optics_FCM.Rmd
6. ACS.Rmd
	rd_data.m
7. Spectral_Decomposition.Rmd
	output files = "Gaussian_Decompositions.csv" (decomposed absorption spectra) & "Pigment_Ratios.csv" (pigment absorption ratios)
8. BB3.Rmd
9. Final_Data.Rmd
	output files = "FINAL_BOTTLE.csv" (data for each culture bottle) & "FINAL_SUMMARY.csv" (summary data for each phytoplankton and experiment)
10. Calculations_for_MS.Rmd
	output files = "Correlations.csv" (Spearman's correlation coefficients)
11. Plots.Rmd

