# README
# R code, functions, and data for the paper
"Covariate-driven factorization by thresholding for multi-blick data" by Xing Gao, Sungwon Lee, Gen Li, and Sungkyu Jung

##### Main functions to apply COBS #####
#In **routines** folder - all functions required to run COBS #
**bisect.dec.R ** - function to find the root for the optimization problem of v
**COBS function in r.R ** - wrapper of functions for simulation results evaluation
**cobs.R ** - wrapper for COBS function with fixed tuning parameter
**cobs.tune.R ** - wrapper for COBS function with BIC
**get.initial.R ** - function to get the initial value v_0 for the iteration
**Sigma.est.R ** - functions to estimate sigma_e^2 and sigma_f^2 for the error terms in the rank-one model
**Sigmas_Full.R ** - function to estimate sigma_0^2 and the diagonal values in Sigma_F for the full model
**summary.cobs.r ** - function to summary the results after applying cobs.R or cobs.tune.R
**update.b.R ** - function to update b in each iteration
**update.v.R ** - function to update v in each iteration

##### Code for simulation comparisons  #####
**Sims_COBS.R ** - code for simulation analysis using COBS
**Sims_OtherMethods1.R ** - code for simulation analysis using SVD, AJIVE, SupSVD
**Sims_OtherMethods2.R ** - code for simulation analysis using SLIDE, GFA, RRR/SRRR, SPCA
**Sims_SIFA.m ** - matlab code to run simulation analysis using SIFA

##### Code for real data applications - using breast cancer data (BRCA) #####
**BRCA_data_testing.R ** - code to run analysis using BRCA data with testing and training data sets for COBS, SLIDE, GFA, AJIVE, SupSVD, SIFA
**BRCA_data.R ** - code to run analysis with full BRCA data sets for COBS
**Read_BRCA.R ** - code to read the ready-to-use BRCA data

#In ** BRCA** folder #
**BRCA_Data_COBS.Rdata** - the ready-to-use BRCA data with n=770, p=600 (200 genes in each group)


##### Other relevant code to the paper and web appendix #####
**Toy Data.R**  - code to generate and to visualize the toy data example
**Data Sim_Small_Sample.R ** - code to generate simulation data for high dimension, small sample cases
**Data Sim_Large_Sample.R ** - code to generate simulation data for low dimension, large sample cases
**Melanoma_data.R ** - code to run image-feature analysis using COBS, SLIDE, GFA, AJIVE, SupSVD, SIFA
**Melanoma ** folder - image-feature data sets 
**OtherMethods** folder -  code/functions for SLIDE and SIFA, downloaded from the authors webpage


