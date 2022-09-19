# GCMI Tutorial Examples

This folder contains some tutorial examples for the use of the GCMI package, reproducing the analyses from:

RAA Ince, BL Giordano, C Kayser, GA Rousselet, J Gross and PG Schyns  
"A statistical framework for neuroimaging data analysis based on mutual information estimated via a Gaussian copula"  
bioRxiv [doi:10.1101/043745](http://dx.doi.org/10.1101/043745)

Each of these scripts will download the required data file if it is not be available, so please take care if you are on a slow or metered internet connection. The data files can be downloaded manually from [https://www.robince.net/data/gcmi](https://www.robince.net/data/gcmi).

## [`discrete_eeg.m`](discrete_eeg.m)

This script covers the analyses in Section 4.1, with a full cap event-related EEG data set with 2 discrete stimulus classes. Topics covered include:

- Calculating GCMI between continuous and discrete variables
- Permutation testing with the method of maximum statistics

## [`continuous_meg.m`](continuous_meg.m)

This script covers the analyses in Section 4.2, with a full helmet continuous design MEG data set with a continuous auditory stimulus feature (speech envelope). Topics covered include:

- Calculating GCMI between two continuous variables.
- Block-permutation testing within a continuous experimental design with the method of maximum statistics.
- Calculating GCMI with multivariate responses (planar gradient magnetic field vector).
- Calculating GCMI in amplitude (eg power) and direction (eg phase) of vector quantities.

## [`eeg_temporal_interaction.m`](eeg_temporal_interaction.m)

This script covers the analyses in Section 4.3, with single sensor event-related EEG data in an event-related design with a continuous stimulus feature (sampled eye visibility). Topics covered include:

- Calculating Conditional Mutual Information (CMI) with correlated stimulus features.
- Calculating cross-temporal interaction information.
- Calculating MI in the multivariate repsonse consisting of raw voltage together with the single-trial instaneous temporal derivative.
- Calculating the emergence of novel MI over time.


Questions / comments : robince@gmail.com
