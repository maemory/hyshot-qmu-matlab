# Hyshot QMU

Personal repo of some old research that I don't want to lose.

---------

## FPVA samples

Paul's initial set of 12*(1+12)=312 FPVA sample values are in the file `PAUL-chemistry-perturbation.txt`. I've collected just the nominal ones in `FPVA-samples.txt`.

## Inflow samples

These were broken into two files. Initially we ran only 14 samples (2*7 dimensions) but reviewers wanted more for confirmation of AS credibility. We then created an additional 36 samples - run only for FFR30. The data for these are in the files `first14_samples.mat` and `second36_samples.mat`, respectively.

## Active Subspace post processing

the file `AS_post_processing.m` handles the inflow uncertainty AS, while the file `FPVA_post_processing.m` takes care of the chemistry.

## tabulated data

Data from the flowfields has been reduced to the excel sheet `activeSubspace_results.xlsx`. Contains both inflow and FPVA results.
