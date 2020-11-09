# R function to estimate adjusted survival probability

## Background
Several statistical programs have been developed to compute direct adjusted survival estimates to consistently represent the results of the Cox proportional hazards model. However, when applied to analyze observational data with large sample sizes or highly stratified treatment groups, the existing programs are either incapable of generating confidence bands and simultaneous p values for comparing adjusted survival outcomes or inefficient to do so. These programs also do not give consideration to the potential occurrence of left truncation in retrospectively collected data. 

This R function was developed for producing direct adjusted survival estimates based on a stratified Cox model with improved computational performance and capability of handling left-truncated and right-censored time-to-event data, to address the increasing needs for analyzing large multi-center observational data sets.

## Function and parameters
``` R
diff.adj.surv.v3 <- function(indata, coxf, seednum=, n.sim=, cb.alpha=, starttime=NULL, endtime=NULL)
```
`indata`: input matrix containing event indicator, time to event and co-variates to be considered in stratified Cox model
`coxf`: Cox model in format of `Survival` object
`seednum`: seed number for replication of simulated results
`n.sim`: number of simulation round for generating survival confidence band
`cb.alpha`: alpha level (typicall 0.05)
`starttime`: start time for confidence band consideration
`endtime`: end time for confidence band consideration
