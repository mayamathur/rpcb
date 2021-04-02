
## Notes on analysis dataset

* There is a [codebook](https://github.com/mayamathur/rpcb/blob/master/Prepped%20data/codebook_for_prepped_data.csv) for the analysis dataset. Variables prefaced with "pw." are various pairwise metrics describing replication success.

* The dataset in "Auxiliary data" is from Olsson-Collentine's public repo and was used to impute the within-pair heterogeneity for sensitivity analyses.

## Notes on analysis methods

### Moderator analyses

We fit a separate model to each pairwise metric, regressing on all candidate moderators for each. Not all pairwise metrics have variances associated with them. We fit different, but comparable, models to pairwise metrics with and without variances. For outcomes that did have variances (pw.diff, pw.FEest), used rma.mv (nested random intercepts of eID within pID) with robust SEs under the CHE working model (multilevel random effects model with constant sampling correlation working model). For outcomes that did not have variances ("pw.ratio", "pw.PIRepInside", "pw.PIRepInside.sens", "pw.Porig", "pw.PorigSens", "pw.PsigAgree1"), used lmer with same specification as above and robust SEs. Robust SEs account are both for clustering and possibility of non-normal residuals (e.g., for ratios). 

## Notes on results files

All results files are in the directory ["Results from R"](https://github.com/mayamathur/rpcb/tree/master/Results%20from%20R). There are three kinds of results files:

### String results for manuscript

Stats that are to appear in the prose of the manuscript are written as strings to the file ["stats_for_paper.csv"](https://github.com/mayamathur/rpcb/blob/master/Results%20from%20R/stats_for_paper.csv).

### Figures

At outcome level, the plots omit 2 pairs with extreme originals (SMD > 80). 


### Tables

The table [pw_metrics_table_exp_level.csv](https://github.com/mayamathur/rpcb/blob/master/Results%20from%20R/Main%20tables/pw_metrics_table_exp_level.csv) summarizes results by aggregating to the experiment level. Numerical and binary variables are summarized as means (or proportions) over outcomes within each experiment. Porig is a harmonic mean p-value over outcomes. 

The table [moderator_regressions_outcome_level.csv](https://github.com/mayamathur/rpcb/blob/master/Results%20from%20R/Main%20tables/moderator_regressions_outcome_level.csv) shows the results of meta-regressions on the moderators.


## To do


