
## Notes on analysis datasets

* The outcome-level dataset `prepped_outcome_level_data.csv` (created by `data_prep.R` and called `do` in `analyze.R`) includes the 136 outcomes for which the estimate direction of the original effect was coded as "Positive" and wasn't `NA`, and hence includes non-quantitative pairs (i.e., those for which we didn't have a quantitative effect size for either the original or the replication). There is a [codebook](https://github.com/mayamathur/rpcb/blob/master/Prepped%20data/codebook_for_prepped_data.csv) for the analysis dataset.
    - Variables prefaced with "pw." are various pairwise metrics describing replication success.
    - The CI limits for ES2 and ES3 were only used as an intermediate step toward approximating the final SE rather than as standalone CIs. we recommend against interpreting these intermediate CIs directly. This is because we calculated these intermediate CI limits by directly transforming CI limits from the native scale, which can yield a CI that is not ideal because of symmetry considerations (see `data_prep.R` for details).

* The dataset `prepped_outcome_level_data_pw_metrics.csv` (created by `analyze.R`) is like `prepped_outcome_level_data.csv`, but includes pairwise metrics of replication success calculated for each pair that had the relevant information for that particular metric, and left `NA` otherwise. This dataset still includes non-quantitative pairs.

* The dataset `prepped_exp_level_data_pw_metrics.csv` (created by `analyze.R` and called `de` there) aggregates the pairwise metrics within each experiment and contains all 48 experiments represented in `do` (i.e., `prepped_outcome_level_data.csv`). For a given experiment, the pairwise metrics are `NA` if any of the outcomes within that experiment did not allow calculation of that pairwise metric.

* The table `pw_metrics_table_exp_level.csv` is like `de` (i.e., `prepped_exp_level_data_pw_metrics.csv`), except that it includes only quantitative pairs (i.e., 33 experiments) and is formatted using character strings for prettiness in the manuscript.  

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



