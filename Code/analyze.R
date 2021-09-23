
# NOTES ---------------------------------------------


# High-level summary of this script:
#  - Calculates pairwise metrics at the outcome level
#  - Then aggregates these at the experiment level to create an experiment-level dataset
#  - Conducts remaining analyses at both levels of analysis (outcome, experiment)

# See READ-ME at https://github.com/mayamathur/rpcb for information about the datasets and code

# Useful fns for interacting with this script (see helper.R):
# - searchBook("p value")
# - vr()
# - wr()


# PRELIMINARIES ---------------------------------------------

rm(list=ls())

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
renv::restore()


library(readxl)
library(dplyr)
library(ggplot2)
library(MetaUtility)
library(robumeta)
library(testthat)
library(Replicate)
library(data.table)
library(metafor)
library(here)
library(ggalt)
library(tidyverse)
library(here)
library(tableone)
library(corrr)
library(clubSandwich)
library(lme4)

# run this only if you want to update the R environment specs:
# renv::snapshot()

# define working directories
root.dir = here()
raw.data.dir = here("Raw data")
prepped.data.dir = here("Prepped data")
code.dir = here("Code")
results.dir = here("Results from R")

# no sci notation
options(scipen=999)

# for plots
colors = c("red", "black")

# source helper fns
setwd(code.dir)
source("helper.R")

# read in prepped data
# this dataset includes all 
setwd(prepped.data.dir)
do = fread("prepped_outcome_level_data.csv")

# read in codebook for easy searching
setwd(raw.data.dir)
cd = fread("codebook_merged.csv")

# Olsson-Collentine data for sensitivity analyses with imputed heterogeneity
# we got this dataset by running their publicly available R code
setwd( here("Auxiliary data") )
dOlsson = fread("het_df_from_their_code.csv")


# WORK ON OUTCOME-LEVEL DATA: CALCULATE PAIRWISE METRICS ---------------------------------------------

# ~ Estimate average heterogeneity in Olsson-Collentine for use in sensitivity analyses


# sanity check
# reported heterogeneity medians: 0.047 for correlations and 0.09 for SMDs
dOlsson %>% group_by(effect_type) %>%
  summarise(median(tau),
            mean( abs(eff_size) ),
            # how big is tau as a ratio vs. the effect size magnitude?
            tauRatio = median(tau / abs(eff_size) ) )
# medians match :)

# variance of tau so I can meta-analyze them
dOlsson$tauSE = ( dOlsson$tau_ci.ub - dOlsson$tau ) / qnorm(.975)
dOlsson$tauVar = dOlsson$tauSE^2
# sanity check
expect_equal( dOlsson$tau[1] + qnorm(.975)*dOlsson$tauSE[1],
              dOlsson$tau_ci.ub[1] )

# average heterogeneity across metas for SMDs only
( m.SMD = robu( tau ~ 1,
                var.eff.size = tauSE^2,
                studynum = rp,
                data = dOlsson[ dOlsson$effect_type == "SMD", ],
                model = "HIER",
                small = TRUE ) )

# heterogeneity value to use in imputed analyses
( t2.imp = m.SMD$mod_info$tau.sq )


# ~ Calculate pairwise metrics ---------------------------------------------

# all the "pw.XXX" metrics are calculating using the ES3's (i.e., SMDs)
# retain non-quant pairs since we will also create plots and metrics for those below
do = do %>% 
  rowwise() %>% 
  mutate( analyze_one_row(origES3,
                          origVar3, 
                          repES3,
                          repVar3,
                          repDirection,
                          t2 = t2.imp) )


# save it
setwd(prepped.data.dir)
fwrite(do, "prepped_outcome_level_data_pw_metrics.csv")

# ~~ Sanity checks for analyze_one_row ---------------------------------------------

# check homogeneous prediction interval
yio = do$origES3
vio = do$origVar3
yir = do$repES3
vir = do$repVar3
pooled.SE = sqrt(vio + vir)

# check which entries should be NA
expect_equal( is.na(yio) | is.na(vio) | is.na(yir) | is.na(vir),
              is.na(do$pw.PIRepInside) )

expect_equal( yio - qnorm(0.975) * pooled.SE,
              do$pw.PILo)
expect_equal( yio + qnorm(0.975) * pooled.SE,
              do$pw.PIHi)
expect_equal( yir >= do$pw.PILo & yir <= do$pw.PIHi,
              do$pw.PIRepInside)

# check heterogeneous prediction interval
pooled.SE = sqrt(vio + vir + c(t2.imp))

expect_equal( yio - qnorm(0.975) * pooled.SE,
              do$pw.PILo.sens)
expect_equal( yio + qnorm(0.975) * pooled.SE,
              do$pw.PIHi.sens)
expect_equal( yir >= do$pw.PILo.sens & yir <= do$pw.PIHi.sens,
              do$pw.PIRepInside.sens)

# check homogeneous Porig
denom = sqrt(vio + vir)
Z = (abs(yio - yir))/denom
pval = as.numeric(2 * (1 - pnorm(Z)))
expect_equal( do$pw.Porig, pval )

# check hetero Porig
denom = sqrt(vio + vir + c(t2.imp))
Z = (abs(yio - yir))/denom
pval = as.numeric(2 * (1 - pnorm(Z)))
expect_equal( do$pw.PorigSens, pval )

# compare homo and hetero results
mean(do$pw.PIRepInside, na.rm = TRUE); mean(do$pw.PIRepInside.sens, na.rm = TRUE)
median(do$pw.Porig, na.rm=TRUE); median(do$pw.PorigSens, na.rm=TRUE)

# check difference
expect_equal( yio - yir, do$pw.diff )
expect_equal( vio + vir, do$pw.diffVar )

# check observed significance agreement
expect_equal( do$repDirection == "Positive",
              do$pw.observedSigAgree )

# check homogeneous PsigAgree (expected significance agreement)
pooled.SE = sqrt( vio + vir )
# because all yio>0:
P = 1 - pnorm( ( qnorm(.975) * sqrt( vir ) - yio ) / pooled.SE )
expect_equal( P,
              do$pw.PsigAgree1 )

# check heterogeneous PsigAgree (expected significance agreement)
pooled.SE = sqrt( 2 * c(t2.imp) + vio + vir )
# because all yio>0:
P2 = 1 - pnorm( ( qnorm(.975) * sqrt( vir ) - yio ) / pooled.SE )
expect_equal( P2,
              do$pw.PsigAgree1.sens )

# compare homo and hetero results
# very close for this one
mean(do$pw.PsigAgree1, na.rm = TRUE); mean(do$pw.PsigAgree1.sens, na.rm = TRUE)

# check FEest
expect_equal( (yio/vio + yir/vir) / (1/vio + 1/vir),
              do$pw.FEest )
expect_equal( 1 / (1/vio + 1/vir),
              do$pw.FEestVar )


# MAKE CODEBOOK --------------------------------------------- 

# this includes only the new variables added to outcome-level data
newVars = names(do)[104: length(names(do))]

# codebook
cb2 = data.frame( variable = newVars,
                  type = NA,
                  description = NA )

# vec: vector of type and description
update_codebook_row = function(var, vec) {
  cb2[ cb2$variable == var, c("type", "description") ] <<- vec
  
}


update_codebook_row( "repSignif", c("bin", "Was replication p<0.05? (This and all other 'rep' variables in this dataset are regarding the FE pooled estimate when there were multiple internal replications.)") )
update_codebook_row( "origSignif", c("bin", "Was original p<0.05?") )

update_codebook_row( "quantPair", c("bin", "Did we have quantitative ES for both original and replication?") )

update_codebook_row( "origES2", c("num", "A meta-analyzable effect size obtained from ES (e.g., log-HR instead of HR), but NOT necessarily an SMD") )
update_codebook_row( "ES2Type", c("char", "Scale of ES2") )


update_codebook_row( "origES3", c("num", "A meta-analyzable effect size on the SMD scale") )

update_codebook_row( "pw.PIRepInside", c("bin", "Was replication inside 95% PI?") )
update_codebook_row( "pw.PIRepInside.sens", c("bin", "Was replication inside 95% PI, allowing for hypothetical heterogeneity?") )

update_codebook_row( "pw.Porig", c("num", "p-value for original inconsistency with replication") )
update_codebook_row( "pw.PorigSens", c("num", "p-value for original inconsistency with replication, allowing for hypothetical heterogeneity") )

update_codebook_row( "pw.ratio", c("num", "origES3 / repES3") )
update_codebook_row( "pw.diff", c("num", "origES3 - repES3") )

update_codebook_row( "pw.observedSigAgree", c("num", "Observed significance agreement") )

update_codebook_row( "pw.PsigAgree1", c("num", "Expected probability of significance agreement under null, assuming no heterogeneity (Mathur & VanderWeele)") )

update_codebook_row( "pw.PsigAgree1.sens", c("num", "Expected probability of significance agreement under null, with imputed heterogeneity (Mathur & VanderWeele)") )

update_codebook_row( "pw.FEest", c("num", "Fixed-effects estimate pooling original and replication") )


# save it
setwd(prepped.data.dir)
fwrite(cb2, "codebook_for_prepped_data.csv")



# META-REGRESSION --------------------------------------------- 

# moderators are at experiment-level, but analysis is at outcome level with 
#  nesting to handle experiment-level and paper-level correlation structure


# ~ Moderator Summary Table and Correlation Matrix --------------------------------------------- 

modVars = c("expAnimal",
            "hasCROLab",
            "hasCoreLab",
            "materialsShared",
            "infoQuality",
            "changes")

CreateTableOne(vars = modVars, data = do)

# moderator correlation matrix
# for categorical mods, simplify by just looking at one level for the corr matrix
modVarsBin = c("expAnimal",
               "hasCROLab",
               "hasCoreLab",
               "materialsShared_b.Yes",
               "infoQuality",
               "changes_b..LT.moderate.success")

corrs = do %>%
  filter(quantPair == TRUE) %>%
  select(modVarsBin) %>%
  correlate( use = "pairwise.complete.obs" ) %>%
  stretch() %>%
  arrange(desc(r)) %>%
  group_by(r) %>%
  filter(row_number()==1)

corrs$r = round(corrs$r, 2)

corrs = corrs[ !is.na(corrs$r), ]
View(corrs)

# save it
setwd(results.dir)
setwd("Additional tables and figures")
write.csv(corrs, "moderator_cormat.csv")


# ~ Analyze Pairwise Metrics That *Do* Have Variances --------------------------------------------- 

# exclude changes variable from meta-regression because it's too homogeneous
#  (all but one outcome are in the same category)
#  so it breaks the meta-regression
table(do$changes)
modVars = c("expAnimal",
            "hasCROLab",
            "hasCoreLab",
            "materialsShared",
            "infoQuality")


outcomesWithVar = c( "pw.diff",
                     "pw.FEest")

# clear the results table to be created
if ( exists("modTable") ) rm(modTable)

for ( i in outcomesWithVar ) {
  
  modTable = safe_analyze_moderators(  .dat = do[ do$quantPair == TRUE, ],
                                       yi.name = i,
                                       # below assumes a standardized naming convention for
                                       #  variances of the pairwise metrics:
                                       vi.name = paste(i, "Var", sep=""),
                                       
                                       # cut out the "pw." part of outcome name
                                       analysis.label = gsub("^.*?\\.","",i),
                                       
                                       modVars = modVars,
                                       
                                       # here's where we control the number of tests counted in Bonferroni
                                       n.tests = length(modVars),
                                       digits = 2 )
  
}


modTable
table(modTable$Problems)


# ~~ Sanity check for one of the outcomes  --------------------------------------------- 
formString = paste( "pw.diff ~ ", paste( modVars, collapse= " + ") )
formString2 = paste( formString, " + (1 | pID / eID)" )

temp = do[ , c("pID", "eID", "pw.diff", "pw.diffVar", modVars ) ]
temp = temp[ complete.cases(temp), ]

V_mat = impute_covariance_matrix(vi = temp$pw.diffVar,
                                 cluster = temp$pID, 
                                 r = 0.6)  # just a working "guess"

model = rma.mv( eval( parse(text=formString2) ),
                V = V_mat,
                random = ~ 1 | pID / eID,
                data = temp,
                sparse = TRUE)

# robust inference
res = conf_int(model, vcov = "CR2")
pvals = coef_test(model, vcov = "CR2")$p_Satt

# visually compare model-based inference to robust inference
# pretty similar
model$pval; pvals

resString = paste( round( model$b, 2 ), 
                   " ",
                   format_CI( res$CI_L, 
                              res$CI_U,
                              2),
                   sep = "" )

# check all stats for this outcome
expect_equal( resString, modTable$Est[ modTable$Analysis == "diff" ] )
expect_equal( as.character( round(pvals, 2) ), modTable$Pval[ modTable$Analysis == "diff" ] ) # Bonferroni p-values
myBonf = pmin(1, pvals * length(modVars) )
expect_equal( format_stat(myBonf, 2),
              modTable$Pval.Bonf[ modTable$Analysis == "diff" ] ) 


# Analyze Pairwise Metrics That *Don't* Have Variances --------------------------------------------- 

outcomesWithoutVar = c("pw.ratio",
                       "pw.PIRepInside",
                       "pw.PIRepInside.sens",
                       "pw.Porig",
                       "pw.PorigSens",
                       "pw.observedSigAgree")

for ( i in outcomesWithoutVar ) {
  modTable = safe_analyze_moderators(  .dat = do[ do$quantPair == TRUE, ],
                                       yi.name = i,
                                       # below assumes a standardized naming convention for
                                       vi.name = NA,
                                       
                                       # cut out the "pw." part of outcome name
                                       analysis.label = gsub("^.*?\\.","",i),
                                       
                                       modVars = modVars,
                                       
                                       n.tests = length(modVars),
                                       digits = 2 )
}


# for some reason, the FE problems get deleted
View(modTable)
table(modTable$Problems)


setwd(results.dir)
setwd("Main tables")
write.csv(modTable, "moderator_regressions_outcome_level.csv")


# ~~ Sanity check for one of the outcomes  --------------------------------------------- 
formString = paste( "pw.Porig ~ ", paste( modVars, collapse= " + ") )
formString2 = paste( formString, " + (1 | pID / eID)" )

temp = do[ , c("pID", "eID", "pw.Porig", modVars ) ]
temp = temp[ complete.cases(temp), ]

model = lmer( eval( parse(text=formString2) ),
              data = temp )

# get robust variances
res = conf_int(model, vcov = "CR2")



# robust inference
res = conf_int(model, vcov = "CR2")
pvals = coef_test(model, vcov = "CR2")$p_Satt

resString = paste( round( fixef(model), 2 ), 
                   " ",
                   format_CI( res$CI_L, 
                              res$CI_U,
                              2),
                   sep = "" )

# check all stats for this outcome
expect_equal( resString, modTable$Est[ modTable$Analysis == "Porig" ] )
# this one is a PITA because of sci notation:
# expect_equal( as.character( round(pvals, 2) ),
#               modTable$Pval[ modTable$Analysis == "Porig" ] )
# Bonferroni p-values
myBonf = pmin(1, pvals * length(modVars) )
expect_equal( format_stat(myBonf, 2),
              modTable$Pval.Bonf[ modTable$Analysis == "Porig" ] ) 



# CREATE EXPT-LEVEL DATASET AND TABLE --------------------------------------------- 

# table of pairwise metrics aggregated at the experiment level (~50 rows)

# look at the outcomes to be aggregated
stringsWith( pattern = "pw", x = names(do) )


# ~ Expt-level dataset --------------------------------------------- 

# includes all outcomes that were in dataset "do"

# entries are numerical and not rounded for plotting, analysis, etc.
# this DOES include the non-quant pairs for plotting purposes
# keep pw.XXX variable names the same to facilitate automated plotting below
de = do %>%
  #filter( quantPair == TRUE ) %>%
  group_by(peID) %>%
  summarise( nOutcomes = n(),
             
             pw.FEest = mean(pw.FEest),
             pw.ratio = mean(pw.ratio),
             
             pw.PIRepInside = mean(pw.PIRepInside),
             pw.PIRepInside.sens = mean(pw.PIRepInside.sens), 
             # note: would have preferred p.hmp here because it's asymptotically exact,
             #  but it was giving error messages
             pw.Porig = harmonic_p(pw.Porig),
             pw.Porig.sens = harmonic_p(pw.PorigSens),
             
             # overall proportion (within this experiment) expected to agree
             # note that we're not requiring repSignif = origSignif because some 
             # originals were interpreted as "positive" even though they weren't p<0.05
             # and this dataset already retains only originals coded as "positive"
             # also no need to condition on repSignif because repDirection already includes this
             pw.SigAgree = 100* mean(repDirection == origDirection),
             pw.PercSigAgree1 = 100 * mean(pw.PsigAgree1),
             pw.PercSigAgree1.sens = 100 * mean(pw.PsigAgree1.sens)
  ) 

View(de)


# save it
setwd(prepped.data.dir)
fwrite(de, "prepped_exp_level_data_pw_metrics.csv") 


# ~ Expt-level table (prettified version of above) -----------------

# table does NOT include non-quant pairs
expTable = do %>%
  filter( quantPair == TRUE ) %>%
  group_by(peID) %>%
  summarise( nOutcomes = n(),
             
             origES3 = round( mean(origES3), 2 ),
             repES3 = round( mean(repES3), 2 ),
             
             FEest = round( mean(pw.FEest), 2 ),
             Ratio = round( mean(pw.ratio), 2 ),
             
             PIRepInside = n_perc_string(pw.PIRepInside),
             PIRepInside.sens = n_perc_string(pw.PIRepInside.sens), 
             
             Porig = format.pval( harmonic_p(pw.Porig), digits = 2, eps = "0.0001" ),
             Porig.sens = format.pval( harmonic_p( pw.PorigSens ), digits = 2, eps = "0.0001" ),
             
             # overall proportion (within this experiment) expected to agree
             SigAgree = n_perc_string( repDirection == origDirection ),
             PercSigAgree1 = paste( round( 100 * mean(pw.PsigAgree1), 0 ), "%", sep ="" ),
             PercSigAgree1.sens = paste( round( 100 * mean(pw.PsigAgree1.sens), 0 ), "%", sep ="" )
  ) 

View(expTable)

# save it
setwd(results.dir)
setwd("Main tables")
fwrite(expTable, "pw_metrics_table_exp_level.csv")



# AT MULTIPLE LEVELS OF ANALYSIS: SUMMARY PLOTS AND STATS -----------------

# read back in
setwd(prepped.data.dir)
de = fread("prepped_exp_level_data_pw_metrics.csv")
do = fread("prepped_outcome_level_data_pw_metrics.csv")


# pick which dataset to use for these analyses
# previously we thought we might do the analyses at both expt and outcome level,
#  but ultimately we are only doing the latter
dat = do


# ~~ Counts ------------------
update_result_csv( name = "n all pairs outcome_level",
                   value = nrow(dat) )

update_result_csv( name = "n (perc) non-quant pairs outcome_level",
                   value = n_perc_string(dat$quantPair == FALSE) )

update_result_csv( name = "n (perc) quant pairs outcome_level",
                   value = n_perc_string(dat$quantPair == TRUE) )

update_result_csv( name = "n (perc) same direction all pairs outcome_level",
                   value = n_perc_string(dat$repDirection == "Positive") )

# sanity check on n_perc_string
myString = paste( sum(dat$repDirection == "Positive"), 
                  " (",
                  round( 100 * mean(dat$repDirection == "Positive") ),
                  "%)", 
                  sep = "")
expect_equal( myString,
              n_perc_string(dat$repDirection == "Positive") )
# end sanity check



# ~~ Prediction intervals ------------------
update_result_csv( name = "prop pw.PIRepInside outcome_level",
                   value = mean_CI( dat$pw.PIRepInside == TRUE,
                                    cluster = dat$pID ) )

update_result_csv( name = "prop pw.PIRepInside.sens outcome_level",
                   value = mean_CI(dat$pw.PIRepInside.sens == TRUE,
                                   cluster = dat$pID) )

update_result_csv( name = "Median Porig outcome_level",
                   value = round( median(dat$pw.Porig, na.rm = TRUE), 3 ) )


# sanity check on mean_CI
mod = lm(dat$pw.PIRepInside.sens ~ 1) 
expect_equal( mean(dat$pw.PIRepInside.sens, na.rm=TRUE),
              as.numeric(mod$coefficients) )
Vmat = vcovCR(mod,
              cluster = dat$pID,
              type = "CR2")
CIs = conf_int(mod, vcov = Vmat)
myString = paste( round( as.numeric(mod$coefficients), 2 ), 
                  " [",
                  round( CIs$CI_L, 2 ),
                  ", ", 
                  round( CIs$CI_U, 2 ),
                  "]",
                  sep = "")
expect_equal( myString,
              mean_CI(dat$pw.PIRepInside.sens == TRUE,
                      cluster = dat$pID) )
# end sanity check


# ~~ Porig ------------------
update_result_csv( name = "Harmonic mean Porig outcome_level",
                   value = round( harmonic_p(dat$pw.Porig), 4 ) )

update_result_csv( name = "prop Porig<0.005 outcome_level",
                   value = mean_CI(dat$pw.Porig < 0.005,
                                   cluster = dat$pID) )

update_result_csv( name = "prop Porig<0.05 outcome_level",
                   value = mean_CI(dat$pw.Porig < 0.05,
                                   cluster = dat$pID) )


update_result_csv( name = "Median PorigSens outcome_level",
                   value = round( median(dat$pw.PorigSens, na.rm = TRUE), 3 ) )

update_result_csv( name = "prop PorigSens<0.005 outcome_level",
                   value = mean_CI(dat$pw.PorigSens < 0.005,
                                   cluster = dat$pID) )

update_result_csv( name = "prop PorigSens<0.05 outcome_level",
                   value = mean_CI(dat$pw.PorigSens < 0.05,
                                   cluster = dat$pID) )


# ~~ Metrics on SMD scale ------------------
# update_result_csv( name = "Median SMD ratio outcome_level",
#                    value = round( median(dat$pw.ratio, na.rm = TRUE), 2 ) )

update_result_csv( name = "Mean SMD diff outcome_level",
                   value = mean_CI(dat$pw.diff,
                                   cluster = dat$pID) )

update_result_csv( name = "Mean SMD FEest outcome_level",
                   value = mean_CI(dat$pw.FEest,
                                   cluster = dat$pID) )

update_result_csv( name = "Mean SMD origES3 outcome_level",
                   value = mean_CI(dat$origES3,
                                   cluster = dat$pID) )

update_result_csv( name = "Mean SMD repES3 outcome_level",
                   value = mean_CI(dat$repES3,
                                   cluster = dat$pID) )


# ~~ Significance agreement ------------------

update_result_csv( name = "prop sigAgree outcome_level",
                   value = mean_CI( dat$repDirection == dat$origDirection,
                                    cluster = dat$pID ) )

update_result_csv( name = "Mean PsigAgree1 outcome_level",
                   value = round( mean(dat$pw.PsigAgree1, na.rm = TRUE), 2 ) )

update_result_csv( name = "Mean PsigAgree1.sens outcome_level",
                   value = round( mean(dat$pw.PsigAgree1.sens, na.rm = TRUE), 2 ) )

