
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
do = fread(" ")

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


update_codebook_row( "repSignif", c("bin", "Was replication p<0.05?") )
update_codebook_row( "origSignif", c("bin", "Was original p<0.05?") )

update_codebook_row( "quantPair", c("bin", "Did we have quantitative ES for both original and replication?") )

update_codebook_row( "origES2", c("num", "A meta-analyzable effect size obtained from ES (e.g., log-HR instead of HR), but NOT necessarily an SMD") )
update_codebook_row( "ES2Type", c("char", "Scale of ES2") )


update_codebook_row( "origES3", c("num", "A meta-analyzable effect size on the SMD scale") )

update_codebook_row( "pw.PIRepInside", c("bin", "Was replication inside 95% PI?") )
update_codebook_row( "pw.PIRepInsid.sens", c("bin", "Was replication inside 95% PI, allowing for hypothetical heterogeneity?") )

update_codebook_row( "pw.Porig", c("num", "p-value for original inconsistency with replication") )
update_codebook_row( "pw.PorigSens", c("num", "p-value for original inconsistency with replication, allowing for hypothetical heterogeneity") )

update_codebook_row( "pw.ratio", c("num", "origES3 / repES3") )
update_codebook_row( "pw.diff", c("num", "origES3 - repES3") )

update_codebook_row( "pw.PsigAgree1", c("num", "Expected probability of significance agreement under null, assuming no heterogeneity (Mathur & VanderWeele)") )

update_codebook_row( "pw.PsigAgree1.sens", c("num", "Expected probability of significance agreement under null, with imputed heterogeneity (Mathur & VanderWeele)") )

update_codebook_row( "pw.FEest", c("num", "Fixed-effects estimate pooling original and replication") )


# save it
setwd(prepped.data.dir)
fwrite(cb2, "codebook_for_prepped_data.csv")



# META-REGRESSION --------------------------------------------- 

# moderators are at experiment-level, but analysis is at outcome level with 
#  nesting to handle experiment-level and paper-level correlation structure

# We will report the above metrics for each pair. Additionally, to summarize the above three metrics across pairs while accounting for their possible non-independence, we will robustly meta-regress each metric in a manner that accounts for clustering within original studies (Hedges et al., 2010). This model provides asymptotically correct inference even when the clustering structure is misspecified, which is important here because of the difficulty of precisely specifying the complex nature of clustering. This will yield average values of Porig, the difference, and the fixed-effects pooled estimate across all pairs. 

# Moderators:
# Animal vs. non-animal
# Type of replication lab (contract research organization [CRO] vs. academic core lab)
# What was requested from original authors and the response? (scored subjectively; Likert scale)
# Was a post hoc modification to protocol needed to complete the experiment? (col W in experiment level)



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
                       "pw.PsigAgree1",
                       "pw.PsigAgree1.sens")

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
             pw.SigAgree = 100* mean(repSignif == origSignif &
                                       repDirection == origDirection),
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
             
             origES3 = round( mean(origES3, 2) ),
             repES3 = round( mean(repES3), 2 ),
             
             FEest = round( mean(pw.FEest), 2 ),
             Ratio = round( mean(pw.ratio), 2 ),
             
             PIRepInside = n_perc_string(pw.PIRepInside),
             PIRepInside.sens = n_perc_string(pw.PIRepInside.sens), 
            
             Porig = format.pval( harmonic_p(pw.Porig), digits = 2, eps = "0.0001" ),
             Porig.sens = format.pval( harmonic_p( pw.PorigSens ), digits = 2, eps = "0.0001" ),
             
             # overall proportion (within this experiment) expected to agree
             SigAgree = n_perc_string( repSignif == origSignif & repDirection == origDirection),
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

#@: All of this is working for the outcome-level data.
#  If we want to run all of this at experiment level as well,
#  will need to add some things to its creation (see notes above when making de dataset)
#@: Could eliminate the extra levels here because I don't think we're going to use them. 

analysisLevels = c("exp_level", "outcome_level")


for ( l in analysisLevels ) {
  
  #@test only
  l = "outcome_level"
  
  if ( l == "outcome_level" ) dat = do
  if ( l == "exp_level" ) dat = de
  
  
 # ~ Calculate one-off summary metrics for manuscript ------------------
  
  # Percent sign agreement: The percentage of replications whose estimates agree in direction with the original study. This could be heuristically compared to the 50% that would be expected by chance if the null holds exactly in every replication (i.e., no effect heterogeneity) and conditional on all originals’ being positive in sign.
  
  #@ can remove after removing the additional analysis levels
  if ( l == "outcome_level" ) {

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
    
  }
  
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
                     value = mean_CI( dat$repSignif == dat$origSignif &
                                        dat$repDirection == dat$origDirection,
                                      cluster = dat$pID ) )
  
  update_result_csv( name = "Mean PsigAgree1 outcome_level",
                     value = round( mean(dat$pw.PsigAgree1, na.rm = TRUE), 2 ) )
  
  update_result_csv( name = "Mean PsigAgree1.sens outcome_level",
                     value = round( mean(dat$pw.PsigAgree1.sens, na.rm = TRUE), 2 ) )
  
}  # end giant loop over analysis levels








  # # ~ RPP-style scatterplot ------------------
  # #@NOT CHECKED BECAUSE I DON'T THINK IT'S STILL IN USE
  # 
  # # exclude 2 really extreme originals because they mess up plot scaling
  # dp = droplevels( dat %>% dplyr::filter(quantPair == TRUE) %>%
  #                    filter(origES3 < 50) )
  # 
  # # # randomly sample for testing purposes
  # # set.seed(2)
  # # dp = dp %>% group_by(ES2type) %>% sample_n( 3, replace = TRUE )
  # # #dp = dat[1:10,]
  # dp$plotID = dp$peoID # with eye toward functionizing
  # 
  # #@move this
  # dp$expType.pretty = NA
  # dp$expType.pretty[ dp$expAnimal == TRUE ] = "Animal"
  # dp$expType.pretty[ dp$expAnimal == FALSE ] = "Not animal"
  # 
  # 
  # min( c(dp$origES3, dp$repES3), na.rm = TRUE )
  # xmin = -4
  # max( c(dp$origES3, dp$repES3), na.rm = TRUE ) 
  # xmax = 30
  # 
  # shapes = c(19, 1)
  # 
  # p = ggplot( data = dp,
  #             aes(x = origES3,
  #                 y = repES3,
  #                 color = expType.pretty,
  #                 shape = repES3 < origES3 ) ) + 
  #   
  #   # null
  #   geom_abline(intercept = 0,
  #               slope = 1,
  #               lty = 2,
  #               color = "gray") +
  #   
  #   geom_hline( yintercept = 0,
  #               lty = 1,
  #               color = "gray" ) +
  #   
  #   geom_vline( xintercept = 0,
  #               lty = 1,
  #               color = "gray" ) +
  #   
  #   geom_point( size = 2.4,
  #               #pch = 1,
  #               alpha = 1) +
  #   
  #   # basic prettifying
  #   theme_bw() +
  #   theme( panel.grid.major=element_blank(),
  #          panel.grid.minor=element_blank() ) +
  #   
  #   scale_color_manual( values = colors ) +
  #   scale_shape_manual( values = shapes ) +
  #   scale_x_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, 5) ) +
  #   scale_y_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, 5) ) +
  #   
  #   labs( color = "Experiment type",
  #         shape = "Original > replication" ) +
  #   xlab("Original study estimate (SMD)") +
  #   ylab("Replication study estimate (SMD)")
  # 
  # 
  # # save it
  # setwd(results.dir)
  # setwd("Main figures")
  # ggsave( paste( "plot_scatter", "_", l, ".pdf", sep = "" ),
  #         width = 6.5,
  #         height = 5)
  # 
  # 
  # 
  # ################# WATERFALL PLOT OF DIFFERENCES ################
  # #@NOT CHECKED BECAUSE I DON'T THINK IT'S STILL IN USE
  # 
  # 
  # dp = dp %>%
  #   arrange( desc(pw.diff) )
  # 
  # dp$ind = 1:nrow(dp)
  # 
  # 
  # 
  # p = ggplot( ) +
  #   
  #   # null
  #   geom_hline(yintercept = 0,
  #              lty = 2,
  #              color = "gray") +
  #   
  #   # color-coded by experiment type
  #   geom_point( data = dp,
  #               aes(x = ind,
  #                   y = pw.diff,
  #                   color = expType.pretty) ) +
  #   
  #   geom_errorbar( data = dp,
  #                  aes(x = ind,
  #                      ymin = pw.diff - qnorm(.975) * sqrt(pw.diffVar),
  #                      ymax = pw.diff + qnorm(.975) * sqrt(pw.diffVar),
  #                      color = expType.pretty),
  #                  alpha = 0.4) +
  #   
  #   
  #   # basic prettifying
  #   theme_bw() +
  #   theme( panel.grid.major=element_blank(),
  #          panel.grid.minor=element_blank() ) +
  #   
  #   scale_color_manual( values = colors ) +
  #   #scale_x_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, 5) ) +
  #   #scale_y_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, 5) ) +
  #   
  #   labs( color = "Experiment type" ) +
  #   xlab("Pair") +
  #   ylab("Replication - original estimate (SMD)")
  # 
  # 
  # # save it
  # setwd(results.dir)
  # setwd("Main figures")
  # ggsave(paste( "plot_waterfall_diffs", "_", l, ".pdf", sep = "" ),
  #        width = 10,
  #        height = 5)
  # 
  # 
  # 
  # 
  # ################# DUMBBELL PLOT OF DIFFERENCES ################
  # #@NOT CHECKED BECAUSE I DON'T THINK IT'S STILL IN USE
  # 
  # 
  # # great info on dumbbell plot in ggplot:
  # # https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a
  # 
  # # dp = droplevels( dat %>% dplyr::filter( !is.na(origES2) & !is.na(repES2) & ES2type != "" ) %>%
  # #                    dplyr::filter(ES2type %in% c("Cohen's d", "Glass' delta", "Log hazard ratio", "Cohen's dz") ) )
  # # # randomly sample for testing purposes
  # # set.seed(2)
  # # dp = dp %>% group_by(ES2type) %>% sample_n( 3, replace = TRUE )
  # # #dp = dat[1:10,]
  # dp$plotID = dp$peoID # with eye toward functionizing
  # 
  # 
  # repColor <- "#0171CE"
  # origColor <- "#DE4433"
  # digits = 2
  # 
  # 
  # dp = dp %>% arrange( desc(origES3) )
  # 
  # # @move this?
  # dp$pw.Porig.cat = ">= 0.05"
  # dp$pw.Porig.cat[ dp$pw.Porig < 0.005 ] = "< 0.005"
  # dp$pw.Porig.cat[ dp$pw.Porig >= 0.005 & dp$pw.Porig < 0.05 ] = "< 0.05"
  # Porig.colors = c("red", "black", "gray")
  # 
  # 
  # p = ggplot() + 
  #   
  #   coord_flip() +
  # 
  #   
  #   # null
  #   geom_vline(xintercept = 0,
  #              lty = 2,
  #              color = "gray") +
  #   
  #   # geom_segment(data = dp,
  #   #              aes(y=plotID,
  #   #                  yend=plotID,
  #   #                  x=0,
  #   #                  xend=.5),
  #   #              color="#b2b2b2",
  #   #              size=0.15) +
  #   
  #   geom_dumbbell(data=dp,
  #                 aes(y=1:nrow(dp),
  #                     x=origES3,
  #                     xend=repES3,
  #                     color = pw.Porig.cat),
  #                 size=1.5,
  #                 
  #                 #color="#b2b2b2",
  #                 size_x=3,
  #                 size_xend = 3,
  #                 colour_x = origColor,
  #                 colour_xend = repColor) +
  #   
  #   scale_color_manual( values = Porig.colors ) +
  #   labs(color = "Porig") +
  #   
  #   
  #   # basic prettifying
  #   theme_bw() +
  #   theme( panel.grid.major=element_blank(),
  #          panel.grid.minor=element_blank() ) +
  #   
  #   xlab("Point estimate") +
  #   ylab("Pair")
  #   
  # 
  # # save it
  # setwd(results.dir)
  # setwd("Main figures")
  # ggsave( paste( "plot_dumbbell", "_", l, ".pdf", sep = "" ),
  #         width = 12,
  #         height = 5)
  
  







# # OLD (BUT DEFINITELY NEED CODE, LIKE FOR WATERFALL):
# ##### Upper Panel: Ordered ratios ######
# 
# dp = dp %>% filter( !is.na(ESgroup) ) %>%
#   arrange(desc(repES2 / origES2))
# 
# dp$ind = 1:nrow(dp)
# 
# 
# # only use CI if variance of ratio is less than 10
# 
# p = ggplot() +
#   
#   # null
#   geom_hline(yintercept = 100,
#              lty = 2,
#              color = "gray") +
#   
#   
#   geom_point( data = dp,
#               aes(x = ind,
#                   y = 100*(repES2 / origES2),
#                   color = ESgroup) ) +
#   
#   # geom_errorbar( data = dp,
#   #                aes(ymin = pw.ratio - qnorm(.975) * pw.ratioVar,
#   #                    ymax = pw.ratio + qnorm(.975) * pw.ratioVar,
#   #                    color = ESgroup) ) +
#   
#   scale_y_log10() +
#   
#   # basic prettifying
#   theme_bw() +
#   theme( panel.grid.major=element_blank(),
#          panel.grid.minor=element_blank() ) +
#   
#   ylab("Replication percent of original (logged axis)")
# 
# 
# 
# #################################### DUMBBELL PLOT################
# 
# 
# # great info on dumbbell plot in ggplot:
# # https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a
# 
# dat = read_interm("intermediate_analysis_dataset_step2.csv")
# 
# 
# library(ggalt)
# library(tidyverse)
# 
# dp = droplevels( dat %>% dplyr::filter( !is.na(origES2) & !is.na(repES2) & ES2type != "" ) %>%
#                    dplyr::filter(ES2type %in% c("Cohen's d", "Glass' delta", "Log hazard ratio", "Cohen's dz") ) )
# # randomly sample for testing purposes
# set.seed(2)
# dp = dp %>% group_by(ES2type) %>% sample_n( 3, replace = TRUE )
# #dp = dat[1:10,]
# dp$plotID = dp$peoID # with eye toward functionizing
# 
# 
# repColor <- "#0171CE"
# origColor <- "#DE4433"
# digits = 2
# 
# stringPanelWidth = 3
# panelSpacer = 1
# 
# max( c(dp$origES2, dp$repES2), na.rm = TRUE )  # use this to inform stringPanel1Start
# stringPanel1Start = 15  # an x-axis coordinate
# stringPanel1End = stringPanel1Start + stringPanelWidth  # x-axis 
# 
# stringPanel2Start = stringPanel1End + panelSpacer
# stringPanel2End = stringPanel2Start + stringPanelWidth
# 
# 
# # for labeling joy, have a df with just the first row to appear in plot
# # which is actually the LAST row in the LAST ES2type
# temp = dp[ dp$ES2type == unique(dp$ES2type)[1], ]  
# # sort in order to match ggplot's ordering
# dpRow = temp[ temp$peoID == rev( sort(temp$peoID) )[1], ] 
# #dpRow = dp[ dp$plotID == dp$plotID[ nrow(dp) ], ]
# #dpRow = dp[ dp$plotID == dp$plotID[ 1 ], ]
# 
# 
# p = ggplot() + 
#   
#   facet_grid(ES2type ~ .,
#              scales = "free",
#              space = "free"
#              #space = "free_y"
#   ) +
#   
#   # null
#   geom_vline(xintercept = 0,
#              lty = 2,
#              color = "gray") +
#   
#   geom_segment(data = dp,
#                aes(y=plotID,
#                    yend=plotID,
#                    x=0,
#                    xend=.5),
#                color="#b2b2b2",
#                size=0.15) +
#   
#   geom_dumbbell(data=dp,
#                 aes(y=plotID,
#                     x=origES2,
#                     xend=repES2),
#                 size=1.5,
#                 color="#b2b2b2",
#                 size_x=3,
#                 size_xend = 3,
#                 colour_x = origColor,
#                 colour_xend = repColor) +
#   
#   # "Original" and "Replication" labels
#   # data arg: only label one of the rows
#   geom_text( data= dpRow,
#              aes(x=repES2,
#                  y=plotID,
#                  label="Replication"),
#              color=repColor,
#              size=3,
#              vjust=-1.5,
#              fontface="bold") +
#   
#   geom_text( data= dpRow,
#              aes(x=origES2,
#                  y=plotID,
#                  label="Original"),
#              color=origColor,
#              size=3,
#              vjust=-1.5,
#              fontface="bold") +
#   
#   # numerical labels
#   geom_text(data=dp,
#             aes(x=repES2,
#                 y=plotID,
#                 label= round(repES2, digits) ),
#             color=repColor,
#             size=2.75,
#             vjust=2.5) +
#   
#   geom_text(data=dp,
#             aes(x=origES2,
#                 y=plotID,
#                 label= round(origES2, digits) ),
#             color=origColor,
#             size=2.75,
#             vjust=2.5) +
#   
#   ##### differences panel
#   geom_rect(data=dp,
#             aes(xmin=stringPanel1Start,  # hard-coded location
#                 xmax=stringPanel1End,
#                 ymin=-Inf,
#                 ymax=Inf),
#             fill="grey") +
#   
#   geom_text(data=dp,
#             aes(label = round( origES2 - repES2, digits ),
#                 y=plotID,
#                 x= mean( c(stringPanel1Start, stringPanel1End) ) ),
#             #fontface="bold",
#             size=3) +
#   # header
#   geom_text(data=dpRow, 
#             aes(x=mean( c(stringPanel1Start, stringPanel1End) ),  # needs to match above
#                 y=plotID,
#                 label="Original - replication"),
#             color="black",
#             size=3.1,
#             vjust=-2,
#             fontface="bold") 
# 
# 
# 
# # bm
# # Porig panel
# p + geom_rect(data=dp,
#               aes(xmin=stringPanel2Start,  # hard-coded location
#                   xmax=stringPanel2End,
#                   ymin=-Inf,
#                   ymax=Inf),
#               fill="grey") +
#   
#   geom_text(data=dp,
#             aes(label = round( pw.Porig, digits ),
#                 y=plotID,
#                 x=mean( c(stringPanel2Start, stringPanel2End) ) ),
#             #fontface="bold",
#             size=3) +
#   
#   geom_text(data=dpRow, 
#             aes(x=mean( c(stringPanel2Start, stringPanel2End) ),  # needs to match above
#                 y=plotID,
#                 #label = TeX("$P_{orig}$")
#                 label="Porig"
#             ),
#             #label = "asdfd",
#             #label = TeX("$P_{orig}$"),
#             #label = expression(paste("DOC (mg ", L^-1,")")),
#             color="black",
#             size=3.1,
#             vjust=-2,
#             fontface="bold") + 
#   
#   # basic prettifying
#   theme_bw() +
#   theme( panel.grid.major=element_blank(),
#          panel.grid.minor=element_blank() ) +
#   
#   xlab("Point estimate") +
#   ylab("Paper, experiment, and outcome")
# 
# 
# 
# 
# 
# 
# 
# # # 15 x 9 works well
# # # save
# # setwd(objects.dir)
# # ggsave( "forest.pdf",
# #         width = 15,
# #         height = 10,
# #         units = "in" )
# # setwd(results.overleaf)
# # ggsave( "forest.pdf",
# #         width = 15,
# #         height = 9,
# #         units = "in" )
# 





################################ GROUP METRICS ################################ 

# Primary: Percentage of replication estimates >0 (but will be an overestimate due to overdispersion)
# @replace with calibrated version?

# FIGURE: Density plot of calibrated estimates (outcome level) for group of replications and group of originals, overlaid, and facetted by EStype






