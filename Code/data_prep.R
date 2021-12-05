
# See READ-ME at https://github.com/mayamathur/rpcb for information about the datasets and code

# Useful fns for interacting with this script (see helper.R):
# - searchBook("p value")
# - vr()
# - wr()

# Important notes:
#   - Direction transformation of CI limits means intermediate CI limits for ES2 and ES3
#    may not be ideal (i.e., asymmetric), but that's because they're only an intermediate
#    step toward approximating the final SE



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                   PRELIMINARIES                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

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
library(data.table)
library(tableone)
library(qdapTools)
library(metafor)
library(fastDummies)
library(here)

# run this only if you want to update the R environment specs:
#renv::snapshot()

# define working directories
root.dir = here()
raw.data.dir = here("Raw data")
prepped.data.dir = here("Prepped data")
code.dir = here("Code")
results.dir = here("Results from R")

setwd(code.dir)
source("helper.R")

# should View2() open tabs?
useView = FALSE

# should we run sanity checks?
run.sanity = TRUE

# read in paper-, experiment-, and outcome-level data
setwd(raw.data.dir)
# we won't actually be using the first of these
dp = read_xlsx("2021_9_19_raw_data.xlsx", sheet = "Paper level data"); nrow(dp)
de = read_xlsx("2021_9_19_raw_data.xlsx", sheet = "Experiment level data"); nrow(de)
do = read_xlsx("2021_9_19_raw_data.xlsx", sheet = "Outcome level data"); nrow(do)


##### Sanity Checks on Hierarchical Data Structure #####
# nesting levels: paper > experiment > outcome
names(de)
# exp-level data have 193 unique paper-exp combos
uni( paste(de$`Paper #`, de$`Experiment #` ) )
# and 53 papers
uni(de$`Paper #`)

# outcome-level data have only 50 unique paper-exp combos
uni( paste(do$`Paper #`, do$`Experiment #` ) ) 
# and 23 papers
uni(do$`Paper #`)

# outcome-level only contains ones with quantitative effect sizes
table( is.na(do$`Original effect size`))

# confirm that outcome-level data have all papers from exp-level data for 
#  which the replication was completed
expect_equal( uni(de$`Paper #`[de$`Replication experiment completed` == "Yes"]),
              uni(do$`Paper #`) )


##### Look at What Info Is in Each Dataset #####
# variables that only appear in one or the other dataset
names(do)[ !names(do) %in% names(de) ]
names(de)[ !names(de) %in% names(do) ]  # moderators

# anything non-overlapping in paper-level data?
names(dp)[ !names(dp) %in% c(names(de), names(do)) ]

# variables in both datasets
names(do)[ names(do) %in% names(de) ]
names(dp)[ names(dp) %in% names(de) ]


##### Merge Datasets #####
d = merge( do, de, by = c("Paper #", "Experiment #"), all=TRUE )
nrow(d)
# because de has a superset of do's papers and experiments:
expect_equal( uni(d$`Paper #`), uni(de$`Paper #`) )
expect_equal( uni(d$`Experiment #`), uni(de$`Experiment #`) )

# this time don't keep papers that don't even give experiment-level data
# i.e., left-join
d = merge( d, dp, by = "Paper #", all.x = TRUE)
nrow(d)



##### Merge Codebooks for Easy Searching #####

# # this only has to be done when the codebooks are updated
# setwd(raw.data.dir)
# for ( i in c("Paper level data dictionary",
#              "Experiment level data dictionar",  # Excel cuts off last char
#              "Outcome level data dictionary") ){
#   cdNew = read_xlsx("2021_5_25_raw_data.xlsx", sheet = i)
#   cdNew$dataset = i
#   if ( i == "Paper level data dictionary" ) cd = cdNew else cd = rbind(cd, cdNew)
# }
# fwrite(cd, "codebook_merged.csv")

# read in codebook for easy searching
setwd(raw.data.dir)
cd = fread("codebook_merged.csv")


################################ 1. RENAME VARIABLES ################################ 

# keep variables required to calculate original's power
# moderators: put "mod" in their name for ease of analysis 

d2 = d %>% rename( pID = "Paper #",
                   eID = "Experiment #",
                   oID = "Effect #",
                   internalID = "Internal replication #",
                   
                   testReported = "Was a statistical test reported in the original paper?.x",
                   
                   origDirection = "Expected difference based on the original paper?",
                   repDirection = "Observed difference in replication? (SMD)",
                   
                   origN = "Original sample size",
                   repN = "Replication sample size",
                   
                   # point estimates that were sometimes available even when 
                   #  couldn't calculate effect sizes
                   origRawDiff = "Original point difference (for representative data)",
                   repRawDiff = "Replication raw difference (if original reported representative data) (SMD)",
                   
                   # main effect sizes for us, but not yet on same scale
                   origES = "Original effect size (SMD)",
                   origESLo = "Original lower CI (SMD)",
                   origESHi = "Original upper CI (SMD)",
                   origdf1 = "Original df1 (SMD)",
                   origdf2 = "Original df2 (SMD)",
                   origPval = "Original p value (SMD)",
                   origSE = "Original standard error (SMD)",
                   
                   # IMPORTANT: all of these are currently defined at the level of the 
                   #  internal replication, NOT the FE-pooled replication that we use in 
                   #  actual analyses
                   repES = "Replication effect size (SMD)",
                   repESLo = "Replication lower CI (SMD)",
                   repESHi = "Replication upper CI (SMD)",
                   repdf1 = "Replication df1 (SMD)",
                   repdf2 = "Replication df2 (SMD)",
                   repPval = "Replication p value (SMD)",
                   repSE = "Replication standard error (SMD)",
                   
                   EStype = "Effect size type (SMD)",
                   
                   # raw moderators
                   expType = "Type of experiment",
                   # has strings that can be used to determine type of lab:
                   labsContracted = "Lab(s) contracted for the experiment",
                   materialsSharedProse = "Key materials offered to be shared",  # "1" suffix because we're going to recode it
                   infoNeeded = "Clarifications asked of original authors",
                   #materialsRequested = "Key materials asked to be shared",
                   infoQuality = "Quality of response from original authors",
                   changesSuccess = "Changes able to be implemented during experimentation?",
                   changesNeededRaw = "Changes needed during experimentation?")

# save intermediate dataset for easy debugging
write_interm(d2, "intermediate_dataset_step1.csv")



################################ 2. RECODE VARIABLES AND MAKE NEW ONES ################################ 


# read back in
d2 = read_interm("intermediate_dataset_step1.csv")

# look at the analysis variables
# easy trick to find analysis variables: their names don't have periods (converted from space upon reading in data)
( analysisVars = names(d2)[ !whichStrings("[.]", names(d2)) ] )
analysisVars = analysisVars[ !analysisVars %in% c("Notes.x", "Notes.y")]

# cast as numeric
( makeNumeric = stringsWith( c("origES", "repES", "Diff", "df", "origN", "repN", "Pval"), names(d2) ) )
# this is being a butt
d2 = d2 %>% mutate_at( .vars = makeNumeric, .funs = as.numeric )
# warning is okay

# repDirection is about both direction and significance
# split it up
d2$repSignif = d2$repPval < 0.05
d2$origSignif = d2$origPval < 0.05
# sanity checks
# agreement with repDirection variable
table(d2$repSignif, d2$repDirection)
table(d2$origSignif, d2$origDirection)

# recode direction variables
d2$origDirection[ d2$origDirection == "" ] = NA
d2$repDirection[ d2$repDirection == "" ] = NA


# - recode "changes able to be implemented" variable
# - changesSuccess (formerly `Changes able to be implemented`) is basically the temporal successor to
#  changesNeededRaw (formerly `Changes needed during experimentation?`)
# - this has too many levels to work well as a moderator, so we are dichotomizing
#  at >=3 (i.e., moderately implemented or better)
# - changesNeededRaw (formerly `Changes needed during experimentation?`) is coded st 
#  0 (*not* NA) means no changes were needed 
# - and changesNeededRaw = 0 should always correspond to changesSuccess = 6 ("not needed")
expect_equal( d2$changesNeededRaw == 0, d2$changesSuccess == 6)

d2$changes = NA
d2$changes[ is.na(d2$changesSuccess) ] = "c. No changes needed"
d2$changes[ d2$changesSuccess < 3 ] = "b. LT moderate success" # "LT" = "less than"
d2$changes[ d2$changesSuccess >= 3 ] = "a. GTE moderate success"  # "GTE" = "greater than or equal"
# sanity check on recoding
table( d2$changes, d2$changesSuccess, useNA = "ifany")
# dummy-coded version of changes
# only used for making the moderator correlation matrix
d2 = dummy_cols(.data = d2, select_columns = "changes")


# recode "clarifications asked" variable
# this is basically the temporal precendent to infoQuality
table(d2$infoNeeded > 0)  # this means that they ALWAYS needed clarifications
# so could just leave infoQuality as continuous variable
# sanity check: if no info needed, info quality should be NA
expect_equal( d2$infoNeeded == FALSE, is.na(d2$infoQuality) )



# were (any) materials actually shared?
# regardless of whether they were requested
d2$materialsShared = NA
d2$materialsShared[ is.na(d2$materialsSharedProse) ] = "c. Not requested"
d2$materialsShared[ !is.na(d2$materialsSharedProse) & d2$materialsSharedProse == "No" ] = "a.No"
d2$materialsShared[ !is.na(d2$materialsSharedProse) & d2$materialsSharedProse != "No" ] = "b.Yes"
# sanity check
table( d2$materialsSharedProse, d2$materialsShared, useNA = "ifany" )
# dummy-coded version of materialsShared
# only used for making the moderator correlation matrix
d2 = dummy_cols(.data = d2, select_columns = "materialsShared")

# recode lab type
d2$hasCROLab = whichStrings( x = d2$labsContracted, pattern = "CRO" )
d2$hasCoreLab = whichStrings( x = d2$labsContracted, pattern = "Core" )


# recode experiment type as animal vs. non-animal
d2$expAnimal = (d2$expType == "Animal")

# indicator for being a completed quantitative pair
# should be 132 per Tim
d2$quantPair = !is.na(d2$origES) & !is.na(d2$repES)
expect_equal( sum(d2$quantPair), 132 )

# ID for paper-exp-outcome combos
d2$peoID = paste( "p", d2$pID, "e", d2$eID, "o", d2$oID, sep = "")
# number of internal replications per outcome
d2 = d2 %>% group_by(peoID) %>%
  mutate(nInternal = n())

# ID for paper-exp combos
# used for making the experiment-level data
d2$peID = paste( "p", d2$pID, "e", d2$eID, sep = "" )


# look at the analysis variables
# easy trick to find analysis variables: their names don't have spaces
( analysisVars = names(d2)[ !whichStrings("[.]", names(d2)) ] )
analysisVars = analysisVars[ !analysisVars %in% c("Notes, Organisms")]

# moderators
modVars = c("expAnimal",
            "hasCROLab",
            "hasCoreLab",
            "materialsShared",
            "infoQuality",
            "changes")

CreateTableOne( dat = d2 %>% select(analysisVars) %>%
                  select( -c("changesNeededRaw", stringsWith("ID", names(d2) ) ) ) )


# save intermediate dataset for easy debugging
write_interm(d2, "intermediate_dataset_step2.csv")



################################ 3. CALCULATE VARIOUS EFFECT SIZES (ES2 AND ES3) ################################ 

# read back in
d2 = read_interm("intermediate_dataset_step2.csv")

##### Sanity checks on effect types #####

# breakdown of ES types
# sanity check:
# Tim said "For awareness, we are at 90 Cohen's d, 15 Cohen's dz, 20 Glass' delta, 6 Pearson's r, and 7 Hazard ratios for effects where there are quantitative pairs (this is counting them as unique paper/experiment/effect/internal replication, so naturally this goes down as they are collapsed on any of those elements for analysis/visualization)."
# that totals 138
d2 %>% filter( quantPair == TRUE ) %>%
  group_by(EStype) %>%
  summarise(n())
# matches :)

# look at which statistical tests yielded which effect size types
t = d2 %>% group_by(EStype, Statistical.test.applied.to.original.data..SMD.) %>%
  summarise(n())
t



##### ES2: Converted to a scale that can be meta-analyzed, but NOT necessarily SMDs #####

# statistics that should be converted to ES2 scale
toConvert = c("origES", "origESLo", "origESHi", "repES", "repESLo", "repESHi")

# convert column-by-column
for (i in toConvert) {
  newName = paste(i, "2", sep="")
  
  temp = convert_to_ES2( d2[[i]], .EStype = d2$EStype )
  
  d2[[newName]] = temp$ES2
  
  # this will be the same for each i, so just do it once
  if (i == toConvert[1]) d2$ES2type = temp$ES2type
}

# sanity checks
data.frame( d2 %>% group_by(EStype, ES2type) %>%
              summarise( meanNA(origES), 
                         meanNA(origES2),
                         meanNA(repES),
                         meanNA(repES2) ) )

table(d2$ES2type)


##### Standard errors #####

# by default, use a simple Z approximation for the SEs
# 1.96 * SE * 2 = full CI width
d2$origSE2 = ( d2$origESHi2 - d2$origESLo2 ) / ( 2 * qnorm(.975) )
d2$repSE2 = ( d2$repESHi2 - d2$repESLo2 ) / ( 2 * qnorm(.975) )
# but for effect sizes that were directly calculated as SMDs (i.e., not 
#  converted from some other scale), we can use the "native" SEs directly:
ind = d2$EStype %in% c("Cohen's d",
                       "Cohen's dz",
                       "Glass' delta" ) 
d2$origSE2[ind] = d2$origSE[ind]

d2$origVar2 = d2$origSE2^2
d2$repVar2 = d2$repSE2^2


##### ES3: Convert all to approximate SMDs #####
# convert point estimates and CI limits
toConvert = c("origES2", "origESLo2", "origESHi2", "repES2", "repESLo2", "repESHi2")

library(stringr)
# convert column-by-column
for (i in toConvert) {
  # name columns with "3", e.g., "origES3", etc.
  newName = str_replace(string = i,
                        pattern = "2",
                        replacement = "3")
  
  temp = convert_to_ES3( d2[[i]], .ES2type = d2$ES2type )
  
  d2[[newName]] = temp$ES3
}



##### Sanity checks on point estimates and CI limits #####
# check each type
table(d2$EStype)

# check all SMDs
# they should have stayed the same upon conversion to ES3
ind = which( d2$EStype %in% c("Cohen's d",
                              "Cohen's dz",
                              "Glass' delta" ) )
expect_equal( d2$origES3[ind], d2$origES[ind] )
expect_equal( d2$repES3[ind], d2$repES[ind] )
expect_equal( d2$origESLo3[ind], d2$origESLo[ind] )

# check all HRs (common-outcome conversions)
ind = which( d2$EStype == "Hazard ratio" )
myRR = ( 1 - 0.5^sqrt(d2$origES[ind]) ) / ( 1 - 0.5^sqrt(1/d2$origES[ind]) )
myOR = sqrt(myRR)
mySMD = log(myOR) * sqrt(3) / pi
expect_equal( d2$origES3[ind], mySMD )
# CI limit
myRRLo = ( 1 - 0.5^sqrt(d2$origESLo[ind]) ) / ( 1 - 0.5^sqrt(1/d2$origESLo[ind]) )
myORLo = sqrt(myRRLo)
mySMDLo = log(myORLo) * sqrt(3) / pi
expect_equal( d2$origESLo3[ind], mySMDLo )

# check all Pearson's r
ind = which( d2$EStype == "Pearson's r" )
mySMD = r_to_d(r = d2$origES[ind],
               sx = 1,
               delta = 1)$d
expect_equal( d2$origES3[ind], mySMD )
# CI limits
mySMDLo = r_to_d(r = d2$origESLo[ind],
                 sx = 1,
                 delta = 1)$d
expect_equal( d2$origESLo3[ind], mySMDLo )


##### Calculate standard errors #####
# Important: For effect sizes that were extracted on a non-SMD
#  scale and then converted, we got SEs by transforming CI limits,
#  then using z approximation to back out the variance
#  when origES had a highly asymmetric CI, so will origES3
#   because we just transformed the CI limits
#  so the approach below approximates the SMD's variance by 
#  using the full CI width
# we did this instead of using delta method on variance because we
#  didn't have original variances
d2$origSE3 = ( d2$origESHi3 - d2$origESLo3 ) / ( 2 * qnorm(.975) )
d2$repSE3 = ( d2$repESHi3 - d2$repESLo3 ) / ( 2 * qnorm(.975) )

d2$origVar3 = d2$origSE3^2
d2$repVar3 = d2$repSE3^2



##### Sanity checks on SEs #####

# compare to delta-method SEs for one EStype (Pearson's R)
ind = which( d2$EStype == "Pearson's r" )
# CI limits
mySMD = r_to_d(r = d2$origES[ind],
               sx = 1,
               delta = 1,
               N = d2$origN[ind] )
# SEs are pretty similar
cbind( sqrt(d2$origVar3[ind]), mySMD$se )
# the CIs themselves are different because of aforementioned asymmetry, which makes sense
cbind( d2$origESLo3[ind], mySMD$lo )

# discrepancies between CI limits and SEs
# will be close but not necessarily exactly equal in the case of highly asymmetric CIs
d2$discrep = abs((d2$origES3 + qnorm(.975) * d2$origSE3) - d2$origESHi3)
summary(d2$discrep)

d2$badCI = ( d2$origES > d2$origESHi ) |  ( d2$origES < d2$origESLo )
expect_equal( TRUE %in% d2$badCI, FALSE )


# look at the initial and fully converted estimates and CI limits
temp = d2 %>% select( peoID, EStype, origES, origESLo, origESHi, origES3, origESLo3, origESHi3, discrep ) %>%
  mutate(loArm = origES3 - origESLo3, hiArm = origESHi3 - origES3 ) %>%
  arrange( desc(EStype, discrep) )

temp = temp %>% mutate_at( names(temp)[ !names(temp) %in% c("peoID", "EStype")],
                           function(x) round(x,2) )

# and for replications
temp = d2 %>% select( peoID, EStype, repES, repESLo, repESHi, repES3, repESLo3, repESHi3, discrep ) %>%
  mutate(loArm = repES3 - repESLo3, hiArm = repESHi3 - repES3 ) %>%
  arrange( desc(EStype, discrep) )

temp = temp %>% mutate_at( names(temp)[ !names(temp) %in% c("peoID", "EStype")],
                           function(x) round(x,2) )

##### Write Intermediate Dataset #####
write_interm(d2, "intermediate_dataset_step3.csv")


################################ 4. POOL INTERNAL REPLICATIONS ################################ 

# only doing this for ES3
# we can use just that scale throughout analyses

# read back in
d2 = read_interm("intermediate_dataset_step3.csv")


# number of internal replications per paper-exp-outcome replication
# used below in sanity checks
t = d2 %>%
  filter(!is.na(repES3) & !is.na(origES3)) %>%
  group_by(peoID) %>%
  summarise(n = n())
table(t$n)

##### Pool internal replication(s) via FE meta-analysis #####
# kept separate from d2 for now to facilitate sanity checks
# temp has same dimensions as d2
temp = d2 %>% group_by(pID, eID, oID) %>%
  mutate( FE = sum(repES3 * 1/repVar3) / sum(1/repVar3),
          FEvar = 1 / sum(1/repVar3),
          # Wald p-value
          FEpval = 2 * ( 1 - pnorm( abs( FE / sqrt(FEvar) ) ) ) )

expect_equal( nrow(temp), nrow(d2) )


##### Sanity checks #####
if ( run.sanity == TRUE ) {
  # manually reproduce for ones with 1 replication
  id = t$peoID[t$n == 1]
  # FE estimate should just be the single estimate
  expect_equal(temp$FE[d2$peoID %in% id], d2$repES3[d2$peoID %in% id])
  expect_equal(temp$FEvar[d2$peoID %in% id], d2$repVar3[d2$peoID %in% id])
  
  # note that the FE p-value may NOT agree with rep p-value even 
  #  when there is only 1 internal replication, presumably because single replications
  #  might not use a Wald test
  cbind( temp$FEpval[d2$peoID %in% id], d2$repPval[d2$peoID %in% id] )
  # compare them visually
  plot( temp$FEpval[d2$peoID %in% id], d2$repPval[d2$peoID %in% id] )
  
  # check one with 2 internal replications
  id = t$peoID[t$n == 2][1]
  # compare to metafor
  mod = rma.uni( yi = d2$repES3[d2$peoID == id],
                 vi = d2$repVar3[d2$peoID == id],
                 method = "FE")
  # point estimate
  expect_equal( unique(temp$FE[d2$peoID %in% id]),
                as.numeric( mod$b ) )
  # variance
  expect_equal( unique(temp$FEvar[d2$peoID %in% id]),
                as.numeric( mod$se^2 ) )
  
  # p-value
  expect_equal( unique(temp$FEpval[d2$peoID %in% id]),
                as.numeric( mod$pval ) )
  
  # check that FE estimate is always equal to replication estimate
  #  when nInternal = 1, and otherwise is not
  temp$agrees = abs(temp$FE - temp$repES3) < 0.001
  table(temp$agrees, temp$nInternal, useNA = "ifany")
}


###### Overwrite the d2 variables with pooled ones #####


# OVERWRITE repPval, etc. (initially defined at the level of internal replications)
#  to be based on the pooled internal replications

# as checked above in sanity check, OK to overwrite the estimate and SE/variance when 
#  there's only 1 internal replication...
d2$repES3 = temp$FE
d2$repVar3 = temp$FEvar
d2$repSE3 = sqrt(temp$FEvar)
# ...but we will NOT overwrite p-value when there's 1 internal replication because 
#  of Wald vs. other p-values issue
idMultiple = t$peoID[t$n > 1]
d2$repPval[ d2$peoID %in% idMultiple ] = temp$FEpval[ temp$peoID %in% idMultiple ]

# sanity checks on p-value overwriting
if ( run.sanity == TRUE ) {
  idSingle = t$peoID[t$n == 1]
  expect_equal( TRUE,
                any( d2$repPval[ d2$peoID %in% idSingle ] != temp$FEpval[ temp$peoID %in% idSingle ] ) )
  
 
  expect_equal( d2$repPval[ d2$peoID %in% idMultiple ],
                    temp$FEpval[ temp$peoID %in% idMultiple ] )
}


d2$repDirection = NA
d2$repDirection[ d2$repES3 > 0 & d2$repPval < 0.05 ] = "Positive"
d2$repDirection[ d2$repES3 > 0 & d2$repPval >= 0.05 ] = "Null-positive"
d2$repDirection[ d2$repES3 < 0 & d2$repPval >= 0.05 ] = "Null-negative"
d2$repDirection[ d2$repES3 < 0 & d2$repPval < 0.05 ] = "Negative"

d2$repSignif = d2$repPval < 0.05

# keep only 1 row per outcome (collapse over internal replications)
d3 = d2[ !duplicated(d2$peoID), ]


write_interm(d3, "intermediate_dataset_step4.csv")


################################ 6. SAVE ANALYSIS DATASETS ################################ 

# key variables for analysis:
#  - ES2 variables: converted to a scale that can be meta-analyzed, but not necessarily SMD
# repPval = "Replication p value",
# 
# EStype = "Effect size type",
# 
# # raw moderators
# expType = "Type of experiment",
# # has strings that can be used to determine type of lab:
# labsContracted = "Lab(s) contracted for the experiment",
# materialsRequested = "Key materials asked to be shared",
# infoQuality

# read back in
d3 = read_interm("intermediate_dataset_step4.csv")

# make exclusions
# remove pairs for which we have no info at all about original
# we also decided to exclude null originals
d3 = d3 %>% filter( !is.na(origDirection ) &
                      origDirection == "Positive" &
                      quantPair == TRUE )
expect_equal( 97, nrow(d3) )



# sanity check against Tim's count of 42 for significance agreement
# in several equivalent ways
expect_equal( sum(d3$repES3 > 0 & d3$repPval < 0.05), 42 )
expect_equal( sum(d3$repES3 > 0 & d3$repSignif == TRUE ), 42 )
expect_equal( sum(d3$repDirection == "Positive" ), 42 )


# outcome-level dataset
setwd(prepped.data.dir)
fwrite(d3, "prepped_outcome_level_data.csv")


