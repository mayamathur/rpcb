

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              PRELIMINARIES                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rm(list=ls())

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


root.dir = "~/Dropbox/Personal computer/Independent studies/2020/RPCB reproducibility cancer biology"
raw.data.dir = paste(root.dir, "Raw data", sep="/")
prepped.data.dir = paste(root.dir, "Prepped data", sep="/")
code.dir = paste(root.dir, "Code (git)", sep="/")

setwd(code.dir)
source("helper.R")

# should View2() open tabs?
useView = FALSE

# should we run sanity checks?
run.sanity = FALSE

# read in paper-, experiment-, and outcome-level data
setwd(raw.data.dir)
# we won't actually be using the first of these
dp = read_xlsx("2020_10_5_raw_data.xlsx", sheet = "Paper level data"); nrow(dp)
de = read_xlsx("2020_10_5_raw_data.xlsx", sheet = "Experiment level data"); nrow(de)
do = read_xlsx("2020_10_5_raw_data.xlsx", sheet = "Outcome level data"); nrow(do)

##### Sanity Checks on Hierarchical Data Structure #####
# nesting levels: paper > experiment > outcome
names(de)
# exp-level data have 196 unique paper-exp combos:
uni( paste(de$`Paper #`, de$`Experiment #` ) )
# and 53 papers
uni(de$`Paper #`)

# outcome-level data have only 50 unique paper-exp combos:
uni( paste(do$`Paper #`, do$`Experiment #` ) ) 
# and 23 papers
uni(do$`Paper #`)

# outcome-level only contains ones with quantitative effect sizes:
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

# # this only had to be done once
# for ( i in c("Paper level data dictionary",
#              "Experiment level data dictionar",  # Excel cuts off last char
#              "Outcome level data dictionary") ){
#   cdNew = read_xlsx("2020_10_5_raw_data.xlsx", sheet = i)
#   cdNew$dataset = i
#   if ( i == "Paper level data dictionary" ) cd = cdNew else cd = rbind(cd, cdNew)
# }
# 
# setwd(raw.data.dir)
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
                   repDirection = "Observed difference in replication?",
                   
                   origN = "Original sample size",
                   repN = "Replication sample size",
                   
                   # point estimates that were sometimes available even when 
                   #  couldn't calculate effect sizes
                   origRawDiff = "Original point difference (for representative data)",
                   repRawDiff = "Replication raw difference (if original reported representative data)",
                   
                   # main effect sizes for us, but not yet on same scale
                   origES = "Original effect size",
                   origESLo = "Original lower CI",
                   origESHi = "Original upper CI",
                   origdf1 = "Original df1",
                   origdf2 = "Original df2",
                   origPval = "Original p value",
                   
                   repES = "Replication effect size",
                   repESLo = "Replication lower CI",
                   repESHi = "Replication upper CI",
                   repdf1 = "Replication df1",
                   repdf2 = "Replication df2",
                   repPval = "Replication p value",
                   
                   EStype = "Effect size type",
                   
                   # raw moderators
                   expType = "Type of experiment",
                   # has strings that can be used to determine type of lab:
                   labsContracted = "Lab(s) contracted for the experiment",
                   materialsRequested = "Key materials asked to be shared",
                   responseQuality = "Quality of response from original authors",
                   # when this is NA, no mods were needed:
                   changesNeededProse = "If modifications were needed for experiment to proceed, what were they?" )

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

# make dummies from variables coded as comma-separated categories
# @IMPORTANT: note that response quality coding exists even if no materials were requested
table( is.na(d2$materialsRequested), d2$responseQuality )  
d2 = recode_checkboxes(.d = d2,
                       var = "materialsRequested")
d2 = d2 %>% rename( "reqAntibodies"  = "Antibodies",
                    "reqCells" = "Cells",
                    "reqPlasmids" = "Plasmids" )


# recode lab type
hasCRO = whichStrings( x = d2$labsContracted, pattern = "CRO" )
hasCore = whichStrings( x = d2$labsContracted, pattern = "Core" )
d2$labType = NA
d2$labType[ hasCRO & !hasCore ] = "c.CRO only"
d2$labType[ hasCRO & hasCore ] = "b.Both"
d2$labType[ !hasCRO & hasCore ] = "a.Core only"
# sanity check
table(d2$labsContracted, d2$labType)

# recode experiment type as animal vs. non-animal
d2$expAnimal = (d2$expType == "Animal")

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

# @@maybe move this?
# @@also, why are "repSignif", "origSignif" here?
# moderators

modVars = c("responseQuality", "expType", "labsContracted",
            #"changesNeededProse", 
         #"repSignif", "origSignif",
         "reqAntibodies", "reqCells", "Organisms", 
         "reqPlasmids", "labType", "expAnimal")

CreateTableOne( dat = d2 %>% select(analysisVars) %>%
                  select( -c("changesNeededProse", stringsWith("ID", names(d2) ) ) ) )



# save intermediate dataset for easy debugging
write_interm(d2, "intermediate_dataset_step2.csv")


################################ 3. CALCULATE SMDS ################################ 

# read back in
d2 = read_interm("intermediate_dataset_step2.csv")

##### Sanity checks on effect types #####

t = d2 %>% group_by(EStype, Statistical.test.applied.to.original.data) %>%
  summarise(n())

# @IMPORTANT: these are not necessarily Pearson's r! 
#  some are Wilcoxon, so aren't comparable
# for the Pearson's r ones, find out if exposure was binary
#  or continuous
t = d2 %>% filter( !is.na(EStype) & EStype == "r" ) %>%
  select(pID,
         eID, 
         oID,
         Statistical.test.applied.to.original.data,
         Original.test.statistic.value,
         What.statistical.test.was.reported.,
         Original.test.statistic.type,
         origES) 
View2(t)
# original statistical test (note how many are Wilcoxon or Spearman)
t %>% group_by(Statistical.test.applied.to.original.data) %>%
  summarise(n())
# paper 12 is coded as EStype = r but has none of the other info
# others are often t-tests, ANOVA, Wilcoxon tests, or Spearman
# does that mean Spearman's were coded as r just like Pearson?

# @save this to show Tim
setwd(prepped.data.dir)
setwd("Intermediate work")
write.csv(t, "rows_with_correlations.csv")


# these make sense
SMDtypes = c("Cliff's delta", "Cohen's d", "Cohen's dz", "Glass' delta")
t = d2 %>% filter( !is.na(EStype) & EStype %in% SMDtypes ) %>%
  select(EStype, 
         pID,
         eID, 
         oID,
         Statistical.test.applied.to.original.data,
         What.statistical.test.was.reported.,
         Original.test.statistic.type,
         origES)
View2(t)
# original statistical test (note: @a lot of NAs here)
t %>% group_by(Statistical.test.applied.to.original.data) %>%
  summarise(n())

# HRs
t = d2 %>% filter( !is.na(EStype) & EStype %in% "Hazard ratio" ) %>%
  select(pID,
         eID, 
         oID,
         Statistical.test.applied.to.original.data,
         What.statistical.test.was.reported.,
         Original.test.statistic.type,
         origES)
# @paper 44: they fit a Cox regression but reported a t-test, but we have a HR?
# Tim confirmed this was a mistake on the original authors' part
# original statistical test (note: @a lot of NAs here)
t %>% group_by(Statistical.test.applied.to.original.data) %>%
  summarise(n())


##### Which ones can eventually be converted to Cohen's d type measure? #####

# Most outcomes were already measured on a standardized mean difference scale (e.g., Cohen’s d, Cohen’s w, Glass’ delta). We approximately converted other effect size measures to standardized mean differences for all analyses (i.e., Pearson’s correlations with continuous independent variables per Mathur & VanderWeele, 2020b; hazard ratios per VanderWeele, 2019 and Hasselblad & Hedges, 1995; and Cohen’s w via XXX).

table(d2$EStype)

# @for Tim: ones for which we want means and SDs
temp = c("between-subjects ANOVA", "Kruskal-Wallis rank sum test", "Wilcoxon signed-rank test", "Spearman's rank correlation")
d2$peoID[ ( d2$EStype %in% c("Cliff's delta", "Cohen's w" ) | d2$Statistical.test.applied.to.original.data %in% temp ) ]


# how many can be converted to Cohen's d type SMD? (of the completed pairs only)
canConvert = c("Cohen's d", "Glass' delta", "Hazard ratio", "r")
table(d2$EStype %in% canConvert & !is.na(d2$origES) & !is.na(d2$repES) )  # `102 of 337 can convert at this point 

# @TEMP: prepare to remove rank-based correlations
d2$canConvert = ( d2$EStype %in% c("Cohen's d", "Glass' delta", "Hazard ratio", "r") ) & 
  ( !d2$Statistical.test.applied.to.original.data %in% c("Spearman's rank correlation", "Wilcoxon signed-rank test") ) & 
  !is.na(d2$origES) & !is.na(d2$repES)
# now only 84 can convert after removing rank correlations
table(d2$canConvert)


# Cohen's d and Glass' delta are comparable

# HR can be converted to the above

# r can be converted to the above IF it's Pearson r

# Cliff's delta is rank-based measure of how often the values in one distribution are larger than
#  the values in a second distribution, so is not comparable to the above
# linear transformation of Mann-Whitney U-stat
# SE of delta: https://www.real-statistics.com/non-parametric-tests/mann-whitney-test/cliffs-delta/
View2(d2 %>% filter(EStype == "Cliff's delta"))

# Cohen's dz is a within-subject measure that standardizes by the SD of the CHANGE scores
#  and so is not comparable

# Cohen's w is a chi-square measure that can be converted IF it was a 2-group design
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0010059
# see also https://stats.stackexchange.com/questions/256261/effect-size-and-standard-error-for-mann-whitney-u-statistic
# @how did they get SEs for this?
View2(d2 %>% filter(EStype == "Cohen's w"))
# paper 9



##### ES2: Converted to a scale that can be meta-analyzed, but NOT necessarily SMDs #####

# Most outcomes were already measured on a standardized mean difference scale (e.g., Cohen’s d, Cohen’s w, Glass’ delta). We approximately converted other effect size measures to standardized mean differences for all analyses (i.e., Pearson’s correlations with continuous independent variables per Mathur & VanderWeele, 2020b; hazard ratios per VanderWeele, 2019 and Hasselblad & Hedges, 1995; and Cohen’s w via XXX).

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

# @ASK TIM FOR THE ORIGINAL SES; THIS IS JUST A TEMPORARY Z APPROXIMATION:
# 1.96 * SE * 2 = full CI width
d2$origSE2 = ( d2$origESHi2 - d2$origESLo2 ) / ( 2 * qnorm(.975) )
d2$repSE2 = ( d2$repESHi2 - d2$repESLo2 ) / ( 2 * qnorm(.975) )

d2$origVar2 = d2$origSE2^2
d2$repVar2 = d2$repSE2^2


##### ES3: SMDs where possible; otherwise NA #####

# do me later

# convert 



# convert column-by-column
for (i in toConvert) {
  newName = paste(i, "2", sep="")
  
  temp = convert_to_ES2( d2[[i]], .EStype = d2$EStype )
  
  d2[[newName]] = temp$ES2
  
  # this will be the same for each i, so just do it once
  if (i == toConvert[1]) d2$ES2type = temp$ES2type
}

# write intermediate data
write_interm(d2, "intermediate_dataset_step3.csv")


################################ 4. POOL INTERNAL REPLICATIONS ################################ 

# read back in
d2 = read_interm("intermediate_dataset_step3.csv")


# number of internal replications per paper-exp-outcome replication
# used below in sanity checks
t = d2 %>%
  filter(!is.na(repES2) & !is.na(origES2)) %>%
  group_by(peoID) %>%
  summarise(n = n())
table(t$n)

##### Pool internal replication(s) via FE meta-analysis #####
# kept separate from d2 for now to facilitate sanity checks
# temp has same dimensions as d2
temp = d2 %>% group_by(pID, eID, oID) %>%
  mutate( FE = sum(repES2 * 1/repVar2) / sum(1/repVar2),
          FEvar = 1 / sum(1/repVar2) )

expect_equal( nrow(temp), nrow(d2) )

# # how come this doesn't work??
# x = d2 %>% group_by(pID, eID, oID) %>%
#   mutate( FE = as.numeric( rma.uni(yi = repES2, vi = repVar2, method = "FE")$b ) )

# # but this works?
# # works
# x = d2 %>% group_by(pID, eID, oID) %>%
#   mutate( FE = mean(repES2) )

##### Sanity checks #####

if ( run.sanity == TRUE ) {
  # manually reproduce for ones with 1 replication
  id = t$peoID[t$n == 1]
  # FE estimate should just be the single estimate
  expect_equal(temp$FE[d2$peoID %in% id], d2$repES2[d2$peoID %in% id])
  expect_equal(temp$FEvar[d2$peoID %in% id], d2$repVar2[d2$peoID %in% id])
  
  # check one with 2 internal replications
  id = t$peoID[t$n == 2][1]
  # point estimate
  expect_equal( unique(temp$FE[d2$peoID %in% id]),
                as.numeric( rma.uni( yi = d2$repES2[d2$peoID == id],
                                     vi = d2$repVar2[d2$peoID == id],
                                     method = "FE")$b ) )
  # variance
  expect_equal( unique(temp$FEvar[d2$peoID %in% id]),
                as.numeric( rma.uni( yi = d2$repES2[d2$peoID == id],
                                     vi = d2$repVar2[d2$peoID == id],
                                     method = "FE")$se^2 ) )
  
  # check that FE estimate is always equal to replication estimate
  #  when nInternal = 1, and otherwise is not
  temp$agrees = abs(temp$FE - temp$repES2) < 0.001
  table(temp$agrees, temp$nInternal, useNA = "ifany")
}



###### Overwrite the d2 variables with pooled ones #####
d2$repES2 = temp$FE
d2$repVar2 = temp$FEvar
d2$repSE2 = sqrt(temp$FEvar)

# keep only 1 row per outcome (collapse over internal replications)
d3 = d2[ !duplicated(d2$peoID), ]

write_interm(d3, "intermediate_dataset_step4.csv")


################################ 5. POOL WITHIN EXPERIMENTS (FOR EXPERIMENT-LEVEL DATASET) ################################ 


# are moderators constant within peID?
# note exclusion of the repSign
nUniques = d3 %>%
  group_by(peID) %>%
  summarise( across(.cols = modVars,
                    .fns = uni ) ) %>%
  select(-peID)

# all unique within peIDs :)
any(nUniques  > 1)

# number of replications per experiment
# used below in sanity checks
t = d3 %>%
  filter(!is.na(repES2) & !is.na(origES2)) %>%
  group_by(peID) %>%
  summarise(n = n())
table(t$n)
# @@max 12 outcomes; sanity-check with Tim

##### Pool outcomes via RE meta-analysis #####
# note that this is RE meta-analysis rather than FE as for the internal replications
# kept separate from d2 for now to facilitate sanity checks
# temp has same dimensions as d2
temp = d3 %>%
  filter(!is.na(repES2) & !is.na(origES2)) %>%
  group_by(peID) %>%
  mutate( robu2(yi = repES2, vi = repVar2) )


##### Sanity checks #####

if ( run.sanity == TRUE ) {
  # manually reproduce the first one
  temp$peID[1]
  
  mod = robu( repES2 ~ 1,
              var.eff.size = repVar2,
              data = d3 %>% filter( peID == temp$peID[1]),
              studynum = 1:nrow(d3 %>% filter( peID == temp$peID[1])),
              small = TRUE )
  
  expect_equal( as.numeric(mod$b.r), unique( temp$est[ temp$peID == temp$peID[1]] ) )
  expect_equal( as.numeric(mod$reg_table$SE^2), unique( temp$var[ temp$peID == temp$peID[1]] ) )        
}            


###### Overwrite the d2 variables with pooled ones #####

# experiment-level version of d3
d3e = d3

d3e$repES2 = temp$est
d3e$repVar2 = temp$var
d3e$repSE2 = sqrt(temp$FEvar)

# keep only 1 row per outcome (collapse over internal replications)
d3e = d3e[ !duplicated(d2$peID), ]

# @@ any other vars that change with pe need to be dealt with here, 
#  such as replication p-value, etc.
# BEWARE: variables like this are currently the FIRST value within that peID, 
#  not the ones from pooling



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
# responseQuality

# make exclusions
# @@check with Tim
# remove pairs for which we have no info at all about original
# we also decided to exclude null originals
# 139 rows total after this
d3 = d3 %>% filter( !is.na(origDirection ) & origDirection == "Positive" )
nrow(d3)


# first level of granularity: outcome-level
setwd(prepped.data.dir)
fwrite(d3, "prepped_outcome_level_data.csv")

# second level of granularity: experiment-level
# @@not done yet


# calculate prediction intervals from original

