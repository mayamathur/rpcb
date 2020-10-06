



# With internal replications, just pool them (fixed-effects meta-analysis) to get a single estimate

# Exclude pairs with null originals for primary analyses

# merge moderators into all datasets used for analysis


# Questions for Tim:
# - Experiment-level data have 196 rows, but outcome-level have only 191 rows
#  - Why does outcome-level data have half as many papers as experiment-level? Is that because of non-quantitative ones?
# - Coded materials requested as dummies. Note that quality of response exists even when no materials were requested

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

root.dir = "~/Dropbox/Personal computer/Independent studies/2020/RPCB reproducibility cancer biology"
raw.data.dir = paste(root.dir, "Raw data", sep="/")
prepped.data.dir = paste(root.dir, "Prepped data", sep="/")
code.dir = paste(root.dir, "Code (git)", sep="/")

setwd(code.dir)
source("helper.R")

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
setwd(prepped.data.dir)
setwd("Intermediate work")
write.csv(d2, "intermediate_dataset_step1.csv")




################################ RECODE VARIABLES AND MAKE NEW ONES ################################ 

# read back in
setwd(prepped.data.dir)
setwd("Intermediate work")
d2 = read.csv("intermediate_dataset_step1.csv")

# look at the analysis variables
# easy trick to find analysis variables: their names don't have periods (converted from space upon reading in data)
( analysisVars = names(d2)[ !whichStrings("[.]", names(d2)) ] )
analysisVars = analysisVars[ !analysisVars %in% c("Notes.x", "Notes.y")]

# cast as numeric
( makeNumeric = stringsWith( c("origES", "repES", "Diff", "df", "origN", "repN", "Pval"), names(d2) ) )
# bm: this is being a butt
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

# look at the analysis variables
# easy trick to find analysis variables: their names don't have spaces
( analysisVars = names(d2)[ !whichStrings("[.]", names(d2)) ] )
analysisVars = analysisVars[ !analysisVars %in% c("Notes")]

CreateTableOne( dat = d2 %>% select(analysisVars) %>%
                  select( -c("changesNeededProse", stringsWith("ID", names(d2) ) ) ) )



# save intermediate dataset for easy debugging
setwd(prepped.data.dir)
setwd("Intermediate work")
write.csv(d2, "intermediate_dataset_step2.csv")


################################ CALCULATE SMDS ################################ 

# read back in
setwd(prepped.data.dir)
setwd("Intermediate work")
d2 = read.csv("intermediate_dataset_step2.csv")

##### Sanity checks on effect types #####


# for the Pearson's r ones, find out if exposure was binary
#  or continuous
t = d2 %>% filter( !is.na(EStype) & EStype == "r" ) %>%
  select(pID,
         eID, 
         oID,
         Statistical.test.applied.to.original.data,
         What.statistical.test.was.reported.,
         Original.test.statistic.type,
         origES) 
View(t)
# paper 12 is coded as EStype = r but has none of the other info
# others are often t-tests, ANOVA, Wilcoxon tests, or Spearman
# does that mean Spearman's were coded as r just like Pearson?

# save this to show Tim
setwd(prepped.data.dir)
setwd("Intermediate work")
write.csv(t, "rows_with_correlations.csv")


# these make sense
SMDtypes = c("Cliff's delta", "Cohen's d", "Cohen's dz", "Glass' delta")
View( d2 %>% filter( !is.na(EStype) & EStype %in% SMDtypes ) %>%
        select(EStype, 
               pID,
               eID, 
               oID,
               Statistical.test.applied.to.original.data,
               What.statistical.test.was.reported.,
               Original.test.statistic.type,
               origES) )


View( d2 %>% filter( !is.na(EStype) & EStype %in% "Hazard ratio" ) %>%
        select(pID,
               eID, 
               oID,
               Statistical.test.applied.to.original.data,
               What.statistical.test.was.reported.,
               Original.test.statistic.type,
               origES) )
# @paper 44: they fit a Cox regression but reported a t-test, and we have a HR?


##### ES2: Converted to a scale that can be meta-analyzed, but NOT necessarily SMDs #####

# Most outcomes were already measured on a standardized mean difference scale (e.g., Cohen’s d, Cohen’s w, Glass’ delta). We approximately converted other effect size measures to standardized mean differences for all analyses (i.e., Pearson’s correlations with continuous independent variables per Mathur & VanderWeele, 2020b; hazard ratios per VanderWeele, 2019 and Hasselblad & Hedges, 1995; and Cohen’s w via XXX).

table(d2$EStype)

# Cohen's d and Glass' delta are comparable

# Cliff's delta is a measure of how often the values in one distribution are larger than the values in a second distribution
# so is not comparable to the above

# Cohen's dz is a within-subject measure that standardizes by the SD of the CHANGE scores
#  and so is not comparable


toConvert = c("origES", "origLo", "origHi", "repES", "repLo", "repHi")

for (i in toConvert) {
  newName = paste(i, "2", sep="")
  
  # bm: need to handle NAs in the correlations
  d2[[newName]] = convert_to_ES2( d2[[i]], .EStype = d2$EStype )
}
  
  
 
  
##### ES3: SMDs where possible; otherwise NA #####





# check agreement of pvalues with repDirection after converting to SMDs
# (because that way null will always be 0)












