

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
library(Replicate)
library(data.table)
library(metafor)

root.dir = "~/Dropbox/Personal computer/Independent studies/2020/RPCB reproducibility cancer biology"
raw.data.dir = paste(root.dir, "Raw data", sep="/")
prepped.data.dir = paste(root.dir, "Prepped data", sep="/")
code.dir = paste(root.dir, "Code (git)", sep="/")
results.dir = paste(root.dir, "Results from R", sep="/")

options(scipen=999)

setwd(code.dir)
source("helper.R")


setwd(prepped.data.dir)
d = fread("prepped_outcome_level_data.csv")

# read in codebook for easy searching
setwd(raw.data.dir)
cd = fread("codebook_merged.csv")

# d = read_xlsx("RP_CB Final Analysis .xlsx")
# 
# names(d)
# 
# 
# dim(d)  # 258
# 
# table(d$`Replication attempted`)  # 233
# table(d$`Experiment completed`) # 190
# 
# d %>% filter( `Experiment completed` == "Yes" ) %>%
#   summarise( length(unique(`Original study title`)),
#              mean(`Number of lab(s) contracted for the entire study`) )
# 
# data.frame( d %>% group_by(`Original study title`) %>%
#               summarise( n(),
#                          comp = mean(`Experiment completed` == "Yes"),
#                          expN0 = max(`Experiment #`, na.rm = TRUE),
#                          max(`Study #`, na.rm = TRUE) ) )


#################################### 

# We will conduct the main analyses at two levels of granularity (outcome-level and experiment-level). 


# distinctions in analysis:
#  - completed pairs (has realized replication outcome) vs. all pairs

# 

# with eye toward functionizing everything later (so we can do for outcome- and exp-level):
dat = d


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                     METRICS FOR ALL PAIRS (INCL NON-COMPLETED)                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Percent sign agreement: The percentage of replications whose estimates agree in direction with the original study. This could be heuristically compared to the 50% that would be expected by chance if the null holds exactly in every replication (i.e., no effect heterogeneity) and conditional on all originals’ being positive in sign.

table(dat$origDirection, dat$repDirection, useNA = "ifany")

mean( dat$repDirection == "Positive" )
sum( dat$repDirection == "Positive" )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#              2. MAIN PAIRWISE METRICS (COMPLETED QUANT PAIRS)                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

################################ CALCULATE PAIRWISE METRICS ################################ 

# Primary: Prediction interval and whether replication falls inside it (assuming t2=0)
# Make forest plot of these, including Porig on the side

# Primary: Porig

# Primary: Ratio of original to replication study estimates

# Secondary: P(Replication p < .05) with 2 expectation benchmarks
#  - from JRSSA paper
#  - the true effect size in each original is equal to the effect size for which it would
#   have had 80% power.

# Porig with assumed zero heterogeneity

# Fixed-effects pooled estimate: A meta-analytic pooled estimate of the original and replication estimates from each pair. 

# As sensitivity analysis:
# Porig and pred interval with imputed heterogeneity: If there is moderate or high within-pair effect heterogeneity, this could make the original studies appear less consistent with the replications than they truly are. As a sensitivity analysis, we will impute the average heterogeneity estimate from Olsson-Collentine, et al. (in press), which was tau=0.13 on the SMD scale, and use this to re-calculate Porig for each pair. These values of Porig will likely be large (i.e., indicating better consistency) than those from main analyses.


# dataset has already been aggregated to either the outcome- or the experiment-level
# if this hits errors, it's likely that you have mismatches in the NA dataframe vs. the filled-in one in analyze_one_row
dat = dat %>% 
  rowwise() %>% 
  mutate( analyze_one_row(origES2,
                          origVar2, 
                          repES2,
                          repVar2,
                          ES2type) )

# # ratio sanity checks:
# # @@note that some ratios and their variances are extremely large:
# inds = which(dat$pw.ratioVar > 100000)
# dat$pw.ratio[inds]
# 
# i = 5
# origES2 = dat$origES2[i]
# origVar2 = dat$origVar2[i]
# repES2 = dat$repES2[i]
# repVar2 = dat$repVar2[i]
# 
# deltamethod( ~ x1/x2,
#              mean = c( origES2, repES2 ), cov = diag( c(origVar2, repVar2) ) )^2
# 
# 
# for ( i in 1:nrow(dat) ) {
#   # debug any rows that are brats
#   chunk = analyze_one_row( dat$origES2[i],
#                    dat$origVar2[i],
#                    dat$repES2[i],
#                    dat$repVar2[i],
#                    dat$ES2type[i])
#   if ( i == 1 ) res = chunk
#   if ( i > 1) res = rbind(res, chunk)
# 
# }

write_interm(dat, "intermediate_analysis_dataset_step2.csv")


################################ META-REGRESSION ################################ 

# moderators are at experiment-level, but analysis is at outcome level with 
#  nesting to handle experiment-level and paper-level correlation structure

# We will report the above metrics for each pair. Additionally, to summarize the above three metrics across pairs while accounting for their possible non-independence, we will robustly meta-regress each metric in a manner that accounts for clustering within original studies (Hedges et al., 2010). This model provides asymptotically correct inference even when the clustering structure is misspecified, which is important here because of the difficulty of precisely specifying the complex nature of clustering. This will yield average values of Porig, the difference, and the fixed-effects pooled estimate across all pairs. 

# Moderators:
# Animal vs. non-animal
# Type of replication lab (contract research organization [CRO] vs. academic core lab)
# What was requested from original authors and the response? (scored subjectively; Likert scale)
# Was a post hoc modification to protocol needed to complete the experiment? (col W in experiment level)

#@@Might be problematic (think about):
# N of original
# ES of original



##### Moderator Summary Table and Correlation Matrix #####

# must use dummies for labType here to get correlations
modVars = c("expAnimal",
            "labType_b.Both",
            "labType_c.CRO.only",
            "reqAntibodies",
            "reqCells",
            "reqPlasmids",
            "responseQuality",
            "changesNeeded")

CreateTableOne(vars = modVars, data = dat)


library(corrr)
corrs = dat %>% select(modVars) %>%
  correlate( use = "pairwise.complete.obs" ) %>%
  stretch() %>%
  arrange(desc(r)) %>%
  group_by(r) %>%
  filter(row_number()==1)


corrs$r = round(corrs$r, 2)

corrs = corrs[ !is.na(corrs$r), ]
View(corrs)

setwd(results.dir)
setwd("Tables to prettify")
write.csv(corrs, "moderator_cormat.csv")


##### Analyze Pairwise Metrics That *Do* Have Variances #####

# not surprisingly given its extreme values and variances, ratio doesn't really work here (V not positive definite)
outcomesWithVar = c( #"pw.ratio", 
  "pw.FEest")

# clear the results table to be created
if ( exists("modTable") ) rm(modTable)

for ( i in outcomesWithVar ) {

  modTable = safe_analyze_moderators(  .dat = dat,
                            yi.name = i,
                            # below assumes a standardized naming convention for
                            #  variances of the pairwise metrics:
                            vi.name = paste(i, "Var", sep=""),
                            
                            # cut out the "pw." part of outcome name
                            analysis.label = strsplit(i, "[.]")[[1]][2],
                            
                            modVars = modVars,
                            
                            n.tests = length(modVars),
                            digits = 2 )
  
}


modTable
table(modTable$Problems)


##### Analyze Pairwise Metrics That *Don't* Have Variances #####

outcomesWithoutVar = c("pw.ratio",
                       "pw.PIRepInside",
                       "pw.PIRepInside.sens",
                       "pw.Porig",
                       "pw.PorigSens",
                       "pw.PsigAgree1")

#if ( exists("modTable") ) rm(modTable)

for ( i in outcomesWithoutVar ) {
  modTable = safe_analyze_moderators(  .dat = dat,
                                       yi.name = i,
                                       # below assumes a standardized naming convention for
                                       vi.name = NA,
                                       
                                       # cut out the "pw." part of outcome name
                                       analysis.label = strsplit(i, "[.]")[[1]][2],
                                       
                                       modVars = modVars,
                                       
                                       n.tests = length(modVars),
                                       digits = 2 )
}


# for some reason, the FE problems get deleted
View(modTable)
table(modTable$Problems)

# @@think about ratio problem: maybe instead use difference, but controlling for original ES?


################################ TABLE: EXPT-LEVEL SUMMARIES OF PAIRWISE METRICS ################################ 


# TABLE of these metrics at the experiment level (~50 rows)

#bm
# quick look at results
# stringsWith( pattern = "pw", x = names(dat) )
# 
takeMean = c("pw.PIRepInside",  # function: 
             "pw.PIRepInside.sens",
             "pw.Porig",
             "pw.PorigSens",
             "pw.ratio",
             "pw.PsigAgree1",
             "pw.FEest")
# this is broken:
# res = dat %>% select(takeMean) %>%
#   mutate( across( .cols = everything(),
#                   .fns = mean) ) 

# x: vector of 0/1s
# @@note that this removes NAs

n_perc_string = function(x, digits = 0) {
  if ( all( is.na(x) ) ) return("All missing")
  x = x[!is.na(x)]
  paste( sum(x), " (", round( 100 * mean(x), digits ), "%)", sep = "" )
}
n_perc_string( c(0,0,0,0,1,1,1,0,0) )

harmonic_p = function(x) {
  library(harmonicmeanp)
  if ( all( is.na(x) ) ) return(NA)
  # @@note: better use p.hmp here becuase it's asymptotically exact,
  #  but it was giving error messages
  hmp.stat( x[ !is.na(x) ] ) 
}

#bm: instead of doing strings, should use numbers throughout and then post-process
#  that way we can use these in the plots
RE_string = function(yi, vi, digits = 2) {
  mod = rma.uni(yi = yi,
                vi = vi,
                method = "REML")
  
  paste( round( mod$b, digits ), 
}

# aggregation fns: 
# plain mean (ratio)
# count and percent (repinside, PsigAGree)
# FE analysis (origES2, repES2, FEest)
# harmonic mean p-value (porig)

expTable = dat %>% group_by(peID) %>%
  summarise( nOutcomes = n(),
             
             FEEst
             Ratio = round( mean(pw.ratio), 2 ),
             
             PIRepInside = n_perc_string(pw.PIRepInside),
             PIRepInside.sens = n_perc_string(pw.PIRepInside.sens), 
             # @@note: better use p.hmp here becuase it's asymptotically exact,
             #  but it was giving error messages
             Porig = harmonic_p( pw.Porig ),
             Porig.sens = harmonic_p( pw.PorigSens ),
             
             # overall proportion (within this experiment) expected to agree
             # @@insert actual sig agreement
             SigAgree = n_perc_string( repSignif == origSignif & repDirection == origDirection),
             PercSigAgree1 = paste( round( 100 * mean(pw.PsigAgree1), 0 ), "%", sep ="" )
  ) 

View(expTable)



#################################### DUMBBELL PLOT ###################################


# great info on dumbbell plot in ggplot:
# https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a

dat = read_interm("intermediate_analysis_dataset_step2.csv")


library(ggalt)
library(tidyverse)

dp = droplevels( dat %>% dplyr::filter( !is.na(origES2) & !is.na(repES2) & ES2type != "" ) %>%
                   dplyr::filter(ES2type %in% c("Cohen's d", "Glass' delta", "Log hazard ratio", "Cohen's dz") ) )
# randomly sample for testing purposes
set.seed(2)
dp = dp %>% group_by(ES2type) %>% sample_n( 3, replace = TRUE )
#dp = dat[1:10,]
dp$plotID = dp$peoID # with eye toward functionizing


repColor <- "#0171CE"
origColor <- "#DE4433"
digits = 2

stringPanelWidth = 3
panelSpacer = 1

max( c(dp$origES2, dp$repES2), na.rm = TRUE )  # use this to inform stringPanel1Start
stringPanel1Start = 15  # an x-axis coordinate
stringPanel1End = stringPanel1Start + stringPanelWidth  # x-axis 

stringPanel2Start = stringPanel1End + panelSpacer
stringPanel2End = stringPanel2Start + stringPanelWidth


# for labeling joy, have a df with just the first row to appear in plot
# which is actually the LAST row in the LAST ES2type
temp = dp[ dp$ES2type == unique(dp$ES2type)[1], ]  
# sort in order to match ggplot's ordering
dpRow = temp[ temp$peoID == rev( sort(temp$peoID) )[1], ] 
#dpRow = dp[ dp$plotID == dp$plotID[ nrow(dp) ], ]
#dpRow = dp[ dp$plotID == dp$plotID[ 1 ], ]


p = ggplot() + 
  
  facet_grid(ES2type ~ .,
             scales = "free",
             space = "free"
             #space = "free_y"
  ) +
  
  # null
  geom_vline(xintercept = 0,
             lty = 2,
             color = "gray") +
  
  geom_segment(data = dp,
               aes(y=plotID,
                   yend=plotID,
                   x=0,
                   xend=.5),
               color="#b2b2b2",
               size=0.15) +
  
  geom_dumbbell(data=dp,
                aes(y=plotID,
                    x=origES2,
                    xend=repES2),
                size=1.5,
                color="#b2b2b2",
                size_x=3,
                size_xend = 3,
                colour_x = origColor,
                colour_xend = repColor) +
  
  # "Original" and "Replication" labels
  # data arg: only label one of the rows
  geom_text( data= dpRow,
             aes(x=repES2,
                 y=plotID,
                 label="Replication"),
             color=repColor,
             size=3,
             vjust=-1.5,
             fontface="bold") +
  
  geom_text( data= dpRow,
             aes(x=origES2,
                 y=plotID,
                 label="Original"),
             color=origColor,
             size=3,
             vjust=-1.5,
             fontface="bold") +
  
  # numerical labels
  geom_text(data=dp,
            aes(x=repES2,
                y=plotID,
                label= round(repES2, digits) ),
            color=repColor,
            size=2.75,
            vjust=2.5) +
  
  geom_text(data=dp,
            aes(x=origES2,
                y=plotID,
                label= round(origES2, digits) ),
            color=origColor,
            size=2.75,
            vjust=2.5) +
  
  ##### differences panel
  geom_rect(data=dp,
            aes(xmin=stringPanel1Start,  # hard-coded location
                xmax=stringPanel1End,
                ymin=-Inf,
                ymax=Inf),
            fill="grey") +
  
  geom_text(data=dp,
            aes(label = round( origES2 - repES2, digits ),
                y=plotID,
                x= mean( c(stringPanel1Start, stringPanel1End) ) ),
            #fontface="bold",
            size=3) +
  # header
  geom_text(data=dpRow, 
            aes(x=mean( c(stringPanel1Start, stringPanel1End) ),  # needs to match above
                y=plotID,
                label="Original - replication"),
            color="black",
            size=3.1,
            vjust=-2,
            fontface="bold") 



# bm
# Porig panel
p + geom_rect(data=dp,
              aes(xmin=stringPanel2Start,  # hard-coded location
                  xmax=stringPanel2End,
                  ymin=-Inf,
                  ymax=Inf),
              fill="grey") +
  
  geom_text(data=dp,
            aes(label = round( pw.Porig, digits ),
                y=plotID,
                x=mean( c(stringPanel2Start, stringPanel2End) ) ),
            #fontface="bold",
            size=3) +
  
  geom_text(data=dpRow, 
            aes(x=mean( c(stringPanel2Start, stringPanel2End) ),  # needs to match above
                y=plotID,
                #label = TeX("$P_{orig}$")
                label="Porig"
            ),
            #label = "asdfd",
            #label = TeX("$P_{orig}$"),
            #label = expression(paste("DOC (mg ", L^-1,")")),
            color="black",
            size=3.1,
            vjust=-2,
            fontface="bold") + 
  
  # basic prettifying
  theme_bw() +
  theme( panel.grid.major=element_blank(),
         panel.grid.minor=element_blank() ) +
  
  xlab("Point estimate") +
  ylab("Paper, experiment, and outcome")







# # 15 x 9 works well
# # save
# setwd(objects.dir)
# ggsave( "forest.pdf",
#         width = 15,
#         height = 10,
#         units = "in" )
# setwd(results.overleaf)
# ggsave( "forest.pdf",
#         width = 15,
#         height = 9,
#         units = "in" )



#################################### RPP-STYLE SCATTERPLOT ###################################

#BM

# great info on dumbbell plot in ggplot:
# https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a

dat = read_interm("intermediate_analysis_dataset_step2.csv")


library(ggalt)
library(tidyverse)

# group certain ES together in plot
dat$ESgroup = NA
dat$ESgroup[ dat$ES2type %in% c("Cohen's d", "Glass' delta", "Cohen's dz") ] = "SMD"
dat$ESgroup[ dat$ES2type == "Log hazard ratio" ] = "Log hazard ratio"
dat$ESgroup[ dat$ES2type == "Fisher's z" ] = "Fisher's z"
table(dat$ES2type, dat$ESgroup)

dp = droplevels( dat %>% dplyr::filter( !is.na(origES2) & !is.na(repES2) & !is.na(ESgroup) ) )
                   
# # randomly sample for testing purposes
# set.seed(2)
# dp = dp %>% group_by(ES2type) %>% sample_n( 3, replace = TRUE )
# #dp = dat[1:10,]
dp$plotID = dp$peoID # with eye toward functionizing



##### Lower Panel: RPP-style scatterplots ######
p = ggplot() + 
  
  facet_grid(. ~ ESgroup,
             scales = "free",
             space = "fixed"
             #space = "free_y"
  ) +

  # null
  geom_abline(intercept = 0,
             slope = 1,
             lty = 2,
             color = "gray") +

  
  geom_point( data = dp,
              aes(x = origES2,
                  y = repES2)) +
  
# basic prettifying
theme_bw() +
  theme( panel.grid.major=element_blank(),
         panel.grid.minor=element_blank() ) +
  
  xlab("Original study estimate") +
  ylab("Replication study estimate")



##### Upper Panel: Ordered differences ######

dp = dp %>% filter( !is.na(ESgroup) ) %>%
  arrange(desc(repES2 / origES2))

dp$ind = 1:nrow(dp)


# only use CI if variance of ratio is less than 10

p = ggplot() +
  
  # null
  geom_hline(yintercept = 100,
              lty = 2,
              color = "gray") +
  
  
  geom_point( data = dp,
              aes(x = ind,
                  y = 100*(repES2 / origES2),
                  color = ESgroup) ) +
  
  # geom_errorbar( data = dp,
  #                aes(ymin = pw.ratio - qnorm(.975) * pw.ratioVar,
  #                    ymax = pw.ratio + qnorm(.975) * pw.ratioVar,
  #                    color = ESgroup) ) +
  
  scale_y_log10() +
  
  # basic prettifying
  theme_bw() +
  theme( panel.grid.major=element_blank(),
         panel.grid.minor=element_blank() ) +
  
  ylab("Replication percent of original (logged axis)")



##### Upper Panel: Ordered ratios ######

dp = dp %>% filter( !is.na(ESgroup) ) %>%
  arrange(desc(repES2 / origES2))

dp$ind = 1:nrow(dp)


# only use CI if variance of ratio is less than 10

p = ggplot() +
  
  # null
  geom_hline(yintercept = 100,
             lty = 2,
             color = "gray") +
  
  
  geom_point( data = dp,
              aes(x = ind,
                  y = 100*(repES2 / origES2),
                  color = ESgroup) ) +
  
  # geom_errorbar( data = dp,
  #                aes(ymin = pw.ratio - qnorm(.975) * pw.ratioVar,
  #                    ymax = pw.ratio + qnorm(.975) * pw.ratioVar,
  #                    color = ESgroup) ) +
  
  scale_y_log10() +
  
  # basic prettifying
  theme_bw() +
  theme( panel.grid.major=element_blank(),
         panel.grid.minor=element_blank() ) +
  
  ylab("Replication percent of original (logged axis)")




################################ GROUP METRICS ################################ 

# Primary: Percentage of replication estimates >0 (but will be an overestimate due to overdispersion)
# @replace with calibrated version?

# FIGURE: Density plot of calibrated estimates (outcome level) for group of replications and group of originals, overlaid, and facetted by EStype






