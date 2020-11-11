

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
#              MAIN PAIRWISE METRICS (COMPLETED QUANT PAIRS)                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

################################ TO CALCULATE FOR EACH PAIR ################################ 

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
dat = dat %>% 
  rowwise() %>% 
  mutate( analyze_one_row(origES2,
                          origVar2, 
                          repES2,
                          repVar2,
                          ES2type) )

# quick look at results
# stringsWith( pattern = "pw", x = names(dat) )
# 
takeMean = c("pw.PIRepInside",
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

colMeans( dat %>% select(takeMean), na.rm = TRUE )


################################ SUMMARIES AFTER THE ABOVE ################################ 

# FOREST PLOT
#bm


# TABLE of these metrics at the experiment level (~50 rows)

# FIGURE: Similar to RPP with original estimates vs. replication estimates at all three levels of granularity (emphasizing experiment level, probably)


#################################### FOREST PLOT ###################################


# great info on dumbbell plot in ggplot:
# https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a

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
    
    
  
  
  
  
  
  
#   +
#   
#   # prettify
#   theme_bw() +
#   theme( panel.grid.major=element_blank(),
#          panel.grid.minor=element_blank() )
# 
# 
# # scale_x_continuous(expand=c(0,0), limits=c(0, .625)) +
# # scale_y_discrete(expand=c(0.2,0))
# 
# p



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



################################ GROUP METRICS ################################ 

# Primary: Percentage of replication estimates >0 (but will be an overestimate due to overdispersion)
# @replace with calibrated version?

# FIGURE: Density plot of calibrated estimates (outcome level) for group of replications and group of originals, overlaid



################################ META-REGRESSION ################################ 

# We will report the above metrics for each pair. Additionally, to summarize the above three metrics across pairs while accounting for their possible non-independence, we will robustly meta-regressed each metric in a manner that accounts for clustering within original studies (Hedges et al., 2010). This model provides asymptotically correct inference even when the clustering structure is misspecified, which is important here because of the difficulty of precisely specifying the complex nature of clustering. This will yield average values of Porig, the difference, and the fixed-effects pooled estimate across all pairs. 

# Moderators:
# Animal vs. non-animal
# Type of replication lab (contract research organization [CRO] vs. academic core lab)
# What was requested from original authors and the response? (scored subjectively; Likert scale)
# Was a post hoc modification to protocol needed to complete the experiment? (col W in experiment level)
# N of original
# ES of original


# table(d$`Original test statistic type`)








