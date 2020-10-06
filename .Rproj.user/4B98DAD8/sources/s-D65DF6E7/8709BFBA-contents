

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

root.dir = "~/Dropbox/Personal computer/Independent studies/2020/RPCB reproducibility cancer biology"
raw.data.dir = paste(root.dir, "Raw data", sep="/")
prepped.data.dir = paste(root.dir, "Prepped data", sep="/")
code.dir = paste(root.dir, "Code (git)", sep="/")
results.dir = paste(root.dir, "Results from R", sep="/")

setwd(code.dir)
source("helper.R")


setwd(prepped.data.dir)

d = read_xlsx("RP_CB Final Analysis .xlsx")

names(d)


dim(d)  # 258

table(d$`Replication attempted`)  # 233
table(d$`Experiment completed`) # 190

d %>% filter( `Experiment completed` == "Yes" ) %>%
  summarise( length(unique(`Original study title`)),
             mean(`Number of lab(s) contracted for the entire study`) )

data.frame( d %>% group_by(`Original study title`) %>%
              summarise( n(),
                         comp = mean(`Experiment completed` == "Yes"),
                         expN0 = max(`Experiment #`, na.rm = TRUE),
                         max(`Study #`, na.rm = TRUE) ) )


#################################### 

# We will conduct the main analyses at two levels of granularity (outcome-level and experiment-level). 


# distinctions in analysis:
#  - completed pairs (has realized replication outcome) vs. all pairs

# 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                     METRICS FOR ALL PAIRS (INCL NON-COMPLETED)                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Percent sign agreement: The percentage of replications whose estimates agree in direction with the original study. This could be heuristically compared to the 50% that would be expected by chance if the null holds exactly in every replication (i.e., no effect heterogeneity) and conditional on all originals’ being positive in sign.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                       MAIN PAIRWISE METRICS (COMPLETED QUANT PAIRS)                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

################################ TO CALCULATE FOR EACH PAIR ################################ 

# Primary: Prediction interval and whether replication falls inside it (assuming t2=0)
# Make forest plot of these, including Porig on the side

# Primary: Porig

# Primary: Ratio of original to replication study estimates

# Secondary: Replication p < .05 with 2 expectation benchmarks
#  - from JRSSA paper
#  - the true effect size in each original is equal to the effect size for which it would
#   have had 80% power.


# Porig with assumed zero heterogeneity

# Fixed-effects pooled estimate: A meta-analytic pooled estimate of the original and replication estimates from each pair. 

# As sensitivity analysis:
# Porig and pred interval with imputed heterogeneity: If there is moderate or high within-pair effect heterogeneity, this could make the original studies appear less consistent with the replications than they truly are. As a sensitivity analysis, we will impute the average heterogeneity estimate from Olsson-Collentine, et al. (in press), which was tau=0.13 on the SMD scale, and use this to re-calculate Porig for each pair. These values of Porig will likely be large (i.e., indicating better consistency) than those from main analyses.


################################ SUMMARIES AFTER THE ABOVE ################################ 


# TABLE of these metrics at the experiment level (~50 rows)

# FIGURE: Similar to RPP with original estimates vs. replication estimates at all three levels of granularity (emphasizing experiment level, probably)



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











