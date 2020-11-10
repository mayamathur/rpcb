
# How bad is the Hasselblad & Hedges conversion from OR to SMD if outcome is common?
# Conclusion: 


setwd("~/Dropbox/Personal computer/Independent studies/2020/RPCB reproducibility cancer biology/Code (git)")

source("helper.R")

library(ggplot2)
library(dplyr)

# simulate binary X, normal Y
# dichotomize Y somewhere

# calculate OR with dichotomized Y and convert to d

# compare to true d with the continuous Y

##### Helper #####
# get true and converted (OR->d) SMD for a given dichotomization point in Y
get_SMDs = function(cutoffQuantile){
  cutoff = quantile( d$y, cutoffQuantile )
  
  d$y2 = d$y > cutoff
  
  ##### Calculate True d with Continuous Y #####
  ( dTrue = cohen_d(x2 = d$x, y = d$y) )
  
  
  ##### Calculate OR and Convert to d #####
  m = glm( y2 ~ x, family = binomial, data = d )
  OR = as.numeric( exp( coef(m)["x"] ) )
  
  dConv = logOR_to_SMD(logOR = log(OR))$SMD
  
  return( data.frame( dTrue, dConv) )
}
# get_SMDs(0.05)


##### Simulate Data (Binary X, Normal Y) #####
seed = c(2014)
n = 10000

#set.seed(seed)
x = rbinom(n = n,
           size = 1, 
           prob = 0.5)

#set.seed(seed+1)
y = rnorm(n = n,
          mean = 2 * x)  # coefficient is population ds
d = data.frame(x, y)


##### Compare True to Converted d #####


res = data.frame( q = seq(0.01, 0.99, 0.05))

res = res %>% rowwise() %>%
  mutate( get_SMDs(q) )



min(res$dConv) 
max(res$dConv)

# true d doesn't depend on cutoff, so there's only 1 unique value
.dTrue = unique(res$dTrue)

ggplot( data = res,
        aes(x = q,
            y = dConv) ) + 
  
  geom_hline( yintercept = .dTrue,
              lty = 2,
              color = "red") +
  geom_line() +
  scale_y_continuous( name = "Converted d (red line = true d)",
                      #limits = c(0, .4),
                      #breaks = seq(0, .4, .05),
                      sec.axis = sec_axis(trans = ~ . / .dTrue,
                                          name = "Converted d / True sample d" ) ) +
  scale_x_continuous( name = "Quantile of Y dichotomization",
                      limits = c(0, 1), 
                      breaks = seq(0, 1, .1)) +
  ggtitle( paste( "True sample d = ", round(.dTrue, 2) ) ) +
  theme_classic()
















