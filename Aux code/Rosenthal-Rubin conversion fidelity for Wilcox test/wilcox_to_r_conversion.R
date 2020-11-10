
# See if the r = Z/sqrt(N) somehow works when r is Wilcoxon instead of Pearson.

###### Helper Fns #####

Z_to_r = function(Z, n) Z/sqrt(n)

library(MASS)
library(expect_equal)
library(clusterPower)
library(testthat)

seed = 451


##### Correct Setting: 2 Binary Variables and Z from Chi-Square Test #####

n = 1000

set.seed(seed)
x = rbinom( n = n,
            size = 1, 
            prob = 0.5 )

prob = rep(0.3, n)
prob[ x == 1 ] = .4

setseed(seed+1)
y = rbinom( n = n, 
            size = 1,
            prob = prob )

d = data.frame(x, y)

cor.test(d$x, d$y)

( rTrue = cor(d$x, d$y) )

# now convert from Z
Z = as.numeric( sqrt( chisq.test(x = d$x, y = d$y)$statistic ) )
( rConv = Z_to_r(Z = Z, n = n) )


# very close, as expected
rTrue; rConv
# relative bias 1.7%
abs(rConv-rTrue)/rTrue



##### Wrong Setting: 2 Continuous Variables and Wilcoxon's Z #####

# marginally std normal
# empirical correlation same as that seen above
# i.e., match Pearson's r to the previous point-biserial
set.seed(seed + 2)
d = mvrnorm( n = n,
             mu = c(0, 0),
             Sigma = matrix( c(1, rTrue, rTrue, 1), nrow = 2 ),
             empirical = TRUE )

d = as.data.frame(d)
names(d) = c("x", "y")

expect_equal( cor(d$x, d$y), rTrue )


### Method 1: get "Z-score" from Wilcoxon p-value and use same conversion
# get Z-score from Wilcoxon
# https://stats.stackexchange.com/questions/133077/effect-size-to-wilcoxon-signed-rank-test
( res = wilcox.test(x = d$x, y = d$y) )
p = res$p.value

( Z = qnorm( 1 - p/2) )
expect_equal( 2 * ( 1 - pnorm(Z) ), p )
rConv2 = Z_to_r(Z, n)

# not close
rTrue; rConv2
# relative bias 95%
abs(rConv2-rTrue)/rTrue


##### Method 2: get "t-score" from Wilcoxon p-value and use r_equivalent #####
# should be very equivalent to the above as n -> infty

# Rosenthal & Rubin, pg 494 does exactly this for Mann-Whitney U
# sanity check using their example on pg 494
t = abs( qt(p = .008, df = 7) )
expect_equal(round(t,2), 3.16)

pOne = p/2  # one-tailed p-value
t = abs( qt(p = pOne, df = n-2) )
expect_equal( 2 * ( 1 - pt(t, df = n-2) ), p )
( rConv3 = sqrt( t^2 / ( t^2 + (n-2) ) ) )

# not close
rTrue; rConv3
# relative bias 95%
abs(rConv3-rTrue)/rTrue




