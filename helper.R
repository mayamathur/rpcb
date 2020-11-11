
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                   ANALYSIS HELPER                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# calculates pairwise metrics for one row of data (whether at outcome- or experiment-level)
analyze_one_row = function(origES2,
                           origVar2, 
                           repES2,
                           repVar2,
                           ES2type){
  
  # test only
  #x = dat[1,]
  
  if( is.na(origES2) | is.na(origVar2) | is.na(repES2) | is.na(repVar2) ) {
    return( data.frame( pw.PILo = NA,
                        pw.PIHi = NA,
                        pw.PIRepInside = NA,
                        
                        pw.PILo.sens = NA,
                        pw.PIHi.sens = NA,
                        pw.PIRepInside.sens = NA,
                        
                        pw.Porig = NA,
                        pw.PorigSens = NA,
                        
                        pw.ratio = NA,
                        
                        pw.PsigAgree1 = NA,
                        
                        pw.FEest = NA,
                        pw.FEvar = NA,
                        pw.FElo = NA,
                        pw.FEhi = NA
    ) )
  }
  
  # prediction interval with t2 = 0
  predInt = pred_int( yio = origES2,
                      vio = origVar2,
                      yir = repES2,
                      vir = repVar2 )
  
  # Porig with t2 = 0
  Porig = p_orig( yio = origES2,
                  vio = origVar2,
                  yr = repES2,
                  vyr = repVar2,
                  t2 = 0 )
  
  # for sensitivity analyses:s
  # prediction interval and Porig with tau = 0.13
  # from Olssen-Collentine
  # only use truly comparable SMDs for this
  if ( !is.na(ES2type) & ES2type %in% c("Cohen's d", "Glass' delta") ){
    #if ( ES2type %in% c("Cohen's d", "Glass' delta") ) {
    Porig.sens = p_orig( yio = origES2,
                         vio = origVar2,
                         yr = repES2,
                         vyr = repVar2,
                         t2 = sqrt(.13) )
    
    # prediction interval
    # as in Replicate::pred_int, but adding in t2 here:
    yio = origES2
    vio = origVar2
    yir = repES2
    vir = repVar2
    pooled.SE = sqrt(vio + vir + .13)
    PILo.sens = yio - qnorm(1 - .05/2) * pooled.SE
    PIHi.sens = yio + qnorm(1 - .05/2) * pooled.SE
    PIinside.sens = (yir > PILo.sens) & (yir < PIHi.sens)
    
    #}
  } else {
    Porig.sens = PILo.sens = PIHi.sens = PIinside.sens = NA
  }
  
  
  # P(replication p < 0.05) - Mathur & VanderWeele
  PsigAgree1 = prob_signif_agree(yio = origES2,
                                 vio = origVar2,
                                 vir = repVar2,
                                 t2 = 0,
                                 null = 0)
  
  
  # @@CALCULATING POWER IS GOING TO BE A PROBLEM UNLESS WE CAN DO IT FOR ONLY 2-GROUP DESIGNS
  #  AND MAYBE CORRELATIONS
  # # P(replication p < 0.05) in counterfactual world in which the true 
  # #  effect size in replication and original is equal to the effect size for which it would
  # #  have had 80% power
  # # i.e., this is just the power to detect that effect size given replication's variance
  # power.t.test( )
  # 
  # # the df are only approximate
  # # @@ could use 
  # if ( !is.na(origN) ) {
  #   df = origN - 2
  #   tcrit = qt(0.975, df = df)
  # } else {
  #   stop("Found a case with origN NA, so should use 1.96 here, I guess (also need to impute df)")
  # }
  # 
  # # mu: the true effect
  # # we will solve for it
  # pwr_at_mu = function(mu){
  #   pt( q = tcrit - mu / sqrt(origVar2),
  #       df = df,
  #       lower.tail = FALSE ) +
  #     pt( q = -tcrit - mu / sqrt(origVar2),
  #         df = df,
  #         lower.tail = TRUE )
  # }
  # 
  # # sanity check
  # expect_equal( 0.05, pwr_at_mu(0) )
  # 
  # 
  # pwr_at_mu(15)
  # power.t.test( n = origN,
  #               delta = 15,
  #               sd = sqrt(origVar2) * sqrt(origN) )$power
  # 
  # # solve for counterfactual mu
  # mu.cfactual = uniroot( function(.mu) eval( pwr_cfactual(.mu) ) = 0.80 )
  
  # FE meta-analysis of original and replications
  #bm
  FEmod = rma.uni( yi = c(origES2, repES2),
                   vi = c(origVar2, repVar2),
                   method = "FE")
  
  
  #return as dataframe for mutate joys
  # "pw" prefix for "pairwise" metrics
  return( data.frame( pw.PILo = predInt$int.lo,
                      pw.PIHi = predInt$int.hi,
                      pw.PIRepInside = predInt$rep.inside,
                      
                      pw.PILo.sens = PILo.sens,
                      pw.PIHi.sens = PIHi.sens,
                      pw.PIRepInside.sens = PIinside.sens,
                      
                      pw.Porig = Porig, 
                      pw.PorigSens = Porig.sens,
                      
                      pw.ratio = origES2 / repES2,
                      
                      pw.PsigAgree1 = PsigAgree1,
                      
                      pw.FEest = as.numeric( FEmod$b ),
                      pw.FEvar = as.numeric( FEmod$se^2 ),
                      pw.FElo = as.numeric( FEmod$ci.lb ), 
                      pw.FEhi = as.numeric( FEmod$ci.ub )
  ) )
}

# x = dat[58,]
# analyze_one_row(x$origES2,
#                 x$origVar2,
#                 x$repES2,
#                 x$repVar2,
#                 x$ES2type)


# analyze a subset or a moderator
# n.tests: for Bonferroni
analyze_moderators = function( .dat,
                               yi.name,  # a pairwise metrics
                               vi.name,  # variances of pairwise metrics (use a single NA if doesn't have one)
                               
                               analysis.label,
                               mod.names,
                               mod.continuous = FALSE,
                               #ql,  # log scale
                               #take.exp,
                               #boot.reps = 2000,
                               n.tests=1,
                               digits = 2
) {
  
  # maybe use MRM metrics???
  
  ##### test only
  analysis.name = "FE"
  .dat = dat
  yi.name = "pw.FEest"
  vi.name = "pw.FEvar"
  #vi = rep(1, nrow(dat))  # variances of pairwise metrics (use 1 for all if no variances)
  analysis.label = "Porig",
  mod.names = c("responseQuality",
                "expType",
                "reqAntibodies",
                "reqCells",
                "reqPlasmids",
                "labType",
                "expAnimal")
  # mod.continuous = FALSE,
  
  n.tests=length(moderators)
  digits = 2
  ##### end test
  
  
  .dat$Yi = .dat[[yi.name]]
  .dat$Vi = .dat[[vi.name]]
  
  # complete cases to avoid "t(x) %*% w : non-conformable arguments" in clubSandwich
  .dat = .dat[ , c("pID", "eID", "Yi", "Vi", mod.names ) ]
  .dat = .dat[ complete.cases(.dat), ]
  
  formString = paste( "Yi ~ ", paste( mod.names, collapse= " + ") )
  
  library(clubSandwich)
  
  
  ##### Set Up Working Correlation Matrix #####
  if ( !is.na(vi.name) ) {
    # from Pustejovsky code "Analyze Tanner-Smith & Lipsey 2015 data.R"
    # CHE: multilevel random effects model with constant sampling correlation working
    # constant sampling correlation working model
    V_mat <- impute_covariance_matrix(vi = .dat$Vi,
                                      cluster = .dat$pID,  #@@ may need to change this?
                                      r = 0.6)
  }
  
  
  if ( is.na(vi.name) ) {
    # no variances, so just set to 1s
    V_mat <- impute_covariance_matrix(vi = rep(1,  nrow(.dat) ),
                                      cluster = .dat$pID,  #@@ may need to change this?
                                      r = 0.6)
  }
  
  
  ##### Fit Random-Effects Model #####
  # fit random effects working model in metafor
  # three-level model: https://stats.stackexchange.com/questions/116659/mixed-effects-meta-regression-with-nested-random-effects-in-metafor-vs-mixed-mod
  # @@the random structure will be only two-level when doing exp-level 
  model <- rma.mv( eval( parse(text=formString) ),
                   V = V_mat,
                   random = ~ 1 | pID / eID,
                   data = .dat,
                   sparse = TRUE)
  
  # above reports model-based (not robust) standard errors
  
  # RVE standard errors
  res <- conf_int(model, vcov = "CR2")
  
  # @@decide on t vs. z (keep consistent with conf_int above)
  #pvals = 2 * pt( abs(res$beta) / res$SE, df = res$df, lower.tail = FALSE )
  pvals = 2 * pnorm( abs(res$beta) / res$SE, lower.tail = FALSE )
  
  # compare to model-based p-values
  #cbind(pvals, model$pval)
  
  ##### Put Results in Dataframe #####
  est.string = paste( round( res$beta, digits ),
                      format_CI( res$CI_L, 
                                 res$CI_U,
                                 digits),
                      sep = " " )
  
  tau.string = round( sqrt(model$tau2), digits)
  
  
  new.row = data.frame( Analysis = analysis.name,
                        Coefficient = row.names(res),
                        n = nrow(.dat),  # complete-cases
                        Est = est.string,
                        Pval = format_stat(pvals, cutoffs = c(.1, .0001) ),
                        Pval.Bonf = format_stat( pmin(pvals*n.tests, 1) ),
                        Tau = tau.string )
  
  
  
  # this should be a global variable
  if ( !exists("resE") ){
    modTable <<- new.row
  } else {
    library(plyr)
    modTable <<- rbind.fill(modTable, new.row)
    detach("package:plyr", unload=TRUE)
  }
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                   DATA-PREP HELPER                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

################################ EFFECT-SIZE CONVERSIONS ################################

# convert to scale that's appropriate for meta-analysis (e.g., log-HR instead of HR),
#  but NOT to a mutually comparable scale across pairs
convert_to_ES2 = function(x,
                          .EStype){  
  
  # # test only
  # x = d2$origES
  # .EStype = d2$EStype
  
  x2 = type = rep(NA, length(x))
  
  # these stay the same
  ind = !is.na(.EStype) & .EStype %in% c("Cliff's delta", "Cohen's d", "Cohen's dz", "Glass' delta")
  x2[ind] = x[ind]
  type[ind] = .EStype[ind]
  
  # @need to handle Cohen's w
  
  # HR -> log-HR
  ind = !is.na(.EStype) & .EStype == "Hazard ratio"
  x2[ind] = log(x[ind])
  type[ind] = "Log hazard ratio"
  
  # @ASSUMES R IS PEARSON CORRELATION; NEED TO CHECK IF TRUE:
  ind = !is.na(.EStype) & .EStype == "r" 
  x2[ind] = r_to_z_NA(x[ind] )
  type[ind] = "Fisher's z"
  
  return( data.frame(ES2 = x2,
                     ES2type = type) )
}

# convert_to_ES2( x = c(-.2, 0.32, -.67),
#                 .EStype = c("Cliff's delta", "Hazard ratio", "r") )



# convert to SMD
convert_to_ES3 = function(x,
                          .ES2type){  
  
  # test only
  x = d2$ES2
  .ES2type = d2$ES2type
  
  x2 = type = rep(NA, length(x))
  
  # these are already SMDs
  # @though note that we're allowing Cohen's dz here
  ind = !is.na(.EStype) & .EStype %in% c("Cohen's d", "Cohen's dz", "Glass' delta")
  x2[ind] = x[ind]
  type[ind] = "SMD"
  
  # convert log-HRs
  ind = !is.na(.EStype) & .EStype == "Log hazard ratio"
  x2[ind] = 
    type[ind] = "Log hazard ratio"
  
  # @ASSUMES R IS PEARSON CORRELATION; NEED TO CHECK IF TRUE:
  ind = !is.na(.EStype) & .EStype == "r" 
  x2[ind] = r_to_z_NA(x[ind] )
  type[ind] = "Fisher's z"
  
  return( data.frame(ES2 = x2,
                     ES2type = type) )
}


# UPDATE THIS IN METAUTILITY
# Pearson's r to Fisher's z
# NOT the Z-score
# handles NAs
r_to_z_NA = function(r){
  z = rep(NA, length(r))
  z[ !is.na(r) ] = r_to_z( r[ !is.na(r) ] )
  return(z)
}
# r_to_z_NA( c(.22, -.9, NA) )



# ALSO PUT IN METAUTILITY
# logHR to logOR by way of logRR
# if rare, it's easy because HR \approx RR
# if common, use TVW's two conversions
# hi: upper CI limit
# goal of the fn is to fill in logOR and varLogRR for all entries
#  the other cols that get added in the interim are just intermediate steps
logHR_to_logOR = function(logHR,
                          rareY = rep(FALSE, length(logOR)),
                          lo = NA,
                          hi = NA){
  
  # # test only
  # logHR = c(NA, NA, log(1.03), log(2), log(2))
  # rareY = c(NA, NA, FALSE, TRUE, FALSE)
  # lo = c(NA, NA, log(.8), NA, log(1.5))
  # hi = c(NA, NA, NA, log(2.5), NA)
  
  d = data.frame( logHR,
                  lo,
                  hi,
                  rareY)
  
  
  ##### Rare Outcome #####
  # if outcome is rare, no transformation needed
  d[ eqNA(d$rareY == TRUE), c("logOR", "loLogOR", "hiLogOR") ] = d[ eqNA(d$rareY == TRUE), c("logHR", "lo", "hi") ]
  
  
  ##### Common Outcome #####
  # if outcome is common, use TVW's two transformations
  # logHR ->(Biometrics Thm2 conversion) logOR ->(sqrt conversion) logRR
  logHR_to_logRR_common = Vectorize( function(logHR){
    
    logRR = rep(NA, length(logHR))
    
    logRR[ !is.na(logHR) ] = log( ( 1 - 0.5^sqrt( exp(logHR[ !is.na(logHR) ]) ) ) / ( 1 - 0.5^sqrt( 1 / exp(logHR[ !is.na(logHR) ]) ) ) )
    
    return(logRR)
  } )
  #logHR_to_logRR_common( c(log(1.4), log(.745), NA) )
  
  
  
  logRR_to_logOR_common = Vectorize( function(logRR){
    
    logOR = rep(NA, length(logRR))
    
    logOR[ !is.na(logRR) ] = log( sqrt( exp( logRR[ !is.na(logRR) ] ) ) )
    
    return(logOR)
  } )
  #logRR_to_logOR_common( c(log(1.4), log(.745), NA) )
  
  
  # first convert logHR -> logRR via 
  #  TVW's Biometrics conversion (Thm 2)
  d[ eqNA(d$rareY == FALSE), c("logRR", "loLogRR", "hiLogRR") ] = logHR_to_logRR_common( d[ eqNA(d$rareY == FALSE), c("logHR", "lo", "hi") ] )
  
  # now convert the RRs to ORs via square-root
  #     @bm: this is hitting an error
  d[ eqNA(d$rareY == FALSE), c("logOR", "loLogOR", "hiLogOR") ] = logRR_to_logOR_common( d[ eqNA(d$rareY == FALSE), c("logRR", "loLogRR", "hiLogRR") ] ) 
  
  ##### Get Variance from CI #####
  # use either lower or upper limit, depending on what's available
  d$lim = d$loLogOR
  d$lim[ is.na(d$loLogOR) ] = d$hiLogOR[ is.na(d$loLogOR) ]
  
  d$varLogOR[ !is.na(d$lim) ] = ci_to_var( est = d$logOR[ !is.na(d$lim) ],
                                           ci.lim = d$lim[ !is.na(d$lim) ] )
  
  return(d)
}

# # sanity checks:
# # test only
# logHR = c(NA, NA, log(1.03), log(2), log(2))
# rareY = c(NA, NA, FALSE, TRUE, FALSE)
# lo = c(NA, NA, log(.8), NA, log(1.5))
# hi = c(NA, NA, NA, log(2.5), NA)
# res = logHR_to_logOR( logHR, rareY, lo, hi )
# 
# # rare
# expect_equal( res$logHR[4], res$logOR[4] )
# expect_equal( ( (res$hi[4] - res$logHR[4]) / qnorm(.975) )^2, res$varLogOR[4] )
# 
# # common
# term = 0.5^sqrt( exp(res$logHR[3]) )
# term2 = 0.5^sqrt( 1/exp(res$logHR[3]) )
# expect_equal( res$logRR[3], log( (1 - term) / (1 - term2) ) )  # check RR
# expect_equal( res$logOR[3], log( sqrt( exp( res$logRR[3] ) ) ) )  # check OR
# 
# term = 0.5^sqrt( exp(res$lo[3]) )
# term2 = 0.5^sqrt( 1/exp(res$lo[3]) )
# expect_equal( res$loLogRR[3], log( (1 - term) / (1 - term2) ) )  # check RR limit
# expect_equal( res$loLogOR[3], log( sqrt( exp( res$loLogRR[3] ) ) ) )  # check OR limit
# expect_equal( ( (res$logOR[3] - res$loLogOR[3]) / qnorm(.975) )^2, res$varLogOR[3] )

# Hasselblad & Hedges conversion that we think only works for common outcomes
logOR_to_SMD = function(logOR,
                        lo = NA,
                        hi = NA){
  
  d = data.frame( logOR,
                  lo,
                  hi)
  
  
  # purpose of the internal fn here is to handle NAs in the above dataframe
  .logOR_to_SMD = Vectorize( function(logOR){
    
    SMD = rep(NA, length(logOR))
    
    SMD[ !is.na(logOR) ] = ( sqrt(3) / pi ) * logOR[ !is.na(logOR) ]
    
    return(SMD)
  } )
  #logHR_to_logRR_common( c(log(1.4), log(.745), NA) )
  
  
  d[ c("SMD", "loSMD", "hiSMD") ] = .logOR_to_SMD( d[ c("logOR", "lo", "hi") ] )
  
  
  ##### Get Variance from CI #####
  # use either lower or upper limit, depending on what's available
  d$lim = d$loSMD
  d$lim[ is.na(d$loSMD) ] = d$hiSMD[ is.na(d$loSMD) ]
  
  d$varSMD[ !is.na(d$lim) ] = ci_to_var( est = d$SMD[ !is.na(d$lim) ],
                                         ci.lim = d$lim[ !is.na(d$lim) ] )
  
  return(d)
}


# # sanity check
# logOR = c(log(1.03), log(.75), NA)
# lo = c(log(.8), NA, NA)
# hi = c(NA, log(.9), NA)
# res = logOR_to_SMD( logOR, lo, hi )
# 
# expect_equal( res$SMD[1], sqrt(3)/pi * res$logOR[1] )
# expect_equal( res$loSMD[1], sqrt(3)/pi * res$lo[1] )
# 
# # instead of what the fn is doing (i.e., convert CI limit itself), try getting var of
# #  OR first from its CI limit
# varLogOR = ( ( res$logOR[1] - res$lo[1] ) / qnorm(.975) )^2
# expect_equal( res$varSMD[1], sqrt(3)/pi * varLogOR )
# # @ WHY SO DIFFERENT?




# ALSO PUT IN METAUTILITY
# tests proposition x
# but returns FALSE instead of NA if x itself is NA
eqNA = function(x){
  !is.na(x) & x == 1
}
eqNA(NA == 5)


# ALSO PUT IN METAUTILITY
# calculates variance (Z- or t-based) given CI limit and point estimate
# ci.lim can be upper or lower
# if df not provided, assumes we're using Z
ci_to_var = Vectorize( function(est, ci.lim, df = NA){
  
  if ( is.na(est) | is.na(ci.lim) ) return(NA)
  # Z or t-stat
  stat = abs(est - ci.lim)
  
  if ( is.na(df) ) crit = qnorm(.975)
  if ( !is.na(df) ) crit = qt(p = 0.975, df = df)
  
  se = abs(est - ci.lim) / crit
  
  return(se^2)
} )
# # sanity check
# res = ci_to_var( est = c(1.05, NA),
#                  ci.lim = c(1.15, NA),
#                  df = c(10, NA) )
# expect_equal( 1.05 + qt(.975, df = 10) * sqrt(res[1]), 1.15 )


# test only
est = 1.05
ci.lim = 1.15
df = 10


# x2 should be 0/1
cohen_d = function( x2, y ) {
  n0 = length(x2[x2==0])
  n1 = length(x2[x2==1])
  m0 = mean( y[x2==0] )
  m1 = mean( y[x2==1] )
  sd0 = sd( y[x2==0] )
  sd1 = sd( y[x2==1] )
  num = (n0 - 1) * sd0^2 + (n1 - 1) * sd1^2
  denom = n0 + n1 - 2
  sig.pool = sqrt( num / denom )
  d = (m1 - m0) / sig.pool
  return(d)
}






################################ MISCELLANEOUS FORMATTING AND CONVENIENCE FNS ################################

# read/write intermediate work
write_interm = function(x, filename){
  setwd(prepped.data.dir)
  setwd("Intermediate work")
  write.csv(x, filename)
}

read_interm = function(filename){
  setwd(prepped.data.dir)
  setwd("Intermediate work")
  read.csv(filename)
}

# like View(), but opens the extra tab if global var useView = TRUE
View2 = function(x){
  if ( useView == TRUE ) View(x) 
}

# quick length(unique) equivalent
uni = function(x){
  length(unique(x))
}

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}

# return strings containing anything in pattern vector
stringsWith = function(pattern, x){
  # make regex expression 
  patterns = paste(pattern, collapse="|")
  x[ grepl(pattern = patterns, x = x)]
}
# stringsWith( pattern = c("dog", "cat"),
#  x = c("dogcat", "horse", "cat", "lion") )


# return indices of strings containing anything in pattern vector
whichStrings = function(pattern, x){
  patterns = paste(pattern, collapse="|")
  grepl(pattern = pattern, x = x)
}


# quickly search codebook
# returns rows of coddebook that have pattern anywhere in current variable name
#  or description
searchBook = function(pattern){
  rows = whichStrings( tolower(pattern), tolower(cd$`Current variable name`) ) |
    whichStrings( tolower(pattern), tolower(cd$`Description of variable`) )
  View(cd[ rows, ])
}
#searchBook("p value")


# fit robumeta with independent effects and return only pooled point estimate and variance
#  handles case where data have only 1 row
robu2 = function(yi, vi) {
  
  dat = data.frame( yi, vi )
  
  # handle case with 1 row in dat
  if ( nrow(dat) == 1 ) {
    return( data.frame( est = yi, 
                        var = vi ) )
  }
  
  if ( nrow(dat) > 1 ) {
    mod = robu( yi ~ 1,
                var.eff.size = vi,
                studynum = 1:length(yi),
                data = dat,
                small = TRUE)
    
    return( data.frame( est = as.numeric(mod$b.r), 
                        var = as.numeric(mod$reg_table$SE)^2 ) )
  }
}
# # test with one row
# robu2( d3[ d3$peID == "p15e1", "repES2" ], d3[ d3$peID == "p15e1", "repVar2" ] )
# # test with multiple rows
# robu2( d3[ d3$peID == "p16e1", "repES2" ], d3[ d3$peID == "p16e1", "repVar2" ] )



# recodes a Qualtrics checkbox question (i.e., a single column with comma-separated options)
#  into its constituent non-mutually-exclusive dummy variables
recode_checkboxes = function( .d, 
                              var ) {
  
  # NAs (character) represent that none of the categories apply
  # don't include these as their own dummy
  .d[[var]][ .d[[var]] == "NA" ] = NA
  
  # split race into dummies
  # https://stackoverflow.com/questions/27630588/split-string-column-to-create-new-binary-columns/27630769#27630769
  library(qdapTools)
  t = mtabulate( strsplit( .d[[var]], ",") )
  
  # remove the old variable and replace with the new one
  .d = .d %>% select(-var)
  .d = cbind(.d, t)
  
  return(.d)
}


# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
update_result_csv = function( name,
                              section = NA,
                              value = NA,
                              print = FALSE ) {
  setwd(results.dir)
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    .res = read.csv( "stats_for_paper.csv",
                     stringsAsFactors = FALSE,
                     colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% .res$name) ) .res[ .res$name %in% name, ] = new.rows
    else .res = rbind(.res, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    .res = new.rows
  }
  
  write.csv( .res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  # also write to Overleaf
  if (exists("overleaf.dir")){
    setwd(overleaf.dir)
    write.csv( .res, 
               "stats_for_paper.csv",
               row.names = FALSE,
               quote = FALSE )
  }
  
  if ( print == TRUE ) {
    View(.res)
  }
}