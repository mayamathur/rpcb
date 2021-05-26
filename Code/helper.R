
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                   ANALYSIS HELPER                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# calculates pairwise metrics for one row of data (whether at outcome- or experiment-level)
# see analyze.R for sanity checks
analyze_one_row = function(origES3,
                           origVar3, 
                           repES3,
                           repVar3,
                           # heterogeneity value to impute in sens analyses
                           t2){
  
  
  # note: need to keep these variables matching with those in the eventual return
  if ( is.na(origES3) | is.na(origVar3) | is.na(repES3) | is.na(repVar3) ) {
    return( data.frame( pw.PILo = NA,
                        pw.PIHi = NA,
                        pw.PIRepInside = NA,
                        
                        pw.PILo.sens = NA,
                        pw.PIHi.sens = NA,
                        pw.PIRepInside.sens = NA,
                        
                        pw.Porig = NA,
                        pw.PorigSens = NA,
                        
                        pw.ratio = NA,
                        # variance must have this naming convention ("Var" suffix)
                        #  to be recognized by analyze_moderators later
                        pw.ratioVar = NA,
                        
                        pw.diff = NA,
                        pw.diffVar = NA,
                        
                        pw.PsigAgree1 = NA,
                        pw.PsigAgree1.sens = NA,
                        
                        pw.FEest = NA,
                        pw.FEestVar = NA,
                        pw.FElo = NA,
                        pw.FEhi = NA
    ) )
  }
  
  # prediction interval with t2 = 0
  predInt = pred_int( yio = origES3,
                      vio = origVar3,
                      yir = repES3,
                      vir = repVar3 )
  
  # Porig with t2 = 0
  Porig = p_orig( yio = origES3,
                  vio = origVar3,
                  yr = repES3,
                  vyr = repVar3,
                  t2 = 0 )
  
  # for sensitivity analyses:
  # prediction interval and Porig with tau from Olssen-Collentine
  # makes sense only for ES3, not ES2, so that we're using SMDs
  Porig.sens = p_orig( yio = origES3,
                       vio = origVar3,
                       yr = repES3,
                       vyr = repVar3,
                       t2 = t2 )
  
  # prediction interval
  # as in Replicate::pred_int, but adding in t2 here:
  yio = origES3
  vio = origVar3
  yir = repES3
  vir = repVar3
  pooled.SE = sqrt(vio + vir + t2)
  PILo.sens = yio - qnorm(0.975) * pooled.SE
  PIHi.sens = yio + qnorm(0.975) * pooled.SE
  PIinside.sens = (yir > PILo.sens) & (yir < PIHi.sens)
  
  
  # P(replication p < 0.05)
  PsigAgree1 = prob_signif_agree(yio = origES3,
                                 vio = origVar3,
                                 vir = repVar3,
                                 t2 = 0,
                                 null = 0)
  
  # P(replication p < 0.05) - with heterogeneity
  PsigAgree1.sens = prob_signif_agree(yio = origES3,
                                 vio = origVar3,
                                 vir = repVar3,
                                 t2 = t2,
                                 null = 0)
  
  # Courtney's counterfactual version
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
  #   pt( q = tcrit - mu / sqrt(origVar3),
  #       df = df,
  #       lower.tail = FALSE ) +
  #     pt( q = -tcrit - mu / sqrt(origVar3),
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
  #               sd = sqrt(origVar3) * sqrt(origN) )$power
  # 
  # # solve for counterfactual mu
  # mu.cfactual = uniroot( function(.mu) eval( pwr_cfactual(.mu) ) = 0.80 )
  
  # FE meta-analysis of original and replications
  FEmod = rma.uni( yi = c(origES3, repES3),
                   vi = c(origVar3, repVar3),
                   method = "FE")
  
  # pw.ratioVar = NA
  # # use delta method to approximate variance of ratio
  # tryCatch({
  library(msm)
  pw.ratioVar = deltamethod( ~ x1/x2,
                             mean = c( origES3, repES3 ), cov = diag( c(origVar3, repVar3) ) )^2
  # }, error = function(err){
  #   browser()
  # })
  
  
  #return as dataframe for mutate joy
  # "pw" prefix for "pairwise" metrics
  return( data.frame( pw.PILo = predInt$int.lo,
                      pw.PIHi = predInt$int.hi,
                      pw.PIRepInside = predInt$rep.inside,
                      
                      pw.PILo.sens = PILo.sens,
                      pw.PIHi.sens = PIHi.sens,
                      pw.PIRepInside.sens = PIinside.sens,
                      
                      pw.Porig = Porig, 
                      pw.PorigSens = Porig.sens,
                      
                      pw.ratio = origES3 / repES3,
                      pw.ratioVar = pw.ratioVar,
                      
                      pw.diff = origES3 - repES3,
                      pw.diffVar = origVar3 + repVar3,
                      
                      pw.PsigAgree1 = PsigAgree1,
                      pw.PsigAgree1.sens = PsigAgree1.sens,
                      
                      pw.FEest = as.numeric( FEmod$b ),
                      pw.FEestVar = as.numeric( FEmod$se^2 ),
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
                               vi.name = NA,  # variances of pairwise metrics (use a single NA if doesn't have one)
                               
                               analysis.label,
                               modVars,
                               
                               n.tests=1,
                               digits = 2
) {
  
  .dat$Yi = .dat[[yi.name]]
  
  if ( !is.na(vi.name) ) .dat$Vi = .dat[[vi.name]]
  if ( is.na(vi.name) ) .dat$Vi = 1
  
  # complete cases to avoid "t(x) %*% w : non-conformable arguments" in clubSandwich
  .dat = .dat[ , c("pID", "eID", "Yi", "Vi", modVars ) ]
  .dat = .dat[ complete.cases(.dat), ]
  
  formString = paste( "Yi ~ ", paste( modVars, collapse= " + ") )

  # if dealing with an outcome that has a variance
  if ( !is.na(vi.name) ) {
    ##### Set Up Working Correlation Matrix #####
    # from Pustejovsky code "Analyze Tanner-Smith & Lipsey 2015 data.R"
    # CHE: multilevel random effects model with constant sampling correlation working model
    # allows between- and within-study heterogeneity as well as correlated errors within clusters
    # latter could reflect, e.g., shared control groups
    # Eq. (3) in Pustejovsky & Tipton (2021)
    V_mat = impute_covariance_matrix(vi = .dat$Vi,
                                     cluster = .dat$pID, 
                                     r = 0.6)  # just a working "guess"
    
    ##### Fit Random-Effects Model #####
    # fit random-effects working model in metafor
    # three-level model: https://stats.stackexchange.com/questions/116659/mixed-effects-meta-regression-with-nested-random-effects-in-metafor-vs-mixed-mod
    # this provides the point estimates, but we will NOT use its SEs because they are parametric
    model = rma.mv( eval( parse(text=formString) ),
                    V = V_mat,
                    #V = Vi,  # if Vi = 1 throughout, makes it similar to lmer
                    random = ~ 1 | pID / eID,
                    data = .dat,
                    sparse = TRUE)
    
    t2 = sqrt(model$tau2)
  }
  
  # if dealing with outcome that doesn't have a variance
  if ( is.na(vi.name) ) {
    ##### Fit Random-Effects Model #####
    #   model-based should match exactly IF variances all equal
    formString2 = paste( formString, " + (1 | pID / eID)" )
    
    model = lmer( eval( parse(text=formString2) ),
                  data = .dat )
    
    t2 = NA
    
  }
  
  # get robust variances
  # use club sandwich estimators throughout instead of plain sandwich
  #  to handle few clusters and/or wrong working model
  # same regardless of rma.mv vs. lmer
  res = conf_int(model, vcov = "CR2")

  
  # compare to model-based p-values
  #cbind(pvals, model$pval)
  
  # p-values with Satterthewaite correction to match conf_int above
  pvals = coef_test(model, vcov = "CR2")$p_Satt
  
  
  ##### Put Results in Dataframe #####
  est.string = paste( round( res$beta, digits ),
                      format_CI( res$CI_L, 
                                 res$CI_U,
                                 digits),
                      sep = " " )
  
  tau.string = round( t2, digits)
  
  
  new.chunk = data.frame( Analysis = analysis.label,
                          Model = ifelse( is.na(vi.name), "lmer", "rma.mv" ),
                          Coefficient = row.names(res),
                          n = nrow(.dat),  # complete-cases
                          Est = est.string,
                          Pval = format_stat(pvals, cutoffs = c(.1, .0001) ),
                          Pval.Bonf = format_stat( pmin(pvals*n.tests, 1) ),
                          Tau = tau.string )
  
  
  # modTable should be a global variable
  if ( !exists("modTable") ){
    modTable = new.chunk
  } else {
    library(plyr)
    modTable = rbind.fill(modTable, new.chunk)
    detach("package:plyr", unload=TRUE)
  }
  
  return(modTable)
} 



# catches warnings and puts the warnings in modTable as a column
safe_analyze_moderators = function(...) {
  
  r = 
    tryCatch(
      withCallingHandlers(
        {
          error_text <- ""
          list(value = analyze_moderators(...), error_text = error_text)
        }, 
        warning = function(e) {
          error_text <<- trimws(paste0("WARNING: ", e))
          invokeRestart("muffleWarning")
        }
      ), 
      error = function(e) {
        return(list(value = NA, error_text = trimws(paste0("ERROR: ", e))))
      }, 
      finally = {
      }
    )
  
  modTable = r$value  # includes any past modTable
  # "rev" part is an ugly, hacky way to only overwrite "Problems" for the most recent analysis that we just ran
  modTable$Problems[ modTable$Analysis == rev(unique(modTable$Analysis))[1] ] = r$error_text
  return(modTable)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                   DATA-PREP HELPER                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


################################ MISC ################################

# ALSO PUT IN METAUTILITY
# tests proposition x
# but returns FALSE instead of NA if x itself is NA
eqNA = function(x){
  !is.na(x) & x == 1
}
#eqNA(NA == 5)

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
  ind = !is.na(.EStype) & .EStype %in% c("Cliff's delta",
                                         "Cohen's d",
                                         "Cohen's dz",
                                         "Glass' delta")
  x2[ind] = x[ind]
  type[ind] = .EStype[ind]
  
  # HR -> log-HR
  ind = !is.na(.EStype) & .EStype == "Hazard ratio"
  x2[ind] = log(x[ind])
  type[ind] = "Log hazard ratio"
  
  # Pearson's r -> Fisher's z
  ind = !is.na(.EStype) & .EStype == "Pearson's r" 
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
  
  # # # test only
  # x = d2$repES2
  # .ES2type = d2$ES2type
  
  # initialize new vector of SMDs
  x2 = rep(NA, length(x))
  
  # these are already SMDs, albeit with different interpretations
  ind = !is.na(.ES2type) & .ES2type %in% c("Cohen's d",
                                           "Cohen's dz",
                                           "Glass' delta")
  # no conversion
  x2[ind] = x[ind]
  
  # convert log-HRss
  ind = !is.na(.ES2type) & .ES2type == "Log hazard ratio"
  x2[ind] = logOR_to_SMD( logHR_to_logOR( x[ind] ) )$SMD
  
  # convert Fisher's z via Mathur & VanderWeele (2020)
  # since SDs are unknown, define contrast of interest as 1 SD change in X
  # @@mention this in Supplement
  ind = !is.na(.ES2type) & .ES2type == "Fisher's z" 

  x2[ind] = r_to_d( z_to_r( x[ind] ), 
                    sx = rep( 1, length( x[ind] ) ),  
                    delta = rep( 1, length( x[ind] ) ) )$d
  
  return( data.frame(ES3 = x2) )
}


# @@MM:update this fn in MetaUtility because this is more general
# Pearson's r to Fisher's z
# NOT the Z-score
# handles NAs
r_to_z_NA = function(r){
  z = rep(NA, length(r))
  z[ !is.na(r) ] = r_to_z( r[ !is.na(r) ] )
  return(z)
}
# r_to_z_NA( c(.22, -.9, NA) )



# @@MM: also put this in MetaUtility
# logHR to logOR by way of logRR
# if rare, it's easy because HR \approx RR
# if common, use TVW's two conversions
# hi: upper CI limit
# goal of the fn is to fill in logOR and varLogRR for all entries
#  the other cols that get added in the interim are just intermediate steps
logHR_to_logOR = function(logHR,
                          rareY = rep(FALSE, length(logHR)),
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
  
  logRR_to_logOR_common = Vectorize( function(logRR){
    
    logOR = rep(NA, length(logRR))
    
    logOR[ !is.na(logRR) ] = log( sqrt( exp( logRR[ !is.na(logRR) ] ) ) )
    
    return(logOR)
  } )

  # first convert logHR -> logRR via 
  #  TVW's Biometrics conversion (Thm 2)
  d[ eqNA(d$rareY == FALSE), c("logRR", "loLogRR", "hiLogRR") ] = logHR_to_logRR_common( d[ eqNA(d$rareY == FALSE), c("logHR", "lo", "hi") ] )
  
  # now convert the RRs to ORs via square-root
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
# # @ VERY DIFFERENT. THINK ABOUT. 



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


# # test only
# est = 1.05
# ci.lim = 1.15
# df = 10


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






################################ FNS FOR AGGREGATING PAIRWISE METRICS WITH STRING OUTPUT ################################

# aggregation fns: 
# plain mean (ratio)
# count and percent (repinside, PsigAgree)
# FE analysis (origES2, repES2, FEest)
# harmonic mean p-value (porig)
n_perc_string = function(x, digits = 0) {
  if ( all( is.na(x) ) ) return("All missing")
  x = x[!is.na(x)]
  paste( sum(x), " (", round( 100 * mean(x), digits ), "%)", sep = "" )
  
  
}
# n_perc_string( c(0,0,0,0,1,1,1,0,0) )


# construct string with sample mean and its CI
# by default, use cluster-robust inference
# also works for binary variables (linear prob model)
# return.n: should it return the ANALYZED sample size (i.e., non-missing on x or cluster var)?
mean_CI = function(x,
                   cluster = NA,
                   robust = TRUE,
                   return.n = FALSE,
                   return.df = FALSE) {
  
  if (all(is.na(x))) return("")
  
  if(robust == TRUE & length(x) != length(cluster)) stop("Cluster var is wrong length")
  
  
  if (robust == TRUE) {
    # remove missing data
    keep = !is.na(x) & !is.na(cluster)
    x = x[keep==TRUE]
    cluster = cluster[keep==TRUE]
    
    if(length(x) != length(cluster)) browser()
    
    # OLS (= linear prob model if binary)
    library(clubSandwich)
    mod = lm(x ~ 1) 
    Vmat = vcovCR(mod,
                  cluster = cluster,
                  type = "CR2")
    CIs = conf_int(mod, vcov = Vmat)
    
    tryCatch({
      
      if ( return.n == FALSE ) {
        
        if ( return.df == FALSE ) {
          return( paste( format_stat( mean(x) ),
                         format_CI( CIs$CI_L, CIs$CI_U ) ) )
        } else {
          return( data.frame( est = format_stat( mean(x) ),
                              lo = format_stat(CIs$CI_L),
                              hi = format_stat(CIs$CI_U) ) )
        }
        
        
      }
      
      if ( return.n == TRUE ) {
        
        if ( return.df == FALSE ) {
          return( paste( format_stat( mean(x) ),
                         format_CI( CIs$CI_L, CIs$CI_U ),
                         ", n =",
                         length(x),
                         sep = " " ) )
        } else {
          return( data.frame( est = format_stat( mean(x) ),
                              lo = format_stat(CIs$CI_L),
                              hi = format_stat(CIs$CI_U),
                              n = length(x) ) )
        }
        
        
      }
      
    }, error = function(err) {
      browser()
    })
    
    # return( paste( format_stat( mean(x) ),
    #                format_CI( CIs$CI_L, CIs$CI_U ) ) )
  }
  
  if (robust == FALSE) {
    x = x[!is.na(x)]
    xbar = mean(x)
    lims = as.numeric( t.test(x)$conf.int )
    
    if ( return.n == FALSE ) {
      
      if ( return.df == FALSE ) {
        return( paste( format_stat( mean(x) ),
                       format_CI( lims[1], lims[2] ) ) )
      }
      
      if (return.df == TRUE) {
        return( data.frame( est = format_stat( mean(x) ),
                            lo = format_stat(lims[1]),
                            hi = format_stat(lims[2]) ) )
      }
      
      
      
    } else {
      
      if ( return.df == FALSE ) {
        return( paste( format_stat( mean(x) ),
                       format_CI( lims[1], lims[2] ),
                       ", n =",
                       length(x),
                       sep = " " ) )
      }
      
      if ( return.df == TRUE ) {
        return( data.frame( est = format_stat( mean(x) ),
                            lo = format_stat(lims[1]),
                            hi = format_stat(lims[2]),
                            n = length(x) ) )
      }
      
    }
    
    
  }
}





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
  
  paste( round( mod$b, digits ), " [", round( mod$ci.lb, digits ), ", ", round( mod$ci.ub, digits ) )
}


################################ MISCELLANEOUS FORMATTING AND CONVENIENCE FNS ################################


# stands for "wipe results"
wr = function(){
  setwd(results.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  # setwd(overleaf.dir)
  # if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
}

# stands for "view results"
vr = function(){
  setwd(results.dir)
  View( read.csv("stats_for_paper.csv") )
}


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
  rows = whichStrings( tolower(pattern), tolower(cd$`Variable name`) ) |
    whichStrings( tolower(pattern), tolower(cd$`Description of variable`) )
  View(cd[ rows, ])
}
#searchBook("p value")


# fit robumeta with independent effects and return only pooled point estimate and variance
#  handles case where data have only 1 row
robu2 = function(yi, vi) {
  
  if ( any( is.na(yi) ) | any( is.na(vi) ) ) {
    return( data.frame( est = NA, 
                        var = NA ) )
  }
  
  dat = data.frame( yi, vi )
  
  # handle case with 1 row in dat
  if ( nrow(dat) == 1 ) {
    return( data.frame( est = yi, 
                        var = vi ) )
  }
  
  if ( nrow(dat) > 1 ) {
    #@make sure to cite:
    # #Tipton, E., & Pustejovsky, J. E. (2015). Small-sample adjustments for tests of moderators and model
    # fit using robust variance estimation in meta-regression. _Journal of Educational and Behavioral
    # Statistics, 40_(6), 604-634. doi: 10.3102/1076998615606099
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