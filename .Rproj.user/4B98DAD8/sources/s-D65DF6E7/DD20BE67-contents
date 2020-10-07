
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



# handles NAs
r_to_z_NA = function(r){
  z = NA
  z[ !is.na(r) ] = r_to_z( r[ !is.na(r) ] )
  return(z)
}
# sr_to_z_NA( c(.22, -.9, NA) )

################################ MISCELLANEOUS ################################

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