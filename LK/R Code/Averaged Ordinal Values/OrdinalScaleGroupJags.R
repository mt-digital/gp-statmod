
# OrdinalScaleGroupJags.R
# John Kruschke, April 2017 and February-March 2018.
#
# This is a lengthy script, but the vast majority of it is data manipulation (for 
# input) and graphics manipulation (for output). The Bayesian and frequentist 
# modeling constitutes only a small portion of the script.

#------------------------------------------------------------------------------
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

source("DBDA2E-utilities.R")
# Define new function for saving multiple formats with one command:
saveGraph2 = function( fileName , typeVec=c("png","pdf","eps") , ... ) {
  for ( typeVal in typeVec ) {
    saveGraph( file=fileName , type=typeVal , ... )
  }
}

#------------------------------------------------------------------------------
# Load/create the data file.
# The ordered-probit JAGS model expects data to be in long format, with one
# rating per row. There are columns for ordinal response, question ID, subject
# ID, group ID. The responses must be integers, from 1 to K. The ID's can be
# non-numeric, and will be converted to integer indices internally.
# N.B.: This version of the program assumes that all items have the same number 
# of ordinal levels.

# Function for generating simulated data:
genGroupData = function( 
  # specify means on underlying scales:
  mQ = c( 1.0 , 1.5 , 2.5 , 3.0 ) , # length is number of questions
  # specify standard deviations of responses on underlying scales:
  sdQ =  c( 1.5 , 2.0 , 2.5 , 3.0 ) , # length is number of questions
  # specify correlations of responses on underlying scales
  rQval = 0.9 ,
  # specify thresholds:
  thresh = matrix( c( 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 ,
                      1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 ,
                      1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 ,
                      1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 ),
                   nrow=length(mQ) , byrow=TRUE ) ,
  # specify number of subjects in the group:
  nS = 100 ,
  showSummary = FALSE
  ) { # generate data...
  require(MASS)
  nQ = length(mQ)
  rQ = matrix( rQval , nrow=nQ , ncol=nQ )
  diag(rQ) = 1.0
  # construct covariance matrix:
  sigmaCov = ( sdQ %o% sdQ ) * rQ
  # construct latent data:
  d = mvrnorm( n=nS , mu=mQ , Sigma=sigmaCov , empirical=TRUE )
  # compute response data:
  y = 0*d
  for ( sIdx in 1:nrow(d) ) {
    for ( qIdx in 1:ncol(d) ) {
      y[sIdx,qIdx] =  max( which( d[sIdx,qIdx] > c(-Inf,thresh[qIdx,]) ) )
    }
  }
  # show data summary on console:
  if ( showSummary ) {
    apply( d , 2 , mean )
    var( d )
    # show response data on console:
    apply( y , 2 , mean )
    var( y )
    cor( y )
  }
  return(y)
}

NperGroup = 500 # 500
rQvalGlobal = 0.8 # 0.8
for ( simIdx in c( "1QeqM" , "1QuneqM" , 
                   "2QeqM" , "2QuneqM" ,
                   "3QeqM" , "3QuneqM" ,
                   "6QeqM" , "6QuneqM" ,
                   "1Q1Gfa" , "4Q1Gfa"  ,"1Q1Gmiss" )[
                     c(1,2,
                       #3,4,
                       5,6,
                       7,8#,
                       #9,10,11
                       )] ) {
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="1Q1Gmiss" ) {
    # specify thresholds for all groups:
    nQ=1 # number of questions (a.k.a. items)
    thresh = matrix( rep( c( 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1sd = 2.75 # group 1 sd
    nullHypVal = 2.6 # null hypothesis value for computing t test and eff size
    # desEffSz = 0.1
    g1m = 2.0 # desEffSz * g1sd + nullHypVal # group 1 mean
    # g2m = NULL
    # g2sd = NULL
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1)) ) #+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } 
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="4Q1Gfa" ) {
    # specify thresholds for all groups:
    nQ=4 # number of questions (a.k.a. items)
    thresh = matrix( rep( c( 1.5 , 2.0 , 2.75 , 3.75 , 5.0 , 6.5 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 4.0 # group 1 mean
    nullHypVal = 4.0 # null hypothesis value for computing t test and eff size
    g1sd = 1.5 # group 1 sd
    # g2m = NULL
    # g2sd = NULL
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1)) ) #+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } 
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="1Q1Gfa" ) {
    # specify thresholds for all groups:
    nQ=1 # number of questions (a.k.a. items)
    thresh = matrix( rep( c( 1.5 , 2.0 , 2.75 , 3.75 , 5.0 , 6.5 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 4.0 # group 1 mean
    nullHypVal = 4.0 # null hypothesis value for computing t test and eff size
    g1sd = 1.5 # group 1 sd
    # g2m = NULL
    # g2sd = NULL
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1)) ) #+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } 
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="1QeqM" ) {
    # specify thresholds for all groups:
    nQ=1 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 4.2 # group 1 mean
    g1sd = 1.0 # group 1 sd
    g2m = 4.2
    g2sd = 4.0
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } 
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="2QeqM" ) {
    # specify thresholds for all groups:
    nQ=2 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 4.2 # group 1 mean
    g1sd = 1.0 # group 1 sd
    g2m = 4.2
    g2sd = 4.0
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data 1
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="3QeqM" ) {
    # specify thresholds for all groups:
    nQ=3 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 4.2 # group 1 mean
    g1sd = 1.0 # group 1 sd
    g2m = 4.2
    g2sd = 4.0
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data 1
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="6QeqM" ) {
    # specify thresholds for all groups:
    nQ=6 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 4.2 # group 1 mean
    g1sd = 1.0 # group 1 sd
    g2m = 4.2
    g2sd = 4.0
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data 1
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="1QuneqM" ) {
    # specify thresholds for all groups:
    nQ=1 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 6.12
    g1sd = 4.0
    g2m = 4.20 # group 2 mean
    g2sd = 1.0 # group 2 sd
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data 
  
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="2QuneqM" ) {
    # specify thresholds for all groups:
    nQ=2 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 6.12
    g1sd = 4.0
    g2m = 4.20 # group 2 mean
    g2sd = 1.0 # group 2 sd
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="3QuneqM" ) {
    # specify thresholds for all groups:
    nQ=3 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 6.12
    g1sd = 4.0
    g2m = 4.20 # group 2 mean
    g2sd = 1.0 # group 2 sd
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data
  
  #-----------------------------------------------------
  # Simulated data:
  if ( simIdx=="6QuneqM" ) {
    # specify thresholds for all groups:
    nQ=6 # number of questions (a.k.a. items)
    thresh = matrix( rep( seq( 1.5 , 4.5 , 1 ) , nQ ),
                     nrow=nQ , byrow=TRUE )  # thresh's of each Q
    g1m = 6.12
    g1sd = 4.0
    g2m = 4.20 # group 2 mean
    g2sd = 1.0 # group 2 sd
    fileNameRoot=paste0("OrdinalScaleGroupJags-Sim-",simIdx)
    likertScaleDensMax = 0.7 # for use in histogram later
    dataScaleDensMax = 0.7 # for use in histogram later
    # Now generate the data
    set.seed(47405)
    y1 = genGroupData(  
      mQ = rep(0.0,nQ)+g1m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g1sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    set.seed(47405)
    y2 = genGroupData(  
      mQ = rep(0.0,nQ)+g2m , # mean of each Q
      sdQ =  rep(1.0,nQ)*g2sd , # sd of each Q
      rQval = rQvalGlobal , # correlation of every pair of Q's
      thresh = thresh , # thresh's of each Q
      nS = NperGroup # number of subjects in this group
    )
    # Build a data frame in long format:
    resp = item = subID = cond = rep( 0 , prod(dim(y1))+prod(dim(y2)) )
    dfIdx = 1
    for ( sIdx in 1:nrow(y1) ) { for ( qIdx in 1:ncol(y1) ) {
      resp[dfIdx] = y1[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 1
      dfIdx = dfIdx+1
    } }
    for ( sIdx in 1:nrow(y2) ) { for ( qIdx in 1:ncol(y2) ) {
      resp[dfIdx] = y2[sIdx,qIdx]
      item[dfIdx] = qIdx
      subID[dfIdx] = sIdx
      cond[dfIdx] = 2
      dfIdx = dfIdx+1
    } }
    datFrm = data.frame( resp=resp , item=item , subID=subID , cond=cond )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
  } # end simulated data
  
  
  #------------------------------------------------------------------------------
  
  # Real data from Pinsof & Haselton (2016) supplied by Torrin:
  if ( !TRUE ) {
    datFrm = read.csv( file="formattedSSMD.csv" )
    yName="resp"
    qName="item"
    sName="subID"
    gName="cond"
    fileNameRoot="OrdinalScaleGroupJags-formattedSSMD"
    likertScaleDensMax = 0.7 # for use in histogram later
  } # end of real data
  
  #==============================================================================
  # ASSEMBLE THE DATA FOR JAGS. 
  #
  # N.B. THE JAGS MODEL ASSUMES THAT ALL ITEMS ARE SCORED WITH POSITIVE
  # INTERCORRELATIONS. 
  
  # Rename and reclass y values for convenience:
  y = as.numeric(datFrm[,yName])
  # Do some checking that data make sense:
  if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
  if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    cat("************************************************\n")
    cat("** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS **\n")
    cat("************************************************\n")
  }
  Ntotal = length(y)
  # Question, Subject, and Group data vectors:
  q = as.numeric(as.factor(datFrm[,qName]))
  s = as.numeric(as.factor(datFrm[,sName])) # not used
  g = as.numeric(as.factor(datFrm[,gName]))
  qLevels = levels(as.factor(datFrm[,qName]))
  sLevels = levels(as.factor(datFrm[,sName])) # not used
  gLevels = levels(as.factor(datFrm[,gName]))
  nQ = max(q)
  nS = max(s) # not used
  nG = max(g)
  # Create threshold matrix. For the first question, threshold 1 and nYlevels-1 are 
  # fixed; other interior thresholds are estimated. For other questions, all thresholds
  # are estimated.
  # ** THIS PRESENTLY ASSUMES THAT ALL ITEMS HAVE THE SAME NUMBER OF LEVELS. **
  # ** IT PRESUMABLY WON'T WORK OTHERWISE WITHOUT MODIFICATION. **
  # Compute number of Y levels for each question:
  nYlevels = aggregate( y , by=list(q) , FUN=max )$x
  thresh = matrix( NA , nrow=nQ , ncol=max(nYlevels)-1 ) # default to NA
  # Fix low thresh of item 1 at 1.5:
  thresh[1,1] = 1 + 0.5 
  # Fix upper thresh of item 1 at K-0.5:
  thresh[1,nYlevels[1]-1] = nYlevels[1] - 0.5
  # Specify the data in a list, for ORDERED-PROBIT model, for later shipment to JAGS:
  ordDataList = list(
    y = y ,
    q = q ,
    g = g ,
    thresh = thresh ,
    nYlevels = nYlevels ,
    nQ = nQ ,
    nG = nG ,
    Ntotal = Ntotal 
  )
  
  #-----------------------------------------------------------------------------
  # Treat data as metric, average Likert scale:
  meanY = aggregate( y , by=list(s=s,g=g) , FUN=mean ) # mean across q, each s
  # Specify the data in a list, for METRIC MODEL, for later shipment to JAGS:
  metDataList = list(
    meanY = meanY[,"x"] ,
    #q = q ,
    g = meanY[,"g"] ,
    #thresh = thresh ,
    nYlevels = nYlevels , # for calibrating the prior
    #nQ = nQ ,
    nG = nG ,
    Ntotal = nrow(meanY) 
  )
  
  #-----------------------------------------------------------------------------
  # FREQUENTIST ANALYSIS OF METRIC MODEL:
  if ( TRUE ) {
  # Open window for plot:
  openGraph(width=3.5*nG,height=2.5*1+0.75)
  layout(matrix(1:nG,nrow=1))
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) , oma=c(0,0,3.5,0) ) # , xpd=NA )
  histBreakInc = 1/nQ
  if ( nQ > 1 ) { xlabText="Mean of Ordinal Resp." } else { xlabText="Ordinal Resp." }
  # make histograms of data:
  for ( gIdx in 1:nG ) {
    ydata = meanY$x[meanY$g==gIdx]
    histInfo = hist( ydata , main=paste0("Group ",gIdx) , 
                     xlab=xlabText , ylab="Probability" ,
                     freq=FALSE , col="pink" , border="white" , cex.lab=1.5 ,
                     breaks=seq(-histBreakInc/2 , 
                                max(nYlevels)+1+histBreakInc/2 , 
                                by=histBreakInc ) , 
                     ylim=c(0,likertScaleDensMax) # NB: must set manually at this time
    )
    #text( max(histInfo$breaks) , likertScaleDensMax , adj=c(1,1) ,
    #      labels=bquote(N==.(length(ydata))) )
    xcomb = seq(0,max(nYlevels)+1,length=201)
    lines( xcomb , dnorm( xcomb-mean(ydata) , 0 , sd(ydata) ) , 
           col="skyblue" , lwd=2 )
  } # end for( gIdx in 1:nG )
  # Annotate with frequentist t-test info:
  if ( nG <= 2 ) {
    sdN = function(x){ sqrt(mean((x-mean(x))^2)) }
    if ( nG == 2 ) {
      tInfo = t.test( meanY$x[meanY$g==1] , meanY$x[meanY$g==2] )
      effSz = ( ( mean(meanY$x[meanY$g==1]) - mean(meanY$x[meanY$g==2]) ) 
                / sqrt(( sdN(meanY$x[meanY$g==1])^2 + sdN(meanY$x[meanY$g==2])^2 )/2) )
    } else {
      tInfo = t.test( meanY$x[meanY$g==1] , mu=nullHypVal )
      effSz = (mean(meanY$x[meanY$g==1])-nullHypVal) / sdN(meanY$x[meanY$g==1])
    }
    show( tInfo )
    show( effSz )
    if ( tInfo$p.value >= 0.001 ) {
      mtext( text=bquote( "d = "* .(round(effSz,2)) 
                               *", t = "* .(round(tInfo$statistic,2)) 
                               *", p = "* .(round(tInfo$p.value,3)) )  , 
             at=c(0.5) , cex=1.5 , side=3 , outer=TRUE , 
             adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
    } else {
      mtext( text=bquote( "d = "* .(round(effSz,2)) 
                               *", t = "* .(round(tInfo$statistic,2)) 
                               *", p < 0.001" ) , 
             at=c(0.5) , cex=1.5 , side=3 , outer=TRUE , 
             adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
    }
  }
  saveGraph2( paste0(fileNameRoot,"-LikertScale-Freq") , typeVec=c("eps","pdf") )
  } # end if frequentist analysis
  
  #-----------------------------------------------------------------------------
  # SPECIFY THE METRIC MODEL FOR AVERAGED MULTIPLE ITEMS FOR JAGS:
  metModelString = "
  model {
    for ( i in 1:Ntotal ) {
      meanY[i] ~ dnorm( mu[g[i]] , 1/sigma[g[i]]^2  )
    }
    for ( gIdx in 1:nG ) { 
      mu[gIdx] ~ dnorm( (1+max(nYlevels))/2 , 1/(max(nYlevels))^2 )
      sigma[gIdx] ~ dunif( 0.01 , max(nYlevels)*10 )
    }
  }
  " # close quote for metModelString
  # Write out metModelString to a text file
  writeLines( metModelString , con="OrdAsOrdAndMet-MetModel.txt" )
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS FOR METRIC MODEL
  parameters = c( "mu" , "sigma" )
  numSavedSteps = 20000 # 20000 
  thinSteps = 5 # 5 
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  saveName=fileNameRoot 
  runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
  nChains=nChainsDefault # from DBDA2E-utilities
  metRunJagsOut <- run.jags( method="parallel" , # runjagsMethod ,
                             model="OrdAsOrdAndMet-MetModel.txt" , 
                             monitor=parameters , 
                             data=metDataList ,  
                             #inits=initsList , 
                             n.chains=nChains ,
                             adapt=adaptSteps ,
                             burnin=burnInSteps , 
                             sample=ceiling(numSavedSteps/nChains) ,
                             thin=thinSteps ,
                             summarise=FALSE ,
                             plots=FALSE )
  metCodaSamples = as.mcmc.list( metRunJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( metCodaSamples , file=paste(saveName,"-Met-Mcmc.Rdata",sep="") )
  }
  #-----------------------------------------------------------------------------
  # Diagnostics of METRIC model:
  if ( !TRUE ) {
    metParameterNames = varnames(metCodaSamples) 
    for ( parName in metParameterNames ) {
      diagMCMC( codaObject=metCodaSamples , parName=parName ,  
                saveName=fileNameRoot , saveType="pdf" )
    }
  }
  #-----------------------------------------------------------------------------
  # Examine posterior distribution of METRIC model:
  metMcmcMat = as.matrix( metCodaSamples )
  metChainLength = nrow(metMcmcMat)
  
  # Open window for plot:
  openGraph(width=3.5*nG,height=2.5*1+0.75)
  layout(matrix(1:nG,nrow=1))
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) , oma=c(0,0,3.5,0) ) # , xpd=NA )
  histBreakInc = 1/nQ
  if ( nQ > 1 ) { xlabText="Mean of Ordinal Resp." } else { xlabText="Ordinal Resp." }
  # make histograms of data:
  for ( gIdx in 1:nG ) {
    ydata = meanY$x[meanY$g==gIdx]
    histInfo = hist( ydata , main=paste0("Group ",gIdx) , 
                     xlab=xlabText , ylab="Probability" ,
                     freq=FALSE , col="pink" , border="white" , cex.lab=1.5 ,
                     breaks=seq(-histBreakInc/2 , 
                                max(nYlevels)+1+histBreakInc/2 , 
                                by=histBreakInc ) , 
                     ylim=c(0,likertScaleDensMax) # NB: must set manually at this time
    )
    #text( max(histInfo$breaks) , likertScaleDensMax , adj=c(1,1) ,
    #      labels=bquote(N==.(length(ydata))) )
    xcomb = seq(0,max(nYlevels)+1,length=201)
    # *** Replace next line with smattering from posterior **
    for ( rIdx in floor(seq(1,nrow(metMcmcMat),length=20)) ) {
      lines( xcomb , dnorm( xcomb , 
                            metMcmcMat[rIdx,paste0("mu[",gIdx,"]")] ,
                            metMcmcMat[rIdx,paste0("sigma[",gIdx,"]")] ) , 
             col="skyblue" )
    }
  } # end for( gIdx in 1:nG )
  # Compute Bayesian effect size and HDI for two groups:
  if ( nG==2 ) {
    m1 = metMcmcMat[,"mu[1]"]
    m2 = metMcmcMat[,"mu[2]"]
    s1 = metMcmcMat[,"sigma[1]"]
    s2 = metMcmcMat[,"sigma[2]"]
    cohenD = (m1-m2)/sqrt((s1^2+s2^2)/2)
    hdiCohenD = HDIofMCMC( cohenD )
    muD = (m1-m2)
    hdiMuD = HDIofMCMC( muD )
    mtext( text=bquote( "d = "* .(round(median(cohenD),2)) 
                        *", 95%HDI = ["* .(round(hdiCohenD[1],2)) 
                        *","* .(round(hdiCohenD[2],2)) *"]" 
                        # *". Diff mu: "* .(round(muD,2))
                        # *"["* .(round(hdiMuD[1],2)) 
                        # *","* .(round(hdiMuD[2],2)) *"]" 
                        ) ,
           at=c(0.5) , cex=1.5 , side=3 , outer=TRUE , 
           adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
  } else {
    mtext( text="Bayesian Metric ... title placeholder" , 
           at=c(0.5) , cex=1.5 , side=3 , outer=TRUE , 
           adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
  }
  saveGraph2( paste0(fileNameRoot,"-LikertScale-Bayes") , typeVec=c("eps","pdf") )
  
  
  
  if ( TRUE ) { # Do ordered-probit analysis
  #=============================================================================
  # THE ORDERED-PROBIT MODEL FOR MULTIPLE ITEMS.
  #
  # N.B.: THIS MODEL ASSUMES THAT ALL ITEMS ARE SCORED WITH
  # POSITIVE CORRELATIONS WITH EACH OTHER.
  ordModelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dcat( pr[i,1:nYlevels[q[i]]] )
      pr[i,1] <- pnorm( thresh[q[i],1] ,
                        mu[g[i]] , 1/sigma[g[i]]^2 )
      for ( k in 2:(nYlevels[q[i]]-1) ) {
        pr[i,k] <- max( 0 ,  pnorm( thresh[q[i], k ]  ,
                                    mu[g[i]] , 1/sigma[g[i]]^2 )
                           - pnorm( thresh[q[i],k-1] ,
                                    mu[g[i]] , 1/sigma[g[i]]^2 ) )
      }
      pr[i,nYlevels[q[i]]] <- (
         1 - pnorm( thresh[q[i],nYlevels[q[i]]-1] ,
                    mu[g[i]] , 1/sigma[g[i]]^2 ) )
    }
    for ( gIdx in 1:nG ) { 
      mu[gIdx] ~ dnorm( (1+max(nYlevels))/2 , 1/(max(nYlevels))^2 )
      sigma[gIdx] ~ dunif( 0.01 , max(nYlevels)*10 )
    }
    # Prior on thresh[q,k]. Stochastic for all except thresh[1,1] and thresh[1,last].
    for ( qIdx in 1 ) {
      for ( kIdx in 2:(nYlevels[qIdx]-2) ) { # 1 and nYlevels-1 are fixed
        thresh[qIdx,kIdx] ~ dnorm( kIdx+0.5 , 1/2^2 )
      }
    }
    for ( qIdx in 2:nQ ) {
      for ( kIdx in 1:(nYlevels[qIdx]-1) ) { 
        thresh[qIdx,kIdx] ~ dnorm( kIdx+0.5 , 1/2^2 )
      }
    }
}
  " # close quote for ordModelString
  # Write out ordModelString to a text file
  writeLines( ordModelString , con="OrdAsOrdAndMet-OrdModel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS FOR ORDERED-PROBIT MODEL.
  parameters = c( "mu" , "sigma" , "thresh" )
  numSavedSteps = 20000 # 20000 
  thinSteps = 5 # 5 
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  saveName=fileNameRoot 
  runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
  nChains=nChainsDefault # from DBDA2E-utilities
  
  ordRunJagsOut <- run.jags( method="parallel" , # runjagsMethod ,
                          model="OrdAsOrdAndMet-OrdModel.txt" , 
                          monitor=parameters , 
                          data=ordDataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  ordCodaSamples = as.mcmc.list( ordRunJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( ordCodaSamples , file=paste(saveName,"-Ord-Mcmc.Rdata",sep="") )
  }
  
  #-----------------------------------------------------------------------------
  # Diagnostics:
  if ( !TRUE ) {
    ordParameterNames = varnames(ordCodaSamples) 
    for ( parName in ordParameterNames ) {
      diagMCMC( codaObject=ordCodaSamples , parName=parName ,  
                saveName=fileNameRoot , saveType="pdf" )
    }
  }
  
  #-----------------------------------------------------------------------------
  # Examine posterior distribution of ordered-probit variant:
  ordMcmcMat = as.matrix( ordCodaSamples )
  ordChainLength = nrow(ordMcmcMat)
  
  # Posterior mu's and sigma's:
  openGraph(height=2.5*nG,width=7)
  par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
  layout( matrix(1:(2*nG),nrow=nG,byrow=TRUE) )
  if ( nG >= 2 ) {
    muLim = range( ordMcmcMat[, grep("^mu\\[",colnames(ordMcmcMat)) ] )
    sigmaLim = range( ordMcmcMat[, grep("^sigma\\[",colnames(ordMcmcMat)) ] )
    for ( gIdx in 1:nG ) {
      plotPost( ordMcmcMat[,paste0("mu[",gIdx,"]")] , xlab=paste0("mu[",gIdx,"]") ,
                xlim=muLim , main=paste0("Group ",gIdx) )
      plotPost( ordMcmcMat[,paste0("sigma[",gIdx,"]")] , xlab=paste0("sigma[",gIdx,"]") ,
                xlim=sigmaLim , main=paste0("Group ",gIdx) )
    }
  } else {
    muLim = range( ordMcmcMat[, grep("^mu",colnames(ordMcmcMat)) ] )
    sigmaLim = range( ordMcmcMat[, grep("^sigma",colnames(ordMcmcMat)) ] )
    for ( gIdx in 1:nG ) {
      plotPost( ordMcmcMat[,paste0("mu")] , xlab=paste0("mu") ,
                xlim=muLim , main=paste0("Group ",gIdx) )
      plotPost( ordMcmcMat[,paste0("sigma")] , xlab=paste0("sigma") ,
                xlim=sigmaLim , main=paste0("Group ",gIdx) )
    }
  }
  saveGraph2( paste0(fileNameRoot,"-MuSigma") , typeVec=c("eps","pdf") )
  
  # Posterior effect size:
  if ( nG == 2 ) {
    postEffSz = ( ( ordMcmcMat[,"mu[1]"] - ordMcmcMat[,"mu[2]"] ) 
                  / sqrt( ( ordMcmcMat[,"sigma[1]"]^2 + ordMcmcMat[,"sigma[2]"]^2 ) / 2 ) )
    openGraph(width=7,height=4)
    plotPost( postEffSz , xlab="Effect Size" , ROPE=c(-0.1,0.1) )
    saveGraph2( paste0(fileNameRoot,"-EffSz") , typeVec=c("eps","pdf") )
  }
  if ( nG == 1 ) {
    postEffSz = ( ordMcmcMat[,"mu"] - nullHypVal ) / ordMcmcMat[,"sigma"]
    openGraph(width=3.5,height=4)
    plotPost( postEffSz , xlab="Effect Size" , ROPE=c(-0.1,0.1) )
    saveGraph2( paste0(fileNameRoot,"-EffSz") , typeVec=c("eps","pdf") )
  }
  
  # Thresholds:
  openGraph(height=min(2.5*nQ,14),width=7)
  par( mar=c(3.5,3.5,2,1) , mgp=c(2.25,0.7,0) )
  layout( matrix(1:nQ,nrow=nQ) )
  for ( qIdx in 1:nQ ) {
    threshCols = paste0("thresh[",qIdx,",",1:(nYlevels[qIdx]-1),"]")
    threshMean = rowMeans( ordMcmcMat[,threshCols] )
    xLim = range(ordMcmcMat[,threshCols])
    nPtToPlot = 1000
    plotIdx = floor(seq(1,nrow(ordMcmcMat),length=nPtToPlot))
    plot( ordMcmcMat[plotIdx,threshCols[1]] , threshMean[plotIdx] , col="skyblue" ,
          xlim=xLim , xlab="Threshold" , ylab="Mean Threshold" , 
          main=paste0("Item ",qIdx) )
    abline(v=mean(ordMcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="skyblue")
    for ( i in 2:length(threshCols) ) {
      points( ordMcmcMat[plotIdx,threshCols[i]] , threshMean[plotIdx] , col="skyblue" )
      abline(v=mean(ordMcmcMat[plotIdx,threshCols[i]]),lty="dashed",col="skyblue")
    }
  }  
  saveGraph2( paste0(fileNameRoot,"-Thresh") , typeVec=c("eps","pdf") )
  
  # Posterior predictive: Histograms of ordinal responses with posterior predicted
  # probabilities superimposed.
  openGraph(height=min(2.5*nQ+0.75,14),width=3.5*nG)
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) , oma=c(0,0,3.5,0) ) # , xpd=NA )
  layout( matrix(1:(nG*nQ),nrow=nQ,ncol=nG,byrow=FALSE) )
  for ( gIdx in 1:nG ) {
    for ( qIdx in 1:nQ ) {
      # Data histogram:
      thisY = y[ g==gIdx & q==qIdx ]
      xLim = c( min(thisY)-0.5 , max(thisY)+0.5 )
      xBreaks = seq( xLim[1] , xLim[2] , 1 )  
      histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                       xlim=xLim , breaks=xBreaks ,
                       xlab="Ordinal Resp." , ylab="p(Resp)" , cex.lab=1.5 , 
                       col="pink" , border="white" , 
                       #yaxt="n" , 
                       main=paste0("Item ",qIdx,", Group ",gIdx) )
      # Posterior predicted probabilities:
      outProb=matrix(0,nrow=ordChainLength,ncol=max(thisY))
      for ( stepIdx in 1:ordChainLength ) {
        if ( nG > 1 ) {
            threshCumProb = pnorm( 
              ordMcmcMat[ stepIdx , paste0("thresh[",qIdx,",",1:(max(thisY)-1),"]") ] ,
              ordMcmcMat[ stepIdx , paste0("mu[",gIdx,"]") ] ,
              ordMcmcMat[ stepIdx , paste0("sigma[",gIdx,"]") ] )
        } else { # if nG == 1
            threshCumProb = pnorm( 
              ordMcmcMat[ stepIdx , paste0("thresh[",qIdx,",",1:(max(thisY)-1),"]") ] ,
              ordMcmcMat[ stepIdx , paste0("mu") ] ,
              ordMcmcMat[ stepIdx , paste0("sigma") ] )
        }
        outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
      }
      outHdi = apply( outProb , 2 , HDIofMCMC )
      outMean = apply( outProb , 2 , median , na.rm=TRUE )
      show(outMean)
      points( x=1:max(thisY) , y=outMean  , pch=19 , cex=1.5 , col="skyblue" )
      segments( x0=1:max(thisY) , y0=outHdi[1,] , 
                x1=1:max(thisY) , y1=outHdi[2,] , lwd=4 , col="skyblue" )
    }
  }
  mtext( text=bquote(list( d==.(round(median(postEffSz),2)) ,
                           "95% HDI = [" * .(round(HDIofMCMC(postEffSz)[1],2))
                           * "," * .(round(HDIofMCMC(postEffSz)[2],2)) * "]" ) ) , 
         at=c(0.5) , cex=ifelse(nG>=2,1.5,1.0) , side=3 , outer=TRUE , 
         adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
  saveGraph2( paste0(fileNameRoot,"-PostPred") , typeVec=c("eps","pdf") )
  
  if ( simIdx=="1Q1Gfa" ) {
    # posterior of difference in predicted probs of nullHypVal+1 vs nullHypVal-1
    outProb=matrix(0,nrow=ordChainLength,ncol=max(thisY))
    for ( stepIdx in 1:ordChainLength ) {
      threshCumProb = pnorm( 
        ordMcmcMat[ stepIdx , paste0("thresh[",qIdx,",",1:(max(thisY)-1),"]") ] ,
        ordMcmcMat[ stepIdx , paste0("mu") ] ,
        ordMcmcMat[ stepIdx , paste0("sigma") ] )
      outProb[stepIdx,] = c(threshCumProb, 1) - c(0, threshCumProb)
    }
    oneUp = round(nullHypVal)+1
    oneDown = round(nullHypVal)-1
    diffPredProb = outProb[,oneUp] - outProb[,oneDown]
    openGraph(height=min(2.5*nQ+0.75,14),width=3.5*nG)
    par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) ) # , xpd=NA )
    plotPost( diffPredProb , xlab="Difference" , 
              main=paste0("p('",oneUp,"') - p('",oneDown,"')") )
    saveGraph2( paste0(fileNameRoot,"-PredDiff") , typeVec=c("eps","pdf") )
  }
  
  } # end if for doing ordered-probit analysis
  
} # end of for(simIdx...
