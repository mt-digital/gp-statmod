####################################################################################
#
# OrderedProbitModel.R
# John K. Kruschke, January-August 2018. 
# Sept 18, 2018: Revised by John K. Kruschke; graphical output revised.
# Sept 28, 2018: Revised by John K. Kruschke; graphical output revised more, and
#                argument caseIDcolName added.
#
### Description:
#
# This R script defines a function, ordinalAndMetricAnalysis(), that analyzes
# groups of ordinal data with an ordered probit model and an inappropriate
# metric model.
#
# The function is used for an example in the article, "Analyzing ordinal data
# with metric models: What could possibly go wrong?" by Torrin M. Liddell and
# John K. Kruschke. Files can be found at https://osf.io/53ce9/
#
### Usage:
#
# source("OrderedProbitModel.R") # source this file to load the function
#
# ordinalResults = ordinalAndMetricAnalysis( dataFileName , yColNames ,
#                                            caseIDcolName=NULL ,
#                                            compareCases=NULL , 
#                                            graphFileType="pdf" ,
#                                            hierarchSD=FALSE )
#
# For a detailed example of using this function, see the script
# OrderedProbitModel-Example.R
#
### Arguments:
#
# dataFileName: The name of the data file. Assumed to be comma separated values.
#               See below for more details.
#
# yColNames: A vector of column names of the frequencies of ordinal levels.
#
# caseIDcolName: Column name of row identifiers or labels. Optional, defaults to
#                NULL.
#
# compareCases: A list of case pairs to compare. Optional, defaults to NULL.
#
# graphFileType: The graphics format of saved graphs. Defaults to "pdf".
#
# hierarchSD: Whether or not to put hierarchical structure on the group standard
#             deviations. Logical TRUE/FALSE, defaults to FALSE.
#
### Output (Value):
#
# A list of two components:
#
# OrdSummaryMat: A matrix, with rows being the parameters in the ordinal model
# and columns being summary statistics of the posterior distribution.
#
# OrdMcmcMat: A matrix, with rows being steps in the MCMC chain and columns
# being the chain number and the parameters in the ordinal model.
#
# N.B.: The function saves the output to the working directory and also produces
# numerous graphs that are saved to the working directory.
#
### Important Dependencies:
#
# This script assumes a particular data structure: The dependent (predicted)
# variable is an ordinal value from a single item. The independent (predictor)
# variable is group. For example, the dependent variable could be a 1 to 5 star
# rating with the independent variable being different movies. (The model makes
# no assumptions about repeated measures, that is, the model does not represent
# any identities of responders that might be the same across groups.) For each
# group, the data are specified as a count of each level of rating scale. The
# data file must be strucured with named columns, one column for the group
# identifier, and K columns for the counts of each level of the 1 through K
# ordinal levels. See OrderedProbitModel-Example.R for more info.
#
# The file DBDA2E-utilities.R is needed for various functions such as
# gammaShRaFromModeSD(), diagMCMC(), saveGraph(), plotPost(), HDIofMCMC(), etc.
# It also installs and loads the R package runjags, so you must be connected to
# the internet if you do not already have runjags installed. The file
# DBDA2E-utilities.R should accompany this file at the article's OSF repository,
# or you can find it in the DBDA2Eprograms.zip file from Step 5 of
# https://sites.google.com/site/doingbayesiandataanalysis/software-installation
# Get the DBDA2Eprograms.zip file, extract it, and copy the file
# DBDA2E-utilities.R into the same folder (working directory) as this file.
#
### Example:
#
# For a detailed example of using this function, see the script
# OrderedProbitModel-Example.R
#
####################################################################################

source("DBDA2E-utilities.R")

ordinalAndMetricAnalysis = function( dataFileName , yColNames , 
                                     caseIDcolName=NULL ,
                                     compareCases=NULL , 
                                     graphFileType=c("pdf","png","eps","jpg")[1] ,
                                     hierarchSD=FALSE ) {
  fileNameRoot = paste0("OrderedProbitModel-",
                        ifelse(hierarchSD,"hierSD-",""),
                        strsplit(dataFileName,"\\.")[[1]][1])
  #-----------------------------------------------------------------------------------
  # Set up data for JAGS
  datFrm = read.csv( dataFileName , header=TRUE )
  # yColNames must be vector of names IN ORDER from 1 to nYlevels
  y = as.matrix(datFrm[,yColNames]) # y is matrix
  z = rowSums(y)
  x = 1:nrow(datFrm)
  Ncases = nrow(datFrm)
  # Threshold 1 and nYlevels-1 are fixed; other thresholds are estimated.
  # This allows all parameters to be interpretable on the response scale.
  nYlevels = length(yColNames)  
  thresh = rep(NA,nYlevels-1)
  thresh[1] = 1 + 0.5
  thresh[nYlevels-1] = nYlevels-1 + 0.5
  # Specify the data in a list, for later shipment to JAGS:
  gammaShRa = unlist( gammaShRaFromModeSD( mode=3.0 , sd=3.0 ) )
  # For the ordered probit model:
  if( hierarchSD ) {
    OrdDataList = list(
      y = y ,
      nYlevels = nYlevels ,
      thresh = thresh ,
      x = x ,
      z = z ,
      Ncases = Ncases ,
      gammaShRa = gammaShRa
    ) 
  } else {
    OrdDataList = list(
      y = y ,
      nYlevels = nYlevels ,
      thresh = thresh ,
      x = x ,
      z = z ,
      Ncases = Ncases 
    )
  }
  # For the metric model:
  if ( hierarchSD ) {
    MetDataList = list(
      y = y ,
      nYlevels = nYlevels ,
      x = x ,
      Ncases = Ncases ,
      gammaShRa = gammaShRa ,
      zeros = rep(0,Ncases) , # for Poisson trick in JAGS
      C = 10 # 100  # for Poisson trick in JAGS
    ) 
  } else {
    MetDataList = list(
      y = y ,
      nYlevels = nYlevels ,
      x = x ,
      Ncases = Ncases ,
      zeros = rep(0,Ncases) , # for Poisson trick in JAGS
      C = 10 # 100  # for Poisson trick in JAGS
    )
  }
  
  #-----------------------------------------------------------------------------------
  # THE *METRIC* MODEL:
  ldsnString = paste( "logdensity.norm( " , 
                      1:nYlevels , 
                      ", mu[x[i]] , 1/sigma[x[i]]^2 )*y[i," ,
                      1:nYlevels , 
                      "]" ,
                      collapse=" + " )
  modelString = paste("
  model {
    for ( i in 1:Ncases ) {
      ldsn[i] <- (", ldsnString , ")
      zeros[i] ~ dpois( -ldsn[i] + C )
    }
    for ( j in 1:Ncases ) { 
      mu[j] ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )
      sigma[j] ~ dgamma( sigmaSh , sigmaRa )
    }
    sigmaSh <- 1 + sigmaMode * sigmaRa
    sigmaRa <- ( ( sigmaMode + sqrt( sigmaMode^2 + 4*sigmaSD^2 ) ) 
                 / ( 2*sigmaSD^2 ) )
    ",
    ifelse( hierarchSD ,
    "sigmaMode ~ dgamma( gammaShRa[1] , gammaShRa[2] ) 
     sigmaSD ~ dgamma( gammaShRa[1] , gammaShRa[2] ) " ,
    "sigmaMode <- 3.0
     sigmaSD <- 3.0" ) ,
  "}") # close quote for modelString paste
  # Write out modelString to a text file
  writeLines( modelString , con=paste0(fileNameRoot,"-MetModel.txt") )
  #-----------------------------------------------------------------------------------
  # Run the metric model:
  if ( hierarchSD ) {
    parameters = c( "mu" , "sigma" , "sigmaMode" , "sigmaSD" ) 
  } else {
    parameters = c( "mu" , "sigma" ) 
  }
  adaptSteps = 500
  burnInSteps = 1000
  numSavedSteps = 24000
  thinSteps = 5
  nChains = 4
  runJagsOut <- run.jags( method="parallel" ,
                          model=paste0(fileNameRoot,"-MetModel.txt") , 
                          monitor=parameters , 
                          data=MetDataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  MetCodaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(fileNameRoot) ) {
    save( MetCodaSamples , file=paste0(fileNameRoot,"-MetModel-Mcmc.Rdata") )
  }
  #-----------------------------------------------------------------------------------
  # Display diagnostics of chains:
  MetParameterNames = varnames(MetCodaSamples) 
  if ( !TRUE ) {
    for ( parName in MetParameterNames ) {
      diagMCMC( codaObject=MetCodaSamples , parName=parName ,
                saveName=NULL , #paste0(fileNameRoot,"-MetModel") , 
                saveType=graphFileType )
    }
  }
  
  #-----------------------------------------------------------------------------------
  # THE *ORDERED PROBIT* MODEL:
  modelString = paste0("
  model {
    for ( i in 1:Ncases ) {
      y[i, ] ~ dmulti( pr[i,1:nYlevels] , z[i] )
      pr[i,1] <- pnorm( thresh[1] , mu[x[i]] , 1/sigma[x[i]]^2 )
      for ( k in 2:(nYlevels-1) ) {
        pr[i,k] <- max( 0 ,  pnorm( thresh[ k ] , mu[x[i]] , 1/sigma[x[i]]^2 )
                            - pnorm( thresh[k-1] , mu[x[i]] , 1/sigma[x[i]]^2 ) )
      }
      pr[i,nYlevels] <- 1 - pnorm( thresh[nYlevels-1] , mu[x[i]] , 1/sigma[x[i]]^2 )
    }
    for ( j in 1:Ncases ) { 
      mu[j] ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )
      sigma[j] ~ dgamma( sigmaSh , sigmaRa )
    }
    sigmaSh <- 1 + sigmaMode * sigmaRa
    sigmaRa <- ( ( sigmaMode + sqrt( sigmaMode^2 + 4*sigmaSD^2 ) ) 
                  / ( 2*sigmaSD^2 ) ) ",
    ifelse( hierarchSD ,
      "sigmaMode ~ dgamma( gammaShRa[1] , gammaShRa[2] ) 
       sigmaSD ~ dgamma( gammaShRa[1] , gammaShRa[2] ) " ,
      "sigmaMode <- 3.0
       sigmaSD <- 3.0" ) , " # open quote for next line
    for ( k in 2:(nYlevels-2) ) {  # 1 and nYlevels-1 are fixed, not stochastic
      thresh[k] ~ dnorm( k+0.5 , 1/2^2 )
    }
  }") # close quote for modelString paste
  # Write out modelString to a text file
  writeLines( modelString , con=paste0(fileNameRoot,"-OrdModel.txt") )
  #-----------------------------------------------------------------------------------
  # Run the ordered-probit model:
  if ( hierarchSD ) {
    parameters = c( "mu" , "sigma" , "thresh" , "sigmaMode" , "sigmaSD" ) 
  } else {
    parameters = c( "mu" , "sigma" , "thresh" ) 
  }
  adaptSteps = 500
  burnInSteps = 1000
  numSavedSteps = 24000
  thinSteps = 5
  nChains = 4
  runJagsOut <- run.jags( method="parallel" ,
                          model=paste0(fileNameRoot,"-OrdModel.txt") , 
                          monitor=parameters , 
                          data=OrdDataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  OrdCodaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(fileNameRoot) ) {
    save( OrdCodaSamples , file=paste0(fileNameRoot,"-OrdModel-Mcmc.Rdata") )
  }
  #-----------------------------------------------------------------------------------
  # For output, consider case names:
  if ( !is.null(caseIDcolName) ) { 
    caseNames = as.character( datFrm[,caseIDcolName]  )
  } else {
    caseNames = paste( "Case" , 1:Ncases )
  }
  titleCaseNameMaxLength = 10 # maximum characters to display from label
  #-----------------------------------------------------------------------------------
  # Display diagnostics of chains:
  OrdParameterNames = varnames(OrdCodaSamples) 
  if ( !TRUE ) {
    for ( parName in OrdParameterNames ) {
      diagMCMC( codaObject=OrdCodaSamples , parName=parName ,
                saveName=NULL , #paste0(fileNameRoot,"-OrdModel") , 
                saveType=graphFileType )
    }
  }
  
  #-----------------------------------------------------------------------------------
  # Display and output posterior information:
  OrdMcmcMat = as.matrix(OrdCodaSamples,chains=TRUE)
  OrdChainLength = NROW( OrdMcmcMat )
  MetMcmcMat = as.matrix(MetCodaSamples,chains=TRUE)
  MetChainLength = NROW( MetMcmcMat )
  
  # Create numerical summary object for returning from this function:
  postStats = c("Mean","Median","Mode","HDImass","HDIlow","HDIhigh","ESS") # +"psrf" below
  OrdSummaryMat = matrix( 0 , nrow=length(OrdParameterNames) , ncol=length(postStats)+1 )
  colnames(OrdSummaryMat) = c(postStats,"psrf")
  rownames(OrdSummaryMat) = OrdParameterNames
  for ( parName in OrdParameterNames ) {
    OrdSummaryMat[parName,postStats] = summarizePost(OrdMcmcMat[,parName],credMass=0.95)[postStats]
    OrdSummaryMat[parName,"psrf"] = coda::gelman.diag( OrdCodaSamples[,parName] )$psrf[1,1]
  }
  #show(OrdSummaryMat)
  
  #-----------------------------------------------
  # Plot thresholds of ordered-probit model:
  threshCols = grep("thresh",colnames(OrdMcmcMat),value=TRUE)
  threshMean = rowMeans( OrdMcmcMat[,threshCols] )
  xLim = range(OrdMcmcMat[,threshCols])
  xLim = c( xLim[1]-0.5 , xLim[2]+0.5 )
  xTickPos = apply( OrdMcmcMat[,threshCols] , 2 , median )
  nPtToPlot = 1000
  plotIdx = floor(seq(1,nrow(OrdMcmcMat),length=nPtToPlot))
  openGraph(width=7.0,height=5.0)
  layout( matrix( c( rep(1,length(xTickPos)-2) , 1+(1:(length(xTickPos)-2)) ) ,
                  nrow=2 , byrow=TRUE ) )
  par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
  # Plot thresh mean x thresh:
  plot( OrdMcmcMat[plotIdx,threshCols[1]] , threshMean[plotIdx] ,
        xlab="Latent Scale" , ylab="Mean Threshold" , cex.lab=1.5 ,
        xlim=xLim , xaxt="n" , #xaxp=xTickPos ,
        main=bquote("Thresholds ("*theta[k]*")") , cex.main=1.75 ,
        col="skyblue" )
  abline(v=mean(OrdMcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="skyblue")
  for ( i in 2:length(threshCols) ) {
    points( OrdMcmcMat[plotIdx,threshCols[i]] , threshMean[plotIdx] , col="skyblue" )
    abline(v=mean(OrdMcmcMat[plotIdx,threshCols[i]]),lty="dashed",col="skyblue")
  }
  axis( 1 , at=xTickPos , labels=round(xTickPos,2) , cex.axis=1.5 )
  for ( k in 1:(length(threshCols)+1) ) {
    text( x=mean(c(1.0,xTickPos,length(xTickPos)+1)[c(k,k+1)]) , 
          y=mean(threshMean) , labels=paste0("'",k,"'") , cex=1.5 )
  }
  # Plot marginals of individual thresh's:
  for ( tIdx in 1+(1:(length(xTickPos)-2)) ) {
    plotPost( OrdMcmcMat[,threshCols[tIdx]] , main=bquote(theta[.(tIdx)]) ,
              xlab="Latent Scale" , cex.main=1.75 )
  }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-OrdModel-Thresh"), type=graphFileType)
  }

  #-----------------------------------------------
  if ( hierarchSD ) {
    # Plot sigmaMode, sigmaSD of ordered-probit model:
    openGraph(width=7,height=3.5)
    layout(matrix(1:2,nrow=1,byrow=TRUE))
    par( mar=c(4,2,4,1) , mgp=c(2.25,0.7,0) )
    plotPost( OrdMcmcMat[,"sigmaMode"] , xlab=bquote(omega[sigma[i]]) ,
              main="Ordered-Probit Modal Sigma" )
    plotPost( OrdMcmcMat[,"sigmaSD"] , xlab=bquote(sigma[sigma[i]]) ,
              main="Ordered-Probit SD of Sigma[i]" )
    if ( !is.null(fileNameRoot) ) {
      saveGraph( file=paste0(fileNameRoot,"-OrdModel-sigmaMode"), type=graphFileType)
    }
    
    # Plot sigmaMode, sigmaSD of metric model:
    openGraph(width=7,height=3.5)
    layout(matrix(1:2,nrow=1,byrow=TRUE))
    par( mar=c(4,2,4,1) , mgp=c(2.25,0.7,0) )
    plotPost( MetMcmcMat[,"sigmaMode"] , xlab=bquote(omega[sigma[i]]) ,
              main="Metric Modal Sigma" )
    plotPost( MetMcmcMat[,"sigmaSD"] , xlab=bquote(sigma[sigma[i]]) ,
              main="Metric SD of Sigma[i]" )
    if ( !is.null(fileNameRoot) ) {
      saveGraph( file=paste0(fileNameRoot,"-MetModel-sigmaMode"), type=graphFileType)
    }
  }
  #-----------------------------------------------
  # Plot data with ORDERED PROBIT posterior predictions and latent mu,sigma annotation:
  nCell = length(z)
  nCol=ceiling(sqrt(nCell))
  nRow=ceiling(nCell/nCol)
  nCol = min(nCol,6)
  nRow = min(nRow,6)
  inchPerPanel = 1.7
  openGraph(width=nCol*inchPerPanel,height=nRow*inchPerPanel)
  layout(matrix(1:(nCol*nRow),nrow=nRow,ncol=nCol,byrow=TRUE))
  par( mar=c(3,3,4,1) , mgp=c(1.8,0.7,0) )
  #sortMuIdx = sort( apply( OrdMcmcMat[,grep("mu",colnames(OrdMcmcMat))] , 2 , median ) ,
  #                  index.return=TRUE )$ix
  sortMuIdx = 1:length(grep("mu",colnames(OrdMcmcMat)))
  for ( caseIdx in sortMuIdx[1:min(length(sortMuIdx),(nRow*nCol))] ) {
    plot( 1:length(yColNames) , y[caseIdx,] , 
          xlab="Rating" , ylab="Frequency" , cex.lab=0.75 ,
          xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(y[caseIdx,])) ,
          type="h" , lend=1 , lwd=10 , col="pink" )
    medMu = median(OrdMcmcMat[,paste0("mu[",caseIdx,"]")])
    medSigma = median(OrdMcmcMat[,paste0("sigma[",caseIdx,"]")])
    if ( nchar(caseNames[caseIdx]) > titleCaseNameMaxLength ) {
      titleCaseName = paste0( caseIdx , ":" ,
                              substr(caseNames[caseIdx],1,titleCaseNameMaxLength) ,
                              "." )
    } else {
      titleCaseName = paste0( caseIdx , ":" , 
                              caseNames[caseIdx] )
    }
    title( main=bquote( atop( .(titleCaseName)
                              *", N="*.(z[caseIdx]) *"  "
                              ,
                              mu*"="*.(round(medMu,2))
                              *", "*sigma*"="*.(round(medSigma,2))
    ) ) , cex.main=1.0 )
    # superimpose posterior predictions:
    predProb = matrix(0,nrow=nrow(OrdMcmcMat),ncol=length(yColNames))
    for ( stepIdx in 1:nrow(OrdMcmcMat) ) {
      threshVec = OrdMcmcMat[stepIdx,grep("thresh",colnames(OrdMcmcMat))]
      thisMu = OrdMcmcMat[stepIdx,paste0("mu[",caseIdx,"]")]
      thisSigma = OrdMcmcMat[stepIdx,paste0("sigma[",caseIdx,"]")]
      predProb[stepIdx,] = ( pnorm( (c(threshVec,Inf)-thisMu)/thisSigma )
                             - pnorm( (c(-Inf,threshVec)-thisMu)/thisSigma ) )
    }
    points( 1:length(yColNames) , 
            apply(predProb,2,median) * sum(y[caseIdx,])  , 
            cex=1 , lwd=2 , col="skyblue" )
    for ( respIdx in 1:length(yColNames) ) {
      lines( rep(respIdx,2) , HDIofMCMC( predProb[,respIdx] ) * sum(y[caseIdx,]) , 
             col="skyblue" , lwd=3 )
    }
  }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-OrdModel-PostPred"), type=graphFileType)
  }
  
  #-----------------------------------------------
  # Plot data with METRIC posterior predictions and latent mu,sigma annotation:
  nCell = length(z)
  nCol=ceiling(sqrt(nCell))
  nRow=ceiling(nCell/nCol)
  nCol = min(nCol,6)
  nRow = min(nRow,6)
  inchPerPanel = 1.7
  openGraph(width=nCol*inchPerPanel,height=nRow*inchPerPanel)
  layout(matrix(1:(nCol*nRow),nrow=nRow,ncol=nCol,byrow=TRUE))
  par( mar=c(3,3,4,1) , mgp=c(1.8,0.7,0) )
  #sortMuIdx = sort( apply( MetMcmcMat[,grep("mu",colnames(MetMcmcMat))] , 2 , median ) ,
  #                  index.return=TRUE )$ix
  sortMuIdx = 1:length(grep("mu",colnames(MetMcmcMat)))
  for ( caseIdx in sortMuIdx[1:min(length(sortMuIdx),(nRow*nCol))] ) {
    plot( 1:length(yColNames) , y[caseIdx,] , 
          xlab="Rating" , ylab="Frequency" , cex.lab=0.75 ,
          xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(y[caseIdx,])) ,
          type="h" , lend=1 , lwd=10 , col="pink" )
    medMu = median(MetMcmcMat[,paste0("mu[",caseIdx,"]")])
    medSigma = median(MetMcmcMat[,paste0("sigma[",caseIdx,"]")])
    if ( nchar(caseNames[caseIdx]) > titleCaseNameMaxLength ) {
      titleCaseName = paste0( caseIdx , ":" ,
                              substr(caseNames[caseIdx],1,titleCaseNameMaxLength) ,
                              "." )
    } else {
      titleCaseName = paste0( caseIdx , ":" , 
                              caseNames[caseIdx] )
    }
    title( main=bquote( atop( .(titleCaseName)
                              *", N="*.(z[caseIdx]) *"  "
                              ,
                              mu*"="*.(round(medMu,2))
                              *", "*sigma*"="*.(round(medSigma,2))
    ) ) , cex.main=1.0 )
    # superimpose posterior predictions:
    xcomb = seq(0,nYlevels+1,length=201)
    for ( stepIdx in seq(1,nrow(MetMcmcMat),length=20) ) {
      thisMu = MetMcmcMat[stepIdx,paste0("mu[",caseIdx,"]")]
      thisSigma = MetMcmcMat[stepIdx,paste0("sigma[",caseIdx,"]")]
      lines( xcomb , dnorm(xcomb,thisMu,thisSigma)*sum(y[caseIdx,]) , 
             col="skyblue"  )
    }
  }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-MetModel-PostPred"), type=graphFileType)
  }
  
  #-----------------------------------------------
  
  # Graphs of comparisons:
  if ( !is.null(compareCases) ) {
    for ( compareIdx in 1:length(compareCases) ) {
      openGraph(width=9,height=5)
      layout(matrix(1:8,nrow=2,byrow=TRUE))
      par( mar=c(4,2,4,1) , mgp=c(2.25,0.7,0) )
      #................................................................
      # ORDINAL model row:
      for ( caseIdx in compareCases[[compareIdx]] ) {
        # plot data:
        plot( 1:length(yColNames) , y[caseIdx,] , 
              xlab="Rating" , ylab="Frequency" , cex.lab=1.5 ,
              xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(y[caseIdx,])) ,
              type="h" , lend=1 , lwd=10 , col="pink" )
        medMu = median(OrdMcmcMat[,paste0("mu[",caseIdx,"]")])
        medSigma = median(OrdMcmcMat[,paste0("sigma[",caseIdx,"]")])
        if ( nchar(caseNames[caseIdx]) > titleCaseNameMaxLength ) {
          titleCaseName = paste0( caseIdx , ":" ,
                                  substr(caseNames[caseIdx],1,titleCaseNameMaxLength) ,
                                  "." )
        } else {
          titleCaseName = paste0( caseIdx , ":" , 
                                  caseNames[caseIdx] )
        }
        title( main=bquote( atop( .(titleCaseName)
                                  *", N="*.(z[caseIdx]) *"  "
                                  ,
                                  mu*"="*.(round(medMu,2))
                                  *", "*sigma*"="*.(round(medSigma,2))
        ) ) , cex.main=1.33 )
        # superimpose ORDINAL-model posterior predictions:
        predProb = matrix(0,nrow=nrow(OrdMcmcMat),ncol=length(yColNames))
        for ( stepIdx in 1:nrow(OrdMcmcMat) ) {
          threshVec = OrdMcmcMat[stepIdx,grep("thresh",colnames(OrdMcmcMat))]
          thisMu = OrdMcmcMat[stepIdx,paste0("mu[",caseIdx,"]")]
          thisSigma = OrdMcmcMat[stepIdx,paste0("sigma[",caseIdx,"]")]
          predProb[stepIdx,] = ( pnorm( (c(threshVec,Inf)-thisMu)/thisSigma )
                                 - pnorm( (c(-Inf,threshVec)-thisMu)/thisSigma ) )
        }
        points( 1:length(yColNames) , 
                apply(predProb,2,median) * sum(y[caseIdx,])  , 
                cex=1.5 , lwd=2 , col="skyblue" )
        for ( respIdx in 1:length(yColNames) ) {
          lines( rep(respIdx,2) , HDIofMCMC( predProb[,respIdx] ) * sum(y[caseIdx,]) , 
                 col="skyblue" , lwd=3 )
        }
      }
      # plot ORDINAL-model posterior of difference of mu's:
      caseIdx1 = compareCases[[compareIdx]][1]
      caseIdx2 = compareCases[[compareIdx]][2]
      postInfo = plotPost( OrdMcmcMat[,paste0("mu[",caseIdx1,"]")] 
                           - OrdMcmcMat[,paste0("mu[",caseIdx2,"]")] ,
                           compVal=0 ,
                           xlab=bquote(mu[.(caseIdx1)]-mu[.(caseIdx2)]) , cex.lab=1.5 ,
                           main=paste0("Ordered Probit") )
      # plot ORDINAL-model posterior of difference of sigma's:
      postInfo = plotPost( OrdMcmcMat[,paste0("sigma[",caseIdx1,"]")] 
                           - OrdMcmcMat[,paste0("sigma[",caseIdx2,"]")] ,
                           compVal=0 ,
                           xlab=bquote(sigma[.(caseIdx1)]-sigma[.(caseIdx2)]) , cex.lab=1.5 ,
                           main=paste0("Ordered Probit") )
      #.......................................................................
      # Now plot METRIC-model row:    
      for ( caseIdx in compareCases[[compareIdx]] ) {
        #yVals = y[caseIdx,]
        #yd = rep(1:5,times=yVals)
        # plot data:
        plot( 1:length(yColNames) , y[caseIdx,] , 
              xlab="Rating" , ylab="Frequency" , cex.lab=1.5 ,
              xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(y[caseIdx,])) ,
              type="h" , lend=1 , lwd=10 , col="pink" )
        medMu = median(MetMcmcMat[,paste0("mu[",caseIdx,"]")])
        medSigma = median(MetMcmcMat[,paste0("sigma[",caseIdx,"]")])
        if ( nchar(caseNames[caseIdx]) > titleCaseNameMaxLength ) {
          titleCaseName = paste0( caseIdx , ":" ,
                                  substr(caseNames[caseIdx],1,titleCaseNameMaxLength) ,
                                  "." )
        } else {
          titleCaseName = paste0( caseIdx , ":" , 
                                  caseNames[caseIdx] )
        }
        title( main=bquote( atop( .(titleCaseName)
                                  *", N="*.(z[caseIdx]) *"  "
                                  ,
                                  mu*"="*.(round(medMu,2))
                                  *", "*sigma*"="*.(round(medSigma,2))
        ) ) , cex.main=1.33 )
        # Superimpose METRIC-model posterior pred normals:
        xcomb = seq(0,nYlevels+1,length=201)
        for ( stepIdx in seq(1,nrow(MetMcmcMat),length=20) ) {
          thisMu = MetMcmcMat[stepIdx,paste0("mu[",caseIdx,"]")]
          thisSigma = MetMcmcMat[stepIdx,paste0("sigma[",caseIdx,"]")]
          lines( xcomb , dnorm(xcomb,thisMu,thisSigma)*sum(y[caseIdx,]) , 
                 col="skyblue"  )
        }
      }
      # plot METRIC-model posterior of difference of mu's:
      caseIdx1 = compareCases[[compareIdx]][1]
      caseIdx2 = compareCases[[compareIdx]][2]
      postInfo = plotPost( MetMcmcMat[,paste0("mu[",caseIdx1,"]")] 
                           - MetMcmcMat[,paste0("mu[",caseIdx2,"]")] ,
                           compVal=0 ,
                           xlab=bquote(mu[.(caseIdx1)]-mu[.(caseIdx2)]) , cex.lab=1.5 ,
                           main=paste0("(Ord. as) Metric") )
      # plot METRIC-model posterior of difference of sigma's:
      postInfo = plotPost( MetMcmcMat[,paste0("sigma[",caseIdx1,"]")] 
                           - MetMcmcMat[,paste0("sigma[",caseIdx2,"]")] ,
                           compVal=0 ,
                           xlab=bquote(sigma[.(caseIdx1)]-sigma[.(caseIdx2)]) , cex.lab=1.5 ,
                           main=paste0("(Ord. as) Metric") )
      # Save graph
      if ( !is.null(fileNameRoot) ) {
        saveGraph( file=paste0(fileNameRoot,"-Compare-",
                               compareCases[[compareIdx]][1],"-",
                               compareCases[[compareIdx]][2]) , 
                   type=graphFileType )
      }
    }
  }

  # Graphs of individual cases within comparisons:
  if ( !is.null(compareCases) ) {
    for ( compareIdx in 1:length(compareCases) ) {
      for ( caseIdx in compareCases[[compareIdx]] ) {
        # Open window for this case:
        openGraph(width=7,height=3)
        layout(matrix(1:3,nrow=1,byrow=TRUE))
        par( mar=c(4,4,4,1) , mgp=c(2.25,0.7,0) )
        # plot data for this case:
        plot( 1:length(yColNames) , y[caseIdx,] , 
              xlab="Rating" , ylab="Frequency" , cex.lab=1.5 ,
              xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(y[caseIdx,])) ,
              type="h" , lend=1 , lwd=10 , col="pink" )
        medMu = median(OrdMcmcMat[,paste0("mu[",caseIdx,"]")])
        medSigma = median(OrdMcmcMat[,paste0("sigma[",caseIdx,"]")])
        if ( nchar(caseNames[caseIdx]) > titleCaseNameMaxLength ) {
          titleCaseName = paste0( caseIdx , ":" ,
                                  substr(caseNames[caseIdx],1,titleCaseNameMaxLength) ,
                                  "." )
        } else {
          titleCaseName = paste0( caseIdx , ":" , 
                                  caseNames[caseIdx] )
        }
        title( main=bquote( atop( .(titleCaseName)
                                  *", N="*.(z[caseIdx]) *"  "
                                  ,
                                  mu*"="*.(round(medMu,2))
                                  *", "*sigma*"="*.(round(medSigma,2))
        ) ) , cex.main=1.33 )
        # superimpose ORDINAL-model posterior predictions for this case:
        predProb = matrix(0,nrow=nrow(OrdMcmcMat),ncol=length(yColNames))
        for ( stepIdx in 1:nrow(OrdMcmcMat) ) {
          threshVec = OrdMcmcMat[stepIdx,grep("thresh",colnames(OrdMcmcMat))]
          thisMu = OrdMcmcMat[stepIdx,paste0("mu[",caseIdx,"]")]
          thisSigma = OrdMcmcMat[stepIdx,paste0("sigma[",caseIdx,"]")]
          predProb[stepIdx,] = ( pnorm( (c(threshVec,Inf)-thisMu)/thisSigma )
                                 - pnorm( (c(-Inf,threshVec)-thisMu)/thisSigma ) )
        }
        points( 1:length(yColNames) , 
                apply(predProb,2,median) * sum(y[caseIdx,])  , 
                cex=1.5 , lwd=2 , col="skyblue" )
        for ( respIdx in 1:length(yColNames) ) {
          lines( rep(respIdx,2) , HDIofMCMC( predProb[,respIdx] ) * sum(y[caseIdx,]) , 
                 col="skyblue" , lwd=3 )
        }
        # Posterior of mean for this case:
        plotPost( OrdMcmcMat[,paste0("mu[",caseIdx,"]")] , cenTend="median" ,
                  xlab="Latent Scale" , main=bquote(mu[.(caseIdx)]) ,
                  cex.main=1.75 )
        # Posterior of sigma for this case:
        plotPost( OrdMcmcMat[,paste0("sigma[",caseIdx,"]")] , cenTend="median" ,
                  xlab="Latent Scale" , main=bquote(sigma[.(caseIdx)]) ,
                  cex.main=1.75 )
        if ( !is.null(fileNameRoot) ) {
          saveGraph( file=paste0(fileNameRoot,"-Case-",caseIdx) , 
                     type=graphFileType )
        }
      }
    }
  }

  #-----------------------------------------------
  # Scatter plot of  METRIC mu vs ORDERED-PROBIT mu:
  medOrdMuVec = apply( OrdMcmcMat[,grep("mu",colnames(OrdMcmcMat))] , 2 , median )
  medMetMuVec = apply( MetMcmcMat[,grep("mu",colnames(MetMcmcMat))] , 2 , median )
  nVec = apply( y , 1 , FUN=function(x){sum(x)} )
  # compute HDIs of the mu's:
  OrdMuHDI = matrix( 0 , nrow=length(grep("mu",colnames(OrdMcmcMat))) , ncol=2 )
  rowIdx = 1
  for ( colName in grep("mu",colnames(OrdMcmcMat)) ) {
    OrdMuHDI[rowIdx,] = HDIofMCMC(OrdMcmcMat[,colName])
    rowIdx = rowIdx+1
  }
  MetMuHDI = matrix( 0 , nrow=length(grep("mu",colnames(MetMcmcMat))) , ncol=2 )
  rowIdx = 1
  for ( colName in grep("mu",colnames(MetMcmcMat)) ) {
    MetMuHDI[rowIdx,] = HDIofMCMC(MetMcmcMat[,colName])
    rowIdx = rowIdx+1
  }
  
  # Plot  METRIC mu w HDI against ORDERED-PROBIT mu w HDI:
  openGraph(width=7,height=7)
  par( mar=c(3.5,3.5,2,1) , mgp=c(2.0,0.7,0) )
  xLim = range(OrdMuHDI)
  plot( x=medOrdMuVec , y=medMetMuVec , 
        xlim=xLim , ylim=range(MetMuHDI) ,  
        type="p" , # col=1:length(medMetMuVec) ,
        xlab="Mu of Ordered Probit model" , ylab="Mu of Metric model" ,
        cex.lab=1.5 )
  # Super-impose S-shaped curves for different ordered-probit sigma's:
  # Define function for computing bin probabilities:
  binProb = function( threshVec , m , s ) {
    phigh = pnorm( c( threshVec,Inf) , mean=m , sd=s )
    plow  = pnorm( c(-Inf,threshVec) , mean=m , sd=s )
    pbin = phigh - plow
    return( pbin )
  }
  medianThresh = apply( OrdMcmcMat[,grep("thresh",colnames(OrdMcmcMat))] , 2 , median )
  medianSigma = apply( OrdMcmcMat[,grep("sigma\\[",colnames(OrdMcmcMat))] , 2 , median )
  muOrdComb = seq( min(xLim) , max(xLim) , length=501 )
  muMet = 0*muOrdComb
  for ( sdVal in seq(min(medianSigma),max(medianSigma),length=2) ) {
    for ( muIdx in 1:length(muOrdComb) ) {
      muOrdVal = muOrdComb[muIdx]
      sigmaVal = sdVal
      pOrd = binProb( medianThresh , muOrdVal , sigmaVal )
      muMet[muIdx] = sum( pOrd * 1:length(pOrd) )
    }
    lines( muOrdComb , muMet , lwd=3 , col="gray" ) 
  }
  # Grid lines:
  abline( h=seq(1,5,by=0.5) , lty="dotted" )
  abline( v=seq(-4,12,by=2.0) , lty="dotted" )
  # Superimpose HDI's:
  for ( muIdx in 1:length(medMetMuVec) ) {
    lines( OrdMuHDI[muIdx,] , rep(medMetMuVec[muIdx],2) , col="skyblue" , lwd=2 ) 
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  for ( muIdx in 1:length(medOrdMuVec) ) {
    lines( rep(medOrdMuVec[muIdx],2) , MetMuHDI[muIdx,] , col="skyblue" , lwd=2 )
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  # replot points with Case Index:
  points( medOrdMuVec , medMetMuVec ) # , col=1:length(medMetMuVec) )
  for ( muIdx in 1:length(medMetMuVec) ) {
    text( medOrdMuVec[muIdx]+0.1 , medMetMuVec[muIdx] , labels=muIdx , 
          adj=c(0.0,0.5) , cex=0.67 )
  }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-MetMuByOrdMu"), type=graphFileType)
  }

  # Plot  METRIC mu w HDI against ORDERED-PROBIT mu w HDI, 
  # *zoomed out* ("dolly back"):
  openGraph(width=7,height=7)
  par( mar=c(3.5,3.5,2,1) , mgp=c(2.0,0.7,0) )
  xLim = c( min( c( 0 , min(OrdMuHDI) ) ) ,
            max( c( length(thresh)+2 , max(OrdMuHDI) ) ) )
  plot( x=medOrdMuVec , y=medMetMuVec , 
        xlim=xLim , ylim=c(1,ncol(y)) , # max(MetMuHDI) ,  
        type="p" , # col=1:length(medMetMuVec) ,
        xlab="Mu of Ordered Probit model" , ylab="Mu of Metric model" ,
        cex.lab=1.5 )
  # Super-impose S-shaped curves for different ordered-probit sigma's:
  # # Define function for computing bin probabilities:
  # binProb = function( threshVec , m , s ) {
  #   phigh = pnorm( c( threshVec,Inf) , mean=m , sd=s )
  #   plow  = pnorm( c(-Inf,threshVec) , mean=m , sd=s )
  #   pbin = phigh - plow
  #   return( pbin )
  # }
  medianThresh = apply( OrdMcmcMat[,grep("thresh",colnames(OrdMcmcMat))] , 2 , median )
  medianSigma = apply( OrdMcmcMat[,grep("sigma\\[",colnames(OrdMcmcMat))] , 2 , median )
  muOrdComb = seq( min(xLim) , max(xLim) , length=501 )
  muMet = 0*muOrdComb
  for ( sdVal in seq(min(medianSigma),max(medianSigma),length=2) ) {
    for ( muIdx in 1:length(muOrdComb) ) {
      muOrdVal = muOrdComb[muIdx]
      sigmaVal = sdVal
      pOrd = binProb( medianThresh , muOrdVal , sigmaVal )
      muMet[muIdx] = sum( pOrd * 1:length(pOrd) )
    }
    lines( muOrdComb , muMet , lwd=3 , col="gray" ) 
  }
  # # Grid lines:
  # abline( h=seq(1,5,by=0.5) , lty="dotted" )
  # abline( v=seq(-4,12,by=2.0) , lty="dotted" )
  # Superimpose HDI's:
  for ( muIdx in 1:length(medMetMuVec) ) {
    lines( OrdMuHDI[muIdx,] , rep(medMetMuVec[muIdx],2) , col="skyblue" , lwd=2 ) 
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  for ( muIdx in 1:length(medOrdMuVec) ) {
    lines( rep(medOrdMuVec[muIdx],2) , MetMuHDI[muIdx,] , col="skyblue" , lwd=2 )
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  # # replot points with Case Index:
  points( medOrdMuVec , medMetMuVec ) # , col=1:length(medMetMuVec) )
  # for ( muIdx in 1:length(medMetMuVec) ) {
  #   text( medOrdMuVec[muIdx]+0.1 , medMetMuVec[muIdx] , labels=muIdx , 
  #         adj=c(0.0,0.5) , cex=0.67 )
  # }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-MetMuByOrdMu-DollyBack"), type=graphFileType)
  }

  #-----------------------------------------------
  # Scatter plot of Metric-model Sigmas and function of ordered-probit mu's, 
  # with curves for constant ordered-probit sigma's
  medMetSigmaVec = apply( MetMcmcMat[,grep("sigma\\[",colnames(MetMcmcMat))] , 2 , median )
  # compute HDIs of the sigma's:
  MetSigmaHDI = matrix( 0 , nrow=length(grep("sigma\\[",colnames(MetMcmcMat))) , ncol=2 )
  rowIdx = 1
  for ( colName in grep("sigma\\[",colnames(MetMcmcMat)) ) {
    MetSigmaHDI[rowIdx,] = HDIofMCMC(MetMcmcMat[,colName])
    rowIdx = rowIdx+1
  }
  
  # Plot  METRIC sigma w HDI against ORDERED-PROBIT mu w HDI:
  openGraph(width=7,height=7)
  par( mar=c(3.5,3.5,2,1) , mgp=c(2.0,0.7,0) )
  xLim = range(OrdMuHDI)
  plot( x=medOrdMuVec , y=medMetSigmaVec , 
        xlim=xLim , ylim=range(MetSigmaHDI) ,  
        type="p" , # col=1:length(medMetMuVec) ,
        xlab="Mu of Ordered Probit model" , ylab="Sigma of Metric model" ,
        cex.lab=1.5 )
  # Super-impose curves for different ordered-probit sigma's.
  # Recall medianSigma, computed earlier, holds ordered-probit sigma's,
  # and muOrdComb is comb for horizontal axis of ordered-probit mu.
  #muOrdComb = seq( -max(OrdMuHDI) , max(OrdMuHDI) , length=501 )
  muOrdComb = seq( min(xLim) , max(xLim) , length=501 )
  sdOrdAsMetric = 0*muOrdComb
  for ( sdVal in seq(min(medianSigma),max(medianSigma),length=2) ) {
    for ( muIdx in 1:length(muOrdComb) ) {
      muOrdVal = muOrdComb[muIdx]
      sigmaVal = sdVal
      pOrd = binProb( medianThresh , muOrdVal , sigmaVal )
      muOrdAsMet = sum( pOrd * 1:length(pOrd) )
      sdOrdAsMetric[muIdx] = sqrt( sum( pOrd*( 1:length(pOrd) - muOrdAsMet )^2 ) )
    }
    lines( muOrdComb , sdOrdAsMetric , lwd=3 , col="gray" ) 
  }
  # Grid lines:
  abline( h=seq(0,5,by=0.5) , lty="dotted" )
  abline( v=seq(-4,12,by=2.0) , lty="dotted" )
  # Superimpose HDI's:
  for ( sdIdx in 1:length(medMetSigmaVec) ) {
    lines( OrdMuHDI[sdIdx,] , rep(medMetSigmaVec[sdIdx],2) , col="skyblue" , lwd=2 ) 
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  for ( muIdx in 1:length(medOrdMuVec) ) {
    lines( rep(medOrdMuVec[muIdx],2) , MetSigmaHDI[muIdx,] , col="skyblue" , lwd=2 )
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  # replot points with Case Index:
  points( medOrdMuVec , medMetSigmaVec ) # , col=1:length(medMetMuVec) )
  for ( muIdx in 1:length(medMetMuVec) ) {
    text( medOrdMuVec[muIdx]+0.1 , medMetSigmaVec[muIdx] , labels=muIdx , 
          adj=c(-0.0,0.5) , cex=0.67 )
  }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-MetSigmaByOrdMu"), type=graphFileType)
  }
  
  # Plot  METRIC sigma w HDI against ORDERED-PROBIT mu w HDI:
  # *zoomed out* ("dolly back"):
  openGraph(width=7,height=7)
  par( mar=c(3.5,3.5,2,1) , mgp=c(2.0,0.7,0) )
  xLim = c( min( c( 0 , min(OrdMuHDI) ) ) ,
            max( c( length(thresh)+2 , max(OrdMuHDI) ) ) )
  plot( x=medOrdMuVec , y=medMetSigmaVec , 
        xlim=xLim , ylim=c(0,1.1*max(MetSigmaHDI)) ,  
        type="p" , # col=1:length(medMetMuVec) ,
        xlab="Mu of Ordered Probit model" , ylab="Sigma of Metric model" ,
        cex.lab=1.5 )
  # Super-impose curves for different ordered-probit sigma's.
  # Recall medianSigma, computed earlier, holds ordered-probit sigma's,
  # and muOrdComb is comb for horizontal axis of ordered-probit mu.
  muOrdComb = seq( min(xLim) , max(xLim) , length=501 )
  sdOrdAsMetric = 0*muOrdComb
  for ( sdVal in seq(min(medianSigma),max(medianSigma),length=2) ) {
    for ( muIdx in 1:length(muOrdComb) ) {
      muOrdVal = muOrdComb[muIdx]
      sigmaVal = sdVal
      pOrd = binProb( medianThresh , muOrdVal , sigmaVal )
      muOrdAsMet = sum( pOrd * 1:length(pOrd) )
      sdOrdAsMetric[muIdx] = sqrt( sum( pOrd*( 1:length(pOrd) - muOrdAsMet )^2 ) )
    }
    lines( muOrdComb , sdOrdAsMetric , lwd=3 , col="gray" ) 
  }
  # # Grid lines:
  # abline( h=seq(0,5,by=0.5) , lty="dotted" )
  # abline( v=seq(-4,12,by=2.0) , lty="dotted" )
  # Superimpose HDI's:
  for ( sdIdx in 1:length(medMetSigmaVec) ) {
    lines( OrdMuHDI[sdIdx,] , rep(medMetSigmaVec[sdIdx],2) , col="skyblue" , lwd=2 ) 
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  for ( muIdx in 1:length(medOrdMuVec) ) {
    lines( rep(medOrdMuVec[muIdx],2) , MetSigmaHDI[muIdx,] , col="skyblue" , lwd=2 )
    #, col=muIdx, lwd=1.5*log10(z[muIdx]) )
  }
  # replot points with Case Index:
  points( medOrdMuVec , medMetSigmaVec ) # , col=1:length(medMetMuVec) )
  # for ( muIdx in 1:length(medMetMuVec) ) {
  #   text( medOrdMuVec[muIdx]+0.1 , medMetSigmaVec[muIdx] , labels=muIdx , 
  #         adj=c(-0.0,0.5) , cex=0.67 )
  # }
  if ( !is.null(fileNameRoot) ) {
    saveGraph( file=paste0(fileNameRoot,"-MetSigmaByOrdMu-DollyBack"), type=graphFileType)
  }

  #------------------------------------------------------------------------------- 
  if ( !is.null(fileNameRoot) ) {
    write.csv( OrdSummaryMat , 
               file=paste0(fileNameRoot,"-OrdModel-ParameterSummary.csv") ,
               row.names=TRUE )
    write.csv( OrdMcmcMat , 
               file=paste0(fileNameRoot,"-OrdModel-McmcMatrix.csv") ,
               row.names=TRUE )
  }
  return( list( OrdSummaryMat=OrdSummaryMat , OrdMcmcMat=OrdMcmcMat ) )
} # end of function definition
