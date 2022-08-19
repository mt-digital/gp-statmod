makeJAGSModelData <- function(N, firstBinValue, nBins, latentMean, latentSd)  # Pre, latentSdPost)
{
    thetaVec = as.vector(matrix(NA, ncol=nBins - 1, nrow=1))

    opinions = simulatedObservation(N, firstBinValue, nBins, latentMean, latentSd)

    modelDataList = list(
        o = opinions,
        thetaVec = thetaVec,
        nBins = nBins
    )
        
    return (modelDataList)
}


# JAGS model definition 
ordSingleModelString = "

  model {

    for (i in 1:N) {

      o[i] ~ dcat(pr[i, 1:nBins])
      pr[i,1] <- pnorm(thetaVec[1], mu, 1/sigma^2)

      for (k in 2:(nBins - 1)) {
        pr[i, k] <- max(
          0, 
          pnorm(thetaVec[k], mu, 1 / (sigma^2)) - pnorm(thetaVec[k - 1], mu, 1 / (sigma^2)) 
        )
      }

      pr[i, nBins <- (1 - pnorm(thetaVec[nBins - 1], mu, 1/sigma^2))
    }

    mu ~ dnorm((1+max(nBins)) / 2, 1 / (max(nBins)^2))
    sigma ~ dunif(0.01 , max(nBins)*10)

    # Prior on thetaVec[q,k]. Stochastic for all except thetaVec[1,1] and thetaVec[1,last].
    for (kIdx in 2:(nBins - 2)) { # 1 and nBins-1 are fixed
      thetaVec[kIdx] ~ dnorm(kIdx + 0.5, 1/(2^2))
    }
}

" # close quote for ordModelString



#   # ASSEMBLE THE DATA FOR JAGS. 
#   #
#   # N.B. THE JAGS MODEL ASSUMES THAT ALL ITEMS ARE SCORED WITH POSITIVE
#   # INTERCORRELATIONS. 
  
# simDataForJAGS <- function(simData)
# {
#   # Rename and reclass y values for convenience:
#   y = as.numeric(simData)
#   # Do some checking that data make sense:
#   if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
#   if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
#   # COMPRESS OUT ANY EMPTY VALUES OF Y:
#   yOrig=y
#   y=as.numeric(factor(y, levels=names(table(y))))
#   if ( any(y != yOrig) ) { 
#     cat("************************************************\n")
#     CAT("** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS **\N")
#     CAT("************************************************\N")
#   }
#   NTOTAL = LENGTH(Y)
#   # QUESTION, SUBJECT, AND GROUP DATA VECTORS:
#   Q = C(1)
#   G = C("PRE", "POST")
#   QLEVELS = LEVELS(AS.FACTOR(Q))
#   # SLEVELS = LEVELS(AS.FACTOR(DATFRM[,SNAME])) # NOT USED
#   GLEVELS = LEVELS(AS.FACTOR(G))
#   NG = MAX(G)
#   # CREATE THRESHOLD MATRIX. FOR THE FIRST QUESTION, THRESHOLD 1 AND NYLEVELS-1 ARE 
#   # FIXED; OTHER INTERIOR THRESHOLDS ARE ESTIMATED. FOR OTHER QUESTIONS, ALL THRESHOLDS
#   # ARE ESTIMATED.
#   # ** THIS PRESENTLY ASSUMES THAT ALL ITEMS HAVE THE SAME NUMBER OF LEVELS. **
#   # ** IT PRESUMABLY WON'T WORK OTHERWISE WITHOUT MODIFICATION. **
#   # COMPUTE NUMBER OF Y LEVELS FOR EACH QUESTION:
#   NYLEVELS = AGGREGATE( Y , BY=LIST(Q) , FUN=MAX )$X
#   THRESH = MATRIX( NA , NROW=NQ , NCOL=MAX(NYLEVELS)-1 ) # DEFAULT TO NA
#   # FIX LOW THRESH OF ITEM 1 AT 1.5:
#   THRESH[1,1] = 1 + 0.5 
#   # FIX UPPER THRESH OF ITEM 1 AT K-0.5:
#   THRESH[1,NYLEVELS[1]-1] = NYLEVELS[1] - 0.5
#   # SPECIFY THE DATA IN A LIST, FOR ORDERED-PROBIT MODEL, FOR LATER SHIPMENT TO JAGS:
#   ORDDATALIST = LIST(
#     Y = Y ,
#     G = G ,
#     THRESH = THRESH ,
#     NYLEVELS = NYLEVELS ,
#     NQ = NQ ,
#     NG = NG ,
#     NTOTAL = NTOTAL 
#   )

#   RETURN (ORDDATALIST)
# }

