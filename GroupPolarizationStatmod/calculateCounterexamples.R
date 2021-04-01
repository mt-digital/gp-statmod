##
#
# Code to calculate the latent and simulated observed pre- and post-
# deliberation standard deviations (4 values per row of input data).
#
# Author: Matthew A. Turner
# Date: 2021-03-29
#


source("experiments.R")
source("model.R")
source("numerical.R")


processData <- function(fileLoc = "data/CaseStudiesMinimal.csv")
{
    # Each case study data row is one condition of one study in the form
    # <FirstAuthorName>-<ConditionDescription>, e.g., Schkade-COSprings-GlobWarm.
    d <- read.csv(fileLoc, row.names="CaseStudyTag")
    # print(d)

    # Initialize numerical results data columns.
    for (newCol in c("LatentSDPre", "LatentSDPost", "SimObsSDPre", "SimObsSDPost"))
    {
        d[, newCol] <- numeric(nrow(d))        
        # print(d[, newCol])
    }

    # Iterate over conditions where our false detection hypothesis is
    # plausible to calculate the pre- and post-deliberation latent 
    # standard deviations, and their simulated observed standard deviations.
    for (n in row.names(d))
    {
        # TODO: adapt Shiny code to calculate what we need

        # Extract the row of interest.
        r <- d[n, ]

        if (as.logical(r$Plausible))
        {
            # Need to create vector of bins which need to be renamed in the
            # solveForLatentSD signatures, in which this is called "kVec".
            opinionBins <- as.numeric(r$MinBinValue):as.numeric(r$MaxBinValue)

            # This SD seemed to give good results -- TODO: need to check about
            # whether an alternative results in a false detection find in
            # the Myers and Bishop study.
            sdGuess <- 1.5 

            # In the cases where 
            latentMeanGuess <- r$ObservedMeanPost

            # Calculate latent and simulated observed pre- and post-discussion
            # standard deviations of opinion.
            #
            # Although solutions exist, the hillclimbing algorithm currently
            # uses randomized step values, and so sometimes fails to find a
            # solution.
            tol <- 1e-3
            converged <- FALSE
            while (!converged)
            {
                latSDPreResult <- solveForLatentSD(opinionBins, 
                                                   latentMeanGuess, 
                                                   r$ObservedMeanPre, 
                                                   1.5)

                latSDPostResult <- solveForLatentSD(opinionBins, 
                                                    latentMeanGuess, 
                                                    r$ObservedMeanPost, 
                                                    1.5)

                converged <- hillclimbSuccess(latSDPreResult, latSDPostResult)
            }

            latSDPre <- d[n, "LatentSDPre"] <- latSDPreResult[1]
            latSDPost <- d[n, "LatentSDPost"] <- latSDPostResult[1]
            
            d[n, "SimObsSDPre"] <- sdObs(opinionBins, 
                                         makeProbVec(opinionBins, 
                                                     latentMeanGuess,
                                                     latSDPre)
                                         )

            d[n, "SimObsSDPost"] <- sdObs(opinionBins, 
                                          makeProbVec(opinionBins, 
                                                     latentMeanGuess,
                                                     latSDPost)
                                         )
        }
    }

    write.csv(d, "data/CaseStudiesMinimalProcessed.csv")


    return (d);
}
