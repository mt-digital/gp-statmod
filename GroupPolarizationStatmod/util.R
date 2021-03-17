##
# The main display element is one row of controls and two rows of plots. 
#
# The top row will be the large N exact calculations of the 
# latent standard deviations given observed pre- and post-deliberation
# means, a hypothesized latent mean, and a reporting model and measurement 
# scheme. The latent mean and observed pre- and post-deliberation means
# can be set using the top row controls, or loaded from preset values
# from published result case studies that contain plausibly
# false demonstrations of the group polarization effect.
#
# The second row displays simulated empirical observations for an N given in 
# the controls at the top of the main display. The numbers and plots 
# are generated using the same latent and observed 
# pre- and post-deliberation means as for the large N exact calculations.
#
# The "reporting model" is operationalized as thresholds that determine
# which opinion scale (e.g. Likert) bin a particpant's opinion is reported as.
# See the code for how this works. The rating scale would be, e.g., a
# 7-point Likert, or the CDQ responses, which rate riskiness from 1 to 10.
#  
# getCaseStudyData <- function() {
#     read.csv('caseStudies.csv') 
#     data.names <- data$CaseStudyName
#     return (data)
# }
# CaseStudyData <- getCaseStudyData();

library(shiny)

source("model.R")
source("numerical.R")

CaseStudyData <- read.csv('caseStudies.csv', row.names='CaseStudyName')


##
# Arguments:
#     caseStudy (string): See CSV for values to use, e.g., 
#         Moscovici1969-DeGaulle.
setCaseStudyValues <- function(caseStudyName, session)
{
    data <- CaseStudyData[caseStudyName, ]
    updateSliderInput(session, "latentMean", 
                      value = data$LatentMean, min = data$MinBinValue - 2,
                      max = data$MaxBinValue + 2)
                      
    updateSliderInput(session, "observedPreDelibMean", 
                      value = data$ObservedPreDelibMean, 
                      min = data$MinBinValue,
                      max = data$MaxBinValue)

    updateSliderInput(session, "observedPostDelibMean", 
                      value = data$ObservedPostDelibMean, 
                      min = data$MinBinValue,
                      max = data$MaxBinValue)
}
getCaseStudyValues <- function(caseStudyName)
{
    return (CaseStudyData[caseStudyName, ])
}


##
# Arguments:
#     input (Shiny input)
#     time (string): Either "pre" or "post" for pre- or post-deliberation.
#
largeNBarplot <- function(input, time)
{
    kVec <- input$minBinValue:input$maxBinValue
    
    if (time == "pre")
        observedMean <- input$observedPreDelibMean
    else
        observedMean <- input$observedPostDelibMean

    # Use a guess of 1.5 for latentSD.
    latentSD = solveForLatentSD(kVec, input$latentMean, observedMean, 1.5)

    barplot(makeProbVec(kVec, input$latentMean, latentSD), 
            # names.arg = kVec, 
            # ylab = "Frequency of response",
            main = paste(
                "Calculated ", time, "-discussion SD = ", latentSD, sep=""
            )
    )
}

