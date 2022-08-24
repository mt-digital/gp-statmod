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
library(grid)
library(gridExtra)
library(lattice)
library(ggplot2)

source("model.R")
source("numerical.R")

# CaseStudyData <- read.csv('caseStudies.csv', row.names='CaseStudyTag')


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
largeNBarplots <- function(input, time)
{
    kVec <- input$minBinValue:input$maxBinValue
    
    observedPreMean <- input$observedPreDelibMean
    observedPostMean <- input$observedPostDelibMean

    # Use a guess of 1.5 for latentSD.
    latentPreSD = solveForLatentSD(kVec, input$latentMean, observedPreMean, 1.5)
    latentPostSD = solveForLatentSD(kVec, input$latentMean, observedPostMean, 1.5)

    preBarplot <- barplot(makeProbVec(kVec, input$latentMean, latentSD), 
                # names.arg = kVec, 
                # ylab = "Frequency of response",
                main = paste(
                    "Calculated ", time, "-discussion SD = ", latentPreSD, sep=""
                ))

    postBarplot <- barplot(makeProbVec(kVec, input$latentMean, latentPostSD), 
                # names.arg = kVec, 
                # ylab = "Frequency of response",
                main = paste(
                    "Calculated ", time, "-discussion SD = ", latentSD, sep=""
                ))

    grid.arrange(grobs=c(preBarplot, postBarplot), ncol = 2) 
}


##
# Friedkin (1999) seems to convert decimal probabilities given in response
# to CDQ questions on Sports, School, and Surgery (in three group-size
# conditions) to percentage probabilities in Table 1 (p. 868). Friedkin says,
# "The subjects' responses (opinions) were restricted to one of 20 probability
# values: .05, .10, .15, ..., 1.00." Yet in Table 1, average choice shifts
# are given that are much greater than 1, and initial and final mean opinions
# have values between 0 and 100. I see no explanation of the difference.
#
# To use Friedkin's reported measurements I need to scale responses apparently
# from the 5 to 100 percentage scale to a 1-to-20 ordinal scale. The results
# from this calculation can be used in the web app for hypothesized latent
# mean and observed pre- and post-deliberation means.
convertFriedkinTo20Bin <- function(friedkinVal)
{
    return ( 20 * (friedkinVal / 100.0) )
}

pooled.sd <- function(x1, x2)
{
  #calculate sample size of each group
  n1 <- length(x1)
  n2 <- length(x2)
  
  #calculate sample variance of each group
  var1 <- var(x1)
  var2 <- var(x2)
  
  #calculate pooled variance between the two groups
  pooled <- ((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2)
  
  #display pooled variance
  return (pooled)
}