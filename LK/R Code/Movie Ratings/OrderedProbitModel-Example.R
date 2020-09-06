# OrderedProbitModel-Example.R
# John K. Kruschke, January-August 2018.
# Sept 28, 2018: Revised by John K. Kruschke; uses new argument caseIDcolName.
#
# This is an example of using the function ordinalAndMetricAnalysis, defined in
# the R script OrderedProbitModel.R. It analyzes ordinal data with ordered probit
# and (inappropriately) metric models.
#
# These R scripts accompany the article, "Analyzing ordinal data with metric
# models: What could possibly go wrong?" by Torrin M. Liddell and John K.
# Kruschke. Files can be found at https://osf.io/53ce9/
# 
# To use this script, you must have in R's working directory the following three
# files:
# (1) OrderedProbitModel.R
# (2) DBDA2E-utilities.R
# (3) the data file, which in this example is MoviesData.csv. The data file
# must be in comma-separated-value (CSV) format, with each row having the
# frequency counts of the levels of a single group or case. The columns must be
# named, but the names can be whatever you want. For example, the top of the
# MovieData.csv file looks like this:
#
#   "ID","Descrip","n1","n2","n3","n4","n5"
#   1,"The Whole Truth",49,70,119,217,245
#   2,"Priceless",67,22,22,60,574
#
# The first row must be the column names. Each subsequent row is for a
# particular movie. This file happens to have movie ID and Description columns,
# but these are not used by the script. The columns named "n1" through "n5"
# contain the frequenices observed at each ordinal level.
#
# See important comments at the top of OrderedProbitModel.R for more information
# about the function.
#
#-----------------------------------------------------------------------------------
# Optional housekeeping:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#-----------------------------------------------------------------------------------

source("OrderedProbitModel.R")

# Example with the movies data.
OrdModelResults = ordinalAndMetricAnalysis( 
  dataFileName = "MoviesData.csv" ,
  yColNames = c("n1","n2","n3","n4","n5") , # column names in data file
  caseIDcolName = "Descrip" ,
  compareCases = list( c(34,4) , # a case of inverted means
                       c(5,6) , # another case of inverted means
                       c(2,25) , # cases of approx equal metric SDs 
                       c(35,29) , # but more power in ordered-probit...
                       c(25,29) , 
                       c(19,27) , 
                       c(34,11) 
  ) , # close parenthesis of compareCases list
  hierarchSD = FALSE
) # close parenthesis of function arguments

#...................................................................................

# Another example, using artificial data in the A,B,C,D configuration of 
# Figure 4 in the article.
OrdModelResults = ordinalAndMetricAnalysis( 
  dataFileName = "ABCDdata.csv" ,
  yColNames = c("n1","n2","n3","n4","n5") , # column names in data file
  caseIDcolName = "GroupID" ,
  compareCases = list( c(2,1) , # false alarm (Type I error)
                       c(4,2) , # miss (Type II error)
                       c(4,3)  # inversion!
  ) , # close parenthesis of compareCases list
  hierarchSD = FALSE
) # close parenthesis of function arguments

#...................................................................................

# The following lines are optional, presented as examples of how to do further
# exploration of the output.

# Display the first few lines of the parameter summary matrix:
OrdModelResults$OrdSummaryMat[1:6,]

# Display the top left of the MCMC chain matrix:
OrdModelResults$OrdMcmcMat[1:6,1:10]

# Plot the posterior distribution of a difference of means:
openGraph() # opens new graphics window
muDiff = OrdModelResults$OrdMcmcMat[,"mu[2]"] - OrdModelResults$OrdMcmcMat[,"mu[1]"]
plotPost( muDiff , xlab=bquote(mu[2]-mu[1]) , main="Difference of Means" )
          
# Plot the posterior distribution of Cohen's d for a difference of means:
openGraph() # opens new graphics window
muCohenD = ( ( OrdModelResults$OrdMcmcMat[,"mu[2]"] 
               - OrdModelResults$OrdMcmcMat[,"mu[1]"] )
             / sqrt( ( OrdModelResults$OrdMcmcMat[,"sigma[2]"]^2 
                       + OrdModelResults$OrdMcmcMat[,"sigma[1]"]^2 )/2 ) )
plotPost( muCohenD , xlab=bquote((mu[2]-mu[1])/sqrt((sigma[1]^2+sigma[2]^2)/2)) ,
          main="Effect Size" )
