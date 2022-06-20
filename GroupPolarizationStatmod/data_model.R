
draftStudiesAnalysisColumns <- 
  c("Article",	"StudyTag",	"ObservedMeanPre",	
    "ObservedMeanPost",	"ObservedShiftSD", "LatentMean",
    "HillclimbStepSize",	"HillclimbSuccessThreshold",	"MinBinValue",	
    "MaxBinValue",	"Plausible",	"PlausibleInt",	"Notes")

studiesAnalysisColumns <- c(draftStudiesAnalysisColumns, c("latPreSD", "latPostSD"))

articleColumns <- c("Authors", "Journal", "Year", "Tag", "StudyTags")

tags <- c("Schkade2010", "Moscovici1969", "Myers1970","Friedkin1999", 
          "Krizan2007", "Abrams1990", "Burnstein1973", "Hogg1990",  "Burnstein1975",
          "Burnstein1973")

authors <- c("Abrams, Wetherell, Cochrane, Hogg, and Turner", 
             "Burnstein and Vinokur", 
             "Burnstein and Vinokur", 
             "Friedkin", 
             "Hogg, Turner, and Davidson", 
             "Krizan and Baron", 
             "Moscovici and Zavalloni", 
             "Myers", 
             "Myers and Bishop", 
             "Schkade, Sunstein, and Hastie")

year <- c(1990, 1973, 1975, 1999, 1990, 2007, 1969, 1975, 1970, 2010)

tags <- c("Abrams1990", "Burnstein1973", "Burnstein1975", "")



journal <- c("")