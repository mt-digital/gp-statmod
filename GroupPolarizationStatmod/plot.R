##
# Some methods for plotting results from statmods.
#

require(fs)
require(tidyverse)

source("experiments.R")


LIKERT_SCALE_DENS_MAX <- 0.7

##
# Plot the frequentist analysis of difference between two groups.
#
# Arguments:
#     inputDf (data.frame): the simulated data to model
#     fittedFreq (list): list of fitted frequency model parameters and summary
#     nBins (int): number of Likert-type bins used for opinion binning
#
plotFreq <- function(inputDf, fittedFreq, nBins)
{

    # Open window for plot:
    openGraph(width=3.5*2, height=2.5*1+0.75)
    layout(matrix(1:2, nrow=1))
    par(mar=c(3.5,3.5,2.5,0.5), mgp=c(2.0,0.7,0), oma=c(0,0,3.5,0)) # , xpd=NA )

    # histBreakInc <- 1/nQ
    # Some that seem to be defaults.
    histBreakInc <- 1
    xlabText="Ordinal Resp." 

    # make histograms of data:
    mainTitle <- c("Pre-discussion", "Post-discussion")
    for ( gIdx in 1:2 ) 
    {
        ydata <- as.numeric(inputDf$resp[inputDf$cond == gIdx])
        histInfo <- hist( 
            ydata, main=mainTitle[gIdx], 
            xlab=xlabText, ylab="Probability",
            freq=FALSE, col="pink", border="white", cex.lab=1.5,
            breaks=seq(-histBreakInc/2, 
                       nBins+1+histBreakInc/2, 
                       by=histBreakInc), 
            ylim=c(0, LIKERT_SCALE_DENS_MAX) # NB: must set manually at this time
        )
        #text( max(histInfo$breaks) , likertScaleDensMax, adj=c(1,1),
        #      labels=bquote(N==.(length(ydata))) )
        # xcomb <- seq(0,max(nYlevels)+1,length=201)
        xcomb <- seq(0, nBins + 1, length=201)
        lines(xcomb, dnorm(xcomb-mean(ydata), 0, sd(ydata)), 
              col="skyblue", lwd=2)
    }
    
    tInfo <- fittedFreq$tInfo
    effSz <- fittedFreq$effSz
    if (tInfo$p.value >= 0.001) {
      mtext(text=bquote( "d = "* .(round(effSz, 2)) 
                               # *", t = "* .(round(tInfo$statistic, 2)) 
                               *", p = "* .(round(tInfo$p.value, 3))) , 
             at=c(0.5), cex=1.5, side=3, outer=TRUE, 
             adj=c(0.5,0.5), padj=c(-0.5,-0.5))
    } else {
      mtext(text=bquote( "d = "* .(round(effSz, 2)) 
                               # *", t = "* .(round(tInfo$statistic, 2)) 
                               *", p < 0.001" ), 
             at=c(0.5), cex=1.5, side=3, outer=TRUE, 
             adj=c(0.5,0.5), padj=c(-0.5,-0.5))
    }
}


plot_metric_cohens <- function(metric_cohens_csv = "data/output/metric_cohens_d.csv",
                               output_file = "figures/metric_cohens.pdf") {
  
  read_csv(metric_cohens_csv) %>% 
    unite(ExperimentID, ArticleTag, TreatmentTag) %>% 
    ggplot(mapping = aes(x = Cohens_d, y = ExperimentID)) + 
      geom_vline(xintercept = 0, color = "red") + 
      geom_vline(xintercept = c(-1, 1), color="red", linetype="dotted") + 
      geom_boxplot() +
      xlim(-2.0, 2.0)
  
  ggsave(output_file)
}


plot_ordinal_cohens <- function(ordinal_data_dir = "data/probit_fits",
                                output_file = "figures/ordinal_cohens.pdf") {
  
  df <- dir_ls(ordinal_data_dir, glob = "*.csv") %>%
    read_csv() %>% unite(ExperimentID, ArticleTag, TreatmentTag)
  
  df$Cohens_d <- cohens_d(df$LatentMeanPrePosteriorMean, df$LatentMeanPostPosteriorMean,
                          df$LatentMeanPrePosteriorSD, df$LatentMeanPostPosteriorSD)
  
  ggplot(df, mapping = aes(x = Cohens_d, y = ExperimentID)) + 
    geom_vline(xintercept = 0, color = "red") + 
    geom_vline(xintercept = c(-1, 1), color="red", linetype="dotted") + 
    geom_boxplot() #+
    # xlim(-2.0, 2.0)
  
  ggsave(output_file)
}