# John K. Kruschke, 2014-Nov-21.
# This function is described at the blog post
# http://doingbayesiandataanalysis.blogspot.com/2014/11/ordinal-probit-regression-transforming.html
# A PDF copy of the blog post is available at
# https://osf.io/fc6zd/
# This function comes with no warranties, guarantees, etc. Use at your own risk.
polrToOrdScale = function( polrObject ) {
  polrThresh = polrObject$zeta
  polrSlopes = polrObject$coefficients
  polrInter = 0.0
  polrSigma = 1.0
  K = length(polrThresh) + 1  # K is number of ordinal levels
  sigmaMult = unname( (K-2)/(polrThresh[K-1]-polrThresh[1]) )
  inter = unname( 1.5 - ( sigmaMult*polrThresh[1] ) )
  respThresh = sigmaMult*polrThresh + inter
  respSigma = sigmaMult*polrSigma
  respB0 = sigmaMult*polrInter + inter
  respSlopes = sigmaMult*polrSlopes
  return( list( sigma=respSigma ,
                b0=respB0 ,
                coefficients=respSlopes ,
                zeta=respThresh ) )
}
