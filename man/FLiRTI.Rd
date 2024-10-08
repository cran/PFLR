% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLiRTI.R
\name{FLiRTI}
\alias{FLiRTI}
\title{FLiRTI Regression Model}
\usage{
FLiRTI(
  Y,
  X,
  d,
  cons,
  domain,
  extra = list(Mf = 6:30, lambda = seq(5e-04, 100, length.out = 50))
)
}
\arguments{
\item{Y}{Vector of length n, centred response.}

\item{X}{Matrix of n x p, covariate matrix, should be dense.}

\item{d}{Integer, degree of the B-spline basis functions.}

\item{cons}{Divide subinterval into how many small ones.}

\item{domain}{The range over which the function X(t) is evaluated and the coefficient function \eqn{\beta}(t) is expanded by the B-spline basis functions.}

\item{extra}{List containing parameters which have default values:
\itemize{
  \item Mf: Mf+1 is the number of knots for the B-spline basis functions that expand \eqn{\beta}(t), default is 6:30.
  \item lambda: Tuning parameter, default is seq(0.0005,100,length.out = 50).
}}
}
\value{
beta: Estimated \eqn{\beta}(t) at discrete points.

extra: List containing other values which may be of use:
       \itemize{
         \item X: Matrix of n x p used for model.
         \item Y: Vector of length n used for model.
         \item domain: The range over which the function X(t) was evaluated and the coefficient function \eqn{\beta}(t) was expanded by the B-spline basis functions.
         \item delta: Estimated cutoff point.
         \item OptM: Optimal number of B-spline knots selected by BIC.
         \item Optlambda: Optimal shrinkage parameter selected by BIC.
       }
}
\description{
Calculates functional regression that's interpretable using the FLiRTI method.
}
\examples{
library(fda)
betaind = 1
snr  = 2
nsim = 200
n    = 50
p    = 21
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)
lambda = seq(0.0005,0.01,length.out = 10)
Mf = 6:13
extra=list(Mf=Mf,lambda=lambda)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=1)
 Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}


fltyfit = FLiRTI(Y=Y[1:n,1],(X[1:n,,1]),d=3,cons=4,domain=domain,extra=extra)


}
