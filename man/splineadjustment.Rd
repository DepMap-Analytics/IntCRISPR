\name{splineadjustment}
\alias{splineadjustment}
\title{Use spline adjustment to correct for differences in screen quality}
\usage{
splineadjustment(data,quantiles=c(0.25,0.5,0.75))
}
\arguments{
  \item{Data}{Data to adjust for screen quality. Matrix, rows genes columns screens.}
  \item{quantile}{Optional vector of quantiles at which to fit knots for spline. Default to 0.2,0.5 and 0.75}


}
\description{
Adjusts for differences in screen quality.
}
\examples{
#Not run
#BroadScreen<-splineadjustment(BroadOverlap,quantiles=c(0.25,0.5,0.75))
}
