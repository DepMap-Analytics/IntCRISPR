\name{ComBatCP}
\alias{ComBatCP}
\title{Use ComBat to get batch correction adjustments.}
\usage{
ComBatCP(dat, batch, mod, par.prior,
                      mean.only, ref.batch, BPPARAM,empBayes=TRUE)
}
\arguments{
  \item{dat}{Matrix of all data to be batch corrected. Rows should be genes columns screens. Recommended data is quantile normalised before passing to this function.}
  \item{batch}{Vector of the batch of origin for each of the screens, should be in the same order as the column names of the Data matrix.}
  \item{par.prior}{See sva package.}
  \item{mod}{Optional. Additional covariate to be fitted by ComBat.}
  \item{mean.only}{Default FALSE. Should ComBat only correct the means across batches (excluding variances)}
  \item{ref.batch}{Optional, default NULL. A level of the batch vector to be used as the reference level of the batch factor in the regression model}
  \item{BBPARAM}{See sva package.}
  \item{empBayes}{Default TRUE. Should an empirical Bayes adjustment be used to calculate the final batch correction vectors}


}
\description{
Calculates the correction vectors based on the ComBat function from the R sva package for removing batch effects. Two datasets of the same models (e.g. screens of cell lines) from two different batches are input and the batch correction vectors returned. Can be used to calculate batch effects from a subset of screens e.g. those cell lines screened in all batches. The correction vectors output from this function can then be applied to non-overlapping screens.
}
\examples{
#Not run
#site<-sapply(colnames(mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])

#CombatOrig<-ComBatCP(as.matrix(mergedData),batch = site[colnames(mergedData)])

}
