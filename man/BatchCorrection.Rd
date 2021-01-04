\name{BatchCorrection}
\alias{BatchCorrection}
\title{Batch correct data sets using a pre-calculated set of correction vectors}
\usage{
BatchCorrection(data1,data2,site1="Broad",site2="Sanger",CombatRes,stdPrior=TRUE,qcThresh=NULL,qcvalues1=NULL,qcvalues2=NULL)
}
\arguments{
  \item{data1}{Matrix of first data set to be batch corrected. Rows should be genes columns screens. Recommended data is quantile normalised before passing to this function.}
    \item{data2}{Matrix of second data set to be batch corrected. Rows should be genes columns screens. Recommended data is quantile normalised before passing to this function.}
  \item{site1}{Default Broad. Character abel for batch of data1}
  \item{site2}{Default Sanger. Character label for batch of data1}
  \item{CombatRes}{Output for the ComBatCP function.}
  \item{stdPrior}{Default TRUE. Should genes be standardised first using mean and variance values from the ComBatRes output.}
  \item{qcThresh}{Optional, default NULL. Minimum quality score to include screen in correction and final dataset}
  \item{qcvalues1}{Default NULL. QC values required when qcThresh is used. The vector should relate to the columns of data1 giving quality scores for each screen.}
  \item{qcvalues2}{Default NULL. QC values required when qcThresh is used. The vector should relate to the columns of data2 giving quality scores for each screen.}


}
\description{
Calculates the correction vectors based on the ComBat function from the R sva package for removing batch effects. Two datasets of the same models (e.g. screens of cell lines) from two different batches are input and the batch correction vectors returned. Can be used to calculate batch effects from a subset of screens e.g. those cell lines screened in all batches. The correction vectors output from this function can then be applied to non-overlapping screens.
}
\examples{
#Not run

#OrigAll<-BatchCorrection(data1=BroadData,data2=SangerData,CombatRes=CombatOrig)
}
