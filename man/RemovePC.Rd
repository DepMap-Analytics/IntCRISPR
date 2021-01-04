\name{RemovePC}
\alias{RemovePC}
\title{Remove principal components from data}
\usage{
RemovePC(data,droppcanumber=1)
}
\arguments{
  \item{data}{Matrix of data to remove principal components from. Rows should be genes columns screens.}
    \item{droppcanumber}{Vector of the indices of the principal components to remove. E.g. to remove the first two principal components droppcanumber=1:2}

}
\description{
Removes the specified number of principal components from the data
}
\examples{

}
