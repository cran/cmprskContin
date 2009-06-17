\name{summaryCmprsk}
\alias{summaryCmprsk}
\title{summary of statistical tests }
\description{
Outputs to the screen a summary of the statistical tests described in Gilbert, McKeague, and Sun (2008). 
}
\usage{
summaryCmprsk(x)
}
\arguments{
  \item{x}{output list from \code{cmprsk} function }
}
\details{
This function accepts the output list from the \code{cmprsk} function and returns to the screen the recommended tests and corresponding test statistics form Gilbert, McKeague and Sun (2008).
}


\references{
Gilbert, P., McKeague, I. and Sun, Y. (2008) \emph{The 2-sample problem for 
failure rates depending on a continuous mark: an application to vaccine efficacy}. 
Biostatistics 9, 2, 263-276.\cr
Andersen, P. K., Borgan, O., Gill, R. D., and Keiding. (1993) \emph{Statistical Models Based on Counting Processes.} New York: Springer.

}

\author{Peter Gilbert \cr
Maintainer: Valerie Obenchain <vobencha@fhcrc.org>}
\seealso{
\code{cmprsk} and \code{plotCmprsk} functions
}
\examples{
## Simulate data
set.seed(100)
gp <- rep(1:2, times=50)
ftime <- 50*runif(100)
ftype <- rep(0:1, each=50)
# uninfected subjects do not have mark value; use 99
infected <- which(ftype==1)
mark <- rep(99, length(gp))
# create mark values for the infected
mark[infected] <- runif(length(infected))
# create some infected subjects with missing mark; use 99 
missing <- sample(infected, 2)
mark[missing] <- 99

## Compute risk
x <- cmprsk(gp,ftime,ftype,mark,nboot=500,T1=2,T2=36,UT1=1,
UT2=42,ttanal=36)

## Summary information
summaryCmprsk(x)

}
\keyword{survival}

