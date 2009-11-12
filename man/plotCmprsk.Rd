\name{plotCmprsk}
\alias{plotCmprsk}
\title{plot estimated vaccine efficacy or relative risk}
\description{
This function plots the cumulative or doubly cumulative estimated vaccine efficacy 
\eqn{VE(ttanal,v)} or relative risk \eqn{RR(ttanal,v)} versus the mark \eqn{v}.
95\% confidence interval estimates versus the mark \eqn{v} are also plotted. 
All plots are output to .ps files.
}
\usage{
plotCmprsk(x, plottype=1, dc=TRUE, p1, p2, filename, main=" ",
ylim=c(-1,1.5), xlim, legloc=1, xlab="Mark", 
ylab="Estimated Risk", ...) 
}
\arguments{
  \item{x}{List output from \code{cmprsk} function}
  \item{plottype}{Binary 0 or 1 to indicate plot type. 1=vaccine efficacy (default), 0=relative risk} 
  \item{dc}{Logical to indicate VE type. TRUE=doubly cumulative VE (default),\cr
 FALSE=cumulative VE} 
  \item{p1}{P-value for testing \eqn{H^0_0:VE(t,v)=0} or \eqn{H^0_0:RR(t,v)=1} for all times \eqn{(t)} and marks \eqn{(v)}. Inputs can be \code{pval11}, \code{pval12}, \code{pval13}, or \code{pval14} from the \code{STATvec} component of the \code{cmprsk} output list. These p-values are for one-sided hypothesis tests. The default is \code{pval13}} 
  \item{p2}{P-value for testing \eqn{H^0:VE(t,v)=VE(t)} or \eqn{H^0:RR(t,v)=RR(t)} for all times \eqn{(t)} and marks \eqn{(v)}. Inputs can be \code{pvalnp1} or \code{pvalnp2} from the \code{STATvec} component of the \code{cmprsk} output list. These p-values are for two-sided hypothesis tests. The default is \code{pvalnp2}.}  
  \item{filename}{String name for the output plot file. All files are output as .ps files, the string does not need to include the .ps extension. This is a required argument.}
  \item{main}{Main title for the plot}
  \item{ylim}{y axis limits for plot}
  \item{xlim}{x axis limits for plot. If no \code{xlim} is supplied, the default used is \cr
 \code{xlim}=(\code{V1},\code{V2}).}
  \item{legloc}{Vector of length 2 specifing the upper right coordinates of the legend; otherwise the legend is placed in the upper right corner of the plot}
  \item{xlab}{x axis label}
  \item{ylab}{y axis label}
  \item{...}{Additional graphical parameters to pass to plot function}
}
\details{
A labeled plot of cumulative or doubly cumulative VE or RR plot with confidence intervals and p-values is produced. Only data with mark values between the smallest and largest mark values used in the estimation are plotted.}
\references{
Gilbert, P., McKeague, I. and Sun, Y. (2008) \emph{The 2-sample problem for failure rates depending on a continuous mark: an application to vaccine efficacy}. Biostatistics 9, 2, 263-276.\cr
Andersen, P. K., Borgan, O., Gill, R. D., and Keiding. (1993) \emph{Statistical Models Based on Counting Processes.} New York: Springer.
}
\author{Peter Gilbert, Ian McKeague, and Yanqing Sun \cr
Maintainer: Valerie Obenchain <vobencha@fhcrc.org>}
\seealso{
\code{summaryCmprsk} and \code{cmprsk} functions
}
\examples{
## Simulate data
set.seed(100)
gp <- rep(0:1, times=50)
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
x <- cmprsk(gp,ftime,ftype,mark,nboot=500,T1=2,T2=36,UT1=1,UT2=42,ttanal=36)

## Plots
V1 <- x[[1]]$V1
V2 <- x[[1]]$V2
ttanal <- x[[1]]$ttanal
xlab <- "Strain Distance v"

# VE doubly cumulative 
plotCmprsk(x,plottype=1,dc=TRUE,filename="plotVEdc",
main=paste("Estimated VE^dc(t,v) vs Distance v, at time-point t=",
ttanal,sep=""),ylim=c(-1,1.5),xlim=c(V1,V2),legloc=1,
xlab=xlab,ylab="Estimated VE^dc(t,v)")

# VE cumulative 
plotCmprsk(x,plottype=1,dc=FALSE,filename="plotVEc",
main=paste("Estimated VE^c(t,v) vs Distance v, at time-point t=",
ttanal,sep=""), ylim=c(-1,1.5), xlim=c(V1,V2), legloc=1,
xlab=xlab,ylab="Estimated VE^c(t,v)")

# RR doubly cumulative 
plotCmprsk(x,plottype=2,dc=TRUE,filename="plotRRdc",
main=paste("Estimated RR^dc(t,v) vs Distance v, at time-point t=",
ttanal,sep=""), ylim=c(0.0,4.5), xlim=c(V1,V2), legloc=1,
xlab=xlab,ylab="Estimated RR^dc(t,v)")

# RR cumulative
plotCmprsk(x,plottype=2,dc=FALSE,filename="plotRRc",
main=paste("Estimated RR^c(t,v) vs Distance v, at time-point t=",
ttanal,sep=""), ylim=c(0.0,4.5), xlim=c(V1,V2), legloc=1,
xlab=xlab,ylab="Estimated RR^c(t,v)")

}
\keyword{survival}

