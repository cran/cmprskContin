\name{cmprsk}
\alias{cmprsk}
\title{continuous mark-specific relative risks for two treatment groups}
\description{
Estimation and testing of continuous mark-specific relative risks in two groups as described in Gilbert, McKeague and Sun (2008, Biostatistics). Motivations for this work were applications in preventive HIV vaccine efficacy trials and therefore results are returned as mark-specific vaccine efficacies (VEs) as well as mark-specific relative risks (RRs) where \eqn{RR=1-VE}.  This function implements the paper's methods for testing of mark-specific hazards ratios (i.e., relative risks), and for estimation 
of mark-specific cumulative incidence ratios (cumulative in time) and of mark-specific doubly cumulative incidence ratios (cumulative in both time and the continuous mark).

}
\usage{
cmprsk(gp, ftime, ftype, mark, nboot=5000, ngrid=25, ngridv=25, 
T1=0, T2=0, UT1=0, UT2=0, ttanal=0, BAND1=0, BAND2=0, TAILSL=1, 
TAILSU=1, V1=0, V2=0, UV1=0, UV2=0, BANDV1=0, BANDV2=0, BANDVLOW=0, 
BANDVUP=0, TAILSV=1)
}
\arguments{
\item{gp}{Vector of 0 and 1 indicating subject group. Missing values are not allowed.}
\item{ftime}{Vector of right-censored failure times. This value should be the earilest
of event time (i.e., HIV infection) and censoring time. Missing values are not allowed.}
\item{ftype}{Vector indicating 1=failure (infected) or 0=right censored (uninfected). Missing values are not allowed.}
\item{mark}{Vector of continuous mark value. This value should be in [0,1] if infected and measured, 99 otherwise. The mark for all subjects with \code{ftype}=0 should be set to 99 (i.e., mark does not exist for these subjects). A mark value of 99 should be used for subjects with \code{ftype}=1 if the mark is missing or was not measured. This function performs a complete-case analysis and the \code{ftype}=1 subjects with a mark value of 99 will be excluded.}
\item{nboot}{Number of bootstrap replicates for computing p-values, \code{nboot} \eqn{\geq 5000} is recommended}
 \item{ngrid}{(optional) Number of failure time gridpoints for computing the test process}
 \item{ngridv}{(optional) Number of mark gridpoints for computing the test process}
 \item{T1}{(optional) Smallest failure time used to compute the test statistics. This value should be just
larger than the smallest observed failure time. The default (indicated by \code{T1}=0) uses the smallest 
failure time such that there are 2 failure times smaller than \code{T1} in each group.}
  \item{T2}{(optional) Largest failure time used to compute the test statistics. This value should be just 
smaller than the largest observed failure time. The default (indicated by \code{T2}=0) is the largest 
failure time such that there are 2 failure times larger than \code{T2} in each group.}
  \item{UT1}{(optional) Lower envelope for smoothing the failure time \eqn{(\code{UT1} \leq \code{T1})}. The default is \eqn{\code{UT1} = \frac{\code{T1}}{2}}.}
  \item{UT2}{(optional) Upper envelope for smoothing the failure time \eqn{(\code{T2} \leq \code{UT2})}. It is recommended to 
set \code{UT2} equal to the largest observed failure time or a bit bigger. The default (indicated by \code{UT2}=0) uses the 
smallest of the largest right-censored failure times for the two groups.}
  \item{ttanal}{(optional) The failure time \eqn{(t)} at which \eqn{VE^c(t,v)} and \eqn{VE^{dc}(t,v)} are evaluated. 
The default (indicated by \code{ttanal}=0) is \code{ttanal}=\code{T2}. 
}
  \item{BAND1}{(optional) Bandwidth for kernel estimation of the group 1 hazard over time. The default (indicated by \code{BAND1}=0) is the optimal bandwidth,
defined as the bandwidth that minimizes an asymptotic approximation of the mean integrated squared error of 
\eqn{\hat \lambda_1(t)} with the second derivative bandwidth set equal to \eqn{\frac{(\code{T2}-\code{T1})}{4}}. Details of this optimization are in Andersen, P. K., Borgan, O., Gill, R. D., and Keiding (1993). }
  \item{BAND2}{Bandwidth for kernel estimation of the group 2 hazard over time. The default (indicated by \code{BAND2}=0) is the optimal bandwidth chosen the same as for group1.}
  \item{TAILSL}{(optional) Left tail-correction for kernel estimation of the hazard functions over time.
1=use correction (default), 2=no correction.  
}
  \item{TAILSU}{(optional) Right tail-correction for kernel estimation of the hazard functions over time.
1=use correction (default), 2=no correction.
}
  \item{V1}{(optional) Left limit for smoothing in the mark \eqn{(V)}. The default (indicated by \code{V1}=0) is the largest of
the smallest observed marks for the 2 groups.
}
  \item{V2}{(optional) Right limit for smoothing in the mark \eqn{(V)}. The default (indicated by \code{V2}=0) is the
smallest of the of the largest observed marks for the 2 groups.
}
  \item{UV1}{(optional) Left envelope for smoothing in the mark; must have \eqn{\code{UV1} \leq \code{V1}}. The default (indicated by \code{UV1}=0) is \eqn{\code{UV1}=\code{V1}}.}
  \item{UV2}{(optional) Right envelope for smoothing in the mark; must have \eqn{\code{V2} \leq \code{UV2}}. The default (indicated by \code{UV2}=0) is \eqn{\code{UV2}=\code{V2}}.}
  \item{BANDV1}{(optional) Bandwidth for smoothing in the group 1 mark. 
The default (indicated by \code{BANDV1}=0) is to use the optimal bandwidth, chosen to minimize the 2-fold cross-validated mean integrated squared error.}
  \item{BANDV2}{(optional) Bandwidth for smoothing in the group 2 mark. 
The default (indicated by \code{BANDV2}=0) is to use the optimal bandwidth, chosen to minimize the 2-fold cross-validated mean integrated squared error.}
  \item{BANDVLOW}{(optional) The lower range of the mark values to search over for computing the optimized \code{BANDV1} and \code{BANDV2} arguments.
The default (indicated by \code{BANDVLOW}=0) is 5\% of the way from the maximum smallest observed mark in the 2 groups to the \code{BANDVUP} value.}
  \item{BANDVUP}{(optional) The upper range of mark values to search over to compute the optimized \code{BANDV1} and \code{BANDV2} arguments.
The default (indicated by \code{BANDVUP}=0) is the maximum mark divided by 2.}
  \item{TAILSV}{(optional) Upper and lower tail-correction for smoothing in the mark. 1=smoothing (default), 2=no smoothing.}
}
\details{

This package performs a complete case analysis, restricting to subjects with complete data. Missing values for \code{gp}, \code{ftime} and \code{ftype} are not allowed. For all subjects with \code{ftype}=0 (uninfected), the mark value should be set to 99. Subjects with \code{ftype}=1 (infected) should have a mark value in [0,1]. If the mark is missing or was not measured for these \code{ftype}=1 subjects, the mark should be set to 99 and they will be excluded from the analysis.

The results are returned as a list. The first list item, \code{STATvec} contains test statistics, p-values and the estimated \eqn{VE(t,v)} or \eqn{RR(t,v)} over the grid of values for \eqn{v}, at the fixed time-point \code{ttanal}. The second and third items, \code{VECmat} and \code{VEDCmat}, are matrices containing point and 95\% confidence interval estimates of the cumulative and doubly cumulative \eqn{VE(t,v)} or \eqn{RR(t,v)} over the grid of values for \eqn{v}, respectively, at time-point \code{ttanal}. The fourth item, \code{timeMark}, is a matrix containing the point and confidence interval estimates of the cumulative and doubly cumulative \eqn{VE(t,v)} or \eqn{RR(t,v)}, over the grid of time-points \eqn{t} between \code{T1} and \code{T2} as well as over the grid of values for \eqn{v}.

}


\value{
A 4-item list is returned, \code{STATvec}, \code{VECmat}, \code{VEDCmat}, and \code{timeMark}. The components of each list are detailed below.\cr
\item{STATvec}{A vector containing test statistics and p-values. \cr
\tabular{ll}{
nsamp1 \tab  Number of group 1 subjects\cr
nsamp2 \tab Number of group 2 subjects\cr
nboot  \tab Number of bootstrap replicates used for computing p-values\cr
T1     \tab Smallest failure time used to compute the test statistics\cr
T2     \tab Largest failure time used to compute the test statistics\cr
ttanal \tab The failure time \eqn{(t)} at which \eqn{VE^c(t,v)} and \eqn{VE^{dc}(t,v)} are evaluated\cr
AvgV1  \tab Average mark among infected group 1 subjects\cr
AvgV2  \tab Average mark among infected group 2 subjects\cr
BAND1  \tab Bandwidth for estimating the group 1  hazard over time\cr
BAND2  \tab Bandwidth for estimating the group 2  hazard over time \cr
BANDV1 \tab Bandwidth for smoothing in the group 1 mark\cr
BANDV2 \tab Bandwidth for smoothing in the group 2 mark\cr 
V1 \tab The smallest mark used for estimation (used for plotting)\cr
V2 \tab The largest mark used for estimation (used for plotting)\cr
VEhatCI \tab Cumulative incidence estimate of overall VE\cr
VEhatPH \tab Proportional hazard estimate of overall VE\cr
LogRankZ \tab Log-rank statistic for comparing overall hazard rates\cr
U11 \tab Nonparametric test statistic for 1-sided 
testing of \eqn{H^0_0:VE(t,v)=0}\cr
U12 \tab Nonparametric test statistic for 1-sided testing 
of \eqn{H^0_0:VE(t,v)=0}\cr
U13 \tab Nonparametric test statistic for 2-sided 
testing of \eqn{H^0_0:VE(t,v)=0}\cr
U14 \tab Nonparametric test statistic for 2-sided
testing of \eqn{H^0_0:VE(t,v)=0} \cr
pval11 \tab p-value corresponding to \code{U11}\cr
pval12 \tab p-value corresponding to \code{U12}\cr
pval13 \tab p-value corresponding to \code{U13}\cr
pval14 \tab p-value corresponding to \code{U14}\cr
Unp1 \tab Nonparametric test statistic for 1-sided testing of \eqn{H_0: VE(t,v) = VE(t)}\cr
Unp2 \tab Nonparametric test statistic for 2-sided testing of \eqn{H_0: VE(t,v) = VE(t)}\cr
pvalnp1 \tab p-value corresponding to \code{Unp1}\cr
pvalnp2 \tab p-value corresponding to \code{Unp2}\cr
Usp1 \tab Semiparametric test statistic for 1-sided testing of \eqn{H_0: VE(t,v) = VE(t)}\cr
Usp2 \tab Semiparametric test statistic for 2-sided testing of \eqn{H_0: VE(t,v) = VE(t)}\cr
pvalsp1 \tab p-value corresponding to \code{Usp1}\cr
pvalsp2 \tab p-value corresponding to \code{Usp2}}
}

\item{VECmat}{A ngridv x 11 matrix containing data for estimating \eqn{VE^c(ttanal,v)} vs \eqn{v}.
\tabular{ll}{
index \tab Index of the analysis time \code{ttanal}\cr
ttanal \tab Failure time at which the cumulative VE is evaluated\cr
mark \tab Mark value \cr
F1 \tab Estimated cumulative incidence function \eqn{F_1(ttanal,v)} for group 1\cr
F2 \tab Estimated cumulative incidence function \eqn{F_2(ttanal,v)} for group 2\cr
varF1 \tab Estimated variance of \code{F1}\cr
varF2 \tab Estimated variance of \code{F2}\cr
VEC \tab Estimated cumulative vaccine efficacy, \eqn{1 - \frac{F_1(ttanal,v)}{F_2(ttanal,v)}}\cr
varVEC \tab Estimated variance of \code{VEC}\cr
CIlow \tab Lower 95\% CI limit for \code{VEC}\cr
CIhigh \tab Upper 95\% CI limit for \code{VEC}}
}

\item{VEDCmat}{An ngridv x 11 matrix containing data for estimating \eqn{VE^{dc}(ttanal,v)} vs \eqn{v}. 
\tabular{ll}{
index \tab Index of the analysis time \eqn{ttanal}\cr
ttanal \tab Failure time at which the doubly cumulative VE is evaluated\cr
mark \tab Mark value\cr
F1dc \tab Estimated doubly cumulative incidence function for group 1\cr
F2dc \tab Estimated doubly cumulative incidence function for group 2\cr
varF1dc \tab Estimated variance of doubly cumulative incidence function for group 1\cr
varF2dc \tab Estimated variance of doubly cumulative incidence function for group 2\cr
VEDC \tab Estimated doubly cumulative vaccine efficacy\cr
varVEDC \tab Estimated variance of \code{VEDC}\cr
CIlow \tab Lower 95\% CI limit for \code{VEDC}\cr
CIhigh \tab Upper 95\% CI limit for \code{VEDC}}
}

\item{timeMark}{A matrix containing data for estimating \eqn{VE^c(t,v)} vs \eqn{(t,v)} at all event times \eqn{t} between \code{T1} and \code{T2} (ie, all times between the smallest and largest failure times used for computing the test statistics).
\tabular{ll}{
nsamp \tab Sample index\cr
eventtime \tab Observed event times\cr
mark \tab Mark value \cr 
F1 \tab Estimated cumulative incidence function \eqn{F_1(t,v)}  for group 1\cr
F2 \tab Estimated cumulative incidence function \eqn{F_2(t,v)} for group 2\cr 
varF1 \tab  Estimated variance of \code{F1} \cr
varF2 \tab Estimated variance of \code{F2} \cr
VEC \tab Estimated cumulative vaccine efficacy \cr
varVEC \tab Estimated variance of cumulative vaccine efficacy} 
}
}

\references{
Gilbert, P., McKeague, I. and Sun, Y. (2008) \emph{The 2-sample problem for failure rates depending on a continuous mark: an application to vaccine efficacy}. Biostatistics 9, 2, 263-276.\cr
Andersen, P. K., Borgan, O., Gill, R. D., and Keiding. (1993) \emph{Statistical Models Based on Counting Processes.} New York: Springer.}


\author{Peter Gilbert, Ian McKeague, and Yanqing Sun \cr
Maintainer: Valerie Obenchain <vobencha@fhcrc.org>}
\seealso{
\code{summaryCmprsk} and \code{plotCmprsk} functions
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
result <- cmprsk(gp,ftime,ftype,mark,nboot=500,T1=2,T2=36,UT1=1,
UT2=42,ttanal=36)

## See cumulative incidence and proportional hazard estimates of overall VE 
result[[1]]$VEhatCI
result[[1]]$VEhatPH

}
\keyword{survival}

