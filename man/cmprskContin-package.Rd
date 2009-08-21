\name{cmprskContin-package}
\alias{cmprskContin-package}
\alias{cmprskContin}
\docType{package}
\title{
continuous mark-specific relative risks for two treatment groups 
}
\description{
Estimation and testing of continuous mark-specific relative risks in two groups as described in Gilbert, McKeague and Sun (2008). Results include point and confidence interval estimates, test statistics, and corresponding p-values.\cr
}

\details{

\tabular{ll}{
Package: \tab cmprskContin\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2009-08-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Motivations for this work were applications in preventive HIV vaccine trials and therefore results are returned as vaccine efficacies (VEs) as well as relative risks (RRs) where \eqn{RR=1-VE}. 
In the following equations, \eqn{\lambda_k} represents the mark-specific hazard rate for group $k$, $k=1,2$. 
Vaccine efficacy is defined as the relative reduction in infection hazard due to vaccination, \eqn{VE(t,v)=1-\frac{\lambda_1(t,v)}{\lambda_2(t,v)}} where \eqn{t} is the time of infection and \eqn{v} is the mark (e.g., the genetic distance of an infecting virus from the vaccine strain).\cr 

Tests are provided for \eqn{H^0_0: \lambda_1(t,v) = \lambda_2(t,v)}, against the following alternative hypotheses: \cr
\eqn{H^0_1: \lambda_1(t,v) \leq \lambda_2(t,v)} \cr
\eqn{H^0_2: \lambda_1(t,v) \neq \lambda_2(t,v)}\cr

If \eqn{H^0_0} is rejected, then it is of interest to asses if the relative risk varies with the mark. For this question, tests are provided for  \eqn{H_0: \frac{\lambda_1(t,v)}{\lambda_2(t,v)}} does not depend on \eqn{v} for \eqn{t \in [0,\tau]}, against the following alternative hypotheses: \cr
\eqn{H_1: \frac{\lambda_1(t,v_1)}{\lambda_2(t,v_1)} \leq \frac{\lambda_1(t,v_2)}{\lambda_2(t,v_2)}} \cr
\eqn{H_2: \frac{\lambda_1(t,v_1)}{\lambda_2(t,v_1)} \neq \frac{\lambda_1(t,v_2)}{\lambda_2(t,v_2)}} \cr

Both integration-type and supremum-type statistics are provided for testing the hypotheses. Test statistics \eqn{U^1_1} and \eqn{U^1_2} provide evidence against \eqn{H^0_0} in the direction of \eqn{H^0_1} and  \eqn{U^1_3} and \eqn{U^1_4} provide evidence against \eqn{H^0_0} in the direction of \eqn{H^0_2}. For \eqn{r} as np (nonparametric) or sp (semiparametric), \eqn{U^r_1} and \eqn{U^r_2} measure departures from \eqn{H_0} in the direction of \eqn{H_1} and \eqn{H_2}, respectively. 

P-values are computed using a Guassian multipliers simulation procedure.

Results provided include the cumulative incidence estimate of overall VE, proportional hazards estimate of overall VE, 
log-rank statistic for comparing overall hazard rates, test statistics for evaluating the hypotheses outlined above, and the corresponding p-values. 

There are 3 functions : \cr
\code{cmprsk}\cr
\code{plotCmprsk}\cr
\code{summaryCmprsk}\cr
}
\author{
Peter Gilbert, Ian McKeague, and Yanqing Sun\cr
Maintainer: Valerie Obenchain <vobencha@fhcrc.org>
}

\references{
Gilbert, P., McKeague, I. and Sun, Y. (2008) \emph{The 2-sample problem for failure rates depending on a continuous mark: an application to vaccine efficacy}. Biostatistics 9, 2, 263-276. \cr
Andersen, P. K., Borgan, O., Gill, R. D., and Keiding. (1993) \emph{Statistical Models Based on Counting Processes.} New York: Springer.
}

\keyword{ package }
\seealso{
\code{cmprsk}, \code{summaryCmprsk} and \code{plotCmprsk} functions
}
