\name{opt.dose}
\alias{opt.dose}
\title{Optimal Interval-valued Dose under the Individualized Interval-valued Dose Rule via (Residual) Jump Q-learning.}
\description{
  This function assigns each individual to one of the subintervals of the entire dosage according to his/her baseline covariates under the estimated I2DR.
}
\usage{
  opt.dose(X,I2DR)
}
\arguments{
  \item{X}{The patient’s baseline covariates, coule be a matrix, including continous or discrete covariates.}
  \item{I2DR}{The Individualized Interval-valued Dose Rule found by the function "JQL" or "RJQL".}
}
\value{

  \item{opt.dose}{The optimal Interval-valued dosage for each individual.}

}
\references{
  Jump Q-learning for Individualized Interval-valued Dose Rule.
}
\examples{
n=50
d=4
x=matrix(runif(n*(d-1),-1,1),nrow=n,ncol=d-1)
a=runif(n,0,1)
y=(1+x[,1])*(a>=0&a<0.35)+(x[,1]-x[,2])*(a>=0.35&a<0.65)+(1-x[,2])*(a>=0.65&a<=1)+rnorm(n,0,1)
rule=find.I2DR(Y=y,A=a,X=x)
n0=10
xnew=matrix(runif(n0*(d-1),-1,1),nrow=n0,ncol=d-1)
opt.dose(X=xnew,I2DR=rule)
}


