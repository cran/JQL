\name{find.I2DR}
\alias{find.I2DR}
\title{Estimating the Individualized Interval-valued Dose Rule via  (Residual) Jump Q-learning.}
\description{
  This function estimates the optimal Individualized Interval-valued Dose Rule (I2DR), and calculates a Wald-type confidence interval for the value function under the estimated optimal I2DR via Bootstrap.
}
\usage{
  find.I2DR(Y,A,X,cm=6,method='JQL',Gamma.list=seq(from=1,to=20,by=2)/5,
  Lambda.list=seq(from=1,to=20,by=2)/5,RF_A.list=c(0,0.25,0.5,0.75,1),
  folds_num=5,alpha=0.95,nboots=500)
  }
\arguments{
  \item{Y}{The patient’s associated response/outcome, the larger the better by convention.}
  \item{A}{The dose level received by each patient, should be continuous.}
  \item{X}{The patient’s baseline covariates, could be a matrix, including continous or discrete covariates.}
  \item{cm}{The constent cm in m=n/cm, where m is the number of total subinterval that diverges with sample size n. The default value is 6.}
  \item{method}{Two methods are available, Jump Q-learning ('JQL') and Residual Jump Q-learning ('RJQL'). The default method is 'JQL'.}
  \item{Gamma.list}{The candidate tuning paramter space for c1 in penalty term gamma=c1 log(n)/n. The default value is seq(from=1,to=20,by=2)/5. If the length of Gamma.list is 1, then the tuning process will be skipped.}
  \item{Lambda.list}{The candidate tuning paramter space for c2 in penalty term lambda=c2 log(n)/n. The default value is seq(from=1,to=20,by=2)/5. If the length of Lambda.list is 1, then the tuning process will be skipped.}
  \item{RF_A.list}{The candidate tuning paramter space for A in fitted E(Y|A=a,X) by Random Forest Regression for method 'RJQL' only. The default value is c(0,0.25,0.5,0.75,1). If the length of RF_A.list is 1, then the tuning process will be skipped.}
  \item{folds_num}{The number of the folds in the cross-validation process. The default value is 5.}
  \item{alpha}{The Confidence level. The default level is 0.95.}
  \item{nboots}{The number of Bootstrap. The default number is 500.}
}
\value{
  An object of class "I2DR" which is a list with
  components:

  \item{Partition}{A partition of the entire dose range.}

  \item{Beta}{The regression coefficients for each partition.}

  \item{Value}{The estimated value function under our proposed I2DR.}

  \item{low_bd}{The lower bound of the confidence interval.}

  \item{up_bd}{The upper bound of the confidence interval.}

  \item{method}{The method used to find the I2DR.}
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
find.I2DR(Y=y,A=a,X=x)
}

