\name{JQL}
\alias{JQL}
\title{Estimating the Individualized Interval-valued Dose Rule via Jump Q-learning.}
\description{
  This function estimates the optimal Individualized Interval-valued Dose Rule (I2DR), and calculates a Wald-type confidence interval for the value function under the estimated optimal I2DR.
}
\usage{
  JQL(Y,A,X,repnum=20,train_pct=0.7,
  c.list=exp(seq(-2,2,length=10)),folds_num=10,alpha=0.95,nboots=500)
}
\arguments{
  \item{Y}{The patient’s associated response/outcome, the larger the better by convention.}
  \item{A}{The dose level received by each patient, should be continuous.}
  \item{X}{The patient’s baseline covariates, coule be a matrix, including continous or discrete covariates.}
  \item{repnum}{The number of the replications to eliminate the bias caused by the heterogeneity of the data, where we randomly split the data into training and testing sets independently. The default value is 20.}
  \item{train_pct}{The percentage of the training data to to train the tuning parameter. The default value is 70\%.}
  \item{c.list}{The candidate tuning paramter space. The default value is exp(seq(-2,2,length=10)).}
  \item{folds_num}{The number of the folds in the cross-validation process. The default value is 10.}
  \item{alpha}{The Confidence level. The default level is 0.95.}
  \item{nboots}{The number of Bootstrap. The default number is 500.}
}
\value{
  An object of class "JQL" which is a list with
  components:

  \item{Partition}{A partition of the entire dose range.}

  \item{Beta}{The regression coefficients for each partition.}

  \item{Value}{The estimated value function under our proposed I2DR.}

  \item{low_bd}{The lower bound of the confidence interval.}

  \item{up_bd}{The upper bound of the confidence interval.}
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
JQL(Y=y,A=a,X=x,folds_num=5,repnum=1)
}

