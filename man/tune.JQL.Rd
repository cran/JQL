\name{tune.JQL}
\alias{tune.JQL}
\title{Tuning function via k-fold cross vaidation for Jump Q-learning.}
\description{
This function uses the cross-validation to train the best tuning parameters lambda_n and gamma_n for Jump Q-learning.
}
\usage{
tune.JQL(sample,cm=6,Gamma.list=seq(from=1,to=20,by=2)/5,
Lambda.list=seq(from=1,to=20,by=2)/5,folds_num=5)
}
\arguments{
  \item{sample}{The training dataset (Y,A,X), where Y is the patient’s associated response/outcome, A is the dose level received by each patient, and X is the patient’s baseline covariates.}
  \item{cm}{The constent cm in m=n/cm, where m is the number of total subinterval that diverges with sample size n. The default value is 6.}
  \item{Gamma.list}{The candidate tuning paramter space for c1 in penalty term gamma=c1 log(n)/n. The default value is seq(from=1,to=20,by=2)/5.}
  \item{Lambda.list}{The candidate tuning paramter space for c2 in penalty term lambda=c2 log(n)/n. The default value is seq(from=1,to=20,by=2)/5.}
  \item{folds_num}{The number of the folds in the cross-validation process. The default value is 5.}
}
\value{
  \item{best_gamma}{The best tuning parameter gamma by minimuming the least square loss function.}
  \item{best_lambda}{The best tuning parameter lambda by minimuming the least square loss function.}
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
sample=data.frame(y=y,a=a,x=x)
tune.JQL(sample)
}
