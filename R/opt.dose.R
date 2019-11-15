#' Optimal Interval-valued Dosage under the Individualized Interval-valued Dose Rule via (Residual) Jump Q-learning.
#'
#' \code{opt.dose} This function assigns each individual to one of the subintervals of the entire dosage according to his/her baseline covariates under the estimated I2DR.
#' @param X The patient’s baseline covariates, could be a matrix, including continous or discrete covariates.
#' @param I2DR The Individualized Interval-based Dose Rule found by the function "JQL" or "RJQL".
#' @export opt.dose
#'
#'
opt.dose <- function(X,I2DR) {

  ############Checking the Formation############
  if(dim(I2DR$Beta)[2]!=(dim(X)[2]+1)){
    stop("The dimension of the covariates is not equal to the decision rule found by JQL!\n")
  }

  n=dim(X)[1]
  d=dim(X)[2]+1
  Xj1=as.matrix(cbind(1,X))
  num_int=length(I2DR$Partition)-1
  Beta=as.matrix(I2DR$Beta)

  opt.dose=matrix(0,nrow=n,ncol=2)
  for(j in 1:n){
    t=as.matrix(Xj1[j,])
    temp=-10000
    for(i in 1:num_int){
      if(t(t)%*%Beta[i,]>temp){
        temp=t(t)%*%Beta[i,]
        flag.dose=i
      }
    }
    opt.dose[j,]=c(I2DR$Partition[flag.dose],I2DR$Partition[flag.dose+1])
  }

  return(list(opt.dose=opt.dose))
}
