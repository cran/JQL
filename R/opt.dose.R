#' Optimal Interval-valued Dosage under the Individualized Interval-valued Dose Rule via Jump Q-learning.
#'
#' \code{opt.dose} This function assigns each individual to one of the subintervals of the entire dosage according to his/her baseline covariates under the estimated I2DR.
#' @param X The patient’s baseline covariates, coule be a matrix, including continous or discrete covariates.
#' @param JQL The Individualized Interval-based Dose Rule found by the function "JQL".
#' @export opt.dose
#'
#'
opt.dose <- function(X,JQL) {

  ############Checking the Formation############
  if(dim(JQL$Beta)[2]!=(dim(X)[2]+1)){
    stop("The dimension of the covariates is not equal to the decision rule found by JQL!\n")
  }

  n=dim(X)[1]
  d=dim(X)[2]+1
  Xj1=as.matrix(cbind(1,X))
  num_int=length(JQL$Partition)-1
  Beta=as.matrix(JQL$Beta)

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
    opt.dose[j,]=c(JQL$Partition[flag.dose],JQL$Partition[flag.dose+1])
  }

  return(list(opt.dose=opt.dose))
}
