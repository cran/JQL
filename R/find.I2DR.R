#' Estimating the Individualized Interval-valued Dose Rule via (Residual) Jump Q-learning.
#'
#' \code{find.I2DR} This function estimates the optimal Individualized Interval-valued Dose Rule (I2DR), and calculates a Wald-type confidence interval for the value function under the estimated optimal I2DR via Bootstrap.
#' @param Y The patient’s associated response/outcome, the larger the better by convention.
#' @param A The dose level received by each patient, should be continuous.
#' @param X The patient’s baseline covariates, could be a matrix, including continous or discrete covariates.
#' @param cm The constent cm in m=n/cm, where m is the number of total subinterval that diverges with sample size n. The default value is 6.
#' @param method Two methods are available, Jump Q-learning ('JQL') and Residual Jump Q-learning ('RJQL'). The default method is 'JQL'.
#' @param Gamma.list The candidate tuning paramter space for c1 in penalty term gamma=c1 log(n)/n. The default value is seq(from=1,to=20,by=2)/5. If the length of Gamma.list is 1, then the tuning process will be skipped.
#' @param Lambda.list The candidate tuning paramter space for c2 in penalty term lambda=c2 log(n)/n. The default value is seq(from=1,to=20,by=2)/5. If the length of Lambda.list is 1, then the tuning process will be skipped.
#' @param RF_A.list The candidate tuning paramter space for A in fitted E(Y|A=a,X) by Random Forest Regression for method 'RJQL' only. The default value is c(0,0.25,0.5,0.75,1). If the length of RF_A.list is 1, then the tuning process will be skipped.
#' @param folds_num The number of the folds in the cross-validation process. The default value is 5.
#' @param alpha The Confidence level. The default level is 0.95.
#' @param nboots The number of Bootstrap. The default number is 500.
#' @importFrom pdist pdist
#' @importFrom stats optim
#' @importFrom caret createFolds
#' @importFrom stats dist
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @export find.I2DR
#'
find.I2DR <- function(Y,A,X,cm=6,method='JQL',Gamma.list=seq(from=1,to=20,by=2)/5,
                      Lambda.list=seq(from=1,to=20,by=2)/5,RF_A.list=c(0,0.25,0.5,0.75,1),
                      folds_num=5,alpha=0.95,nboots=500) {

  ############Checking the Formation############
  if((length(A)!=length(Y))|(length(A)!=dim(X)[1])){
    stop("The length of each data is not equal!\n")
  }

  data=data.frame(Y,A,X)
  a=(A-min(A))/(max(A)-min(A)) # convert the dosage domain to [0,1]
  n=length(Y) # sample size
  d=dim(X)[2]+1

  sample=data.frame(y=Y,a=a,x=X)

  # Use the cross-validation to train the best tuning parameters lambda_n and gamma_n

  # Set the range of the best tuning parameters
  Gamma.list=Gamma.list*log(n)/n
  Lambda.list=Lambda.list*log(n)/n

  if(method=='JQL'){
    if(length(Gamma.list)==1&length(Lambda.list)==1){
      gamma_n=Gamma.list
      lambda_n=Lambda.list
    }else{
      best.paras=tune.JQL(sample=sample,cm=cm,Gamma.list=Gamma.list,Lambda.list=Lambda.list,folds_num=folds_num)
      # Recalculate the optimal I2DR under the best tuning parameters
      gamma_n=best.paras$best_gamma
      lambda_n=best.paras$best_lambda
    }
  }

  if(method=='RJQL'){
    if(length(Gamma.list)==1&length(Lambda.list)==1&length(RF_A.list)==1){
      best_a=RF_A.list
      gamma_n=Gamma.list
      lambda_n=Lambda.list
    }else{
      best.paras=tune.RJQL(sample=sample,cm=cm,Gamma.list=Gamma.list,Lambda.list=Lambda.list,RF_A.list=RF_A.list,folds_num=folds_num)
      # Recalculate the optimal I2DR under the best tuning parameters
      best_a=best.paras$best_a
      gamma_n=best.paras$best_gamma
      lambda_n=best.paras$best_lambda
    }
    # Refit the random forest with best prediction
    rf <- randomForest(as.matrix(cbind(a,X)), Y)
    best_trainXA=cbind(best_a,X)
    colnames(best_trainXA)=NULL
    Y_fit=predict(rf, as.matrix(best_trainXA))
    sample=data.frame(y=Y-Y_fit,a=a,x=X)
  }


  # Generate the optimal I2DR under the selected best tuning parameters.
  n=nrow(sample)
  m=n/cm

  # Find best partition
  tao=Bel=rep(0,m)

  for(r in 1:m){
      Bel[r]=1000000

      if(r<m){
        phi_data=sample[sample$a<(r/m),]
      }else{
        phi_data=sample[sample$a<=1,]
      }

      phi_X=matrix(1,nrow=length(phi_data$y),ncol=d)
      phi_X[,2:d]=as.matrix(phi_data[,3:(d+1)])
      phi=t(phi_X)%*%phi_data$y

      cp_phi=t(phi_X)%*%phi_X+n*lambda_n*r/m*diag(d)

      for(l in 1:r){

        if(l==1){
          Bel_l_1=-gamma_n
        }else{
          Bel_l_1=Bel[l-1]
        }

        # Calculating the distance
        l_data=sample[sample$a<((l-1)/m),]
        l_X=matrix(1,nrow=length(l_data$y),ncol=d)
        l_X[,2:d]=as.matrix(l_data[,3:(d+1)])

        betaI=solve(cp_phi-t(l_X)%*%l_X-n*lambda_n*(l-1)/m*diag(d))%*%(phi-t(l_X)%*%as.matrix(l_data$y))

        ##############################################

        if(r<m){
          int_data=sample[sample$a>=((l-1)/m)&sample$a<(r/m),]
        }else{
          int_data=sample[sample$a>=((l-1)/m)&sample$a<=1,]
        }
        int_X=matrix(1,nrow=length(int_data$y),ncol=d)
        int_X[,2:d]=as.matrix(int_data[,3:(d+1)])

        c=Bel_l_1+gamma_n+sum((as.matrix(int_data$y)-int_X%*%betaI)^2)/n+((r-l+1)/m)*lambda_n*sum(betaI^2)
        if(c<=Bel[r]){
          Bel[r]=c
          tao[r]=l-1
        }
      }
    }

    ##############################################

  # Segmentation from partition

  r=m
  l=tao[r]
  intvl_ep=NULL
  seg_rule=NULL
  while(r>0){
      intvl_ep=c(r/m,intvl_ep)
      # Calculating the linear rule
      if(r<m){
        rdgdata=sample[sample$a>=(l/m)&sample$a<(r/m),]
      }else{
        rdgdata=sample[sample$a>=(l/m)&sample$a<=1,]
      }

      designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
      designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
      betaI=solve(t(designX)%*%designX+n*lambda_n*(r-l)/m*diag(d))%*%(t(designX)%*%as.matrix(rdgdata$y))

      ##############################################
      seg_rule=c(betaI,seg_rule)
      r=l
      l=tao[r]
    }

  # Calculate the value function V(d)
  Xj1=cbind(rep(1,n),X)
  value=rep(0,n)
  seg_rule=as.matrix(seg_rule)
  for(j in 1:n){
    t=as.matrix(Xj1[j,])
    temp=-10000
    for(i in 1:length(intvl_ep)){
      if(t(t)%*%seg_rule[(d*i-(d-1)):(d*i)]>temp){
        temp=t(t)%*%seg_rule[(d*i-(d-1)):(d*i)]
      }
    }
    value[j]=temp
  }
  M_n=mean(value) # estimated value

  if(method=='JQL'){
    cat('Method: Jump Q-learning \n')
  }
  if(method=='RJQL'){
    cat('Method: Residual Jump Q-learning \n')
  }

  # Output the partition and the corresponding regression coefficients
  num_int=length(intvl_ep)
  intvl_ep=intvl_ep*(max(A)-min(A))+min(A)
  seg_rule=as.matrix(seg_rule)
  Beta=matrix(0,nrow=num_int,ncol=d)
  cat('Interval:[',min(A),',',intvl_ep[1],']\n')
  Beta[1,]=seg_rule[1:d]
  cat('Best rule:',Beta[1,],'\n')
  if(num_int>1){
    for(i in 2:num_int){
      cat('Interval: (',intvl_ep[i-1],',',intvl_ep[i],']\n')
      Beta[i,]=seg_rule[((i-1)*d+1):(i*d)]
      cat('Best rule:',Beta[i,],'\n')
    }
  }

  # Output the value function

  cat('Value function:',M_n,'\n')

  # Use Bootstrap to contrust the Wald-type confidence interval for value function by fixing the partition found above
  sample=as.matrix(sample)
  data_resam=array(0,dim=c(n,ncol(sample),nboots))
  beta_Boots=array(0,dim=c(nboots,length(intvl_ep),d))
  value_Boots=rep(0,nboots)
  intvls=c(0,intvl_ep)
  for(k in 1:nboots){
    index=sample(c(1:n),n,replace=TRUE)
    data_resam[,,k]=k_data=sample[as.vector(index),]
    i=1
    while(i<length(intvls)){
      if(i==(length(intvls)-1)){
        rdgdata=k_data[k_data[,2]>=intvls[i]&k_data[,2]<=intvls[i+1],]
        if(is.null(dim(rdgdata)[1])==F){
          designX=matrix(1,nrow=length(rdgdata[,1]),ncol=d)
          designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
          beta_Boots[k,i,]=solve(t(designX)%*%designX+(intvls[i+1]-intvls[i])*diag(d))%*%t(designX)%*%rdgdata[,1]
        }
        if(is.null(dim(rdgdata)[1])==T){
          designX=matrix(1,nrow=length(rdgdata[1]),ncol=d)
          designX[,2:d]=as.matrix(rdgdata[3:(d+1)])
          beta_Boots[k,i,]=solve(t(designX)%*%designX+(j-i)/m*diag(d))%*%t(designX)*rdgdata[1]
        }
      }else{
        rdgdata=k_data[k_data[,2]>=intvls[i]&k_data[,2]<intvls[i+1],]
        if(is.null(dim(rdgdata)[1])==F){
          designX=matrix(1,nrow=length(rdgdata[,1]),ncol=d)
          designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
          beta_Boots[k,i,]=solve(t(designX)%*%designX+(intvls[i+1]-intvls[i])*diag(d))%*%t(designX)%*%rdgdata[,1]
        }
        if(is.null(dim(rdgdata)[1])==T){
          designX=matrix(1,nrow=length(rdgdata[1]),ncol=d)
          designX[,2:d]=as.matrix(rdgdata[3:(d+1)])
          beta_Boots[k,i,]=solve(t(designX)%*%designX+(j-i)/m*diag(d))%*%t(designX)*rdgdata[1]
        }
      }
      i=i+1
    }

    Xj1=cbind(rep(1,n),data_resam[,3:(d+1),k])
    value=rep(0,n)
    for(j in 1:n){
      t=as.matrix(Xj1[j,])
      temp=-10000
      for(i in 1:length(intvl_ep)){
        b=as.matrix(beta_Boots[k,i,])
        if(t(t)%*%b>temp){
          temp=t(t)%*%b
        }
      }
      value[j]=temp
    }
    value_Boots[k]=mean(value)
  }

  # the 95% Wald-type CI for the value

  low_bd=M_n-qnorm((1+alpha)/2,0,1)*sd(value_Boots)
  up_bd=M_n+qnorm((1+alpha)/2,0,1)*sd(value_Boots)

  cat(alpha*100,'% Wald-type CI for the value: (',low_bd,',',up_bd,')\n')

  return(list(Partition=c(min(A),intvl_ep),Beta=Beta,Value=M_n,low_bd=low_bd,up_bd=up_bd, method=method))

}
