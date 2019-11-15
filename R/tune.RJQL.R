#' Tuning function via k-fold cross vaidation for Residual Jump Q-learning.
#'
#' \code{tune.RJQL} This function uses the cross-validation to train the best tuning parameters lambda_n and gamma_n for Residual Jump Q-learning.
#' @param sample The training dataset that includes {Y,A,X}, where Y is the patient’s associated response/outcome, A is the dose level received by each patient, and X is the patient’s baseline covariates.
#' @param cm The constent cm in m=n/cm, where m is the number of total subinterval that diverges with sample size n. The default value is 6.
#' @param Gamma.list The candidate tuning paramter space for c1 in penalty term gamma=c1 log(n)/n. The default value is seq(from=1,to=20,by=2)/5.
#' @param Lambda.list The candidate tuning paramter space for c2 in penalty term lambda=c2 log(n)/n. The default value is seq(from=1,to=20,by=2)/5.
#' @param RF_A.list The candidate tuning paramter space for A in fitted E(Y|A=a,X) by Random Forest Regression for method 'RJQL' only. The default value is c(0,0.25,0.5,0.75,1).
#' @param folds_num The number of the folds in the cross-validation process. The default value is 5.
#' @importFrom pdist pdist
#' @importFrom stats optim
#' @importFrom caret createFolds
#' @importFrom stats dist
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @export tune.RJQL
#'
tune.RJQL=function(sample,cm=6,Gamma.list=seq(from=1,to=20,by=2)/5,Lambda.list=seq(from=1,to=20,by=2)/5,RF_A.list=c(0,0.25,0.5,0.75,1),folds_num=5){

  n=nrow(sample)
  d=ncol(sample)-1
  # create the cross-validation folds, 10 folds default.
  info=matrix(0,nrow=folds_num,ncol=length(RF_A.list)*length(Gamma.list)*length(Lambda.list))
  folds=createFolds(y=c(1:n), k = folds_num, list = TRUE, returnTrain = T)

  tune_x=sample[,3:(d+1)]
  tune_a=sample[,2]
  tune_y=sample[,1]
  rf <- randomForest(as.matrix(cbind(tune_a,tune_x)), tune_y)


  for(k in 1:folds_num){
    # For kth folds
    for(a_index in 1:length(RF_A.list)){
      # Demeaned data
      cv_fit_xa=cbind(RF_A.list[a_index],tune_x)
      colnames(cv_fit_xa)=NULL
      Y_fit=predict(rf, as.matrix(cv_fit_xa))
      sample=data.frame(y=tune_y-Y_fit,a=tune_a,x=tune_x)
      k_data=sample[folds[[k]],]
      nk=nrow(k_data)
      mk=nk/cm

      # Find best partition

      tao=Bel=matrix(0,nrow=mk,ncol=length(RF_A.list)*length(Gamma.list)*length(Lambda.list))

      for(r in 1:mk){

        Bel[r,]=rep(1000000,length(RF_A.list)*length(Gamma.list)*length(Lambda.list))

        if(r<mk){
          phi_data=k_data[k_data$a<(r/mk),]
        }else{
          phi_data=k_data[k_data$a<=1,]
        }

        phi_X=matrix(1,nrow=length(phi_data$y),ncol=d)
        phi_X[,2:d]=as.matrix(phi_data[,3:(d+1)])
        phi=t(phi_X)%*%phi_data$y
        XTX_phi=t(phi_X)%*%phi_X

        for(l in 1:r){

          # Calculating the distance
          l_data=k_data[k_data$a<((l-1)/mk),]
          l_X=matrix(1,nrow=length(l_data$y),ncol=d)
          l_X[,2:d]=as.matrix(l_data[,3:(d+1)])

          cp_phi=XTX_phi-t(l_X)%*%l_X
          decomp=eigen(cp_phi,symmetric=TRUE)
          eigen_value=decomp$values
          U=decomp$vectors
          phi_l_Xl_Y=phi-t(l_X)%*%as.matrix(l_data$y)

          ##############################################

          if(r<mk){
            int_data=k_data[k_data$a>=((l-1)/mk)&k_data$a<(r/mk),]
          }else{
            int_data=k_data[k_data$a>=((l-1)/mk)&k_data$a<=1,]
          }
          int_X=matrix(1,nrow=length(int_data$y),ncol=d)
          int_X[,2:d]=as.matrix(int_data[,3:(d+1)])

          for(lam_index in 1:length(Lambda.list))
          {
            lambda_n=Lambda.list[lam_index]

            betaI=(U%*%((1/(eigen_value+nk*lambda_n*(r-l+1)/mk))*t(U)))%*%phi_l_Xl_Y
            add_error=sum((as.matrix(int_data$y)-int_X%*%betaI)^2)/nk+((r-l+1)/mk)*lambda_n*sum(betaI^2)
            for(gam_index in 1:length(Gamma.list))
            {
              gamma_n=Gamma.list[gam_index]

              if(l==1){
                Bel_l_1=-gamma_n
              }else{
                Bel_l_1=Bel[l-1,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]
              }

              c=Bel_l_1+gamma_n+add_error

              if(c<=Bel[r,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]){
                Bel[r,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]=c
                tao[r,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]=l-1
              }
            }

          }
        }
        ##############################################
      }

      # Segmentation from partition
      for(lam_index in 1:length(Lambda.list))
      {
        lambda_n=Lambda.list[lam_index]

        for(gam_index in 1:length(Gamma.list))
        {
          gamma_n=Gamma.list[gam_index]
          r=mk
          l=tao[r,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]
          intvl_ep=NULL
          seg_rule=NULL
          while(r>0){
            intvl_ep=c(r/mk,intvl_ep)

            # Calculating the linear rule
            if(r<mk){
              rdgdata=k_data[k_data$a>=(l/mk)&k_data$a<(r/mk),]
            }else{
              rdgdata=k_data[k_data$a>=(l/mk)&k_data$a<=1,]
            }

            designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
            designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
            betaI=solve(t(designX)%*%designX+nk*lambda_n*(r-l)/mk*diag(d))%*%(t(designX)%*%as.matrix(rdgdata$y))

            ##############################################
            seg_rule=c(betaI,seg_rule)
            r=l
            l=tao[r,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]
          }

          # Use the left k-fold to calculate the least square loss function.
          rk_data=sample[-folds[[k]],]

          for(i in 1:length(intvl_ep)){
            if(i==1){
              int_data=rk_data[rk_data$a>=0&rk_data$a<intvl_ep[1],]
            }else{
              if(i==length(intvl_ep)){
                int_data=rk_data[rk_data$a>=intvl_ep[i-1]&rk_data$a<=intvl_ep[i],]
              }else{
                int_data=rk_data[rk_data$a>=intvl_ep[i-1]&rk_data$a<intvl_ep[i],]
              }
            }
            designX=matrix(1,nrow=length(int_data$y),ncol=d)
            designX[,2:d]=as.matrix(int_data[,3:(d+1)])
            info[k,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]=info[k,(a_index-1)*(length(Lambda.list)*length(Gamma.list))+(lam_index-1)*length(Gamma.list)+gam_index]+sum((int_data$y-designX%*%as.matrix(seg_rule[(d*i-(d-1)):(d*i)]))^2)
          }
        }
      }
    }
  }

  # Select the best tuning parameters by minimuming the least square loss function
  loc=min(which(apply(info,2,mean)==min(apply(info,2,mean))))
  loc_gam=loc%%length(Gamma.list)
  if(loc_gam==0){
    loc_gam=length(Gamma.list)
  }
  loc_low=(loc-loc_gam)/length(Gamma.list)+1
  loc_lam=loc_low%%length(Lambda.list)
  if(loc_lam==0){
    loc_lam=length(Lambda.list)
  }
  loc_a=(loc_low-loc_lam)/length(Lambda.list)+1

  best_a=RF_A.list[loc_a]
  best_gamma=Gamma.list[loc_gam]
  best_lambda=Lambda.list[loc_lam]
  return(list(best_a=best_a,best_gamma=best_gamma,best_lambda=best_lambda))
}

