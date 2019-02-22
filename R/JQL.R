#' Estimating the Individualized Interval-valued Dose Rule via Jump Q-learning.
#'
#' \code{JQL} This function estimates the optimal Individualized Interval-valued Dose Rule (I2DR), and calculates a Wald-type confidence interval for the value function under the estimated optimal I2DR.
#' @param Y The patient’s associated response/outcome, the larger the better by convention.
#' @param A The dose level received by each patient, should be continuous.
#' @param X The patient’s baseline covariates, coule be a matrix, including continous or discrete covariates.
#' @param repnum The number of the replications to eliminate the bias caused by the heterogeneity of the data, where we randomly split the data into training and testing sets independently. The default value is 20.
#' @param train_pct The percentage of the training data to to train the tuning parameter. The default value is 70%.
#' @param c.list The candidate tuning paramter space. The default value is exp(seq(-2,2,length=10)).
#' @param folds_num The number of the folds in the cross-validation process. The default value is 10.
#' @param alpha The Confidence level. The default level is 0.95.
#' @param nboots The number of Bootstrap. The default number is 500.
#' @importFrom pdist pdist
#' @importFrom stats optim
#' @importFrom caret createFolds
#' @importFrom stats dist
#' @importFrom stats runif
#' @importFrom stats sd
#' @export JQL
#'
JQL <- function(Y,A,X,repnum=20,train_pct=0.7,c.list=exp(seq(-2,2,length=10)),folds_num=10,alpha=0.95,nboots=500) {

  ############Checking the Formation############
  if((length(A)!=length(Y))|(length(A)!=dim(X)[1])){
    stop("The length of each data is not equal!\n")
  }


  A=(A-min(A))/(max(A)-min(A)) # convert the dosage domain to [0,1]
  n=length(Y) # sample size
  data=data.frame(Y,A,X)
  d=dim(X)[2]+1

  # To eliminate the bias caused by the heterogeneity of the data,
  # we randomly split the data into training and testing sets 20 times independently (or specified by the user).

  resul=matrix(0,nrow=repnum,ncol=200)

  set.seed(2333)
  seeds <- ceiling(runif(repnum, 10000, 1e+09))

  for(iii in 1:repnum)
  {
    set.seed(seeds[iii])
    cat("Rep",iii,"\n")

    # pre-settings

    n_train=floor(n*train_pct) # number of the traning data (70% of sample size)
    train_list=sample(nrow(data), n_train) # randomly select the training data index
    train_data=data[train_list,]

    # get the best tuning parameter C_n; default candidate list: exp(seq(-2,2,length=10))

    x=train_data[,3:(d+1)]
    a=train_data[,2]
    y=train_data[,1]
    sample=data.frame(y,a,x)
    row.names(sample)=c(1:n_train)

    info=matrix(0,nrow=10,ncol=length(c.list))
    # create the cross-validation folds, 10 folds default.
    folds=createFolds(y=c(1:n_train), k = folds_num, list = TRUE, returnTrain = T)

    for(k in 1:folds_num){

      # for ith folds
      k_data=sample[folds[[k]],]
      n=nrow(k_data)
      m=n/d

      for(c in 1:length(c.list))
      {
        gamma_n=c.list[c]*log(n)/n

        #find best partition
        p=Bel=rep(0,m)

        for(r in 1:m){
          Bel[r]=1000000
          for(l in 1:r){
            if(l==1){
              Bel_l_1=-gamma_n
            }else{
              Bel_l_1=Bel[l-1]
            }
            #calculating the distance
            if(r<m){
              rdgdata=k_data[k_data$a>=((l-1)/m)&k_data$a<(r/m),]
            }else{
              rdgdata=k_data[k_data$a>=((l-1)/m)&k_data$a<=1,]
            }

            if(length(rdgdata$y)>1){
              designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
              designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
              betaI=solve(t(designX)%*%designX+(r-l+1)/m*diag(d))%*%t(designX)%*%rdgdata$y
              dis=sum((rdgdata$y-designX%*%betaI)^2)/n+
                ((r-l+1)/m)*(1/n)*sum(betaI^2)
            }
            if(length(rdgdata$y)==1){
              designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
              designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
              betaI=solve(t(designX)%*%designX+(r-l+1)/m*diag(d))%*%t(designX)*rdgdata$y
              dis=sum((rdgdata$y-designX%*%betaI)^2)/n+
                ((r-l+1)/m)*(1/n)*sum(betaI^2)
            }
            ##############################################
            b=Bel_l_1+gamma_n+dis
            if(b<=Bel[r]){
              Bel[r]=b
              p[r]=l-1
            }
          }
        }

        ##############################################

        #segmentation from partition

        r=m
        l=p[r]
        intvl_ep=NULL
        seg_rule=NULL
        while(r>0){
          intvl_ep=c(r/m,intvl_ep)
          #calculating the linear rule
          if(r<m){
            rdgdata=k_data[k_data$a>=(l/m)&k_data$a<(r/m),]
          }else{
            rdgdata=k_data[k_data$a>=(l/m)&k_data$a<=1,]
          }

          if(length(rdgdata$y)>1){
            designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
            designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
            betaI=solve(t(designX)%*%designX+(r-l)/m*diag(d))%*%t(designX)%*%rdgdata$y
          }
          if(length(rdgdata$y)==1){
            designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
            designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
            betaI=solve(t(designX)%*%designX+(r-l)/m*diag(d))%*%t(designX)*rdgdata$y
          }
          ##############################################
          seg_rule=c(betaI,seg_rule)
          r=l
          l=p[r]
        }

        # use the left fold to calculate the LSE of the objective function for choosing the best tuning parameter.

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
          info[k,c]=info[k,c]+sum((int_data$y-designX%*%as.matrix(seg_rule[(d*i-(d-1)):(d*i)]))^2)
        }
      }
    }

    # get the optimal tuning parameter with the minimum LSE.
    cmax=c.list[min(which(apply(info,2,mean)==min(apply(info,2,mean))))]


    # use the whole training data to generate the optimal I2DR under the selected best tuning parameter.
    x=train_data[,3:(d+1)]
    a=train_data[,2]
    y=train_data[,1]
    tr_data=data.frame(y,a,x)
    n=nrow(tr_data)
    m=n/d

    gamma_n=cmax*log(n)/n

    #find best partition
    p=Bel=rep(0,m)

    for(r in 1:m){
      Bel[r]=1000000
      for(l in 1:r){
        if(l==1){
          Bel_l_1=-gamma_n
        }else{
          Bel_l_1=Bel[l-1]
        }
        #calculating the distance
        if(r<m){
          rdgdata=tr_data[tr_data$a>=((l-1)/m)&tr_data$a<(r/m),]
        }else{
          rdgdata=tr_data[tr_data$a>=((l-1)/m)&tr_data$a<=1,]
        }

        if(length(rdgdata$y)>1){
          designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
          designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
          betaI=solve(t(designX)%*%designX+(r-l+1)/m*diag(d))%*%t(designX)%*%rdgdata$y
          dis=sum((rdgdata$y-designX%*%betaI)^2)/n+
            ((r-l+1)/m)*(1/n)*sum(betaI^2)
        }
        if(length(rdgdata$y)==1){
          designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
          designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
          betaI=solve(t(designX)%*%designX+(r-l+1)/m*diag(d))%*%t(designX)*rdgdata$y
          dis=sum((rdgdata$y-designX%*%betaI)^2)/n+
            ((r-l+1)/m)*(1/n)*sum(betaI^2)
        }
        ##############################################
        b=Bel_l_1+gamma_n+dis
        if(b<=Bel[r]){
          Bel[r]=b
          p[r]=l-1
        }
      }
    }

    ##############################################

    #segmentation from partition
    r=m
    l=p[r]
    intvl_ep=NULL
    seg_rule=NULL
    while(r>0){
      intvl_ep=c(r/m,intvl_ep) # store the change point locations
      #calculating the linear rule
      if(r<m){
        rdgdata=tr_data[tr_data$a>=(l/m)&tr_data$a<(r/m),]
      }else{
        rdgdata=tr_data[tr_data$a>=(l/m)&tr_data$a<=1,]
      }

      if(length(rdgdata$y)>1){
        designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
        designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
        betaI=solve(t(designX)%*%designX+(r-l)/m*diag(d))%*%t(designX)%*%rdgdata$y
      }
      if(length(rdgdata$y)==1){
        designX=matrix(1,nrow=length(rdgdata$y),ncol=d)
        designX[,2:d]=as.matrix(rdgdata[,3:(d+1)])
        betaI=solve(t(designX)%*%designX+(r-l)/m*diag(d))%*%t(designX)*rdgdata$y
      }
      ##############################################
      seg_rule=c(betaI,seg_rule) # store the regression coefficients for each partition
      r=l
      l=p[r]
    }

    #calculate the value function V(d)
    n=dim(X)[1]
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
    M_n=mean(value)

    resu=c(cmax,M_n,length(intvl_ep),intvl_ep,seg_rule)
    resul[iii,]=c(resu,rep(0,length=200-length(resu)))
}


  ############ Output the optimal I2DR by setting the seed with highest value function ###########
  A=data[,2]
  val=resul[which(resul[,2]==max(resul[,2])),2]
  num_int=resul[which(resul[,2]==max(resul[,2])),3]
  intvl_ep=resul[which(resul[,2]==max(resul[,2])),4:(num_int+3)]*(max(A)-min(A)) +min(A)
  seg_rule=as.matrix(resul[which(resul[,2]==max(resul[,2])),(num_int+3):200])
  Beta=matrix(0,nrow=num_int,ncol=d)
  # output the partition and the corresponding regression coefficients
  cat('Interval:',min(A),intvl_ep[1],'\n')
  Beta[1,]=seg_rule[1:d]
  cat('Best rule:',Beta[1,],'\n')
  if(num_int>1){
    for(i in 2:num_int){
      cat('Interval:',intvl_ep[i-1],intvl_ep[i],r/m,'\n')
      Beta[i,]=seg_rule[((i-1)*d+1):(i*d)]
      cat('Best rule:',Beta[i,],'\n')
    }
  }

  # Output the value function

  cat('Value function:',val,'\n')

  # Use Bootstrap to contrust the Wald-type confidence interval for value function by fixing the partition found above

  intvl_ep=resul[which(resul[,2]==max(resul[,2])),4:(num_int+3)]
  data=as.matrix(data)
  data_resam=array(0,dim=c(n,ncol(data),nboots))
  beta_Boots=array(0,dim=c(nboots,length(intvl_ep),d))
  value_Boots=rep(0,nboots)
  intvls=c(0,intvl_ep)
  for(k in 1:nboots){
    index=sample(c(1:n),n,replace=TRUE)
    data_resam[,,k]=k_data=data[as.vector(index),]
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

  low_bd=val-qnorm((1+alpha)/2,0,1)*sd(value_Boots)
  up_bd=val+qnorm((1+alpha)/2,0,1)*sd(value_Boots)

  cat(alpha*100,'% Wald-type CI for the value: (',low_bd,',',up_bd,')\n')


  return(list(Partition=c(min(A),intvl_ep),Beta=Beta,Value=val,low_bd=low_bd,up_bd=up_bd))
}
