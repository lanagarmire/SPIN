########################################
####### get name of each edge in connectivity matrix
########################################

## n.rows: number of rows in a matrix
## n.cols: number of cols in a matrix

edge_name = function(n.rows,n.cols){
  rows = 1:n.rows
  cols = 1:n.cols
  return(outer(rows,cols,FUN=function(i,j){paste("edge_",i,"--",j,sep="")}))
}

########################################
####### SVM for scalar-on-network regression
########################################

## train_dat: list of training data
##   $Y_ls: list of brain connectivity matrix
##   $C_df: matrix of scalar outcome
##   $M: number of training samples
## test_dat: list of testing data
##   $Y_ls: list of brain connectivity matrix
##   $C_df: matrix of scalar outcome
##   $M: number of training samples
## N: number of region of interests in brain connectivity network
## kern: kernel to use for SVR
## edge_names: name of each edge

SVR_scalar_on_network = function(train_dat,test_dat=NULL,edge_features=TRUE,N,kern="linear",edge_names=NULL,cachesize=40){
  library(e1071)
  
  ## prepare training data
  train_dat_df = cbind(train_dat$C_df,train_dat$cov)
  colnames(train_dat_df) = c("outcome",colnames(train_dat$cov))
  
  if(edge_features){
    train_edge_df = matrix(NA,train_dat$M,N*(N-1)/2)
    for(i in 1:train_dat$M){
      train_edge_df[i,] = train_dat$Y_ls[[i]][lower.tri(train_dat$Y_ls[[i]],diag=FALSE)]
    }
    colnames(train_edge_df) = edge_names
    train_dat_df = cbind(train_dat_df, train_edge_df)
  }
  
  ## fit SVR
  model = svm(outcome ~., data = train_dat_df, kernel=kern, scale=FALSE, cachesize=cachesize)
  
  ## test
  if(!is.null(test_dat)){
    test_dat_df = test_dat$cov
    
    if(edge_features){
      test_edge_df = matrix(NA,test_dat$M,N*(N-1)/2)
      for(i in 1:test_dat$M){
        test_edge_df[i,] = test_dat$Y_ls[[i]][lower.tri(test_dat$Y_ls[[i]],diag=FALSE)]
      }
      colnames(test_edge_df) = edge_names
      test_dat_df = cbind(test_dat_df, test_edge_df)
    }
    
    C_pred = predict(model,test_dat_df)
    pred_MSE = mean((C_pred - test_dat$C_df)^2)
    pred = data.frame(pred = C_pred)
    pred_RMSE = list(pred=sqrt(pred_MSE),train=sqrt(mean((predict(model,train_dat_df) - train_dat$C_df)^2)))
    pred_RMSE_norm = list(pred = sqrt(pred_MSE/var(test_dat$C_df)*(test_dat$M-1)/test_dat$M), 
                          train = sqrt(pred_RMSE$train)/sd(train_dat$C_df)*sqrt((nrow(train_dat_df)-1)/nrow(train_dat_df)))
  }else{
    C_pred = predict(model,train_dat_df)
    pred_MSE = mean((C_pred - train_dat$C_df)^2)
    pred = data.frame(pred = C_pred)
    pred_RMSE = sqrt(pred_MSE)
    pred_RMSE_norm = sqrt(pred_MSE/var(train_dat$C_df)*(nrow(train_dat_df)-1)/nrow(train_dat_df))
  }
  
  return(list(train_mod = model, pred = pred, pred_RMSE = pred_RMSE, pred_RMSE_norm = pred_RMSE_norm))
}

########################################
####### lasso for scalar-on-network regression
########################################

## train_dat: list of training data
##   $Y_ls: list of brain connectivity matrix
##   $C_df: matrix of scalar outcome
##   $M: number of training samples
## test_dat: list of testing data
##   $Y_ls: list of brain connectivity matrix
##   $C_df: matrix of scalar outcome
##   $M: number of training samples
## N: number of region of interests in brain connectivity network
## edge_names: name of each edge

lasso_scalar_on_network = function(train_dat,test_dat=NULL,edge_features=TRUE, N, edge_names=NULL, fix_clinical=TRUE,nested_tune=FALSE,lambda=NULL,cov_threshold=0.001,n.folds){
  library(glmnet)
  
  ## prepare training data
  train_dat_df = cbind(train_dat$C_df,train_dat$cov)
  colnames(train_dat_df) = c("outcome",colnames(train_dat$cov))
  num_cov = ncol(train_dat$cov)
  
  if(edge_features){
    train_edge_df = matrix(NA,train_dat$M,N*(N-1)/2)
    for(i in 1:train_dat$M){
      train_edge_df[i,] = train_dat$Y_ls[[i]][lower.tri(train_dat$Y_ls[[i]],diag=FALSE)]
    }
    colnames(train_edge_df) = edge_names
    train_dat_df = cbind(train_dat_df, train_edge_df)
  }
  
  if(fix_clinical){
    if(nested_tune){
      cv_las_fit = cv.glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=1,family='gaussian',intercept=TRUE, nfolds=n.folds, penalty.factor=c(rep(0,num_cov),rep(1,ncol(train_edge_df))))
      las_fit = glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=1,family='gaussian',intercept=TRUE,lambda=cv_las_fit$lambda.min, penalty.factor=c(rep(0,num_cov),rep(1,ncol(train_edge_df))))
    }else{
      las_fit = glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=1,family='gaussian',intercept=TRUE,lambda=lambda, penalty.factor=c(rep(0,num_cov),rep(1,ncol(train_edge_df))))
    }
  }else{
    if(nested_tune){
      cv_las_fit = cv.glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=1,family='gaussian',intercept=TRUE, nfolds=n.folds)
      las_fit = glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=1,family='gaussian',intercept=TRUE,lambda=cv_las_fit$lambda.min)
    }else{
      las_fit = glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=1,family='gaussian',intercept=TRUE,lambda=lambda)
    }
  }
  lasso_coefs = coef(las_fit)
  cov_sds = apply(train_dat_df[,-1],2,sd)
  if(fix_clinical){
    effect_sizes = cov_sds[(num_cov+1):length(cov_sds)]*abs(lasso_coefs[(num_cov+2):length(lasso_coefs)])
    non_zero_idx = which(effect_sizes > cov_threshold*(max(train_dat_df[,1])-min(train_dat_df[,1])))
    non_zero_idx = c(1:(num_cov+1),non_zero_idx+num_cov+1)
  }else{
    effect_sizes = cov_sds*abs(lasso_coefs[-1])
    non_zero_idx = which(effect_sizes > cov_threshold*(max(train_dat_df[,1])-min(train_dat_df[,1])))
    non_zero_idx = c(1,non_zero_idx+1)
  }
  ## lasso_nonzero = data.frame(features = lasso_coefs@Dimnames[[1]][ which(lasso_coefs != 0 ) ],
  ##                            coefs    = lasso_coefs              [ which(lasso_coefs != 0 ) ])
  lasso_nonzero = data.frame(features = lasso_coefs@Dimnames[[1]][ non_zero_idx ],
                             coefs    = lasso_coefs              [ non_zero_idx ])
  
  ## test
  if(!is.null(test_dat)){
    test_dat_df = test_dat$cov
    
    if(edge_features){
      test_edge_df = matrix(NA,test_dat$M,N*(N-1)/2)
      for(i in 1:test_dat$M){
        test_edge_df[i,] = test_dat$Y_ls[[i]][lower.tri(test_dat$Y_ls[[i]],diag=FALSE)]
      }
      colnames(test_edge_df) = edge_names
      test_dat_df = cbind(test_dat_df, test_edge_df)
    }
    
    ## C_pred = predict(las_fit, newx = test_dat_df)
    C_pred = cbind(1,test_dat_df[,non_zero_idx[-1]-1])%*%lasso_nonzero$coef
    pred_MSE = mean((C_pred - test_dat$C_df)^2)
    pred = data.frame(pred=C_pred)
    pred_RMSE = list(pred=sqrt(pred_MSE),train=sqrt(mean((predict(las_fit, newx = train_dat_df[,2:ncol(train_dat_df)]) - train_dat$C_df)^2)))
    pred_RMSE_norm = list(pred = sqrt(pred_MSE/var(test_dat$C_df)*(test_dat$M-1)/test_dat$M), 
                          train = sqrt(pred_RMSE$train)/sd(train_dat$C_df)*sqrt((nrow(train_dat_df)-1)/nrow(train_dat_df)))
  }else{
    C_pred = predict(las_fit, newx = train_dat_df[,2:ncol(train_dat_df)])
    pred_MSE = mean((C_pred - train_dat$C_df)^2)
    pred = data.frame(pred = C_pred)
    pred_RMSE = sqrt(pred_MSE)
    pred_RMSE_norm = sqrt(pred_MSE/var(train_dat$C_df)*(nrow(train_dat_df)-1)/nrow(train_dat_df))
  }
  
  if(nested_tune){
    return(list(train_mod = las_fit, lasso_nonzero = lasso_nonzero, cvfit = cv_las_fit, lambda=lambda,
                pred = pred, pred_RMSE = pred_RMSE, pred_RMSE_norm = pred_RMSE_norm))
  }else{
    return(list(train_mod = las_fit, lasso_nonzero = lasso_nonzero, lambda=lambda,
                pred = pred, pred_RMSE = pred_RMSE, pred_RMSE_norm = pred_RMSE_norm))
  }
}

########################################
####### elastic net for scalar-on-network regression
########################################

## train_dat: list of training data
##   $Y_ls: list of brain connectivity matrix
##   $C_df: matrix of scalar outcome
##   $M: number of training samples
## test_dat: list of testing data
##   $Y_ls: list of brain connectivity matrix
##   $C_df: matrix of scalar outcome
##   $M: number of training samples
## N: number of region of interests in brain connectivity network
## edge_names: name of each edge
## alpha_vec: alpha values to use for CV

elastic_net_scalar_on_network = function(train_dat,test_dat=NULL,edge_features=TRUE, N, edge_names, fix_clinical=TRUE, nested_tune=FALSE,lambda=NULL,alpha=NULL,cov_threshold=0.001){
  library(glmnet)
  library(doParallel)
  
  ## train
  train_dat_df = cbind(train_dat$C_df,train_dat$cov)
  colnames(train_dat_df) = c("outcome",colnames(train_dat$cov))
  num_cov = ncol(train_dat$cov)
  
  if(edge_features){
    train_edge_df = matrix(NA,train_dat$M,N*(N-1)/2)
    for(i in 1:train_dat$M){
      train_edge_df[i,] = train_dat$Y_ls[[i]][lower.tri(train_dat$Y_ls[[i]],diag=FALSE)]
    }
    colnames(train_edge_df) = edge_names
    train_dat_df = cbind(train_dat_df, train_edge_df)
  }
  
  if(fix_clinical){
    if(nested_tune){
      alpha_vec=seq(0.05, 0.95, 0.05)
      search = foreach(i = alpha_vec, .combine = rbind) %dopar% {
        cv = cv.glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],family='gaussian', parallel = TRUE, intercept=TRUE, alpha = i,penalty.factor=c(rep(0,num_cov),rep(1,ncol(train_edge_df))))
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
      }
      
      cv_ela = search[search$cvm == min(search$cvm), ]
      ela_fit=glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=cv_ela$alpha,lambda=cv_ela$lambda.min,family='gaussian',intercept=TRUE,penalty.factor=c(rep(0,num_cov),rep(1,ncol(train_edge_df))))
    }else{
      ela_fit=glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=alpha,lambda=lambda,family='gaussian',intercept=TRUE,penalty.factor=c(rep(0,num_cov),rep(1,ncol(train_edge_df))))
    }
  }else{
    if(nested_tune){
      alpha_vec=seq(0.05, 0.95, 0.05)
      search = foreach(i = alpha_vec, .combine = rbind) %dopar% {
        cv = cv.glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],family='gaussian', parallel = TRUE, intercept=TRUE, alpha = i)
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
      }
      
      cv_ela = search[search$cvm == min(search$cvm), ]
      ela_fit=glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=cv_ela$alpha,lambda=cv_ela$lambda.min,family='gaussian',intercept=TRUE)
    }else{
      ela_fit=glmnet(train_dat_df[,2:ncol(train_dat_df)],train_dat_df[,1],alpha=alpha,lambda=lambda,family='gaussian',intercept=TRUE)
    }
  }
  ela_coefs = coef(ela_fit)
  cov_sds = apply(train_dat_df[,-1],2,sd)
  if(fix_clinical){
    effect_sizes = cov_sds[(num_cov+1):length(cov_sds)]*abs(ela_coefs[(num_cov+2):length(ela_coefs)])
    non_zero_idx = which(effect_sizes > cov_threshold*(max(train_dat_df[,1])-min(train_dat_df[,1])))
    non_zero_idx = c(1:(num_cov+1),non_zero_idx+num_cov+1)
  }else{
    effect_sizes = cov_sds*abs(ela_coefs[-1])
    non_zero_idx = which(effect_sizes > cov_threshold*(max(train_dat_df[,1])-min(train_dat_df[,1])))
    non_zero_idx = c(1,non_zero_idx+1)
  }
  
  ## ela_nonzero = data.frame(features = ela_coefs@Dimnames[[1]][ which(ela_coefs != 0 ) ],
  ##                          coefs    = ela_coefs              [ which(ela_coefs != 0 ) ])
  ela_nonzero = data.frame(features = ela_coefs@Dimnames[[1]][ non_zero_idx ],
                             coefs    = ela_coefs              [ non_zero_idx ])
  
  ## test
  if(!is.null(test_dat)){
    test_dat_df = test_dat$cov
    
    if(edge_features){
      test_edge_df = matrix(NA,test_dat$M,N*(N-1)/2)
      for(i in 1:test_dat$M){
        test_edge_df[i,] = test_dat$Y_ls[[i]][lower.tri(test_dat$Y_ls[[i]],diag=FALSE)]
      }
      colnames(test_edge_df) = edge_names
      test_dat_df = cbind(test_dat_df, test_edge_df)
    }
    
    ##C_pred = predict(ela_fit,newx=test_dat_df)
    C_pred = cbind(1,test_dat_df[,non_zero_idx[-1]-1])%*%ela_nonzero$coef
    pred_MSE = mean((C_pred - test_dat$C_df)^2)
    pred = data.frame(pred=C_pred)
    pred_RMSE = list(pred=sqrt(pred_MSE),train=sqrt(mean((predict(ela_fit,newx=train_dat_df[,2:ncol(train_dat_df)]) - train_dat$C_df)^2)))
    pred_RMSE_norm = list(pred = sqrt(pred_MSE/var(test_dat$C_df)*(test_dat$M-1)/test_dat$M), 
                          train = sqrt(pred_RMSE$train)/sd(train_dat$C_df)*sqrt((nrow(train_dat_df)-1)/nrow(train_dat_df)))
  }else{
    C_pred = predict(ela_fit,newx=train_dat_df[,2:ncol(train_dat_df)])
    pred_MSE = mean((C_pred - train_dat$C_df)^2)
    pred = data.frame(pred = C_pred)
    pred_RMSE = sqrt(pred_MSE)
    pred_RMSE_norm = sqrt(pred_MSE/var(train_dat$C_df)*(nrow(train_dat_df)-1)/nrow(train_dat_df))
  }
  
  if(nested_tune){
    return(list(train_mod = ela_fit, ela_nonzero = ela_nonzero, cvfit = search, pred = pred, pred_RMSE = pred_RMSE, pred_RMSE_norm = pred_RMSE_norm))
  }else{
    return(list(train_mod = ela_fit, ela_nonzero = ela_nonzero, pred = pred, pred_RMSE = pred_RMSE, pred_RMSE_norm = pred_RMSE_norm))
  }
}

