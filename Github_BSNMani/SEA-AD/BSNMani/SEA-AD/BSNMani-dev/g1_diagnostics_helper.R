library(coda)

procrustes_dist = function(X,Y,type="full"){
  return(procdist(X,Y,type))
}

angular_dist=function(X,Y,clip=TRUE){
  k = ncol(X)
  cos_angle = 1/k*sum(diag(t(X)%*%Y))
  cos_angle = pmin(pmax(cos_angle, -1), 1)
  return(acos(cos_angle))
}

geodesic_dist = function(X,Y,clip=FALSE){
  cos_theta_matrix = t(X) %*% Y
  
  # Singular value decomposition to find cosines of the principal angles
  svd_result = svd(cos_theta_matrix)
  cos_thetas = svd_result$d
  #print(cos_thetas)
  
  # Ensure the values are within the correct range before applying acos
  if(clip){
    cos_thetas = pmin(pmax(cos_thetas, -1), 1)
  }
  
  # Calculate the principal angles
  thetas = acos(cos_thetas)
  
  # Calculate the geodesic distance
  geodesic_distance = sqrt(sum(thetas^2))
  
  res = list(svd=svd_result,principal_angles=thetas,geodesic_distance=geodesic_distance)
  return(res)
}

gelman.diag.v2 = function (x, autoburnin = TRUE){
  x <- as.mcmc.list(x)
  
  if (nchain(x) < 2) 
    stop("You need at least two chains")
  if (autoburnin && start(x) < end(x)/2) 
    x <- window(x, start = end(x)/2 + 1)
  
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)
  x <- lapply(x, as.matrix)
  #S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, Nvar, Nchain))
  #W <- apply(S2, c(1, 2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), nrow = Nvar, ncol = Nchain)
  #B <- Niter * var(t(xbar))
  
  #w <- diag(W)
  #b <- diag(B)
  b = Niter * apply(t(xbar),2,var)
  s2 = sapply(x,FUN=function(y){apply(y,2,var)},simplify=TRUE)
  w = rowMeans(s2)
  # s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/Nchain
  var.b <- (2 * b^2)/(Nchain - 1)
  var_1 = sapply(1:Nvar,FUN=function(i){var(t(s2)[,i],t(xbar^2)[,i])})
  var_2 = sapply(1:Nvar,FUN=function(i){var(t(s2)[,i],t(xbar)[,i])})
  cov.wb <- (Niter/Nchain) * (var_1 - 2 * muhat * var_2)
  #cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * muhat * var(t(s2), t(xbar)))
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate))
  dimnames(psrf) <- list(xnames, c("Point est."))
  out <- list(psrf = psrf)
  class(out) <- "gelman.diag"
  out
}


## results from one chain
## res
#  --$ mcmc
#    --$ par 1
#    --$ par 2 ...
#  --$ stepsize

combine_epochs = function(mcmc_res_ls_1, mcmc_res_ls_2){ ## mcmc_res_ls_2 must come after mcmc_res_ls_1
  res = list()
  res$stepsize = mcmc_res_ls_2$stepsize
  res$runtime = mcmc_res_ls_1$runtime + mcmc_res_ls_2$runtime

  mcmc_ls = list()
  mcmc_ls$Lambda_flat = abind(mcmc_res_ls_1$mcmc$Lambda_flat, mcmc_res_ls_2$mcmc$Lambda_flat, along=3)
  mcmc_ls$X = abind(mcmc_res_ls_1$mcmc$X, mcmc_res_ls_2$mcmc$X, along=3)
  mcmc_ls$U = abind(mcmc_res_ls_1$mcmc$U, mcmc_res_ls_2$mcmc$U, along=3)
  mcmc_ls$s2 = rbind(mcmc_res_ls_1$mcmc$s2, mcmc_res_ls_2$mcmc$s2)
  mcmc_ls$t2_lambda = rbind(mcmc_res_ls_1$mcmc$t2_lambda, mcmc_res_ls_2$mcmc$t2_lambda)
  mcmc_ls$llk = rbind(mcmc_res_ls_1$mcmc$llk, mcmc_res_ls_2$mcmc$llk)
  mcmc_ls$dat_llk = rbind(mcmc_res_ls_1$mcmc$dat_llk, mcmc_res_ls_2$mcmc$dat_llk)
  mcmc_ls$llk_grad = cbind(mcmc_res_ls_1$mcmc$llk_grad, mcmc_res_ls_2$mcmc$llk_grad)
  mcmc_ls$subnet_recon = cbind(mcmc_res_ls_1$mcmc$subnet_recon, mcmc_res_ls_2$mcmc$subnet_recon)
  mcmc_ls$acceptance_rate = rbind(mcmc_res_ls_1$mcmc$acceptance_rate, mcmc_res_ls_2$mcmc$acceptance_rate)
  res$mcmc = mcmc_ls

  return(res)
}

# save as a list of parameters values, each element is a list of length(#chains)
# dat: M, N, q, Y_flat_full, n_samps, n_burnin, n_chains
gather_samples = function(mcmc_res_ls,dat,pars,var_df,sign_flip=TRUE){
  samps_ls = list()
  q = dat$q
  n_samps = dat$n_samps
  N = dat$N

  for(i in 1:nrow(var_df)){
    var = var_df$var[i]
    if(var == "acceptance_rate"){
      samps_ls[["acpt_ratio"]] = lapply(mcmc_res_ls, FUN=function(x){mcmc(x$mcmc[[var]][,2])})
      samps_ls[["stepsize"]] = lapply(mcmc_res_ls, FUN=function(x){mcmc(x$mcmc[[var]][,3])})
    }else if(var_df$dim[i]==1){
      samps_ls[[var]] = lapply(mcmc_res_ls, FUN=function(x){mcmc(x$mcmc[[var]])})
    }else if(var_df$dim[i]==2){
      samps_ls[[var]] = lapply(mcmc_res_ls, FUN=function(x){mcmc(t(x$mcmc[[var]]))})
    }else if(var_df$dim[i]==3){
      if(var == "U" & sign_flip){
        U_ls = lapply(mcmc_res_ls,FUN=function(x){x$mcmc$U})
        U_ls_flipped = U_ls
        ones_mat = matrix(1,dat$N,1)
        for(i in 2:dat$n_chains){
          flip = apply(-1.0*U_ls_flipped[[i]] - U_ls[[1]],2,FUN=function(x){mean(x^2)})
          orig = apply(U_ls[[i]] - U_ls[[1]],2,FUN=function(x){mean(x^2)})
          sign_flip = as.numeric(flip<orig)
          for(k in 1:q){
            sign_flip[k] = ifelse(sign_flip[k]==0,1,-1)
          }
          U_ls_flipped[[i]] = array(ones_mat %*% matrix(sign_flip,1,q),dim=c(N,q,n_samps)) * U_ls[[i]]
          print(as.numeric(flip<orig))
        }
        samps_ls[[var]] = lapply(U_ls_flipped, FUN=function(x){mcmc(t(apply(x,3,c)))})
      }else{
        samps_ls[[var]] = lapply(mcmc_res_ls, FUN=function(x){mcmc(t(apply(x$mcmc[[var]],3,c)))})
      }
    }
  }
  return(samps_ls)
}

subset_samps = function(samps_ls, n_samps, n_burnin){
  samps_subset = lapply(samps_ls,FUN=function(x){
    lapply(x,FUN=function(y){if(is.null(dim(y))){
      mcmc(y[(n_burnin+1):n_samps])
    }else{
      mcmc(y[(n_burnin+1):n_samps,])
    }})
  })
  return(samps_subset)
}

posterior_mean = function(samps_subset){
  mcmc_summary_ls = lapply(samps_subset,FUN=function(x){
    summary(mcmc.list(x))$statistics})
  mean_ls = lapply(mcmc_summary_ls,FUN=function(x){if(is.null(dim(x))){return(x["Mean"])}else{return(x[,"Mean"])}})
  return(mean_ls)
}

DIC_net = function(mean_ls, Y_flat, M, N, q){
  U_mean = polar_expansion(matrix(mean_ls$X,N,q))
  D_bar = -2*Y_llk_new(t(Y_flat), U_mean, t(U_mean), t(matrix(mean_ls$Lambda,M,q)), mean_ls$s2, M, N)
  bar_D = -2*mean_ls$dat_llk
  pDIC = bar_D - D_bar
  DIC = D_bar +2*pDIC

  return(DIC)
}

DIC = function(mean_ls, Y_flat, C, Z, M, N, q,r){
  U_mean = polar_expansion(matrix(mean_ls$X,N,q))
  D_bar = -2*(C_llk_new(C, c(mean_ls$t2), matrix(mean_ls$Lambda,M,q)%*%mean_ls$d[1:q], Z%*%mean_ls$d[(q+1):(q+r)], M) + Y_llk_new(t(Y_flat), U_mean, t(U_mean), t(matrix(mean_ls$Lambda,M,q)), mean_ls$s2, M, N))
  bar_D = -2*mean_ls$dat_llk
  pDIC = bar_D - D_bar
  DIC = D_bar +2*pDIC

  return(DIC)
}

DIC_joint = function(mean_ls, Y_flat, C, Z, M, N, q,r){
  U_mean = polar_expansion(matrix(mean_ls$X,N,q))
  D_bar = -2*(Y_C_llk_new(t(Y_flat), U_mean, t(U_mean), t(matrix(mean_ls$Lambda,M,q)), c(mean_ls$s2),
                            C, mean_ls$d[1:q], mean_ls$d[(q+1):(q+r)], Z, c(mean_ls$t2), 
                            matrix(mean_ls$Lambda,M,q)%*%mean_ls$d[1:q], Z%*%mean_ls$d[(q+1):(q+r)], M, N))
  bar_D = -2*mean_ls$dat_llk
  pDIC = bar_D - D_bar
  DIC = D_bar +2*pDIC

  return(DIC)
}

Rhat = function(samps_subset, auto.burn=FALSE, ret_ls=FALSE){
  Rhat_ls = lapply(samps_subset,FUN=function(x){gelman.diag(mcmc.list(x),autoburnin=auto.burn,multivariate=FALSE)})
  Rhat_conv = sapply(Rhat_ls, FUN=function(x){mean(x$psrf[,1])})
  Rhat_df = data.frame(par = names(samps_subset), mean = Rhat_conv)
    if(ret_ls){
        return(Rhat_ls)
     }else{return(Rhat_df)}  
}

Rhat.v2 = function(samps_subset, auto.burn=FALSE){
  Rhat_ls = lapply(samps_subset,FUN=function(x){gelman.diag.v2(mcmc.list(x),autoburnin=auto.burn)})
  return(Rhat_ls)
  #Rhat_conv = sapply(Rhat_ls, FUN=function(x){mean(x$psrf[,1])})
  #Rhat_df = data.frame(par = names(samps_subset), mean = Rhat_conv)
  #if(ret_ls){
  #  return(Rhat_ls)
  #}else{return(Rhat_df)}  
}

ess = function(samps_subset){
  ess_ls = lapply(samps_subset, FUN=function(x){effectiveSize(mcmc.list(x))})
  ess_sum = sapply(ess_ls, mean)
  ess_df = data.frame(par = names(samps_subset), mean=ess_sum)

  return(ess_df)
}

geweke = function(samps_subset,chain_idx){
  gg_ls = lapply(samps_subset, FUN=function(x){geweke.diag(x[[chain_idx]])$z})
  gg_sum = sapply(gg_ls, mean)
  gg_df = data.frame(par = names(samps_subset), mean=gg_sum)

  return(gg_df)
}

RMSE = function(samps,truth){
  Nvar = length(samps)
  RMSE_ls = lapply(1:Nvar,FUN=function(x){
    sqrt(mean(sapply(samps[[x]],FUN=function(y){
      if(!is.null(dim(y))){
        apply(y,1,FUN=function(z){mean((z-truth[[x]])^2)})
      }else{
        mean((y-truth[[x]])^2)
      }
    })))
  })
  names(RMSE_ls) = names(samps)
  return(RMSE_ls)
}

RMSE_norm = function(samps,truth){
  Nvar = length(samps)
  RMSE_ls = lapply(1:Nvar,FUN=function(x){
    sqrt(mean(sapply(samps[[x]],FUN=function(y){
      if(!is.null(dim(y))){
        apply(y,1,FUN=function(z){mean((z-truth[[x]])^2)/mean(truth[[x]]^2)})
      }else{
        mean((y-truth[[x]])^2)/mean(truth[[x]]^2)
      }
    })))
  })
  names(RMSE_ls) = names(samps)
  return(RMSE_ls)
}

SE = function(samps,truth){
  Nvar = length(samps)
  SE_ls = lapply(1:Nvar,FUN=function(x){
    sapply(samps[[x]],FUN=function(y){
      if(!is.null(dim(y))){
        apply(y,1,FUN=function(z){mean((z-truth[[x]])^2)/mean(truth[[x]]^2)})
      }else{
        ((y-truth[[x]])^2)/mean(truth[[x]]^2)
      }
    })
  })
  names(SE_ls) = names(samps)
  return(SE_ls)
}

subnetwork_heatmap = function(idx_ref, U, N, q, col.side, col.heat){
  modules = names(idx_ref)
  modules_sz = sapply(idx_ref,FUN=function(x){length(x)})
  idx_vec = unname(unlist(idx_ref))
  U_arrange = U[idx_vec,]
  borders = cumsum(modules_sz)
  borders = borders[-length(borders)]
  borders = borders
  print(paste("subnetwork_heatmap_v2: ", dim(U_arrange)[1], "x", dim(U_arrange)[2]))   ## Modified by Tong Liu on 3.1.2025

  #par(mfrow=c(ceiling(q/2),2),mai=c(0.3,0.3,0.3,0.3))
  for(i in 1:q){
    subnet_q = U_arrange[,i,drop=FALSE] %*% t(U_arrange[,i,drop=FALSE])
    heatmap.2(subnet_q, col=col.heat, scale="none",dendrogram='none',main=paste(i),lhei=c(2,5),
      lwid=c(2,5), Rowv=FALSE, Colv=FALSE, trace='none', colsep=borders, rowsep=borders, sepcolor="black",
      RowSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},col.side,modules_sz))),
      ColSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},col.side,modules_sz))))
    legend("left",legend = modules, col = col.side, lty= 1, lwd = 5, cex=.7)
  }
}

subnetwork_heatmap_v2 = function(idx_ref, U_sub, col.side, col.heat, ignore.diag=TRUE, title, breaks_vec, na.color, legend_position="left",lhei_ratio=c(2,5),lwid_ratio=c(2,5)){
  modules = names(idx_ref)
  modules_sz = sapply(idx_ref,FUN=function(x){length(x)})
  idx_vec = unname(unlist(idx_ref))
  U_arrange = U_sub[idx_vec,idx_vec]
  borders = cumsum(modules_sz)
  borders = borders[-length(borders)]
  borders = borders
  print(paste("subnetwork_heatmap_v2: ", dim(U_arrange)[1], "x", dim(U_arrange)[2]))  ## Modified by Tong Liu on 3.1.2025

  if(ignore.diag){
    diag(U_arrange) = NA
  }

  heatmap.2(U_arrange, col=col.heat, scale="none",dendrogram='none',main=title,lhei=lhei_ratio,breaks = breaks_vec, na.color=na.color,
    lwid=lhei_ratio, Rowv=FALSE, Colv=FALSE, trace='none', colsep=borders, rowsep=borders, sepcolor="black",
    RowSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},col.side,modules_sz))),
    ColSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},col.side,modules_sz))))
  legend(legend_position,legend = modules, col = col.side, lty= 1, lwd = 5, cex=.7, border=NA)
}

make_perm=function(q){
  library(combinat)
  res = Reduce(rbind,permn(1:q))
  return(res)
}
