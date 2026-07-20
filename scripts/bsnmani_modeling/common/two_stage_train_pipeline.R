########################################################
##### two-stage scalar-on-network regression
########################################################
resamp = function(mean, sd, len){
  a = rep(-1,len)
  while(sum(a < 0)!=0){
    a = rnorm(n=len, mean=mean, sd=sd)
  }
  return(a)
}

fname_sanity = function(fname,dir_tags,basename_tags){
    dir_path = dirname(fname)
    err_msg = ""
    if(prod(sapply(dir_tags,FUN=function(x){grepl(x,dir_path)}))!=1){
        err_msg = sprintf("%sdirectory path: %s doesn't match experiment conditions\n", err_msg, dir_path)
    }
    
    base_fname = basename(fname)
    if(prod(sapply(basename_tags,FUN=function(x){grepl(x,base_fname)}))!=1){
        err_msg = sprintf("%sbase file name: %s doesn't match experiment conditions\n", err_msg, base_fname)
    }
    
    if(nchar(err_msg)!=0){
      cat(err_msg)
      stop()
    }
}

two_stage_single_chain_train = function(n_roi, q_val,q_val_0, n_pol, FC_dat, idx_ls, mask = TRUE, clinical_df, cov_df,   r_val, ## data
                                        noise,  t2_lambda_0, ## g1-initialize
                                        nu0, s20, eta0, t20, SGI_g1, k1, stepsize, acpt_step, target_acpt, tune, fix_ls, ## g1
                                        rho0, psi20, gamma0, kappa20, omega0, phi20, ## g2
                                        k2, SGI_g2, g2_weighted, ## MH
                                        seed, burn_in, mcmc_sample, chain_idx,
                                        save_path_g1, save_path_g2, save_path_MH,
                                        fname_suffix_g1, fname_suffix_g2, fname_suffix_MH,
                                        run_g1,run_g2,run_MH,load_g1,load_g2,load_MH,save_g1,save_g2,save_MH,
                                        dir_tags, basename_tags){## others
  g1_path = fs::path(save_path_g1,"g1",paste("g1_res",SGI_g1,"chain",chain_idx,fname_suffix_g1,sep="_"),ext="RDS")
  g2_path = fs::path(save_path_g2,"g2",paste("g2_res","chain",chain_idx,fname_suffix_g2,sep="_"),ext="RDS")
  MH_path = fs::path(save_path_MH,"MH",paste("MH_res",SGI_g1,SGI_g2,"chain",chain_idx,fname_suffix_MH,sep="_"),ext="RDS")

  message("checking g1 path...")
  fname_sanity(g1_path, c(dir_tags,"g1"), c(basename_tags,"g1"))

  message("checking g2 path...")
  fname_sanity(g2_path, c(dir_tags,"g2"), c(basename_tags,"g2"))

  message("checking MH path...")
  fname_sanity(MH_path, c(dir_tags,"MH"), c(basename_tags,"MH"))

  print("preprocessing data")
  FC_arr_mask = FC_dat
  for(i in 1:(dim(FC_arr_mask)[3])){
    diag(FC_arr_mask[,,i]) = 0
  }
  
  if(mask){
    FC_dat = FC_arr_mask
  }
  FC_flat_full = matrix(NA, n_pol, n_roi*n_roi)
  for(i in 1:n_pol){
    mat1 = FC_dat[,,i]
    FC_flat_full[i,] = c(mat1)
  }
  
  set.seed(12345)
  region_idx = sort(sapply(idx_ls,FUN=function(x){sample(x,1)}))
  n_region_pairs = length(region_idx) * (length(region_idx) - 1) /2
  
  print("initialize parameters for network model")
  FC_dat_mean = apply(FC_dat, c(1,2), mean)

  # edit by TONG !!!
  saveRDS(FC_dat_mean, file = "/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/TTTest/FC_dat_mean_Astrocyte.RDS")
  print("saved mean FC matrix for debuging")
  U_mean = eigen(FC_dat_mean)$vectors[,1:q_val]
  
  set.seed(12345)
  optim_try_1 = optim(par = c(U_mean),
                      fn = fn_X_optim_prof, gr = gr_X_optim_prof, method = "BFGS", control=list(fnscale=-1, maxit = 1e+4, trace=TRUE, reltol=1e-10),
                      t2_lambda=t2_lambda_0, nu0=nu0, s20=s20,
                      conn_arr = FC_arr_mask, conn_flat_full = FC_flat_full,
                      M = n_pol, N = n_roi, q = q_val)
  
  BFGS_init_1 = list(X = matrix(optim_try_1$par[1:(n_roi*q_val)],n_roi,q_val),
                     U = polar_expansion(matrix(optim_try_1$par[1:(n_roi*q_val)],n_roi,q_val)))
  BFGS_init_1$Lambda = Lambda_cls(BFGS_init_1$U, FC_arr_mask, n_pol, q_val)
  BFGS_init_1$s2 = s2_cls(BFGS_init_1$U, BFGS_init_1$Lambda, FC_arr_mask, n_pol, n_roi)
  
  ### initialize
  print("make sparse grid")
  sgrid_g1 = createSparseGrid(type=SGI_g1, dimension = 1, k=k1)
  sgrid_g2 = createSparseGrid(type=SGI_g2, dimension = 3, k=k2) ## test how big this value can get
  
  print("run two-stage sampling")
  
  dir.create(fs::path(save_path_g1,"g1"),recursive = TRUE, showWarnings = FALSE)
  dir.create(fs::path(save_path_g2,"g2"),recursive = TRUE, showWarnings = FALSE)
  dir.create(fs::path(save_path_MH,"MH"),recursive = TRUE, showWarnings = FALSE)
  
  set.seed(seed)
  init_ls=list()
  init_ls$X = BFGS_init_1$X + matrix(rnorm(n_roi*q_val,0,noise),n_roi,q_val)
  init_ls$Lambda = matrix(rnorm(n_pol*q_val,0,noise),n_pol,q_val)+BFGS_init_1$Lambda
  init_ls$s2 = resamp(BFGS_init_1$s2, noise, 1)[1]
  init_ls$t2_lambda = t2_lambda_0
  
  ## g1  
  if(run_g1){
    message("running stage-one sampling")
    set.seed(seed)
    start = Sys.time()
    g1_res = hybrid_MALA_g1(X_0 = init_ls$X, Lambda_0_flat = init_ls$Lambda, s2_0 = init_ls$s2, t2_lambda_0 = init_ls$t2_lambda,
                            nu0=nu0, s20=s20, eta0=eta0, t20=t20,
                            Y=FC_dat, M=n_pol, N=n_roi, q=q_val, SGI_nodes=sgrid_g1$nodes, SGI_wts=sgrid_g1$weights, weighted=(SGI_g1 == "GQN"), region_select = region_idx,
                            mcmc_sample=mcmc_sample, stepsize=stepsize, acpt_step=acpt_step, Gibbs_step=1, MALA_step=1, target_acpt=target_acpt, tune=tune, fixed=fix_ls)
    end = Sys.time()
    runtime = difftime(end,start,units="secs")
    g1_res$runtime = runtime
    g1_res$init = init_ls
    
    #print(str(g1_res))

    if(save_g1){
      message("saving stage-one sampling results")
      saveRDS(g1_res, file = fs::path(save_path_g1,"g1",paste("g1_res",SGI_g1,"chain",chain_idx,fname_suffix_g1,sep="_"),ext="RDS"))      
    }
  }else if(load_g1){
    message("loading stage-one sampling results")
    g1_res = readRDS(file = fs::path(save_path_g1,"g1",paste("g1_res",SGI_g1,"chain",chain_idx,fname_suffix_g1,sep="_"),ext="RDS"))
  }
  #g1_fname = fs::path(save_path_g1,"g1",paste("g1_res",SGI_g1,"q",q_val,"chain",chain_idx,fname_suffix_g1,sep="_"),ext="RDS")
  #cat(sprintf("export G1_FNAME=%s\n",g1_fname))

  ## g2
  set.seed(seed)
  if(run_g2){
    message("running stage-two sampling")
    Lambda_pos_mean = apply(g1_res$mcmc$Lambda_flat[,,(burn_in+1):mcmc_sample],c(1,2),mean)
    LM_dat = Reduce(cbind,list(clinical_df, cov_df, Lambda_pos_mean))
    #LM_dat = Reduce(cbind,list(clinical_df, cov_df, init_ls$Lambda))
    colnames(LM_dat) = c("Y", paste("var",1:(r_val-1),sep="_"), paste("lambda",1:q_val,sep="_"))
    LM_dat = as.data.frame(LM_dat)
    model = lm(Y ~ ., data=LM_dat)
  
    init_ls=list()
    init_ls$beta = unname(model$coefficients[(r_val+1):(r_val+q_val)])
    init_ls$alpha = unname(model$coefficients[1:r_val])
    init_ls$t2 = var(model$residuals)
    init_ls$t2_beta = var(init_ls$beta)
    init_ls$t2_alpha = var(init_ls$alpha)
  
    start = Sys.time()
    g2_res = g2(beta_0 = init_ls$beta, alpha_0 = init_ls$alpha, t2_0 = init_ls$t2, t2_beta_0 = init_ls$t2_beta, t2_alpha_0 = init_ls$t2_alpha, 
                Lambda_flat0 = g1_res$mcmc$Lambda_flat[,,1], ##initial values
                rho0=rho0, psi20=psi20, gamma0=gamma0, kappa20=kappa20, omega0=omega0, phi20=phi20, ##paras
                Lambda_flat_cube=g1_res$mcmc$Lambda_flat, C=matrix(clinical_df), Z=as.matrix(cbind(1,cov_df)), M=n_pol, N=n_roi, q=q_val, r=r_val, ###data
                mcmc_sample = mcmc_sample, fixed=fix_ls) ## additional parameters
    end = Sys.time()
    runtime = difftime(end,start,units="secs")
    g2_res$runtime = runtime
    g2_res$init = init_ls
    
    if(save_g2){    
      message("saving stage-two sampling results")  
      saveRDS(g2_res, file = fs::path(save_path_g2,"g2",paste("g2_res","chain",chain_idx,fname_suffix_g2,sep="_"),ext="RDS"))
    }
  }else if(load_g2){
    message("loading stage-two sampling results")
    g2_res = readRDS(fs::path(save_path_g2,"g2",paste("g2_res","chain",chain_idx,fname_suffix_g2,sep="_"),ext="RDS"))
  }
  #g2_fname = fs::path(save_path_g2,"g2",paste("g2_res","q",q_val,"chain",chain_idx,fname_suffix_g2,sep="_"),ext="RDS")
  #cat(sprintf("export G2_FNAME=%s\n",g2_fname))

  if(run_MH){
    message("running Metropolis Hastings sampling")  
    start = Sys.time()
    MH_sampling = A_lambda_MH(Lambda_flat=g1_res$mcmc$Lambda_flat, U=g1_res$mcmc$U, X=g1_res$mcmc$X, s2=g1_res$mcmc$s2, t2_lambda=g1_res$mcmc$t2_lambda, 
                              d=g2_res$mcmc$d, t2=g2_res$mcmc$t2, t2_beta=g2_res$mcmc$t2_beta, t2_alpha=g2_res$mcmc$t2_alpha, 
                              Y_llk=g1_res$mcmc$dat_llk, C_llk=g2_res$mcmc$dat_llk,
                              rho0=rho0, psi20=psi20, gamma0=gamma0, kappa20=kappa20, omega0=omega0, phi20=phi20,
                              SGI_nodes=sgrid_g2$nodes, SGI_wts=sgrid_g2$weights, g1_samples=mcmc_sample, 
                              C=matrix(clinical_df), Z=as.matrix(cbind(1,cov_df)), 
                              M=n_pol, N=n_roi, q=q_val,q_0 = q_val_0, r=r_val, scheme=SGI_g2, weighted=g2_weighted)
    end = Sys.time()
    runtime = difftime(end,start,units="secs")
    MH_sampling$runtime = runtime 
    MH_sampling$total_runtime = runtime + g2_res$runtime + g1_res$runtime
    MH_sampling$init_g1 = g1_res$init
    MH_sampling$init_g2 = g2_res$init
    
    if(save_MH){
      message("saving Metropolis Hastings sampling results")
      saveRDS(MH_sampling, file = fs::path(save_path_MH,"MH",paste("MH_res",SGI_g1,SGI_g2,"chain",chain_idx,fname_suffix_MH,sep="_"),ext="RDS"))      
    }
  }else if(load_MH){
    message("loading existing Metropolis Hastings sampling results")
    MH_sampling =  readRDS(fs::path(save_path_MH,"MH",paste("MH_res",SGI_g1,SGI_g2,"chain",chain_idx,fname_suffix_MH,sep="_"),ext="RDS"))
  }
  #MH_fname = fs::path(save_path_MH,"MH",paste("MH_res",SGI_g1,SGI_g2,"q",q_val,"chain",chain_idx,fname_suffix_MH,sep="_"),ext="RDS")
  #cat(sprintf("export MH_FNAME=%s\n",MH_fname))
  
  #return(list(mcmc = MH_sampling))
}