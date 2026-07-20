########### BFGS functions #############
trace = function(x){sum(diag(x))}

prior_lambda_i = function(lambda_i,t2_lambda, q){
  return((-0.5/t2_lambda*sum(lambda_i*lambda_i))-0.5*q*log(t2_lambda))
}
prior_X = function(X){
  return(-0.5*trace(t(X)%*%X))
}
prior_s2 = function(s2,nu0,s20){
  return((-1-nu0/2)*log(s2) - 0.5*nu0*s20/s2) # s2
}
## X, U
### param (vec(X: nxq),vec(Lambda: Mxq),alpha,beta,s2,t2)
joint_pos = function(param, t2_lambda, nu0, s20, conn_flat_full, # M x (n*n)
  M, N, q){

  X = matrix(param[1:(q*N)],N,q)
  U = polar_expansion(X)
  vec_lambda = param[(q*N+1):(q*(M+N))]
  Lambda_mat = matrix(vec_lambda, M, q)
  s2 = param[1+q*(M+N)]

  ## add data llk
  dat_llk = Y_llk_new(t(conn_flat_full), U, t(U), t(Lambda_mat), s2, M, N)
  p_lambda = 0
  for(i in 1:M){
    lambda_i = Lambda_mat[i,]
    p_lambda = p_lambda + prior_lambda_i(lambda_i, t2_lambda, q)
  }
  p_X = prior_X(X)
  p_s2 = prior_s2(s2,nu0,s20) # s2

  return(list(dat_llk = dat_llk, prior_lambda = p_lambda, prior_X = p_X, prior_s2 = p_s2))
}

### param (vec(X: nxq),vec(Lambda: Mxq),alpha,beta,s2,t2)
fn_X_optim = function(param, t2_lambda, nu0, s20, conn_flat_full, M, N, q){
  res_ls = joint_pos(param, t2_lambda, nu0, s20, conn_flat_full, M, N, q)
  return(sum(unlist(res_ls)))
}

################# try profile llk
Lambda_cls = function(U, conn_arr, M, q){  ###MLE of lambda wrt llk(Y|~)
  res = matrix(NA,M,q)
  for(i in 1:M){
    tmp = t(U) %*% conn_arr[,,i] %*% U
    res[i,] = diag(tmp)
  }
  return(res)
}

s2_cls = function(U,Lambda,conn_arr,M,N){
  E = matrix(0, M, N*N)
  for(i in 1:M){
    Y_mean_i = U %*% diag(Lambda[i,]) %*% t(U)
    E[i,] = c(conn_arr[,,i]-Y_mean_i)
  }
  res = sum(E^2)*2/(M*N*(N+1))
  return(res)
}

joint_pos_prof = function(param, t2_lambda, nu0, s20, conn_arr, conn_flat_full, M, N, q){ #param: vec(X)
  X = matrix(param,N,q)
  U = polar_expansion(X)

  Lambda_mat = Lambda_cls(U, conn_arr, M, q)

  s2 = s2_cls(U,Lambda_mat,conn_arr,M,N)

  ## add data llk
  dat_llk = Y_llk_new(t(conn_flat_full), U, t(U), t(Lambda_mat), s2, M, N)
  p_lambda = 0
  for(i in 1:M){
    #lambda_i = vec_lambda[((i-1)*q+1):(i*q)]
    lambda_i = Lambda_mat[i,]
    p_lambda = p_lambda + prior_lambda_i(lambda_i, t2_lambda, q)
  }
  p_X = prior_X(X)
  p_s2 = prior_s2(s2,nu0,s20) # s2

  return(list(dat_llk = dat_llk, prior_lambda = p_lambda, prior_X = p_X, prior_s2 = p_s2))
}

fn_X_optim_prof = function(param, t2_lambda, nu0, s20, conn_arr, conn_flat_full, M, N, q){
  res_ls = joint_pos_prof(param, t2_lambda, nu0, s20, conn_arr, conn_flat_full, M, N, q)
  return(sum(unlist(res_ls)))
}

gr_prior_lambda_i = function(lambda_i,t2_lambda){
  return(-1.0/t2_lambda*lambda_i)
}

gr_prior_s2 = function(s2,nu0,s20){
  return((-1-nu0/2)/s2 + nu0*s20/(2*(s2)^2))
}

gr_dat_llk = function(vec_lambda, X, s2, conn_arr, M, N, q){
  U = polar_expansion(X)
  Lambda_mat = matrix(vec_lambda, M, q)

  d_lambda_mat = matrix(NA, M, q) #M x q
  d_s2 = -0.25*M*N*(N+1)/s2
  for(i in 1:M){
    lambda_i = Lambda_mat[i,]
    d_lambda_i = -1.0/s2*(lambda_i - diag(t(U)%*%conn_arr[,,i]%*%U))
    d_lambda_mat[i,] = d_lambda_i

    d_s2 = d_s2 + 0.5/(s2^2)*norm(conn_arr[,,i]-U%*%diag(lambda_i)%*%t(U),"F")^2
  }

  return(list(d_lambda = d_lambda_mat, d_s2 = d_s2))
}

gr_optim = function(vec_lambda, X, s2, t2_lambda, nu0, s20, conn_arr, M, N, q){
  U = polar_expansion(X)
  gr_dat = gr_dat_llk(vec_lambda, X, s2, conn_arr, M, N, q)
  Lambda_mat = matrix(vec_lambda, M, q)

  gr_s2 = gr_prior_s2(s2,nu0,s20) + gr_dat$d_s2
  gr_X = eval_dlog_mat_norm_2_new(X, U = U, Ut = t(U), I_q = diag(q), dSdX = box_prod_2(diag(N), diag(q), N, N, q, q),
    I_Nq = diag(N*q), I_N = diag(N), Lambda_flat = Lambda_mat, Lambda_flat_t = t(Lambda_mat),
    Y = conn_arr, q2 = q^2, Nq = N*q, box_idx = box_prod_3(q, N), I_q_flat = matrix(c(diag(q)),1,q^2),
    s2 = s2, N = N, q = q, M = M)

  gr_lambda = matrix(NA, M, q)
  for(i in 1:M){
    lambda_i = Lambda_mat[i,]
    #print(lambda_i)
    gr_lambda[i,] = gr_prior_lambda_i(lambda_i,t2_lambda) + gr_dat$d_lambda[i,]
    #print(gr_lambda[i,])
  }
  return(list(gr_X = gr_X, gr_lambda = gr_lambda, gr_s2 = gr_s2))
}

gr_X_optim_prof = function(param, t2_lambda, nu0, s20, conn_arr, conn_flat_full, M, N, q){
  X = matrix(param,N,q)
  U = polar_expansion(X)

  Lambda_mat = Lambda_cls(U, conn_arr, M, q)
  s2 = s2_cls(U,Lambda_mat,conn_arr,M,N)

  gr = gr_optim(c(Lambda_mat), X, s2, t2_lambda, nu0, s20, conn_arr, M, N, q)
  return(c(c(gr$gr_X)))
}
