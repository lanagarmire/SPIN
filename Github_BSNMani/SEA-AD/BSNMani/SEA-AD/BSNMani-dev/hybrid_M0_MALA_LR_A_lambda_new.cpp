// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]

#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double log_A_lambda_integrand(const double& t2, const double& t2_alpha, const double& t2_beta, const double& gamma0, const double& kappa20, const int& M, const int& q, const double& rho0, const double& psi20, const int& r, const double omega0, const double& phi20, const mat& Lambda_flat, const mat& Z, const mat& C, const mat& C_mat){
  double C2_sum = accu(C%C);
  mat W_r(q+r, q+r, fill::zeros);
  vec x2_r_vec(r, fill::zeros);
  vec x3_q_vec(q, fill::zeros);
  mat X = join_rows(Lambda_flat, Z);
  mat Xt = X.t();
  mat XC = X%C_mat;
  vec Qx = sum(XC.t(),1);
  rowvec Qx_t = Qx.t();
  mat Rx = Xt*X;
  mat Rx_W;

  x2_r_vec.fill(1/t2_alpha);
  x3_q_vec.fill(1/t2_beta);
  W_r.diag() = join_cols(x3_q_vec, x2_r_vec);
  Rx_W = W_r + Rx/t2;
  
  double res = -(1+0.5*(M+rho0))*log(t2) - (1+0.5*(r+omega0))*log(t2_alpha) - (1+0.5*(q+gamma0))*log(t2_beta);
  res += -0.5*(rho0*psi20/t2 + omega0*psi20/t2_alpha + gamma0*kappa20/t2_beta + C2_sum/t2);
  res += 0.5/pow(t2,2)*as_scalar(Qx_t*inv_sympd(Rx_W)*Qx);
  res += -0.5*log_det_sympd(Rx_W);
  //double res = pow(t2,-1-0.5*(M+rho0))*pow(t2_alpha,-1-0.5*(r+omega0))*pow(t2_beta, -1-0.5*(q+gamma0))*exp(-0.5*rho0*psi20/t2)*exp(-0.5*omega0*psi20/t2_alpha)*exp(-0.5*gamma0*kappa20/t2_beta)*exp(-0.5*C2_sum/t2)*exp(0.5/pow(t2,2)*as_scalar(Qx_t*inv_sympd(Rx_W)*Qx))/sqrt(det(Rx_W));
  //std::cout << -(1+0.5*(M+rho0))*log(t2) - (1+0.5*(r+omega0))*log(t2_alpha) - (1+0.5*(q+gamma0))*log(t2_beta) << std::endl;
  //std::cout << -0.5*(rho0*psi20/t2 + omega0*psi20/t2_alpha + gamma0*kappa20/t2_beta + C2_sum/t2) << std::endl;
  //std::cout << 0.5/pow(t2,2)*as_scalar(Qx_t*inv_sympd(Rx_W)*Qx) << std::endl;
  //std::cout << -0.5*log_det_sympd(Rx_W) << std::endl;
  return res;
};

// [[Rcpp::export]]
vec compute_log_A_lambda(const bool& weighted, const double& gamma0, const double& kappa20, const int& M, const int& q, const double& rho0, const double& psi20, const int& r, const double omega0, const double& phi20, const mat& Lambda_flat, const mat& Z, const mat& C, const mat& C_mat, const mat& SGI_nodes){
  double C2_sum = accu(C%C);
  mat W_r(q+r, q+r, fill::zeros);
  vec x2_r_vec(r, fill::zeros);
  vec x3_q_vec(q, fill::zeros);
  int n_nodes = SGI_nodes.n_rows;
  mat X = join_rows(Lambda_flat, Z);
  mat Xt = X.t();
  mat XC = X%C_mat;
  vec Qx = sum(XC.t(),1);
  rowvec Qx_t = Qx.t();
  mat Rx = Xt*X;
  mat Rx_W;
  vec log_h_p2_x1(n_nodes, fill::zeros);
  vec log_h_p2_x2(n_nodes, fill::zeros);
  vec log_h_p2_x3(n_nodes, fill::zeros);
  
  if(weighted){
    mat SGI_2_nodes = SGI_nodes % SGI_nodes;
    mat exp_neg_SGI_nodes = exp(-1.0*SGI_nodes);
    
    log_h_p2_x1 = -0.5*(M+rho0)*SGI_nodes.col(0) - 0.5*(r+omega0)*SGI_nodes.col(1) - 0.5*(q+gamma0)*SGI_nodes.col(2) + 0.5*sum(SGI_2_nodes,1);
    //std::cout << "log_h_p2_x1" << log_h_p2_x1 << std::endl;
    log_h_p2_x2 = -0.5*(rho0*psi20*exp_neg_SGI_nodes.col(0) + omega0*psi20*exp_neg_SGI_nodes.col(1) + gamma0*kappa20*exp_neg_SGI_nodes.col(2) + C2_sum*exp_neg_SGI_nodes.col(0));
    //std::cout << "log_h_p2_x2" << log_h_p2_x2 << std::endl;
    
    for(int i=0; i<n_nodes; i++){
      x2_r_vec.fill(exp_neg_SGI_nodes(i,1));
      x3_q_vec.fill(exp_neg_SGI_nodes(i,2));
      W_r.diag() = join_cols(x3_q_vec, x2_r_vec);
      Rx_W = W_r + exp_neg_SGI_nodes(i,0)*Rx;
      
      log_h_p2_x3[i] = 0.5*pow(exp_neg_SGI_nodes(i,0),2)*as_scalar(Qx_t*inv_sympd(Rx_W)*Qx)-0.5*log_det_sympd(Rx_W);
    }
    //std::cout << "log_h_p2_x3" << log_h_p2_x3 << std::endl;
    //vec log_h_p2_x1 = -0.5*exp(-1.0*SGI_nodes.col(0))*(C2_sum+rho0*psi20) - 0.5*(M+rho0)*SGI_nodes.col(0) + 0.5*SGI_2_nodes.col(0);
    //vec log_h_p2_x2 = -0.5*(r+omega0)*SGI_nodes.col(1) - 0.5*omega0*phi20*exp(-1.0*SGI_nodes.col(1)) + 0.5*SGI_2_nodes.col(1);
    //vec log_h_p2_x3 = -0.5*(q+gamma0)*SGI_nodes.col(2) - 0.5*gamma0*kappa20*exp(-1.0*SGI_nodes.col(2)) + 0.5*SGI_2_nodes.col(2);
    //for(int r=0; r<n_nodes; r++){
    //x2_r_vec.fill(exp(-1.0*SGI_nodes(r,1)));
    //x3_q_vec.fill(exp(-1.0*SGI_nodes(r,2)));
    //W_r.diag() = join_cols(x3_q_vec, x2_r_vec);
    //Rx_W = exp(-1.0*SGI_nodes(r,0))*Rx+W_r;
    //log_h_p2_x1[r] += -0.5*log_det_sympd(Rx_W) + 0.5*exp(-2.0*SGI_nodes(r,0))*as_scalar(Qx.t()*solve(Rx_W,I_qr)*Qx);
    //}
  }else{
    mat log_SGI_nodes_inv = log(1/SGI_nodes);
    mat log_log_SGI_nodes_inv = log(log_SGI_nodes_inv);
    log_h_p2_x1 = (0.5*(M+rho0)-1)*log_log_SGI_nodes_inv.col(0) + (0.5*(r+omega0)-1)*log_log_SGI_nodes_inv.col(1) + (0.5*(q+gamma0)-1)*log_log_SGI_nodes_inv.col(2);
    //std::cout << "log_h_p2_x1" << log_h_p2_x1 << std::endl;
    log_h_p2_x2 = -((0.5*(rho0*psi20+C2_sum)-1)*log_SGI_nodes_inv.col(0) + (0.5*(omega0*phi20)-1)*log_SGI_nodes_inv.col(1) + (0.5*(gamma0*kappa20)-1)*log_SGI_nodes_inv.col(2));
    //std::cout << "log_h_p2_x2" << log_h_p2_x2 << std::endl;
    for(int i=0; i<n_nodes; i++){
      x2_r_vec.fill(log_SGI_nodes_inv(i,1));
      x3_q_vec.fill(log_SGI_nodes_inv(i,2));
      W_r.diag() = join_cols(x3_q_vec, x2_r_vec);
      Rx_W = W_r+log_SGI_nodes_inv(i,0)*Rx;
      
      log_h_p2_x3[i] = 0.5*pow(log_SGI_nodes_inv(i,0),2)*as_scalar(Qx_t*inv_sympd(Rx_W)*Qx)-0.5*log_det_sympd(Rx_W);
      //std::cout<< pow(log_SGI_nodes_inv(i,0),2) <<std::endl;
      //std::cout<< as_scalar(Qx_t*inv_sympd(Rx_W)*Qx) <<std::endl;
      //std::cout<< log_det_sympd(Rx_W) <<std::endl;
    }
    //std::cout << "log_h_p2_x3" << log_h_p2_x3 << std::endl;
    //std::cout << (0.5*(M+rho0)+1)*log_log_SGI_nodes_inv.col(0) + (0.5*(r+omega0)+1)*log_log_SGI_nodes_inv.col(1) + (0.5*(q+gamma0)+1)*log_log_SGI_nodes_inv.col(2) << std::endl;
    //std::cout << -((0.5*(rho0*psi20+C2_sum))*log_SGI_nodes_inv.col(0) + (0.5*(omega0*phi20))*log_SGI_nodes_inv.col(1) + (0.5*(gamma0*kappa20))*log_SGI_nodes_inv.col(2)) << std::endl;
    //std::cout << 0.5*pow(log_SGI_nodes_inv(0,0),2)*as_scalar(Qx_t*inv_sympd(Rx_W)*Qx) << std::endl;
    //std::cout << -0.5*log_det_sympd(Rx_W) << std::endl;
  }

  //vec wh = log_h_p2_x1 + log_h_p2_x2 + log_h_p2_x3 + log(SGI_wts);
  vec wh = log_h_p2_x1 + log_h_p2_x2 + log_h_p2_x3;
  return wh;
};

// data structure
struct IS_dat{
  // data; hyperparameters
  mat C; //M x 1
  mat Z; //M x r
  int M; // population size
  int N;
  int q; // # reduced dimensions
  int q0;
  int r; // # covariates
  double rho0; //t2
  double psi20; //t2
  double gamma0; //t2_beta
  double kappa20; //t2_beta
  double omega0; //t2_alpha
  double phi20; //t2_alpha
  mat SGI_nodes; //n_nodes x 3
  vec SGI_wts;
  int g1_samples;
  char scheme;
  bool weighted;
  
  cube Lambda_flat; // g1 results; M by q by iters
  cube U;
  cube X;
  vec s2;
  vec t2_lambda;
  mat d;
  vec t2;
  vec t2_beta;
  vec t2_alpha;
  vec Y_llk;
  vec C_llk;
};

struct IS_comp_par{
  //mat X;
  mat C_mat; //repmat(C,1,W_dim)
};

vec eval_log_A_lambda(const IS_dat& dat, const IS_comp_par& comp_par, const mat& Lambda_flat){
  return compute_log_A_lambda(dat.weighted, dat.gamma0, dat.kappa20, dat.M, dat.q, dat.rho0, dat.psi20, dat.r, dat.omega0, dat.phi20, Lambda_flat, dat.Z, dat.C, comp_par.C_mat, dat.SGI_nodes);
}

typedef vec(*funcPtr_A_lambda)(const IS_dat& dat, const IS_comp_par& comp_par, const mat& Lambda_flat);

struct IS_paras{
  double T; //maximum of l_lambda
  //use for MH
  mat Lambda;
  double l_lambda;
  mat U;
  mat X;
  double s2;
  double t2_lambda;
  vec d;
  double t2;
  double t2_beta;
  double t2_alpha;
  double C_llk;
  double Y_llk;
};

List log_A_lambda_updates_rej(IS_paras& paras, IS_dat& dat, IS_comp_par& comp_par){
  uvec acpt_rej(dat.g1_samples);
  uvec acpt_idx;
  vec l_lambda_vec(dat.g1_samples);
  double log_unif;
  vec l_lambda_SGI;
  mat Lam_tmp;
  
  for(int iter=0; iter<dat.g1_samples; iter++){
    Lam_tmp = dat.Lambda_flat.slice(iter);
    l_lambda_SGI = eval_log_A_lambda(dat, comp_par, Lam_tmp);
    l_lambda_vec[iter] = log(abs(sum(dat.SGI_wts % exp(l_lambda_SGI - max(l_lambda_SGI))))) + max(l_lambda_SGI);
  }
  paras.T = max(l_lambda_vec);
  
  for(int iter=0; iter<dat.g1_samples; iter++){
    log_unif = log(R::runif(0,1));
    //Lam_tmp = dat.Lambda_flat.slice(iter);
    //l_lambda_SGI = eval_log_A_lambda(dat, comp_par, Lam_tmp);
    //l_lambda_vec[iter] = log(abs(sum(dat.SGI_wts % exp(l_lambda_SGI - max(l_lambda_SGI))))) + max(l_lambda_SGI);
    if(log_unif <= l_lambda_vec[iter] - paras.T){
      acpt_rej[iter] = 1;
    }else{
      acpt_rej[iter] = 0;
    }
    std::cout << iter << std::endl;
  }
  
  acpt_idx = arma::find(acpt_rej==1);
  vec t2_lambda_list = dat.t2_lambda.elem(acpt_idx);
  vec s2_list = dat.s2.elem(acpt_idx);
  vec t2_alpha_list = dat.t2_alpha.elem(acpt_idx);
  vec t2_beta_list = dat.t2_beta.elem(acpt_idx);
  vec t2_list = dat.t2.elem(acpt_idx);
  vec Y_llk_list = dat.Y_llk.elem(acpt_idx);
  vec C_llk_list = dat.C_llk.elem(acpt_idx);
  vec Y_C_llk_list = Y_llk_list + C_llk_list;
  mat d_list = dat.d.cols(acpt_idx);
  cube Lambda_flat_list = dat.Lambda_flat.slices(acpt_idx);
  cube X_list = dat.X.slices(acpt_idx);
  cube U_list = dat.U.slices(acpt_idx);
  
  return List::create(Named("l_lambda") = l_lambda_vec,
                      Named("max_l_lambda") = paras.T,
                      Named("keep_idx") = acpt_rej,
                      Named("acpt_idx") = acpt_idx+1,
                      Named("Lambda_flat") = Lambda_flat_list,
                      Named("t2_lambda") = t2_lambda_list,
                      Named("X") = X_list,
                      Named("U") = U_list,
                      Named("s2") = s2_list,
                      Named("d") = d_list,
                      Named("t2_alpha") = t2_alpha_list,
                      Named("t2_beta") = t2_beta_list,
                      Named("t2") = t2_list,
                      Named("Y_llk") = Y_llk_list,
                      Named("C_llk") = C_llk_list,
                      Named("Y_C_llk") = Y_C_llk_list,
                      Named("acpt_ratio") = accu(acpt_rej)/dat.g1_samples,
                      Named("scheme") = dat.scheme);
};

void update_l_lambda(IS_paras& paras, IS_dat& dat, IS_comp_par& comp_par){
  vec l_lambda_SGI;
  l_lambda_SGI = eval_log_A_lambda(dat, comp_par, paras.Lambda);
  paras.l_lambda = log(abs(sum(dat.SGI_wts % exp(l_lambda_SGI - max(l_lambda_SGI))))) + max(l_lambda_SGI);
};

// [[Rcpp::export]]
List A_lambda_rej(cube& Lambda_flat, cube& U, cube& X, vec& s2, vec& t2_lambda, mat& d, vec& t2, vec& t2_beta, vec& t2_alpha, vec& Y_llk, vec& C_llk,
                  double& rho0, double& psi20, double& gamma0, double& kappa20, double& omega0, double& phi20,
                  mat& SGI_nodes, vec& SGI_wts, int& g1_samples, mat& C, mat& Z, int& M, int& N, int& q, int& r, char& scheme, bool& weighted){
  IS_dat dat;
  dat.C = C; //M x 1
  dat.Z = Z; //M x r
  dat.M = M; // population size
  dat.N = N;
  dat.q = q; // # reduced dimensions
  dat.r = r; // # covariates
  dat.rho0 = rho0; //t2
  dat.psi20 = psi20; //t2
  dat.gamma0 = gamma0; //t2_beta
  dat.kappa20 = kappa20; //t2_beta
  dat.omega0 = omega0; //t2_alpha
  dat.phi20 = phi20; //t2_alpha
  dat.Lambda_flat = Lambda_flat; // g1 results; M by q by iters
  dat.U = U;
  dat.X = X;
  dat.s2 = s2;
  dat.t2_lambda = t2_lambda;
  dat.d = d;
  dat.t2 = t2;
  dat.t2_beta = t2_beta;
  dat.t2_alpha = t2_alpha;
  dat.Y_llk = Y_llk;
  dat.C_llk = C_llk;
  dat.SGI_nodes = SGI_nodes; //n_nodes x 3
  dat.SGI_wts = SGI_wts;
  dat.g1_samples = g1_samples;
  dat.scheme = scheme;
  dat.weighted = weighted;
  
  IS_comp_par comp;
  comp.C_mat = repmat(C,1,q+r);
  
  IS_paras paras;
  paras.Lambda = Lambda_flat.slice(0);
  paras.U = U.slice(0);
  paras.X = X.slice(0);
  paras.s2 = s2[0];
  paras.t2_lambda = t2_lambda[0];
  paras.d = d.col(0);
  paras.t2 = t2[0];
  paras.t2_beta = t2_beta[0];
  paras.t2_alpha = t2_alpha[0];
  paras.C_llk = C_llk[0];
  paras.Y_llk = Y_llk[0];
  update_l_lambda(paras, dat, comp);;
  
  List A_lambda_res = log_A_lambda_updates_rej(paras, dat, comp);
  return A_lambda_res;
};


void update_Lambda(IS_paras& paras, IS_dat& dat, int& iter){
  paras.Lambda = dat.Lambda_flat.slice(iter);
};

void update_U(IS_paras& paras, IS_dat& dat, int& iter){
  paras.U = dat.U.slice(iter);
};

void update_X(IS_paras& paras, IS_dat& dat, int& iter){
  paras.X = dat.X.slice(iter);
};

void update_s2(IS_paras& paras, IS_dat& dat, int& iter){
  paras.s2 = dat.s2[iter];
};

void update_t2_lambda(IS_paras& paras, IS_dat& dat, int& iter){
  paras.t2_lambda = dat.t2_lambda[iter];
};

void update_d(IS_paras& paras, IS_dat& dat, int& iter){
  paras.d = dat.d.col(iter);
};

void update_t2(IS_paras& paras, IS_dat& dat, int& iter){
  paras.t2 = dat.t2[iter];
};

void update_t2_beta(IS_paras& paras, IS_dat& dat, int& iter){
  paras.t2_beta = dat.t2_beta[iter];
};

void update_t2_alpha(IS_paras& paras, IS_dat& dat, int& iter){
  paras.t2_alpha = dat.t2_alpha[iter];
};

void update_Y_llk(IS_paras& paras, IS_dat& dat, int& iter){
  paras.Y_llk = dat.Y_llk[iter];
};

void update_C_llk(IS_paras& paras, IS_dat& dat, int& iter){
  paras.C_llk = dat.C_llk[iter];
};

bool MH_one_step(IS_paras& paras, IS_dat& dat, IS_comp_par& comp_par, int& iter){
  double l_lambda_n = paras.l_lambda;
  update_Lambda(paras, dat, iter);
  update_l_lambda(paras, dat, comp_par);
  double log_unif = log(R::runif(0,1));
  
  vec diff_vec(2);
  diff_vec[0] = 0;
  diff_vec[1] = paras.l_lambda - l_lambda_n;
  double thres = min(diff_vec);
  if(log_unif<thres){
    return true;
  }else{
    return false;
  }
};

List log_A_lambda_updates_MH(IS_paras& paras, IS_dat& dat, IS_comp_par& comp_par){
  // Tracking samples
  cube Lambda_flat_list=zeros<cube>(dat.M, dat.q, dat.g1_samples);
  vec t2_lambda_list = zeros<vec>(dat.g1_samples);
  cube U_list=zeros<cube>(dat.N, dat.q0, dat.g1_samples);
  cube X_list = zeros<cube>(dat.N, dat.q0, dat.g1_samples);
  vec s2_list = zeros<vec>(dat.g1_samples);
  mat d_list = zeros<mat>(dat.q+dat.r,dat.g1_samples);
  vec t2_alpha_list = zeros<vec>(dat.g1_samples);
  vec t2_beta_list = zeros<vec>(dat.g1_samples);
  vec t2_list = zeros<vec>(dat.g1_samples);
  vec Y_llk_list = zeros<vec>(dat.g1_samples);
  vec C_llk_list = zeros<vec>(dat.g1_samples);
  
  vec Y_C_llk_list = zeros<vec>(dat.g1_samples);
  vec acpt_rej = zeros<vec>(dat.g1_samples);
  vec l_lambda_vec = zeros<vec>(dat.g1_samples);
  
  // initialize
  // collect samples
  Lambda_flat_list.slice(0) = paras.Lambda;
  t2_lambda_list[0] = paras.t2_lambda;
  U_list.slice(0) = paras.U;
  X_list.slice(0) = paras.X;
  s2_list[0] = paras.s2;
  d_list.col(0) = paras.d;
  t2_alpha_list[0] = paras.t2_alpha;
  t2_beta_list[0] = paras.t2_beta;
  t2_list[0] = paras.t2;
  Y_llk_list[0] = paras.Y_llk;
  C_llk_list[0] = paras.C_llk;
  l_lambda_vec[0] = paras.l_lambda;
  
  // collect samples  [[bookmark]]
  mat Lambda_current;
  double l_lambda_current;
  
  for(int iter=1; iter<dat.g1_samples; iter++){
    Lambda_current = paras.Lambda;
    l_lambda_current = paras.l_lambda;
    
    if(MH_one_step(paras, dat, comp_par, iter)){
      update_U(paras,dat,iter);
      update_X(paras,dat,iter);
      update_s2(paras,dat,iter);
      update_t2_lambda(paras,dat,iter);
      update_d(paras,dat,iter);
      update_t2(paras,dat,iter);
      update_t2_beta(paras,dat,iter);
      update_t2_alpha(paras,dat,iter);
      update_Y_llk(paras,dat,iter);
      update_C_llk(paras,dat,iter);
      
      acpt_rej[iter] = 1;
    }else{
      paras.Lambda = Lambda_current;
      paras.l_lambda = l_lambda_current;
    }
    
    //save samples 
    Lambda_flat_list.slice(iter)=paras.Lambda;
    t2_lambda_list[iter] = paras.t2_lambda;
    U_list.slice(iter)=paras.U;
    X_list.slice(iter) = paras.X;
    s2_list[iter] = paras.s2;
    d_list.col(iter) = paras.d;
    t2_alpha_list[iter] = paras.t2_alpha;
    t2_beta_list[iter] = paras.t2_beta;
    t2_list[iter] = paras.t2;
    Y_llk_list[iter] = paras.Y_llk;
    C_llk_list[iter] = paras.C_llk;
    Y_C_llk_list[iter] = Y_llk_list[iter] + C_llk_list[iter];
    l_lambda_vec[iter] = paras.l_lambda;
  }
  
  // return results 
  return List::create(Named("Lambda_flat") = Lambda_flat_list,
                      Named("t2_lambda") = t2_lambda_list,
                      Named("X")=X_list,
                      Named("U")=U_list,
                      Named("s2") = s2_list,
                      Named("d") = d_list,
                      Named("t2_alpha") = t2_alpha_list,
                      Named("t2_beta") = t2_beta_list,
                      Named("t2") = t2_list,
                      Named("Y_llk") = Y_llk_list,
                      Named("C_llk") = C_llk_list,
                      Named("Y_C_llk") = Y_C_llk_list,
                      Named("l_lambda") = l_lambda_vec,
                      Named("acpt")=acpt_rej,
                      Named("scheme") = dat.scheme);
};

// [[Rcpp::export]]
List A_lambda_MH(cube& Lambda_flat, cube& U, cube& X, vec& s2, vec& t2_lambda, mat& d, vec& t2, vec& t2_beta, vec& t2_alpha, vec& Y_llk, vec& C_llk,
                 double& rho0, double& psi20, double& gamma0, double& kappa20, double& omega0, double& phi20,
                 mat& SGI_nodes, vec& SGI_wts, int& g1_samples, mat& C, mat& Z, int& M, int& N, int& q, int& q_0, int& r, char& scheme, bool& weighted){
  IS_dat dat;
  dat.C = C; //M x 1
  dat.Z = Z; //M x r
  dat.M = M; // population size
  dat.N = N;
  dat.q0 = q_0;
  dat.q = q; // # reduced dimensions
  dat.r = r; // # covariates
  dat.rho0 = rho0; //t2
  dat.psi20 = psi20; //t2
  dat.gamma0 = gamma0; //t2_beta
  dat.kappa20 = kappa20; //t2_beta
  dat.omega0 = omega0; //t2_alpha
  dat.phi20 = phi20; //t2_alpha
  dat.Lambda_flat = Lambda_flat; // g1 results; M by q by iters
  dat.U = U;
  dat.X = X;
  dat.s2 = s2;
  dat.t2_lambda = t2_lambda;
  dat.d = d;
  dat.t2 = t2;
  dat.t2_beta = t2_beta;
  dat.t2_alpha = t2_alpha;
  dat.Y_llk = Y_llk;
  dat.C_llk = C_llk;
  dat.SGI_nodes = SGI_nodes; //n_nodes x 3
  dat.SGI_wts = SGI_wts;
  dat.g1_samples = g1_samples;
  dat.scheme = scheme;
  dat.weighted = weighted;
  
  IS_comp_par comp;
  comp.C_mat = repmat(C,1,q+r);
  
  IS_paras paras;
  paras.Lambda = Lambda_flat.slice(0);
  paras.U = U.slice(0);
  paras.X = X.slice(0);
  paras.s2 = s2[0];
  paras.t2_lambda = t2_lambda[0];
  paras.d = d.col(0);
  paras.t2 = t2[0];
  paras.t2_beta = t2_beta[0];
  paras.t2_alpha = t2_alpha[0];
  paras.C_llk = C_llk[0];
  paras.Y_llk = Y_llk[0];
  update_l_lambda(paras, dat, comp);;
  
  List A_lambda_res = log_A_lambda_updates_MH(paras, dat, comp);
  return A_lambda_res;
};
