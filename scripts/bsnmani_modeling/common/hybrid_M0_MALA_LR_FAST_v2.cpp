// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]

#include <iostream>
//#include "000.macros.h"
//#include "000.templates-types.h"
//#include "000.types.h"
//#include "000.utils.h"
//#include "logSumExp_lowlevel.h"
//#include "logSumExp_lowlevel_template.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat polar_expansion(const mat& X){
  mat XtX = X.t()*X;
  mat res =   X*inv(sqrtmat_sympd(XtX));
  return res; //cost of saving inv(); X(XtX)^(-1/2)
};

// [[Rcpp::export]]
double adjust_acceptance(double accept,double sgm,double target = 0.5){
  double y = 1. + 1000.*(accept-target)*(accept-target)*(accept-target);
  if (y < .9)
    y = .9;
  if (y > 1.1)
    y = 1.1;
  sgm *= y;
  return sgm;
};

//[[Rcpp::export]]
mat kron_r(const mat& A, const mat& B){
  return kron(A,B);
};

// [[Rcpp::export]]
mat box_prod_2(const mat& A, const mat& B, const int& m1, const int& n1, const int& m2, const int& n2){
  umat idx_mat = linspace<umat>(0,(n1*n2-1),(n1*n2));
  idx_mat.reshape(n2,n1);
  uvec idx_vec = vectorise(idx_mat.t());
  mat res = kron(A,B);
  //std::cout<<size(res)<<std::endl;
  return res.cols(idx_vec);
};

// [[Rcpp::export]]
mat compute_d_UtYU(const mat& U, const mat& Ut, const cube& Y, const int& M, const int& q){
  mat res(q,M);
  for(int i=0; i<M; i++){
    res.col(i) = diagvec(Ut*Y.slice(i)*U);
  }
  return res;
};

// compute log(w_r * h(x_r)); compatible with logsumexp
// [[Rcpp::export]]
vec compute_wh(const int& M, const int& N, const int& q, const double& s2, const mat& d_UtYU, const double& eta0, const double& t20, const vec& SGI_nodes, const vec& SGI_wts, const bool& weighted){
  vec h_Theta_Y_p1(SGI_wts.n_elem);
  vec h_Theta_Y_p2(SGI_wts.n_elem);
  if(weighted){
    h_Theta_Y_p1.fill(log(sqrt(2*datum::pi)));
    h_Theta_Y_p1 += -0.5*(q*M+eta0)*SGI_nodes+0.5*pow(SGI_nodes,2)-0.5*eta0*t20*exp(-1.0*SGI_nodes) - 0.5*q*M*log(1/s2+exp(-1.0*SGI_nodes));
    h_Theta_Y_p2 = 0.5*pow(s2,-2)/(1/s2+exp(-1.0*SGI_nodes));
  }else{
    h_Theta_Y_p1 = (0.5*(q*M+eta0)-1)*log(log(1/SGI_nodes)) - 0.5*q*M*log(1/s2+log(1/SGI_nodes)) - (0.5*eta0*t20-1)*log(SGI_nodes);
    h_Theta_Y_p2 = 0.5*pow(s2,-2)/(1/s2+log(1/SGI_nodes));
  }
  vec wh = h_Theta_Y_p1 + h_Theta_Y_p2*accu(d_UtYU % d_UtYU) + log(SGI_wts);
  return wh;
};

// [[Rcpp::export]]
double eval_log_l_Theta_Y(const vec& wh_vec, const int& M, const int& N, const double& s2, const double& accu_YY, const mat& X, const double& nu0, const double& s20){
  // log_l_Theta_Y
  //double res = -0.25*M*N*(N+1)*log(s2)-0.5/s2*accu_YY+log(sum(exp(wh_vec - max(wh_vec))))+max(wh_vec)-0.5*accu(X%X)+(1-nu0/2.0)*log(s2)-0.5*nu0*s20/s2;
  double res = -0.25*M*N*(N+1)*log(s2)-0.5/s2*accu_YY+log(sum(exp(wh_vec - max(wh_vec))))+max(wh_vec)-0.5*accu(X%X)-(1+nu0/2.0)*log(s2)-0.5*nu0*s20/s2;
  return res;
};

// [[Rcpp::export]]
vec eval_d_log_h_s2(const int& M, const int& q, const double& s2, const mat& d_UtYU, const vec& SGI_nodes, const bool& weighted){
  // d_log_h_Theta_Y_s2
  vec nodes_t2_lambda = SGI_nodes;
  if(weighted){
    nodes_t2_lambda = exp(nodes_t2_lambda);
  }else{
    nodes_t2_lambda = 1.0/(log(1.0/nodes_t2_lambda));
  }
  vec eval_d_log_h_Theta_Y_s2_p1 = 0.5*q*M/(pow(s2,2)/nodes_t2_lambda+s2);
  vec eval_d_log_h_Theta_Y_s2_p2 = 0.5*(2*s2/nodes_t2_lambda+1)/(pow(s2,2)*pow(s2/nodes_t2_lambda+1,2));

  vec res = eval_d_log_h_Theta_Y_s2_p1 - eval_d_log_h_Theta_Y_s2_p2*accu(d_UtYU%d_UtYU);
  //vec res = eval_d_log_h_Theta_Y_s2_p1;
  return res;
};

// [[Rcpp::export]]
double eval_d_log_l_s2(const vec& wh_vec, const double& nu0, const double& s20, const int& M, const int& N, const int& q, const double& s2, const double& accu_YY, const mat& d_UtYU, const vec& SGI_nodes, const bool& weighted){
  // d_log_h_Theta_Y_s2
  vec nodes_t2_lambda = SGI_nodes;
  if(weighted){
    nodes_t2_lambda = exp(nodes_t2_lambda);
  }else{
    nodes_t2_lambda = 1.0/(log(1.0/nodes_t2_lambda));
  }
  vec eval_d_log_h_Theta_Y_s2_p1 = 0.5*q*M/(pow(s2,2)/nodes_t2_lambda+s2);
  vec eval_d_log_h_Theta_Y_s2_p2 = 0.5*(2*s2/nodes_t2_lambda+1)/(pow(s2,2)*pow(s2/nodes_t2_lambda+1,2));

  // d_log_l_Theta_Y_s2
  double res = (-1-nu0/2.0-0.25*M*N*(N+1))/s2 + 0.5/(pow(s2,2))*(nu0*s20+accu_YY);
  //std::cout<<res<<std::endl;
  vec wh_vec_a = wh_vec + log(eval_d_log_h_Theta_Y_s2_p1);
  //std::cout<<wh_vec_a<<std::endl;
  vec wh_vec_b = wh_vec + log(eval_d_log_h_Theta_Y_s2_p2);
  //std::cout<<wh_vec_b<<std::endl;
  res += exp(log(sum(exp(wh_vec_a - max(wh_vec_a))))+max(wh_vec_a) - log(sum(exp(wh_vec - max(wh_vec))))-max(wh_vec));
  res -= exp(log(sum(exp(wh_vec_b - max(wh_vec_b))))+max(wh_vec_b) - log(sum(exp(wh_vec - max(wh_vec))))-max(wh_vec)) * accu(d_UtYU%d_UtYU);
  return res;
};

//[[Rcpp::export]]
mat eval_dU_dX(const mat& I_N, const mat& I_q, const mat& X, const mat& dSdX, const int& N, const int& q, const int& q2){
  mat S = X.t();
  mat P = solve(S*X,I_q);
  mat R = sqrtmat_sympd(P);
  mat dU_X = kron(R,I_N) + kron(I_q,X)*solve(kron(R,I_q)+kron(I_q,R),eye(q2,q2))*kron(P,-1.0*P)*(kron(S,I_q)*dSdX+kron(I_q,S)); //changed here
  return dU_X;
};

//[[Rcpp::export]]
mat eval_dftf_U(const mat& Y, const mat& U, const int& N, const int& q){
  mat Ut = U.t();
  mat UtYi;
  mat UtYU = Ut*Y*U;
  vec d_UtYU = UtYU.diag();
  UtYi = Ut*Y;
  UtYi = UtYi.each_col() % d_UtYU;
  return UtYi;
};

//[[Rcpp::export]]
rowvec eval_d_tr_XtX(const mat& X){
  rowvec res = vectorise(-1.0*X).as_row();
  return res;
};

//[[Rcpp::export]]
rowvec eval_d_log_h_X_1(const mat& U, const int& M, const int& q, const int& N, const cube& Y, const mat& d_UtYU){
  mat Ut = U.t();
  mat UtYi;
  mat UY_diag_sum(N*q,M);
  for(int i=0; i<M; i++){
    UtYi = Ut*Y.slice(i);
    UtYi = UtYi.each_col() % d_UtYU.col(i);
    UY_diag_sum.col(i) = vectorise(UtYi.t());
  }

  rowvec res = sum(UY_diag_sum.t());
  return res;
};

//[[Rcpp::export]]
vec eval_d_log_l_X(const vec& wh_vec, const vec& SGI_nodes, const bool& weighted, const double& s2, const mat& d_UtYU, const cube& Y, const int& M, const int& N, const int& q, const mat& U, const mat& X, const mat& I_N, const mat& I_q, const mat& dSdX, const int& q2){
  vec nodes_t2_lambda = SGI_nodes;
  if(weighted){
    nodes_t2_lambda = exp(nodes_t2_lambda);
  }else{
    nodes_t2_lambda = 1.0/(log(1.0/nodes_t2_lambda));
  }

  // compute dU_x
  mat S = X.t();
  mat P = solve(S*X,I_q);
  mat R = sqrtmat_sympd(P);
  mat dU_X = kron(R,I_N) + kron(I_q,X)*solve(kron(R,I_q)+kron(I_q,R),eye(q2,q2))*kron(P,-1.0*P)*(kron(S,I_q)*dSdX+kron(I_q,S)); //changed here

  mat Ut = U.t();
  mat UtYi;
  mat UY_diag_sum(N*q,M);
  for(int i=0; i<M; i++){
    UtYi = Ut*Y.slice(i);
    UtYi = UtYi.each_col() % d_UtYU.col(i);
    UY_diag_sum.col(i) = vectorise(UtYi.t());
  }

  rowvec res = vectorise(-1.0*X).as_row();
  //res += sum(wh_vec % (2/(pow(s2,2))/(1.0/s2+1.0/nodes_t2_lambda)))*reshape(sum(UY_diag_sum.t())*dU_X,N,q)/sum(wh_vec);
  vec wh_vec_wt = wh_vec + log((2/(pow(s2,2))/(1.0/s2+1.0/nodes_t2_lambda)));
  res += exp(log(sum(exp(wh_vec_wt - max(wh_vec_wt))))+max(wh_vec_wt) - log(sum(exp(wh_vec - max(wh_vec))))-max(wh_vec)) * sum(UY_diag_sum.t())*dU_X;
  return res.as_col();
};

// [[Rcpp::export]]
mat Lambda_new_t_Y(const int& q, const int& M, const mat& d_UtYU, const cube& Y, const double& t2_lambda, const double& s2){
  mat var_mat(q,q,fill::zeros);
  var_mat.diag() += 1.0/(1.0/t2_lambda+1.0/s2);

  mat mu_mat = d_UtYU/s2;

  mat Z_0 = randn(q,M);
  mat A = chol(var_mat, "lower");
  mat Lambda_flat_new = A*Z_0 + var_mat*mu_mat;

  return Lambda_flat_new;
};

// [[Rcpp::export]]
double t2_lambda_new(const double& t2_lambda_A, const double& t2_lambda_B, const double& accu_L2){
  double t2_lambda = 1.0/randg<double>(distr_param(t2_lambda_A, 2.0/(t2_lambda_B + accu_L2)));
  return t2_lambda;
};

// [[Rcpp::export]]
double Y_llk_new(const mat& Y_flat, const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2, const int& M, const int& N){
  mat Y_mean;
  mat Y_mean_flat(N*N,M);

  double s = sqrt(s2);
  mat Y_smat(N*N,M, fill::value(s));
  for(int i = 0; i < M; i++){
    Y_mean = U*diagmat(Lambda_flat_t.col(i))*Ut;
    Y_mean_flat.col(i) = vectorise(Y_mean);
  }

  double dat_llk = accu(log_normpdf(Y_flat, Y_mean_flat, Y_smat));
  dat_llk += M*(0.5*N*N*log(2*datum::pi)-0.25*N*(N+1)*log(2*datum::pi)+0.25*N*(N-1)*log(2)+0.5*N*N*log(s2)-0.25*N*(N+1)*log(s2));

  return dat_llk;
};

// [[Rcpp::export]]
double Y_llk_new_RP(const mat& Y_flat, const mat& U, const mat& Ut, const mat& Lambda_flat_t, const vec& s2, const int& M, const int& N){
  mat Y_mean;
  mat Y_mean_flat(N*(N+1)/2,M);
  
  vec s = sqrt(s2);
  mat Y_smat(N*(N+1)/2,M);
  Y_smat.each_col() = s;
  for(int i = 0; i < M; i++){
    Y_mean = U*diagmat(Lambda_flat_t.col(i))*Ut;
    Y_mean_flat.col(i) = Y_mean(trimatl_ind(size(Y_mean)));
  }
  
  double dat_llk = accu(log_normpdf(Y_flat, Y_mean_flat, Y_smat));

  return dat_llk;
};

// [[Rcpp::export]]
vec test_tril(const mat& X){
  vec res = X(trimatl_ind(size(X),-1));
  return res;
}

// [[Rcpp::export]]
uvec box_prod_3(const int& n1, const int& n2){
  umat idx_mat = linspace<umat>(0,(n1*n2-1),(n1*n2));
  idx_mat.reshape(n2,n1);
  uvec idx_vec = vectorise(idx_mat.t());
  return idx_vec;
};

// [[Rcpp::export]]
mat eval_dlog_mat_norm_2_new(const mat& X, const mat& U, const mat& Ut, const mat& I_q, const mat& dSdX, const mat& I_Nq, const mat& I_N, const mat& Lambda_flat, const mat& Lambda_flat_t, const cube& Y, const int& q2, const int& Nq,
                         const uvec& box_idx, const mat& I_q_flat, const double& s2, const int& N, const int& q, const int& M){
  //mat Ut = U.t();
  mat S = X.t();
  mat P = solve(S*X,I_q);
  mat R = sqrtmat_sympd(P);
  //mat Lambda_flat_t = Lambda_flat.t();

  mat UtkI_q = kron(Ut,I_q);
  mat Lambda_i, Y_i, UtYi;
  mat res_temp_1(Nq,Nq,fill::zeros);
  mat res_temp_2(q,N,fill::zeros);
  for(int epoch=0;epoch<M;epoch++){
    Lambda_i = diagmat(Lambda_flat_t.col(epoch));
    Y_i = Y.slice(epoch);
    UtYi = Ut*Y_i;
    res_temp_1 += box_prod_2(Y_i,Lambda_i,N,N,q,q);
    res_temp_2 += Lambda_i*UtYi;
  }

  mat res_temp = UtkI_q*res_temp_1 + kron(I_q,res_temp_2);
  mat res = -1.0*S + trans(reshape(I_q_flat*1.0/s2*res_temp*(kron(R,I_N) + kron(I_q,X)*solve(kron(R,I_q)+kron(I_q,R),eye(q2,q2))*kron(P,-1.0*P)*(kron(S,I_q)*dSdX+kron(I_q,S))), N, q));
  //rowvec res = vectorise(-1.0*X).as_row() + I_q_flat*1.0/s2*res_temp*(kron(R,I_N) + kron(I_q,X)*solve(kron(R,I_q)+kron(I_q,R),eye(q2,q2))*kron(P,-1.0*P)*(kron(S,I_q)*dSdX+kron(I_q,S)));

  return res;
};

//[[Rcpp::export]]
vec subnetwork_region(const mat& U, const mat& Ut, const uvec& region_select, const uvec& upper_idx){
  mat subnet = U.rows(region_select)*Ut.cols(region_select);
  //std::cout<<subnet<<std::endl;
  vec res = subnet(upper_idx);
  return res;
}

// data structure
struct hybrid_dat{
  // data; hyperparameters
  cube Y; // N x N x M
  int M; // population size
  int N; // # roi
  int q; // # reduced dimensions
  double nu0; //s2
  double s20; //s2
  double eta0; //t2_lambda
  double t20; //t2_lambda
  vec SGI_nodes;
  vec SGI_wts;
  bool weighted;
  uvec region_select;
};

struct hybrid_comp_par{
  mat I_q;
  mat dSdX;
  mat I_N;
  int q2;
  double accu_YY; //accu(Y%Y)
  mat Y_flat; //(N*N x M)
  double t2_lambda_A;
  double t2_lambda_B;
  uvec region_upper_idx;
};

mat gibbs_Lambda_t_Y(const hybrid_dat& dat, const mat& d_UtYU, const double& t2_lambda, const double& s2){
  return Lambda_new_t_Y(dat.q, dat.M, d_UtYU, dat.Y, t2_lambda, s2);
}

double gibbs_t2_lambda(const hybrid_comp_par& comp_par, const double& accu_L2){
  return t2_lambda_new(comp_par.t2_lambda_A, comp_par.t2_lambda_B, accu_L2);
}

double eval_dat_llk(const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2, const hybrid_dat& dat, const hybrid_comp_par& comp_par){
  return Y_llk_new(comp_par.Y_flat, U, Ut, Lambda_flat_t, s2, dat.M, dat.N);
}

mat eval_d_UtYU(const mat& U, const mat& Ut, const hybrid_dat& dat){
  return compute_d_UtYU(U, Ut, dat.Y, dat.M, dat.q);
}

vec eval_wh(const double& s2, const mat& d_UtYU, const hybrid_dat& dat){
  return compute_wh(dat.M, dat.N, dat.q, s2, d_UtYU, dat.eta0, dat.t20, dat.SGI_nodes, dat.SGI_wts, dat.weighted);
}

double eval_pi(const vec& wh_vec, const double& s2, const mat& X, const hybrid_dat& dat, const hybrid_comp_par& comp_par){
  return eval_log_l_Theta_Y(wh_vec, dat.M, dat.N, s2, comp_par.accu_YY, X, dat.nu0, dat.s20);
}

// [[Rcpp::export]]
vec eval_dpi_check(const vec& wh_vec, const double& s2, const mat& d_UtYU, const mat& U, const mat& X, const vec& SGI_nodes, const bool& weighted, const cube& Y, const int& M, const int& N, const int& q, const mat& I_N, const mat& I_q, const mat& dSdX, const int& q2, const double& nu0, const double& s20, const double& accu_YY){
  vec d_log_l_X = eval_d_log_l_X(wh_vec, SGI_nodes, weighted, s2, d_UtYU, Y, M, N, q, U, X, I_N, I_q, dSdX, q2);
  int Nq = N*q;
  vec res(Nq+1,fill::zeros);
  res.subvec(0,Nq-1) = d_log_l_X;
  res[Nq] = eval_d_log_l_s2(wh_vec, nu0, s20, M, N, q, s2, accu_YY, d_UtYU, SGI_nodes, weighted);
  return res;
}

vec eval_dpi(const vec& wh_vec, const double& s2, const mat& d_UtYU, const mat& U, const mat& X, const hybrid_dat& dat, const hybrid_comp_par& comp_par){
  vec d_log_l_X = eval_d_log_l_X(wh_vec, dat.SGI_nodes, dat.weighted, s2, d_UtYU, dat.Y, dat.M, dat.N, dat.q, U, X, comp_par.I_N, comp_par.I_q, comp_par.dSdX, comp_par.q2);
  int Nq = dat.N*dat.q;
  vec res(Nq+1,fill::zeros);
  res.subvec(0,Nq-1) = d_log_l_X;
  res[Nq] = eval_d_log_l_s2(wh_vec, dat.nu0, dat.s20, dat.M, dat.N, dat.q, s2, comp_par.accu_YY, d_UtYU, dat.SGI_nodes, dat.weighted);
  return res;
}

vec eval_subnetwork(const mat& U, const mat& Ut, const hybrid_dat& dat, const hybrid_comp_par& comp_par){
  vec res = subnetwork_region(U, Ut, dat.region_select, comp_par.region_upper_idx);
  return res;
}

// pointers to gibbs parameter update function using data struct as input
typedef mat(*funcPtr_Lambda_t_Y)(const hybrid_dat& dat, const mat& d_UtYU, const double& t2_lambda, const double& s2);
typedef double(*funcPtr_t2_lambda)(const hybrid_comp_par& comp_par, const double& accu_L2);
typedef double(*funcPtr_dat_llk)(const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2, const hybrid_dat& dat, const hybrid_comp_par& comp_par);
typedef mat(*funcPtr_d_UtYU)(const mat& U, const mat& Ut, const hybrid_dat& dat);
typedef vec(*funcPtr_wh)(const double& s2, const mat& d_UtYU, const hybrid_dat& dat);
typedef double(*funcPtr_pi)(const vec& wh_vec, const double& s2, const mat& X, const hybrid_dat& dat, const hybrid_comp_par& comp_par);
typedef vec(*funcPtr_dpi)(const vec& wh_vec, const double& s2, const mat& d_UtYU, const mat& U, const mat& X, const hybrid_dat& dat, const hybrid_comp_par& comp_par);
typedef vec(*funcPtr_subnet)(const mat& U, const mat& Ut, const hybrid_dat& dat, const hybrid_comp_par& comp_par);

// gibbs parameter struct
struct Lambda_paras{
  mat Lambda_flat; // M x q
  mat Lambda_flat_t; //q x M
  double accu_L2; //accu(Lambda%Lambda)
};

struct MALA_paras{
  mat U; //N x q
  mat X; //N x q
  mat Ut;
  mat d_UtYU;//q x M
  double s2;
  vec wh_vec; //w_r*h(x_r)
  vec subnet;
};

struct hybrid_paras{
  Lambda_paras Lambda;
  double t2_lambda;
  MALA_paras MALA_par; // parameter that needs to be updated
  //MALA stepsize
  double w; //stepsize
  double pi;
  vec dpi;
  double llk; // value of data llk
  funcPtr_Lambda_t_Y gibbs_Lambda_t_Y;
  funcPtr_t2_lambda gibbs_t2_lambda;
  funcPtr_wh eval_wh_func;
  funcPtr_pi pi_func;
  funcPtr_dpi dpi_func;
  funcPtr_dat_llk dat_llk_func; // pointer to data llk function
  funcPtr_d_UtYU d_UtYU_func;
  funcPtr_subnet subnet_func;
  bool Lambda_fix;
  bool t2_lambda_fix;
  bool U_fix;
  bool s2_fix;
};

//[[Rcpp::export]]
mat test_reshape(const vec& X_vec, const int& n_rows, const int& n_cols){
  mat res = reshape(X_vec,n_rows,n_cols);
  return res;
};

//[[Rcpp::export]]
mat test_X_update(const mat& current_X, const vec& dpi, const vec& T, const double& w){
  vec update = w*w/2.0*dpi + w*T;
  mat res(current_X.n_rows, current_X.n_cols);
  //res = trans(current_X.t() + trans(reshape(update.subvec(0,current_X.n_elem-1),current_X.n_rows,current_X.n_cols))); //match numerator layout
  res = current_X + reshape(update.subvec(0,current_X.n_elem-1),current_X.n_rows,current_X.n_cols);
  return res;
}

//[[Rcpp::export]]
double test_s2_update(const double& current_s2, const double& dpi, const double& T, const double& w){
  double res = current_s2 + w*w/2.0*dpi + w*T;
  return res;
}

void update_MALA(hybrid_paras& paras){// if  h < min{1,xxx}
  // current value of X (X_l)
  mat X_tmp = paras.MALA_par.X;
  double s2_tmp = paras.MALA_par.s2;

  // random normal vector t
  vec T(X_tmp.n_elem+1, fill::randn);

  // update
  vec update = paras.w*paras.w/2.0*paras.dpi + paras.w*T;
  if(!paras.U_fix){
    //paras.MALA_par.X = trans(X_tmp.t() + trans(reshape(update.subvec(0,X_tmp.n_elem-1),X_tmp.n_rows,X_tmp.n_cols))); //match numerator layout
    paras.MALA_par.X = X_tmp + reshape(update.subvec(0,X_tmp.n_elem-1),X_tmp.n_rows,X_tmp.n_cols); //match numerator layout
    paras.MALA_par.U = polar_expansion(paras.MALA_par.X);
    paras.MALA_par.Ut = paras.MALA_par.U.t();
  }
  if(!paras.s2_fix){
    paras.MALA_par.s2 = s2_tmp + update[X_tmp.n_elem];
  }
};

void update_pi(hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp){
  paras.pi = paras.pi_func(paras.MALA_par.wh_vec, paras.MALA_par.s2, paras.MALA_par.X, dat, comp);
};

void update_dpi(hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp){
  paras.dpi = paras.dpi_func(paras.MALA_par.wh_vec, paras.MALA_par.s2, paras.MALA_par.d_UtYU, paras.MALA_par.U, paras.MALA_par.X, dat, comp);
};

void update_Lambda_Y(hybrid_paras& paras, hybrid_dat& dat){
  paras.Lambda.Lambda_flat_t = paras.gibbs_Lambda_t_Y(dat, paras.MALA_par.d_UtYU, paras.t2_lambda, paras.MALA_par.s2);
  paras.Lambda.Lambda_flat = paras.Lambda.Lambda_flat_t.t();
  paras.Lambda.accu_L2 = accu(paras.Lambda.Lambda_flat%paras.Lambda.Lambda_flat);
};

void update_t2_lambda(hybrid_paras& paras, hybrid_comp_par& comp_par){
  paras.t2_lambda = paras.gibbs_t2_lambda(comp_par, paras.Lambda.accu_L2);
};

void update_dat_llk(hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp_par){
  paras.llk = paras.dat_llk_func(paras.MALA_par.U, paras.MALA_par.Ut, paras.Lambda.Lambda_flat_t, paras.MALA_par.s2, dat, comp_par);
};

void update_d_UtYU(hybrid_paras& paras, hybrid_dat& dat){
  paras.MALA_par.d_UtYU = paras.d_UtYU_func(paras.MALA_par.U, paras.MALA_par.Ut, dat);
};

void update_wh(hybrid_paras& paras, hybrid_dat& dat){
  paras.MALA_par.wh_vec = paras.eval_wh_func(paras.MALA_par.s2, paras.MALA_par.d_UtYU, dat);
};

void update_subnet(hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp_par){
  paras.MALA_par.subnet = paras.subnet_func(paras.MALA_par.U, paras.MALA_par.Ut, dat, comp_par);
};

//[[Rcpp::export]]
List test_MALA_one_step(const mat& X_n, const double& s2_n, const vec& log_dpi_n, const double& log_pi_n, const double& w, const int& N, const int& q, const int& M, const mat& X, const mat& U, const double& s2, const cube& Y, const vec& SGI_nodes, const vec& SGI_wts, const bool& weighted, const double& nu0, const double& s20, const double& eta0, const double& t20, const mat& dSdX){
  int Nq = N*q;
  mat I_N = eye(N,N);
  mat I_q = eye(q,q);
  mat d_UtYU = compute_d_UtYU(U, U.t(), Y, M, q);
  vec wh_vec = compute_wh(M, N, q, s2, d_UtYU, eta0, t20, SGI_nodes, SGI_wts, weighted);
  double pi = eval_log_l_Theta_Y(wh_vec, M, N, s2, accu(Y%Y), X, nu0, s20);
  vec dpi = eval_dpi_check(wh_vec, s2, d_UtYU, U, X, SGI_nodes, weighted, Y, M, N, q, I_N, I_q, dSdX, q*q, nu0, s20, accu(Y%Y));

  vec diff(Nq+1, fill::zeros);
  //diff.subvec(0,Nq-1) = vectorise(trans(X - X_n));
  diff.subvec(0,Nq-1) = vectorise(X - X_n);
  diff[Nq] = s2 - s2_n;

  // compute −|(|Z−X_l−w^2 ∇logπ_X (X_l )/2|)|^2/2w^2  (p_Xn)
  vec p_Xn_vec = diff - w*w*log_dpi_n/2.0;
  double p_Xn = -1.0*accu(p_Xn_vec%p_Xn_vec)/(2.0*w*w);

  // compute −|(|X_l−Z−w^2 ∇logπ_X (Z)/2|)|^2/2w^2  (p_Z)
  vec p_Z_vec = -1.0*diff - w*w*dpi/2.0;
  double p_Z = -1.0*accu(p_Z_vec%p_Z_vec)/(2.0*w*w);

  // decide whether or not to update
  double thres = pi - log_pi_n - p_Xn + p_Z;

  return List::create(Named("p_Xn_vec") = p_Xn_vec,
                      Named("p_Xn") = p_Xn,
                      Named("p_Z_vec") = p_Z_vec,
                      Named("p_Z")=p_Z,
                      Named("thres")=thres);
}

bool MALA_one_step(hybrid_paras& paras,hybrid_dat& dat, hybrid_comp_par& comp){
  // current values (X_l, s2_l, log_pi(), log_dpi())
  mat X_n = paras.MALA_par.X; //(X_l)
  double s2_n = paras.MALA_par.s2; //(s2_l)
  vec log_dpi_n = paras.dpi; //∇logπ_X ()
  double log_pi_n = paras.pi; //logπ_X ()
  int Nq = paras.MALA_par.X.n_elem;

  // compute V (Z, log_pi(Z), log_dpi(Z))
  update_MALA(paras);// Z
  update_d_UtYU(paras, dat);
  update_wh(paras,dat);
  update_pi(paras, dat, comp);
  update_dpi(paras, dat, comp);

  vec diff(Nq+1, fill::zeros);
  //diff.subvec(0,Nq-1) = vectorise(trans(paras.MALA_par.X - X_n));
  diff.subvec(0,Nq-1) = vectorise(paras.MALA_par.X - X_n);
  diff[Nq] = paras.MALA_par.s2 - s2_n;

  // compute −|(|Z−X_l−w^2 ∇logπ_X (X_l )/2|)|^2/2w^2  (p_Xn)
  vec p_Xn_vec = diff - paras.w*paras.w*log_dpi_n/2.0;
  double p_Xn = -1.0*accu(p_Xn_vec%p_Xn_vec)/(2.0*paras.w*paras.w);

  // compute −|(|X_l−Z−w^2 ∇logπ_X (Z)/2|)|^2/2w^2  (p_Z)
  vec p_Z_vec = -1.0*diff - paras.w*paras.w*paras.dpi/2.0;
  double p_Z = -1.0*accu(p_Z_vec%p_Z_vec)/(2.0*paras.w*paras.w);

  // decide whether or not to update
  double log_unif = log(R::runif(0,1));
  //double thres = paras.pi - log_pi_Xn + p_Xn - p_Z;
  double thres = paras.pi - log_pi_n - p_Xn + p_Z;
  if(log_unif<thres){
    return true;
  }else{
    return false;
  }
};

// hybrid MCMC updates with paras and data struct
List hybrid_MALA_g1_updates(MALA_paras& MALA0, Lambda_paras& Lambda0, double& t2_lambda_0, //initial values
                    hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp,  //data and parameters
                    int mcmc_sample, double stepsize, int acpt_step, int Gibbs_step, int MALA_step, double target_acpt, bool tune){ // other MCMC parameters
  paras.w = stepsize;
  int acpt_iters = mcmc_sample/(acpt_step*MALA_step); //number of times we check acpt rate
  mat acpt_rate = zeros<mat>(acpt_iters,3); //track acceptance rate
  int acpt_ct = 0;
  int k = 0;
  double stepsize0 = stepsize;
  int k_MALA=0; //track MALA step
  int k_Gibbs=0; //track Gibbs step

  //cout << "initialize parameters \n";
  // Tracking Gibbs Updates
  cube Lambda_flat_list=zeros<cube>(dat.M, dat.q, mcmc_sample/Gibbs_step);
  vec t2_lambda_list = zeros<vec>(mcmc_sample/Gibbs_step);

  // track MALA Updates
  cube U_list=zeros<cube>(dat.N, dat.q, mcmc_sample/MALA_step);
  cube X_list = zeros<cube>(dat.N, dat.q, mcmc_sample/MALA_step);
  vec s2_list = zeros<vec>(mcmc_sample/MALA_step);
  vec llk_list = zeros<vec>(mcmc_sample/MALA_step);
  mat llk_grad_list = zeros<mat>(dat.N*dat.q+1, mcmc_sample/MALA_step);
  vec dat_llk_list = zeros<vec>(mcmc_sample);
  mat subnet_list = zeros<mat>((dat.region_select.n_elem)*(dat.region_select.n_elem-1)/2,mcmc_sample);
  std::cout << size(subnet_list) << std::endl;

  if(mcmc_sample>0){
    // initialize
    paras.Lambda.Lambda_flat = Lambda0.Lambda_flat;
    paras.Lambda.Lambda_flat_t = Lambda0.Lambda_flat_t;
    paras.Lambda.accu_L2 = Lambda0.accu_L2;

    paras.t2_lambda = t2_lambda_0;

    paras.MALA_par.X = MALA0.X;
    paras.MALA_par.U = MALA0.U;
    paras.MALA_par.Ut = MALA0.Ut;
    paras.MALA_par.d_UtYU = MALA0.d_UtYU;
    paras.MALA_par.s2 = MALA0.s2;

    update_wh(paras, dat);
    update_pi(paras, dat, comp); //log_pi()
    update_dpi(paras, dat, comp);//dlog_pi()
    update_dat_llk(paras, dat, comp); //update data llk
    std::cout << "here 1" << std::endl;
    update_subnet(paras, dat, comp); //update select subnetwork
    std::cout << "here 2" << std::endl;

    //double current_pi=paras.pi;
    //mat current_dpi=paras.dpi_2;
    double current_pi;
    mat current_dpi;
    mat current_d_UtYU;
    mat current_X;
    mat current_U;
    mat current_Ut;
    double current_s2;
    vec current_wh_vec;
    mat U_neg; // for sign flipping comparison
    mat U_flipped; // sign flipped U
    umat flip_sign_mat(dat.N, dat.q);
    rowvec current_diff, neg_diff;
    urowvec which_flip;
    urowvec ones_vec(dat.q, fill::value(1));

    // collect samples
    Lambda_flat_list.slice(0) = paras.Lambda.Lambda_flat;
    s2_list[0] = paras.MALA_par.s2;
    t2_lambda_list[0] = paras.t2_lambda;
    U_list.slice(0) = paras.MALA_par.U;
    X_list.slice(0) = paras.MALA_par.X;

    //std::cout<< "initialized" << std::endl;

    for(int iter=0;iter<mcmc_sample;iter++){
      // update parameter
      if((iter+1)%Gibbs_step==0){
        if(!paras.Lambda_fix){
          update_Lambda_Y(paras, dat);
        }
        if(!paras.t2_lambda_fix){
          update_t2_lambda(paras, comp); //t2_lambda[iter]
        }
        Lambda_flat_list.slice(k_Gibbs) = paras.Lambda.Lambda_flat;
        t2_lambda_list[k_Gibbs] = paras.t2_lambda;
        k_Gibbs++;
      }

      //std::cout << "g1-2 done" << std::endl;
      current_U=paras.MALA_par.U; //U[iter-1]
      current_Ut=paras.MALA_par.Ut; //U[iter-1]
      current_X=paras.MALA_par.X; //X[iter-1]
      current_s2=paras.MALA_par.s2; //s2[iter-1]
      current_pi=paras.pi; //pi(X[iter-1], U[iter-1], s2[iter-1])
      current_dpi=paras.dpi; //dpi(X[iter-1], U[iter-1], s2[iter-1])
      current_d_UtYU=paras.MALA_par.d_UtYU;
      current_wh_vec = paras.MALA_par.wh_vec;

      // adjust acceptance rate
      if((iter+1)%(MALA_step*acpt_step)==0 & iter > 0){
        if(!paras.U_fix | !paras.s2_fix){
          if(k < acpt_iters){
            acpt_rate(k,0) = iter+1;
            acpt_rate(k,1) = acpt_ct*1.0/acpt_step;///MALA_step;
            if(tune){
              paras.w = adjust_acceptance(acpt_rate(k,1),paras.w,target_acpt); //improper update 03/21/22
             }
            acpt_rate(k,2) = paras.w;
          }
          k++;
          acpt_ct = 0;
        }
      }

      // update U with MALA
      if((iter+1)%MALA_step==0){
        if(!paras.U_fix | !paras.s2_fix){
          //MALA = MALA_one_step(paras,dat,comp);
          //std::cout << "MALA pre" << std::endl;
          if(MALA_one_step(paras,dat,comp)){ // if accept the new X value
            //std::cout << "MALA accept" << std::endl;
            // 05/17/22 sign flip
            U_neg = -1.0*paras.MALA_par.U;
            current_diff = sum(square(current_U - paras.MALA_par.U),0);
            neg_diff = sum(square(current_U - U_neg),0);
            which_flip = -1*(neg_diff < current_diff);
            which_flip += (neg_diff >= current_diff);

            if(sum(which_flip != ones_vec) > 0){
              std::cout << "flip sign at iter: " << iter+1 << std::endl;
              std::cout << "sign_flips = " << which_flip << std::endl;
              flip_sign_mat.each_row() = which_flip;
              U_flipped = paras.MALA_par.U % flip_sign_mat;
              paras.MALA_par.U = U_flipped;
              paras.MALA_par.Ut = U_flipped.t();
              update_d_UtYU(paras, dat);
              update_wh(paras,dat);
              update_pi(paras, dat, comp); //pi(X[iter], U_flipped[iter], s2[iter])
              update_dpi(paras, dat, comp); //dpi(X[iter], U_flipped[iter], s2[iter])
            }
            acpt_ct++;
            //cout << paras.pi << endl;
          }else{// if not accept
            //std::cout<< "MALA reject" << std::endl;
            paras.MALA_par.U=current_U;
            paras.MALA_par.Ut=current_Ut;
            paras.MALA_par.X=current_X;
            paras.MALA_par.s2=current_s2;
            paras.pi = current_pi;
            paras.dpi = current_dpi;
            paras.MALA_par.d_UtYU = current_d_UtYU;
            paras.MALA_par.wh_vec = current_wh_vec;
          }
        }
        U_list.slice(k_MALA) = paras.MALA_par.U;
        X_list.slice(k_MALA) = paras.MALA_par.X;
        s2_list[k_MALA] = paras.MALA_par.s2;
        llk_list[k_MALA] = paras.pi;
        llk_grad_list.col(k_MALA) = paras.dpi;
        k_MALA++;
      }
      update_dat_llk(paras, dat, comp); // data llk
      update_subnet(paras, dat, comp); //update select subnetwork
      dat_llk_list[iter] = paras.llk;
      subnet_list.col(iter) = paras.MALA_par.subnet;
    }
  }
  return List::create(Named("Lambda_flat") = Lambda_flat_list,
                      Named("s2") = s2_list,
                      Named("t2_lambda") = t2_lambda_list,
                      Named("X")=X_list,
                      Named("U")=U_list,
                      Named("llk")=llk_list,
                      Named("llk_grad")=llk_grad_list,
                      Named("dat_llk")=dat_llk_list,
                      Named("subnet_recon")=subnet_list,
                      Named("acceptance_rate")=acpt_rate);
};

// [[Rcpp::export]]
List hybrid_MALA_g1(mat& X_0, mat& Lambda_0_flat, double& s2_0, double& t2_lambda_0, //initial values
            double& nu0, double& s20, double& eta0, double& t20, //paras
            cube& Y, int& M, int& N, int& q, vec& SGI_nodes, vec& SGI_wts, bool& weighted, uvec& region_select,//data
            int mcmc_sample, double stepsize, int acpt_step, int Gibbs_step, int MALA_step, double target_acpt, bool tune, CharacterVector fixed){ // other MCMC parameters
  // build hybrid data structure
  hybrid_dat dat;
  dat.Y = Y;
  dat.M = M;
  dat.N = N;
  dat.q = q;
  dat.nu0 = nu0;
  dat.s20 = s20;
  dat.eta0 = eta0;
  dat.t20 = t20;
  dat.SGI_nodes = SGI_nodes;
  dat.SGI_wts = SGI_wts;
  dat.weighted = weighted;
  dat.region_select = region_select-1;

  // build hybrid computational paramerters
  hybrid_comp_par comp;
  comp.I_q = eye(q,q);
  comp.I_N = eye(N,N);
  comp.dSdX = box_prod_2(comp.I_N, comp.I_q, N, N, q, q);
  comp.q2 = q*q;
  comp.accu_YY = accu(Y%Y);
  mat Y_flat_mat(N*N, M);
  mat Y_i;
  for(int i=0; i<M; i++){
    Y_i = Y.slice(i);
    //Y_flat_mat.col(i) = Y_i(comp.Y_tril_ind);
    Y_flat_mat.col(i) = vectorise(Y_i);
  }
  comp.Y_flat = Y_flat_mat;
  comp.t2_lambda_A = (nu0+M*q)/2.0;
  comp.t2_lambda_B = nu0*t20;
  mat A(region_select.n_elem,region_select.n_elem,fill::zeros);
  comp.region_upper_idx = trimatu_ind(size(A), 1);

  // build hybrid parameters
  Lambda_paras Lambda_par0;
  Lambda_par0.Lambda_flat = Lambda_0_flat;
  Lambda_par0.Lambda_flat_t = Lambda_0_flat.t();
  Lambda_par0.accu_L2 = accu(Lambda_0_flat % Lambda_0_flat);

  MALA_paras MALA0_st;
  MALA0_st.X=X_0;
  MALA0_st.U=polar_expansion(X_0);
  MALA0_st.Ut = MALA0_st.U.t();
  MALA0_st.d_UtYU = compute_d_UtYU(MALA0_st.U, MALA0_st.Ut, Y, M, q);
  MALA0_st.s2 = s2_0;

  hybrid_paras paras;
  paras.gibbs_Lambda_t_Y = gibbs_Lambda_t_Y;
  paras.gibbs_t2_lambda = gibbs_t2_lambda;
  paras.dpi_func = eval_dpi;
  paras.pi_func = eval_pi;
  paras.dat_llk_func = eval_dat_llk;
  paras.d_UtYU_func = eval_d_UtYU;
  paras.eval_wh_func = eval_wh;
  paras.subnet_func = eval_subnetwork;

  CharacterVector lambda_vec = {"lambda"};
  CharacterVector s2_vec = {"s2"};
  CharacterVector t2_lambda_vec = {"t2_lambda"};
  CharacterVector U_vec = {"U"};

  LogicalVector lambda_fixed = in(lambda_vec,fixed);
  LogicalVector s2_fixed = in(s2_vec,fixed);
  LogicalVector t2_lambda_fixed = in(t2_lambda_vec,fixed);
  LogicalVector U_fixed = in(U_vec,fixed);

  paras.U_fix = U_fixed[0];
  paras.s2_fix = s2_fixed[0];
  paras.Lambda_fix = lambda_fixed[0];
  paras.t2_lambda_fix = t2_lambda_fixed[0];
  std::cout<< "U is fixed " << paras.U_fix <<std::endl;
  std::cout<< "s2 is fixed " <<paras.s2_fix <<std::endl;
  std::cout<< "Lambda is fixed " << paras.Lambda_fix <<std::endl;
  std::cout<< "t2_lambda is fixed " << paras.t2_lambda_fix <<std::endl;

  //run hybrid
  List hybrid_update_list = hybrid_MALA_g1_updates(MALA0_st, Lambda_par0, t2_lambda_0, paras, dat, comp, mcmc_sample, stepsize, acpt_step,
    Gibbs_step, MALA_step, target_acpt, tune);
  return List::create(Named("mcmc")=hybrid_update_list,
                      Named("stepsize")=paras.w);
};
