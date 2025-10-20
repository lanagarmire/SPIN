// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]

#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double s2_new(const double& s2_A, const double& s2_B, const double& s_L_f_d_UtYU, const double& accu_L2, const double& accu_YY){
  double new_s2 = 1.0/randg<double>(distr_param(s2_A, 2.0/(s2_B+accu_YY+accu_L2-2*s_L_f_d_UtYU)));
  return new_s2;
};

// [[Rcpp::export]]
double t2_new(const double& t2_A, const double& t2_B, const mat& C, const mat& Z, const mat& Lambda_flat_beta, const mat& Z_alpha){
  mat C_res = C-Lambda_flat_beta-Z_alpha;
  //cout << accu(C_res%C_res) << endl;
  double new_t2 = 1.0/randg<double>(distr_param(t2_A, 2.0/(t2_B+accu(C_res%C_res))));
  return new_t2;
};

// [[Rcpp::export]]
double log_prior_Lambda(const mat& Lambda_flat, const double& t2_lambda, const int& M, const int& q){
  mat lambda_mean(M,q,fill::zeros);
  mat lambda_sd(M,q,fill::value(sqrt(t2_lambda)));
  return accu(log_normpdf(Lambda_flat, lambda_mean, lambda_sd));
};

// [[Rcpp::export]]
mat Lambda_new_t(const mat& I_q, const int& q, const int& M, const mat& d_UtYU, const mat& Z_alpha, const mat& C, const cube& Y, const double& t2_lambda, const vec& beta, const double& t2, const double& s2){
  mat beta_s_mat = reshape(beta,q,1);
  mat var_mat = beta_s_mat*beta_s_mat.t()/t2;
  var_mat.diag() +=1.0/t2_lambda+1.0/s2;
  var_mat = solve(var_mat, I_q);

  //mat Ut = U.t();
  //mat mu_mat = d_UtYU;
  //for(int i=0; i<M; i++){
  //  mu_mat.col(i) = diagvec(Ut*Y.slice(i)*U);
  //}
  //mu_mat /= s2;
  //mu_mat += 1.0/t2*beta*trans(C-Z*alpha);
  mat mu_mat = d_UtYU/s2 + 1.0/t2*beta*trans(C-Z_alpha);

  mat Z_0 = randn(q,M);
  mat A = chol(var_mat, "lower");
  mat Lambda_flat_new = A*Z_0 + var_mat*mu_mat;

  //return Lambda_flat_new;
  return Lambda_flat_new;
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
vec d_new(const int& r, const int& q, const int& W_dim, const mat& I_W, const mat& Z, const mat& C_mat, const double& t2_beta, const double& t2_alpha, const double& t2, const mat& Lambda_flat){
  mat X = join_rows(Lambda_flat, Z);
  mat Xt = X.t();

  vec t2_beta_vec(q, fill::value(t2_beta));
  vec t2_alpha_vec(r, fill::value(t2_alpha));

  mat var_mat = 1.0/t2*Xt*X;
  var_mat.diag() += 1.0/join_cols(t2_beta_vec, t2_alpha_vec);
  var_mat = solve(var_mat, I_W);
  mat XC = X%C_mat;
  vec mu = var_mat*1.0/t2*sum(XC.t(),1);

  mat d_mat = mvnrnd(mu, var_mat);
  vec d_res = d_mat.col(0);
  return d_res;
};

// [[Rcpp::export]]
double t2_lambda_new(const double& t2_lambda_A, const double& t2_lambda_B, const double& accu_L2){
  double t2_lambda = 1.0/randg<double>(distr_param(t2_lambda_A, 2.0/(t2_lambda_B + accu_L2)));
  return t2_lambda;
};

// [[Rcpp::export]]
double t2_beta_new(double const& t2_beta_A, double const& t2_beta_B, const vec& beta){
  double t2_beta = 1.0/randg<double>(distr_param(t2_beta_A,2.0/(t2_beta_B+sum(beta%beta))));
  return t2_beta;
};

// [[Rcpp::export]]
double t2_alpha_new(const double& t2_alpha_A, const double& t2_alpha_B, const vec& alpha){
  double t2_alpha = 1.0/randg<double>(distr_param(t2_alpha_A,2.0/(t2_alpha_B+sum(alpha%alpha))));
  return t2_alpha;
};

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

  mat Lambda_i, Y_i, UtYi;
  mat res_temp_1(q2,Nq,fill::zeros);
  mat res_temp_2(q,N,fill::zeros);
  for(int epoch=0;epoch<M;epoch++){
    Lambda_i = diagmat(Lambda_flat_t.col(epoch));
    Y_i = Y.slice(epoch);
    UtYi = Ut*Y_i;
    res_temp_1 += kron(Lambda_i, UtYi);
    res_temp_2 += Lambda_i*UtYi;
  }

  mat res_temp = res_temp_1.cols(box_idx) + kron(res_temp_2,I_q);
  mat res = -1.0*X + reshape(I_q_flat*1.0/s2*res_temp*(kron(I_N,R) + kron(X,I_q)*(solve(kron(I_q,R)+kron(R,I_q),eye(q2,q2))*(-1.0*kron(P,P))*(kron(I_q,S)*dSdX+kron(S,I_q)))), N, q);
  return res;
};

// [[Rcpp::export]]
double Y_C_llk_new_old(const mat& Y_flat, const uvec& Y_tril_ind, const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2,
                   const mat& C, const vec& beta, const vec& alpha, const mat& Z, const double& t2, const mat& Lambda_flat_beta,
                   const mat& Z_alpha, const int& M, const int& N){
  double dat_llk = 0.0;
  //mat Ut= U.t();
  //mat Lambda_flat_t=Lambda_flat.t();
  mat Y_mean;
  mat Y_mean_flat(N*(N+1)/2,M);

  mat C_mean = Lambda_flat_beta + Z_alpha;
  double t = sqrt(t2);
  mat C_smat(M, 1, fill::value(t));
  dat_llk = accu(log_normpdf(C, C_mean, C_smat));

  double s = sqrt(s2);
  mat Y_smat(N*(N+1)/2,M, fill::value(s));
  for(int i = 0; i < M; i++){
    Y_mean = U*diagmat(Lambda_flat_t.col(i))*Ut;
    Y_mean_flat.col(i) = Y_mean(Y_tril_ind);
  }

  dat_llk += accu(log_normpdf(Y_flat, Y_mean_flat, Y_smat));
  return dat_llk;
};

// [[Rcpp::export]]
double Y_C_llk_new(const mat& Y_flat, const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2,
                   const mat& C, const vec& beta, const vec& alpha, const mat& Z, const double& t2, const mat& Lambda_flat_beta,
                   const mat& Z_alpha, const int& M, const int& N){
  mat C_mean = Lambda_flat_beta + Z_alpha;
  double t = sqrt(t2);
  mat C_smat(M, 1, fill::value(t));
  double dat_llk = accu(log_normpdf(C, C_mean, C_smat));

  mat Y_mean;
  mat Y_mean_flat(N*N,M);
  double s = sqrt(s2);
  mat Y_smat(N*N,M, fill::value(s));
  for(int i = 0; i < M; i++){
    Y_mean = U*diagmat(Lambda_flat_t.col(i))*Ut;
    Y_mean_flat.col(i) = vectorise(Y_mean);
  }

  dat_llk += accu(log_normpdf(Y_flat, Y_mean_flat, Y_smat));
  dat_llk += M*(0.5*N*N*log(2*datum::pi)-0.25*N*(N+1)*log(2*datum::pi)+0.25*N*(N-1)*log(2)+0.5*N*N*log(s2)-0.25*N*(N+1)*log(s2));
  return dat_llk;
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
double C_llk_new(const mat& C, const double& t2, const mat& Lambda_flat_beta, const mat& Z_alpha, const int& M, bool scale){
  mat C_mean = Lambda_flat_beta + Z_alpha;
  double dat_llk;
  if(scale){
    double t = sqrt(t2);
    mat C_smat(M, 1, fill::value(t));
    dat_llk = accu(log_normpdf(C, C_mean, C_smat));
  }else{
    dat_llk = -0.5*M*log(t2) - 0.5/t2*accu(pow(C-C_mean,2));
  }

  return dat_llk;
}

// [[Rcpp::export]]
double Y_llk_i_new(const mat& Y, const mat& U, const mat& Ut, const mat& Lambda, const double& s2, const int& N){
  double s = sqrt(s2);
  vec Y_svec(N*N, fill::value(s));
  double dat_llk = accu(log_normpdf(vectorise(Y), vectorise(U*Lambda*Ut), Y_svec));
  return dat_llk;
};

// [[Rcpp::export]]
mat compute_d_UtYU(const mat& U, const mat& Ut, const cube& Y, const int& M, const int& q){
  mat res(q,M);
  for(int i=0; i<M; i++){
    res.col(i) = diagvec(Ut*Y.slice(i)*U);
  }
  return res;
};

// [[Rcpp::export]]
double sum_L_f_d_UtYU(const mat& Lambda_flat, const mat& d_UtYU, const int& M){
  mat tmp(1,1,fill::zeros);

  for(int i=0; i<M; i++){
    tmp += Lambda_flat.row(i)*d_UtYU.col(i);
  }

  return as_scalar(tmp);
};

// [[Rcpp::export]]
double eval_log_mat_norm_new(const mat& X, const double& s_L_f_d_UtYU, const double& s2){
  //double res = s_L_f_d_UtYU/s2 - 0.5*accu(X%X);
  //return res;
  return s_L_f_d_UtYU/s2 - 0.5*accu(X%X);
};

// [[Rcpp::export]]
double grad_C_t2(const double& t2, const int& q, const mat& C, const int& M, const mat& Lambda_flat_beta, const mat& Z_alpha){
  return -0.5*M/t2+0.5/(pow(t2,2))*accu(pow(C-Lambda_flat_beta-Z_alpha,2));
};

// [[Rcpp::export]]
vec grad_C_beta(const mat& Lambda_flat, const mat& C, const mat& Lambda_flat_beta, const mat& Z_alpha, const double& t2, const int& q){
  mat scale = repmat(C-Lambda_flat_beta-Z_alpha,1,q);
  return sum(Lambda_flat.t() % scale.t(),1)/t2;
};

// [[Rcpp::export]]
vec grad_C_alpha(const mat& Z, const mat& C, const mat& Lambda_flat_beta, const mat& Z_alpha, const double& t2, const int& r){
  mat scale = repmat(C-Lambda_flat_beta-Z_alpha,1,r);
  return sum(Z.t() % scale.t(),1)/t2;
};

// [[Rcpp::export]]
vec grad_C_d(const mat& Lambda_flat, const mat& Z, const mat& C, const mat& Lambda_flat_beta, const mat& Z_alpha, const int& q, const int& r, const double& t2){
  vec grad_d(q+r);
  grad_d.subvec(0,q-1) = grad_C_beta(Lambda_flat, C, Lambda_flat_beta, Z_alpha, t2, q);
  grad_d.subvec(q,q+r-1) = grad_C_alpha(Z, C, Lambda_flat_beta, Z_alpha, t2, r);
  return grad_d;
};

// [[Rcpp::export]]
vec grad_C_lambda(const vec& beta, const mat& Z, const mat& C, const mat& Lambda_flat_beta, const mat& Z_alpha, const double& t2, const int& M, const int& q){
  mat grad_lambda_mat(q,M);
  mat beta_mat(beta);
  grad_lambda_mat = beta_mat*(C-Lambda_flat_beta-Z_alpha).t()/t2;
  return vectorise(grad_lambda_mat);
}

// data structure
struct hybrid_dat{
  // data; hyperparameters
  cube Y; // N x N x M
  mat C; //M x 1
  mat Z; //M x r
  int M; // population size
  int N; // # roi
  int q; // # reduced dimensions
  int r; // # covariates
  double nu0; //s2
  double s20; //s2
  double rho0; //t2
  double psi20; //t2
  double eta0; //t2_lambda
  double t20; //t2_lambda
  double gamma0; //t2_beta
  double kappa20; //t2_beta
  double omega0; //t2_alpha
  double phi20; //t2_alpha
};

struct hybrid_comp_par{
  double s2_A;
  double s2_B;
  double t2_A;
  double t2_B;
  mat I_q;
  mat dSdX;
  mat I_Nq;
  mat I_N;
  int q2;
  int Nq;
  uvec box_idx;
  mat I_q_flat;
  mat I_W;
  int W_dim;
  double t2_lambda_A;
  double t2_lambda_B;
  double t2_beta_A;
  double t2_beta_B;
  double t2_alpha_A;
  double t2_alpha_B;
  double accu_YY; //accu(Y%Y)
  mat C_mat; //repmat(C,1,W_dim)
  //mat Y_flat; //(N*(N+1)/2 x M)
  mat Y_flat; //(N*N x M)
  uvec Y_tril_ind;
};

// gibbs parameter update function using data struct as input
double gibbs_s2(const hybrid_comp_par& comp_par, const double& s_L_f_d_UtYU, const double& accu_L2){
  return s2_new(comp_par.s2_A, comp_par.s2_B, s_L_f_d_UtYU, accu_L2, comp_par.accu_YY);
}

double gibbs_t2(const hybrid_comp_par& comp_par, const hybrid_dat& dat, const vec& beta, const mat& Lambda_flat_beta, const mat& Z_alpha){
  return t2_new(comp_par.t2_A, comp_par.t2_B, dat.C, dat.Z, Lambda_flat_beta, Z_alpha);
}

mat gibbs_Lambda_t(const hybrid_comp_par& comp_par, const hybrid_dat& dat, const mat& d_UtYU, const mat& Z_alpha, const double& t2_lambda, const vec& beta, const double& t2, const double& s2){
  return Lambda_new_t(comp_par.I_q, dat.q, dat.M, d_UtYU, Z_alpha, dat.C, dat.Y, t2_lambda, beta, t2, s2);
}

mat gibbs_Lambda_t_Y(const hybrid_dat& dat, const mat& d_UtYU, const double& t2_lambda, const double& s2){
  return Lambda_new_t_Y(dat.q, dat.M, d_UtYU, dat.Y, t2_lambda, s2);
}

vec gibbs_d(const hybrid_dat& dat, const hybrid_comp_par& comp_par, const double& t2_beta, const double& t2_alpha, const double& t2, const mat& Lambda_flat){
  return d_new(dat.r, dat.q, comp_par.W_dim, comp_par.I_W, dat.Z, comp_par.C_mat, t2_beta, t2_alpha, t2, Lambda_flat);
}

double gibbs_t2_lambda(const hybrid_comp_par& comp_par, const double& accu_L2){
  return t2_lambda_new(comp_par.t2_lambda_A, comp_par.t2_lambda_B, accu_L2);
}

double gibbs_t2_beta(const hybrid_comp_par& comp_par, const vec& beta){
  return t2_beta_new(comp_par.t2_beta_A, comp_par.t2_beta_B, beta);
}

double gibbs_t2_alpha(const hybrid_comp_par& comp_par, const vec& alpha){
  return t2_alpha_new(comp_par.t2_alpha_A, comp_par.t2_alpha_B, alpha);
}

mat eval_dpi_2(const mat& X, const mat& U, const mat& Ut, const mat& Lambda_flat, const mat& Lambda_flat_t, const double& s2, const hybrid_dat& dat, const hybrid_comp_par& comp_par){
  return eval_dlog_mat_norm_2_new(X, U, Ut, comp_par.I_q, comp_par.dSdX, comp_par.I_Nq, comp_par.I_N, Lambda_flat, Lambda_flat_t,
                                  dat.Y, comp_par.q2, comp_par.Nq, comp_par.box_idx, comp_par.I_q_flat, s2, dat.N, dat.q, dat.M);
}

double eval_dat_llk(const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2, const vec& beta, const vec& alpha,
                    const double& t2, const mat& Lambda_flat_beta, const mat& Z_alpha, const hybrid_dat& dat, const hybrid_comp_par& comp_par){
  return Y_C_llk_new(comp_par.Y_flat, U, Ut, Lambda_flat_t, s2, dat.C, beta, alpha, dat.Z, t2, Lambda_flat_beta,
                     Z_alpha, dat.M, dat.N);
}

mat eval_d_UtYU(const mat& U, const mat& Ut, const hybrid_dat& dat){
  return compute_d_UtYU(U, Ut, dat.Y, dat.M, dat.q);
}

double eval_sum_L_f_d_UtYU(const mat& Lambda_flat, const mat& d_UtYU, const hybrid_dat& dat){
  return sum_L_f_d_UtYU(Lambda_flat, d_UtYU, dat.M);
}

double eval_pi(const mat& X, const double& s_L_f_d_UtYU, const double& s2){
  return eval_log_mat_norm_new(X, s_L_f_d_UtYU, s2);
}

// pointers to gibbs parameter update function using data struct as input
typedef double(*funcPtr_s2)(const hybrid_comp_par& comp_par, const double& s_L_f_d_UtYU, const double& accu_L2);
typedef double(*funcPtr_t2)(const hybrid_comp_par& comp_par, const hybrid_dat& dat, const vec& beta, const mat& Lambda_flat_beta, const mat& Z_alpha);
typedef mat(*funcPtr_Lambda_t)(const hybrid_comp_par& comp_par, const hybrid_dat& dat, const mat& d_UtYU, const mat& Z_alpha, const double& t2_lambda, const vec& beta, const double& t2, const double& s2);
typedef mat(*funcPtr_Lambda_t_Y)(const hybrid_dat& dat, const mat& d_UtYU, const double& t2_lambda, const double& s2);
typedef vec(*funcPtr_d)(const hybrid_dat& dat, const hybrid_comp_par& comp_par, const double& t2_beta, const double& t2_alpha, const double& t2, const mat& Lambda_flat);
typedef double(*funcPtr_t2_lambda)(const hybrid_comp_par& comp_par, const double& accu_L2);
typedef double(*funcPtr_t2_beta)(const hybrid_comp_par& comp_par, const vec& beta);
typedef double(*funcPtr_t2_alpha)(const hybrid_comp_par& comp_par, const vec& alpha);
typedef mat(*funcPtr_dpi_2)(const mat& X, const mat& U, const mat& Ut, const mat& Lambda_flat, const mat& Lambda_flat_t, const double& s2, const hybrid_dat& dat, const hybrid_comp_par& comp_par);
typedef double(*funcPtr_dat_llk)(const mat& U, const mat& Ut, const mat& Lambda_flat_t, const double& s2, const vec& beta, const vec& alpha, const double& t2, const mat& Lambda_flat_beta, const mat& Z_alpha, const hybrid_dat& dat, const hybrid_comp_par& comp_par);
typedef mat(*funcPtr_d_UtYU)(const mat& U, const mat& Ut, const hybrid_dat& dat);
typedef double(*funcPtr_s_L_f_d_UtYU)(const mat& Lambda_flat, const mat& d_UtYU, const hybrid_dat& dat);
typedef double(*funcPtr_pi)(const mat& X, const double& s_L_f_d_UtYU, const double& s2);

// gibbs parameter struct
struct Lambda_paras{
  mat Lambda_flat; // M x q
  mat Lambda_flat_t; //q x M
  double accu_L2; //accu(Lambda%Lambda)
  double log_m;
};

struct d_paras{
  vec alpha; // r
  vec beta; // q
  vec d; // r+q
  mat Z_alpha; // M x 1
};

struct X_paras{
  mat U; //N x q
  mat X; //N x q
  mat Ut;
  mat d_UtYU;//q x M
};

struct hybrid_paras{
  double s2;
  double t2;
  Lambda_paras Lambda;
  d_paras d;
  double t2_lambda;
  double t2_beta;
  double t2_alpha;
  X_paras X_par; // parameter that needs to be updated
  mat Lambda_flat_beta; // M x 1
  double s_L_f_d_UtYU;
  //MALA stepsize
  double w; //stepsize
  double pi;
  mat dpi_2;
  double llk; // value of data llk
  funcPtr_s2 gibbs_s2;
  funcPtr_t2 gibbs_t2;
  funcPtr_Lambda_t gibbs_Lambda_t;
  funcPtr_Lambda_t_Y gibbs_Lambda_t_Y;
  funcPtr_d gibbs_d;
  funcPtr_t2_lambda gibbs_t2_lambda;
  funcPtr_t2_beta gibbs_t2_beta;
  funcPtr_t2_alpha gibbs_t2_alpha;
  funcPtr_pi pi_func;
  funcPtr_dpi_2 dpi_fun_2;
  funcPtr_dat_llk dat_llk_func; // pointer to data llk function
  funcPtr_d_UtYU d_UtYU_func;
  funcPtr_s_L_f_d_UtYU s_L_f_d_UtYU_func;
};

void update_X(hybrid_paras& paras){// if  h < min{1,xxx}
  // current value of X (X_l)
  mat X_tmp = paras.X_par.X;

  // random normal matrix T
  mat T(size(X_tmp), fill::randn);

  // update
  X_tmp = X_tmp + paras.w*paras.w/2.0*paras.dpi_2 + paras.w*T;
  paras.X_par.X = X_tmp;
  paras.X_par.U = polar_expansion(X_tmp);
  paras.X_par.Ut = paras.X_par.U.t();
};

void update_pi(hybrid_paras& paras){
  paras.pi = paras.pi_func(paras.X_par.X, paras.s_L_f_d_UtYU, paras.s2);
};

void update_dpi_2(hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp){
  paras.dpi_2 = paras.dpi_fun_2(paras.X_par.X, paras.X_par.U, paras.X_par.Ut, paras.Lambda.Lambda_flat, paras.Lambda.Lambda_flat_t, paras.s2, dat, comp);
};

void update_s2(hybrid_paras& paras, hybrid_comp_par& comp_par){
  paras.s2 = paras.gibbs_s2(comp_par, paras.s_L_f_d_UtYU, paras.Lambda.accu_L2);
};

void update_t2(hybrid_comp_par& comp_par, hybrid_dat& dat, hybrid_paras& paras){
  paras.t2 = paras.gibbs_t2(comp_par, dat, paras.d.beta, paras.Lambda_flat_beta, paras.d.Z_alpha);
};

void update_Lambda(hybrid_comp_par& comp_par, hybrid_paras& paras, hybrid_dat& dat){
  paras.Lambda.Lambda_flat_t = paras.gibbs_Lambda_t(comp_par, dat, paras.X_par.d_UtYU, paras.d.Z_alpha, paras.t2_lambda, paras.d.beta, paras.t2, paras.s2);
  paras.Lambda.Lambda_flat = paras.Lambda.Lambda_flat_t.t();
  paras.Lambda.accu_L2 = accu(paras.Lambda.Lambda_flat%paras.Lambda.Lambda_flat);
};

void update_Lambda_Y(hybrid_paras& paras, hybrid_dat& dat){
  paras.Lambda.Lambda_flat_t = paras.gibbs_Lambda_t_Y(dat, paras.X_par.d_UtYU, paras.t2_lambda, paras.s2);
  paras.Lambda.Lambda_flat = paras.Lambda.Lambda_flat_t.t();
  paras.Lambda.accu_L2 = accu(paras.Lambda.Lambda_flat%paras.Lambda.Lambda_flat);
};

void update_d(hybrid_paras& paras, hybrid_comp_par& comp_par, hybrid_dat& dat){
  paras.d.d = paras.gibbs_d(dat, comp_par, paras.t2_beta, paras.t2_alpha, paras.t2, paras.Lambda.Lambda_flat);
  paras.d.beta = paras.d.d.subvec(0,dat.q-1);
  paras.d.alpha = paras.d.d.subvec(dat.q, comp_par.W_dim-1);
  paras.d.Z_alpha = dat.Z*paras.d.alpha;
};

void update_t2_lambda(hybrid_paras& paras, hybrid_comp_par& comp_par){
  paras.t2_lambda = paras.gibbs_t2_lambda(comp_par, paras.Lambda.accu_L2);
};

void update_t2_beta(hybrid_paras& paras, hybrid_comp_par& comp_par){
  paras.t2_beta = paras.gibbs_t2_beta(comp_par, paras.d.beta);
};

void update_t2_alpha(hybrid_paras& paras, hybrid_comp_par& comp_par){
  paras.t2_alpha = paras.gibbs_t2_alpha(comp_par, paras.d.alpha);
};

void update_dat_llk(hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp_par){
  paras.llk = paras.dat_llk_func(paras.X_par.U, paras.X_par.Ut, paras.Lambda.Lambda_flat_t, paras.s2, paras.d.beta, paras.d.alpha, paras.t2, paras.Lambda_flat_beta, paras.d.Z_alpha, dat, comp_par);
};

void update_d_UtYU(hybrid_paras& paras, hybrid_dat& dat){
  paras.X_par.d_UtYU = paras.d_UtYU_func(paras.X_par.U, paras.X_par.Ut, dat);
};

void update_s_L_f_d_UtYU(hybrid_paras& paras, hybrid_dat& dat){
  paras.s_L_f_d_UtYU = paras.s_L_f_d_UtYU_func(paras.Lambda.Lambda_flat, paras.X_par.d_UtYU, dat);
};

bool MALA_one_step(hybrid_paras& paras,hybrid_dat& dat, hybrid_comp_par& comp){
  // current values (X_l, log_pi(X_l), log_dpi(X_l))
  mat X_n = paras.X_par.X; //(X_l)
  mat log_dpi_Xn = paras.dpi_2; //∇logπ_X (X_l )
  double log_pi_Xn = paras.pi; //logπ_X (X_l )

  // compute V (Z, log_pi(Z), log_dpi(Z))
  update_X(paras);// Z
  update_d_UtYU(paras, dat);
  update_s_L_f_d_UtYU(paras, dat);
  update_pi(paras);
  update_dpi_2(paras, dat, comp);

  // compute −|(|Z−X_l−w^2 ∇logπ_X (X_l )/2|)|^2/2w^2  (p_Xn)
  mat p_Xn_mat = paras.X_par.X - X_n - paras.w*paras.w*log_dpi_Xn/2.0;
  double p_Xn = -1.0*accu(p_Xn_mat%p_Xn_mat)/(2.0*paras.w*paras.w);

  // compute −|(|X_l−Z−w^2 ∇logπ_X (Z)/2|)|^2/2w^2  (p_Z)
  mat p_Z_mat = X_n - paras.X_par.X - paras.w*paras.w*paras.dpi_2/2.0;
  double p_Z = -1.0*accu(p_Z_mat%p_Z_mat)/(2.0*paras.w*paras.w);

  // decide whether or not to update
  double log_unif = log(R::runif(0,1));
  //double thres = paras.pi - log_pi_Xn + p_Xn - p_Z;
  double thres = paras.pi - log_pi_Xn - p_Xn + p_Z;
  if(log_unif<thres){
    return true;
  }else{
    return false;
  }
};

bool Lambda_rejection_one_step(hybrid_paras& paras,hybrid_dat& dat){
  update_Lambda_Y(paras, dat);
  // decide whether or not to reject
  //double log_unif = sum(log(R::runif(dat.M,0,1)));
  double log_unif = sum(log(randu(dat.M)));
  double thres = C_llk_new(dat.C, paras.t2, paras.Lambda.Lambda_flat*paras.d.beta, paras.d.Z_alpha, dat.M, false);
  std::cout << "thres = " << thres << endl;
  double current_log_m = paras.Lambda.log_m;
  std::cout << "log_unif = " << log_unif << endl;
  if(current_log_m <= thres){
    paras.Lambda.log_m = thres;
  }
  if(log_unif<=(thres-current_log_m)){
    return true;
  }else{
    return false;
  }
};

// hybrid MCMC updates with paras and data struct
List hybrid_MALA_updates(X_paras& X0, Lambda_paras& Lambda0, d_paras& d0, double& s2_0, double& t2_0, //initial values
                    double& t2_lambda_0, double& t2_beta_0, double& t2_alpha_0, //initial values
                    hybrid_paras& paras, hybrid_dat& dat, hybrid_comp_par& comp,  //data and parameters
                    int mcmc_sample, double stepsize, int acpt_step, int Gibbs_step, int MALA_step, double target_acpt, bool tune, bool two_step, CharacterVector fixed){ // other MCMC parameters
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
  mat d_list = zeros<mat>(dat.q+dat.r, mcmc_sample/Gibbs_step);
  vec s2_list = zeros<vec>(mcmc_sample/Gibbs_step);
  vec t2_list = zeros<vec>(mcmc_sample/Gibbs_step);
  vec t2_lambda_list = zeros<vec>(mcmc_sample/Gibbs_step);
  vec t2_beta_list = zeros<vec>(mcmc_sample/Gibbs_step);
  vec t2_alpha_list = zeros<vec>(mcmc_sample/Gibbs_step);

  // track MALA Updates
  cube U_list=zeros<cube>(dat.N, dat.q, mcmc_sample/MALA_step);
  cube X_list = zeros<cube>(dat.N, dat.q, mcmc_sample/MALA_step);
  vec llk_list = zeros<vec>(mcmc_sample/MALA_step);
  cube llk_grad_list = zeros<cube>(dat.N, dat.q, mcmc_sample/MALA_step);
  vec dat_llk_list = zeros<vec>(mcmc_sample);
  vec rejection_bound_list = zeros<vec>(mcmc_sample);

  //character vectors
  CharacterVector lambda_vec = {"lambda"};
  CharacterVector d_vec = {"d"};
  CharacterVector beta_vec = {"beta"};
  CharacterVector alpha_vec = {"alpha"};
  CharacterVector s2_vec = {"s2"};
  CharacterVector t2_vec = {"t2"};
  CharacterVector t2_lambda_vec = {"t2_lambda"};
  CharacterVector t2_beta_vec = {"t2_beta"};
  CharacterVector t2_alpha_vec = {"t2_alpha"};
  CharacterVector X_vec = {"X"};
  CharacterVector U_vec = {"U"};

  LogicalVector lambda_fixed = in(lambda_vec,fixed);
  LogicalVector d_fixed = in(d_vec,fixed);
  LogicalVector beta_fixed = in(beta_vec,fixed);
  LogicalVector alpha_fixed = in(alpha_vec,fixed);
  LogicalVector s2_fixed = in(s2_vec,fixed);
  LogicalVector t2_fixed = in(t2_vec,fixed);
  LogicalVector t2_lambda_fixed = in(t2_lambda_vec,fixed);
  LogicalVector t2_beta_fixed = in(t2_beta_vec,fixed);
  LogicalVector t2_alpha_fixed = in(t2_alpha_vec,fixed);
  LogicalVector X_fixed = in(X_vec,fixed);
  LogicalVector U_fixed = in(U_vec,fixed);

  if(mcmc_sample>0){
    // initialize
    paras.Lambda.Lambda_flat = Lambda0.Lambda_flat;
    paras.Lambda.Lambda_flat_t = Lambda0.Lambda_flat_t;
    paras.Lambda.accu_L2 = Lambda0.accu_L2;
    paras.Lambda.log_m = Lambda0.log_m;

    paras.d.d = d0.d;
    paras.d.alpha = d0.alpha;
    paras.d.beta = d0.beta;
    paras.d.Z_alpha = d0.Z_alpha;

    paras.s2 = s2_0;
    paras.t2 = t2_0;
    paras.t2_lambda = t2_lambda_0;
    paras.t2_beta = t2_beta_0;
    paras.t2_alpha = t2_alpha_0;

    paras.X_par.X = X0.X;
    paras.X_par.U = X0.U;
    paras.X_par.Ut = X0.Ut;
    paras.X_par.d_UtYU = X0.d_UtYU;

    paras.Lambda_flat_beta = Lambda0.Lambda_flat*d0.beta;
    paras.s_L_f_d_UtYU = sum_L_f_d_UtYU(Lambda0.Lambda_flat, X0.d_UtYU, dat.M);

    update_pi(paras); //log_pi(X_0)
    update_dpi_2(paras, dat, comp);//dlog_pi(X_0)
    update_dat_llk(paras, dat, comp); //update data llk

    //double current_pi=paras.pi;
    //mat current_dpi=paras.dpi_2;
    double current_pi;
    mat current_dpi;
    mat current_d_UtYU;
    double current_s_L_f_d_UtYU;
    mat current_X;
    mat current_U;
    mat current_Ut;
    mat U_neg; // for sign flipping comparison
    mat U_flipped; // sign flipped U
    umat flip_sign_mat(dat.N, dat.q);
    rowvec current_diff, neg_diff;
    urowvec which_flip;
    urowvec ones_vec(dat.q, fill::value(1));
    mat current_Lambda;

    // collect samples
    Lambda_flat_list.slice(0) = paras.Lambda.Lambda_flat;
    d_list.col(0) = paras.d.d;
    s2_list[0] = paras.s2;
    t2_list[0] = paras.t2;
    t2_lambda_list[0] = paras.t2_lambda;
    t2_beta_list[0] = paras.t2_beta;
    t2_alpha_list[0] = paras.t2_alpha;
    U_list.slice(0) = paras.X_par.U;
    X_list.slice(0) = paras.X_par.X;

    for(int iter=0;iter<mcmc_sample;iter++){
      // update parameter
      if((iter+1)%Gibbs_step==0){
        if(!lambda_fixed[0]){
          if(two_step){
            current_Lambda = paras.Lambda.Lambda_flat;
            if(Lambda_rejection_one_step(paras, dat)){
              rejection_bound_list[iter] = paras.Lambda.log_m;
            }else{
              paras.Lambda.Lambda_flat = current_Lambda;
              paras.Lambda.Lambda_flat_t = current_Lambda.t();
              paras.Lambda.accu_L2 = accu(current_Lambda % current_Lambda);
              rejection_bound_list[iter] = paras.Lambda.log_m;
            }
          }else{
            update_Lambda(comp, paras, dat);
          }
        }
        if(!d_fixed[0]){
          update_d(paras, comp, dat); //d[iter]
        }
        if(!lambda_fixed[0] | !d_fixed[0]){
          paras.Lambda_flat_beta = paras.Lambda.Lambda_flat*paras.d.beta;
        }
        if(!s2_fixed[0]){
          update_s2(paras, comp); //s2[iter]
        }
        if(!t2_fixed[0]){
          update_t2(comp, dat, paras); //t2[iter]
        }
        if(!t2_lambda_fixed[0]){
          update_t2_lambda(paras, comp); //t2_lambda[iter]
        }
        if(!t2_beta_fixed[0]){
          update_t2_beta(paras, comp); //t2_beta[iter]
        }
        if(!t2_alpha_fixed[0]){
          update_t2_alpha(paras, comp); //t2_alpha[iter]
        }
        //update pi and dpi (since s2 and Lambda have been updated)
        //cout << paras.pi << endl;
        update_s_L_f_d_UtYU(paras, dat); //(Lambda)
        update_pi(paras); //log_pi(X_0)
        update_dpi_2(paras, dat, comp);//dlog_pi(X_0)

        Lambda_flat_list.slice(k_Gibbs) = paras.Lambda.Lambda_flat;
        d_list.col(k_Gibbs) = paras.d.d;
        s2_list[k_Gibbs] = paras.s2;
        t2_list[k_Gibbs] = paras.t2;
        t2_lambda_list[k_Gibbs] = paras.t2_lambda;
        t2_beta_list[k_Gibbs] = paras.t2_beta;
        t2_alpha_list[k_Gibbs] = paras.t2_alpha;
        k_Gibbs++;
      }
      //cout << paras.pi << endl;

      current_pi=paras.pi; //pi(X[iter-1], U[iter-1], Lambda[iter], s2[iter])
      current_dpi=paras.dpi_2; //dpi(X[iter-1], U[iter-1], Lambda[iter], s2[iter])
      current_U=paras.X_par.U; //U[iter-1]
      current_Ut=paras.X_par.Ut; //U[iter-1]
      current_X=paras.X_par.X; //X[iter-1]
      current_d_UtYU=paras.X_par.d_UtYU;
      current_s_L_f_d_UtYU=paras.s_L_f_d_UtYU;

      // adjust acceptance rate
      if((iter+1)%(MALA_step*acpt_step)==0 & iter > 0){
        if(!U_fixed[0] & !X_fixed[0]){
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
        if(!U_fixed[0] & !X_fixed[0]){
          //MALA = MALA_one_step(paras,dat,comp);
          //cout << MALA << endl;
          if(MALA_one_step(paras,dat,comp)){ // if accept the new X value
            //cout << "here" << endl;
            // 05/17/22 sign flip
            U_neg = -1.0*paras.X_par.U;
            current_diff = sum(square(current_U - paras.X_par.U),0);
            neg_diff = sum(square(current_U - U_neg),0);
            which_flip = -1*(neg_diff < current_diff);
            which_flip += (neg_diff >= current_diff);

            if(sum(which_flip != ones_vec) > 0){
              std::cout << "flip sign at iter: " << iter+1 << std::endl;
              std::cout << "sign_flips = " << which_flip << std::endl;
              flip_sign_mat.each_row() = which_flip;
              U_flipped = paras.X_par.U % flip_sign_mat;
              paras.X_par.U = U_flipped;
              paras.X_par.Ut = U_flipped.t();
              update_d_UtYU(paras, dat);
              update_s_L_f_d_UtYU(paras, dat);
              update_pi(paras); //pi(X[iter], U_flipped[iter], Lambda[iter], s2[iter])
              update_dpi_2(paras, dat, comp); //dpi(X[iter], U_flipped[iter], Lambda[iter], s2[iter])
            }
            //current_pi=paras.pi;
            //current_dpi=paras.dpi_2;
            acpt_ct++;
            //cout << paras.pi << endl;
          }else{// if not accept
            paras.X_par.U=current_U;
            paras.X_par.Ut=current_Ut;
            paras.X_par.X=current_X;
            paras.pi = current_pi;
            paras.dpi_2 = current_dpi;
            paras.X_par.d_UtYU = current_d_UtYU;
            paras.s_L_f_d_UtYU = current_s_L_f_d_UtYU;
          }
        }
        U_list.slice(k_MALA) = paras.X_par.U;
        X_list.slice(k_MALA) = paras.X_par.X;
        llk_list[k_MALA] = paras.pi;
        llk_grad_list.slice(k_MALA) = paras.dpi_2;
        k_MALA++;
      }

      //update_d_UtYU(paras, dat);
      //update_s_L_f_d_UtYU(paras, dat);
      update_dat_llk(paras, dat, comp); // data llk
      dat_llk_list[iter] = paras.llk;
    }
  }
  return List::create(Named("Lambda_flat") = Lambda_flat_list,
                      Named("d") = d_list,
                      Named("s2") = s2_list,
                      Named("t2") = t2_list,
                      Named("t2_lambda") = t2_lambda_list,
                      Named("t2_beta") = t2_beta_list,
                      Named("t2_alpha") = t2_alpha_list,
                      Named("X")=X_list,
                      Named("U")=U_list,
                      Named("llk")=llk_list,
                      Named("llk_grad")=llk_grad_list,
                      Named("dat_llk")=dat_llk_list,
                      Named("acceptance_rate")=acpt_rate,
                      Named("rejection_bound")=rejection_bound_list);
};

// [[Rcpp::export]]
List hybrid_MALA(mat& X_0, mat& Lambda_0_flat, vec& beta_0, vec& alpha_0, double& s2_0, double& t2_0, //initial values
            double& t2_lambda_0, double& t2_beta_0, double& t2_alpha_0, //initial values
            double& nu0, double& s20, double& rho0, double& psi20, double& eta0, double& t20, //paras
            double& gamma0, double& kappa20, double& omega0, double& phi20, double& log_m0,//paras
            cube& Y, mat& C, mat& Z, int& M, int& N, int& q, int& r, //data
            int mcmc_sample, double stepsize, int acpt_step, int Gibbs_step, int MALA_step, double target_acpt, bool tune, bool two_step, CharacterVector fixed){ // other MCMC parameters
  // build hybrid data structure
  hybrid_dat dat;
  dat.Y = Y;
  dat.C = C;
  dat.Z = Z;
  dat.M = M;
  dat.N = N;
  dat.q = q;
  dat.r = r;
  dat.nu0 = nu0;
  dat.s20 = s20;
  dat.rho0 = rho0;
  dat.psi20 = psi20;
  dat.eta0 = eta0;
  dat.t20 = t20;
  dat.gamma0 = gamma0;
  dat.kappa20 = kappa20;
  dat.omega0 = omega0;
  dat.phi20 = phi20;

  // build hybrid computational paramerters
  hybrid_comp_par comp;
  comp.I_q = eye(q,q);
  comp.I_N = eye(N,N);
  comp.Nq = N*q;
  comp.I_Nq = eye(comp.Nq,comp.Nq);
  comp.dSdX = box_prod_2(comp.I_q, comp.I_N, q, q, N, N);
  comp.q2 = q*q;
  comp.box_idx = box_prod_3(q, N);
  comp.I_q_flat = reshape(vectorise(comp.I_q),1,comp.q2);
  //comp.s2_A = (nu0 + M*N*N)/2.0;
  comp.s2_A = 0.5*nu0 + 0.25*M*N*(N+1);
  comp.s2_B = nu0*s20;
  comp.t2_A = (rho0+M)/2.0;
  comp.t2_B = rho0*psi20;
  comp.W_dim = q+r;
  comp.I_W = eye(comp.W_dim, comp.W_dim);
  comp.t2_lambda_A = (nu0+M*q)/2.0;
  comp.t2_lambda_B = nu0*t20;
  comp.t2_beta_A = (gamma0 + q)/2.0;
  comp.t2_beta_B = gamma0 * kappa20;
  comp.t2_alpha_A = (omega0 + r)/2.0;
  comp.t2_alpha_B = omega0*phi20;
  comp.accu_YY = accu(Y%Y);
  comp.C_mat = repmat(C,1,comp.W_dim);
  comp.Y_tril_ind = trimatl_ind(size(comp.I_N));
  //mat Y_flat_mat(N*(N+1)/2, M);
  mat Y_flat_mat(N*N, M);
  mat Y_i;
  for(int i=0; i<M; i++){
    Y_i = Y.slice(i);
    //Y_flat_mat.col(i) = Y_i(comp.Y_tril_ind);
    Y_flat_mat.col(i) = vectorise(Y_i);
  }
  comp.Y_flat = Y_flat_mat;

  // build hybrid parameters
  d_paras d0;
  vec d0_d(q+r);
  d0.d = d0_d;
  d0.d.subvec(0,q-1) = beta_0;
  d0.d.subvec(q,q+r-1) = alpha_0;
  d0.beta = beta_0;
  d0.alpha = alpha_0;
  d0.Z_alpha = Z*alpha_0;

  Lambda_paras Lambda_par0;
  Lambda_par0.Lambda_flat = Lambda_0_flat;
  Lambda_par0.Lambda_flat_t = Lambda_0_flat.t();
  Lambda_par0.accu_L2 = accu(Lambda_0_flat % Lambda_0_flat);
  Lambda_par0.log_m = log_m0;

  X_paras X0_st;
  X0_st.X=X_0;
  X0_st.U=polar_expansion(X_0);
  X0_st.Ut = X0_st.U.t();
  X0_st.d_UtYU = compute_d_UtYU(X0_st.U, X0_st.Ut, Y, M, q);

  hybrid_paras paras;
  paras.gibbs_s2 = gibbs_s2;
  paras.gibbs_t2 = gibbs_t2;
  paras.gibbs_Lambda_t = gibbs_Lambda_t;
  paras.gibbs_Lambda_t_Y = gibbs_Lambda_t_Y;
  paras.gibbs_d = gibbs_d;
  paras.gibbs_t2_lambda = gibbs_t2_lambda;
  paras.gibbs_t2_beta = gibbs_t2_beta;
  paras.gibbs_t2_alpha = gibbs_t2_alpha;
  paras.dpi_fun_2 = eval_dpi_2;
  paras.pi_func = eval_pi;
  paras.dat_llk_func = eval_dat_llk;
  paras.d_UtYU_func = eval_d_UtYU;
  paras.s_L_f_d_UtYU_func = eval_sum_L_f_d_UtYU;

  //run hybrid
  List hybrid_update_list = hybrid_MALA_updates(X0_st, Lambda_par0, d0, s2_0, t2_0, t2_lambda_0, t2_beta_0, t2_alpha_0, paras, dat, comp, mcmc_sample, stepsize, acpt_step, Gibbs_step, MALA_step,
    target_acpt, tune, two_step, fixed);
  return List::create(Named("mcmc")=hybrid_update_list,
                      Named("stepsize")=paras.w);
};
