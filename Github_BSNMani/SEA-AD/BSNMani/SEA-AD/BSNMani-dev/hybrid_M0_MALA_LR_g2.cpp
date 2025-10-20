// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]

#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double t2_new(const double& t2_A, const double& t2_B, const mat& C, const mat& Z, const mat& Lambda_flat_beta, const mat& Z_alpha){
  mat C_res = C-Lambda_flat_beta-Z_alpha;
  //cout << accu(C_res%C_res) << endl;
  double new_t2 = 1.0/randg<double>(distr_param(t2_A, 2.0/(t2_B+accu(C_res%C_res))));
  return new_t2;
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
double C_llk_new(const mat& C, const double& t2, const mat& Lambda_flat_beta, const mat& Z_alpha, const int& M){
  mat C_mean = Lambda_flat_beta + Z_alpha;
  double dat_llk;
  double t = sqrt(t2);
  mat C_smat(M, 1, fill::value(t));
  dat_llk = accu(log_normpdf(C, C_mean, C_smat));

  return dat_llk;
}

// data structure
struct g2_dat{
  // data; hyperparameters
  mat C; //M x 1
  mat Z; //M x r
  int M; // population size
  int N; // # roi
  int q; // # reduced dimensions
  int r; // # covariates
  cube Lambda_flat_cube; // M x q x mcmc_samples
  double rho0; //t2
  double psi20; //t2
  double gamma0; //t2_beta
  double kappa20; //t2_beta
  double omega0; //t2_alpha
  double phi20; //t2_alpha
  int mcmc_samples;
};

struct g2_comp_par{
  double t2_A;
  double t2_B;
  mat I_W;
  int W_dim;
  double t2_beta_A;
  double t2_beta_B;
  double t2_alpha_A;
  double t2_alpha_B;
  mat C_mat; //repmat(C,1,W_dim)
};

double gibbs_t2(const g2_comp_par& comp_par, const g2_dat& dat, const vec& beta, const mat& Lambda_flat_beta, const mat& Z_alpha){
  return t2_new(comp_par.t2_A, comp_par.t2_B, dat.C, dat.Z, Lambda_flat_beta, Z_alpha);
}

vec gibbs_d(const g2_dat& dat, const g2_comp_par& comp_par, const double& t2_beta, const double& t2_alpha, const double& t2, const mat& Lambda_flat){
  return d_new(dat.r, dat.q, comp_par.W_dim, comp_par.I_W, dat.Z, comp_par.C_mat, t2_beta, t2_alpha, t2, Lambda_flat);
}

double gibbs_t2_beta(const g2_comp_par& comp_par, const vec& beta){
  return t2_beta_new(comp_par.t2_beta_A, comp_par.t2_beta_B, beta);
}

double gibbs_t2_alpha(const g2_comp_par& comp_par, const vec& alpha){
  return t2_alpha_new(comp_par.t2_alpha_A, comp_par.t2_alpha_B, alpha);
}

double eval_dat_llk(const double& t2, const mat& Lambda_flat_beta, const mat& Z_alpha, const g2_dat& dat, const g2_comp_par& comp_par){
  return C_llk_new(dat.C, t2, Lambda_flat_beta, Z_alpha, dat.M);
}

typedef double(*funcPtr_t2)(const g2_comp_par& comp_par, const g2_dat& dat, const vec& beta, const mat& Lambda_flat_beta, const mat& Z_alpha);
typedef vec(*funcPtr_d)(const g2_dat& dat, const g2_comp_par& comp_par, const double& t2_beta, const double& t2_alpha, const double& t2, const mat& Lambda_flat);
typedef double(*funcPtr_t2_beta)(const g2_comp_par& comp_par, const vec& beta);
typedef double(*funcPtr_t2_alpha)(const g2_comp_par& comp_par, const vec& alpha);
typedef double(*funcPtr_dat_llk)(const double& t2, const mat& Lambda_flat_beta, const mat& Z_alpha, const g2_dat& dat, const g2_comp_par& comp_par);

struct d_paras{
  vec alpha; // r
  vec beta; // q
  vec d; // r+q
  mat Z_alpha; // M x 1
  mat Lambda_flat_beta; // M x 1
};

struct g2_paras{
  double t2;
  d_paras d;
  mat Lambda; //M x q
  double t2_beta;
  double t2_alpha;
  double llk; // value of data llk
  funcPtr_t2 gibbs_t2;
  funcPtr_d gibbs_d;
  funcPtr_t2_beta gibbs_t2_beta;
  funcPtr_t2_alpha gibbs_t2_alpha;
  funcPtr_dat_llk dat_llk_func; // pointer to data llk function
};

void update_t2(g2_comp_par& comp_par, g2_dat& dat, g2_paras& paras){
  paras.t2 = paras.gibbs_t2(comp_par, dat, paras.d.beta, paras.d.Lambda_flat_beta, paras.d.Z_alpha);
};

void update_d(g2_paras& paras, g2_comp_par& comp_par, g2_dat& dat){
  paras.d.d = paras.gibbs_d(dat, comp_par, paras.t2_beta, paras.t2_alpha, paras.t2, paras.Lambda);
  paras.d.beta = paras.d.d.subvec(0,dat.q-1);
  paras.d.alpha = paras.d.d.subvec(dat.q, comp_par.W_dim-1);
  paras.d.Z_alpha = dat.Z*paras.d.alpha;
  paras.d.Lambda_flat_beta = paras.Lambda*paras.d.beta;
};

void update_t2_beta(g2_paras& paras, g2_comp_par& comp_par){
  paras.t2_beta = paras.gibbs_t2_beta(comp_par, paras.d.beta);
};

void update_t2_alpha(g2_paras& paras, g2_comp_par& comp_par){
  paras.t2_alpha = paras.gibbs_t2_alpha(comp_par, paras.d.alpha);
};

void update_dat_llk(g2_paras& paras, g2_dat& dat, g2_comp_par& comp_par){
  paras.llk = paras.dat_llk_func(paras.t2, paras.d.Lambda_flat_beta, paras.d.Z_alpha, dat, comp_par);
};

void update_Lambda(g2_paras& paras, g2_dat& dat, int& iter){
  paras.Lambda = dat.Lambda_flat_cube.slice(iter);
}

// hybrid MCMC updates with paras and data struct
List g2_updates(d_paras& d0, double& t2_0, double& t2_beta_0, double& t2_alpha_0, //initial values
                    g2_paras& paras, g2_dat& dat, g2_comp_par& comp,  //data and parameters
                    int& mcmc_sample, CharacterVector fixed){ // other MCMC parameters
  // Tracking Gibbs Updates
  mat d_list = zeros<mat>(dat.q+dat.r, mcmc_sample);
  vec t2_list = zeros<vec>(mcmc_sample);
  vec t2_beta_list = zeros<vec>(mcmc_sample);
  vec t2_alpha_list = zeros<vec>(mcmc_sample);
  cube Lambda_flat_list=zeros<cube>(dat.M, dat.q, mcmc_sample);

  vec dat_llk_list = zeros<vec>(mcmc_sample);

  //character vectors
  CharacterVector d_vec = {"d"};
  CharacterVector beta_vec = {"beta"};
  CharacterVector alpha_vec = {"alpha"};
  CharacterVector t2_vec = {"t2"};
  CharacterVector t2_beta_vec = {"t2_beta"};
  CharacterVector t2_alpha_vec = {"t2_alpha"};

  LogicalVector d_fixed = in(d_vec,fixed);
  LogicalVector beta_fixed = in(beta_vec,fixed);
  LogicalVector alpha_fixed = in(alpha_vec,fixed);
  LogicalVector t2_fixed = in(t2_vec,fixed);
  LogicalVector t2_beta_fixed = in(t2_beta_vec,fixed);
  LogicalVector t2_alpha_fixed = in(t2_alpha_vec,fixed);

  if(mcmc_sample>0){
    // initialize
    paras.d.d = d0.d;
    paras.d.alpha = d0.alpha;
    paras.d.beta = d0.beta;
    paras.d.Z_alpha = d0.Z_alpha;
    paras.d.Lambda_flat_beta = d0.Lambda_flat_beta;

    paras.t2 = t2_0;
    paras.t2_beta = t2_beta_0;
    paras.t2_alpha = t2_alpha_0;

    update_dat_llk(paras, dat, comp); //update data llk

    // collect samples
    d_list.col(0) = paras.d.d;
    t2_list[0] = paras.t2;
    t2_beta_list[0] = paras.t2_beta;
    t2_alpha_list[0] = paras.t2_alpha;
    Lambda_flat_list.slice(0) = paras.Lambda;

    for(int iter=0;iter<mcmc_sample;iter++){
      // update parameter
      update_Lambda(paras, dat, iter);
      if(!d_fixed[0]){
        update_d(paras, comp, dat); //d[iter]
      }
      if(!t2_fixed[0]){
        update_t2(comp, dat, paras); //t2[iter]
      }
      if(!t2_beta_fixed[0]){
        update_t2_beta(paras, comp); //t2_beta[iter]
      }
      if(!t2_alpha_fixed[0]){
        update_t2_alpha(paras, comp); //t2_alpha[iter]
      }

      d_list.col(iter) = paras.d.d;
      t2_list[iter] = paras.t2;
      t2_beta_list[iter] = paras.t2_beta;
      t2_alpha_list[iter] = paras.t2_alpha;
      Lambda_flat_list.slice(iter) = paras.Lambda;

      update_dat_llk(paras, dat, comp); // data llk
      dat_llk_list[iter] = paras.llk;
    }
  }
  return List::create(Named("d") = d_list,
                      Named("t2") = t2_list,
                      Named("t2_beta") = t2_beta_list,
                      Named("t2_alpha") = t2_alpha_list,
                      Named("Lambda_flat") = Lambda_flat_list,
                      Named("dat_llk")=dat_llk_list);
};

// [[Rcpp::export]]
List g2(vec& beta_0, vec& alpha_0, double& t2_0, double& t2_beta_0, double& t2_alpha_0, mat& Lambda_flat0, //initial values
            double& rho0, double& psi20, double& gamma0, double& kappa20, double& omega0, double& phi20, //paras
            cube& Lambda_flat_cube, mat& C, mat& Z, int& M, int& N, int& q, int& r, //data
            int& mcmc_sample, CharacterVector fixed){ // other MCMC parameters
  // build hybrid data structure
  g2_dat dat;
  dat.C = C;
  dat.Z = Z;
  dat.M = M;
  dat.N = N;
  dat.q = q;
  dat.r = r;
  dat.Lambda_flat_cube = Lambda_flat_cube;
  std::cout<< "break 1" << std::endl;
  dat.rho0 = rho0;
  std::cout<< "break 2" << std::endl;
  dat.psi20 = psi20;
  dat.gamma0 = gamma0;
  dat.kappa20 = kappa20;
  dat.omega0 = omega0;
  dat.phi20 = phi20;

  // build hybrid computational paramerters
  g2_comp_par comp;
  comp.t2_A = (rho0+M)/2.0;
  comp.t2_B = rho0*psi20;
  comp.W_dim = q+r;
  comp.I_W = eye(comp.W_dim, comp.W_dim);
  comp.t2_beta_A = (gamma0 + q)/2.0;
  comp.t2_beta_B = gamma0 * kappa20;
  comp.t2_alpha_A = (omega0 + r)/2.0;
  comp.t2_alpha_B = omega0*phi20;
  comp.C_mat = repmat(C,1,comp.W_dim);

  // build hybrid parameters
  d_paras d0;
  vec d0_d(q+r);
  d0.d = d0_d;
  d0.d.subvec(0,q-1) = beta_0;
  std::cout<< "beta" << std::endl;
  d0.d.subvec(q,q+r-1) = alpha_0;
  std::cout<< "alpha" << std::endl;
  d0.beta = beta_0;
  d0.alpha = alpha_0;
  d0.Z_alpha = Z*alpha_0;
  d0.Lambda_flat_beta = Lambda_flat0*beta_0;

  g2_paras paras;
  paras.gibbs_t2 = gibbs_t2;
  paras.gibbs_d = gibbs_d;
  paras.gibbs_t2_beta = gibbs_t2_beta;
  paras.gibbs_t2_alpha = gibbs_t2_alpha;
  paras.dat_llk_func = eval_dat_llk;
  paras.Lambda = Lambda_flat0;

  //run hybrid
  List hybrid_update_list = g2_updates(d0, t2_0, t2_beta_0, t2_alpha_0, paras, dat, comp, mcmc_sample, fixed);
  return List::create(Named("mcmc")=hybrid_update_list);
};
