#include <RcppEigen.h>
#include <R.h>
#include <cmath>
#include "rand_functions.h"

// [[Rcpp::depends(RcppEigen)]]

double update_sigma(const int &n, const Eigen::MatrixXd &y, const Eigen::MatrixXd &x, const double &alpha, const Eigen::MatrixXd &beta){
  double c = n/2.0 + 1.0;
  double d = ((y.array()-alpha) - (x*beta).array()).square().sum() / 2;
  return 1/std::sqrt(R::rgamma(c, 1/d));
}

double update_gamma(const int &p, const double &M, const double &lambda, const Eigen::ArrayXd &psi){
  double e = 2 + p*lambda;
  double f = M/(2*lambda) + 0.5*psi.sum();
  return R::rgamma(e, 1/f);
}

Eigen::ArrayXd update_psi(int const &p, const Eigen::MatrixXd &beta, const double &gamma_, const double &lambda){
  double m = lambda - 0.5;
  Eigen::ArrayXd out(p);
  for(int j=0; j<p; j++){
    out(j) = rgig1(m, beta(j, 0)*beta(j, 0), gamma_);
  }
  return out;
}

void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
  R::R_FlushConsole();
#endif
}

Eigen::MatrixXd update_alpha_beta(const Eigen::MatrixXd &X_t_X, 
                                  const Eigen::MatrixXd &X_t_y, 
                                  const double &sigma, 
                                  const Eigen::ArrayXd &psi, 
                                  const Eigen::MatrixXd &I){
  int n = X_t_X.rows();
  Eigen::ArrayXd psi_zero(n);
  psi_zero(0) = 0.0;
  psi_zero.tail(X_t_X.cols()-1) =  1/psi;
  Eigen::MatrixXd diag_ = psi_zero.matrix().asDiagonal();
  double var = sigma*sigma;
  // can use Cholesky decomposition to invert here as (X_t_X + var*diag_)
  // is semi positive definite
  Eigen::MatrixXd temp = (X_t_X + var*diag_).llt().solve(I);
  return rmvnorm(1, temp*X_t_y, var*temp);
}

double update_lambda(const int &p, 
                     const Eigen::ArrayXd &psi, 
                     const double &gamma, 
                     const double &lambda, 
                     const double &M,
                     const double &proposal_sd){
  
  double gamma_ = gamma/2;
  double temp = lambda / gamma_;
  double lambda_prop = lambda * std::exp(proposal_sd * R::rnorm(0.0, 1.0));
  double gamma_prop = lambda_prop / temp;
  
  double logaccept = std::log(lambda_prop) - std::log(lambda) - 0.5*(lambda_prop - lambda)
                     + p * lambda_prop * std::log(gamma_prop) - p * std::lgamma(lambda_prop)
                     - p * lambda * std::log(gamma_) + p * std::lgamma(lambda)
                     + (lambda_prop - lambda) * psi.log().sum()
                     + (gamma_ - gamma_prop)*psi.sum();
  
  // we are using a log-normal proposal so need to adjust probability (change of variables)
  logaccept = logaccept - std::log(lambda) + std::log(lambda_prop);
  
  if(logaccept > std::log(R::runif(0.0, 1.0))){
    return lambda_prop;
  } else {
    return lambda;
  }
  
}

// [[Rcpp::export]]
Eigen::MatrixXd get_chains(const Eigen::MatrixXd &y, 
                           const Eigen::MatrixXd &x, 
                           const Eigen::MatrixXd &X_t_X,
                           const Eigen::MatrixXd &X_t_y,
                           const int &n_samp,
                           const int &n_thin,
                           const int &n,
                           const int &p,
                           const double &M,
                           const double &init_alpha, 
                           const Eigen::MatrixXd &init_beta,
                           const Eigen::ArrayXd &init_psi, 
                           const double &init_sigma, 
                           const double &init_lambda,
                           const double &init_gamma, 
                           const double &init_prop_sd,
                           const int &tuning_time, 
                           const int &tuning_freq, 
                           const int &verbose){
  
  double proposal_sd = init_prop_sd;
  Eigen::MatrixXd chains((n_samp+1), 2*p + 4);
  
  double curr_alpha = init_alpha;
  Eigen::MatrixXd curr_beta = init_beta;
  Eigen::ArrayXd curr_psi = init_psi;
  double curr_sigma = init_sigma;
  double curr_lambda = init_lambda;
  double curr_gamma = init_gamma;
  Eigen::ArrayXd new_alpha_beta;
  
  // some indices and counters
  int sigma_ind = 2*p + 1;
  int lambda_ind = 2*p + 2;
  int gamma_ind = 2*p + 3;
  int verb_freq_count = 100;
  int tuning_count = 0;
  int temp = 0;
  int interupt_count = 0;
  
  // store initial chain states
  chains(0, 0) = curr_alpha;
  chains.block(0, 1, 1, p) = curr_beta.transpose();
  chains.block(0, p+1, 1, p ) = curr_psi.transpose();
  chains(0, sigma_ind) = curr_sigma;
  chains(0, lambda_ind) = curr_lambda;
  chains(0, gamma_ind) = curr_gamma;
  
  // get random number state
  GetRNGstate();
  
  Eigen::ArrayXd log_lambdas(tuning_time);
  Eigen::ArrayXd K(tuning_freq);
  Eigen::ArrayXd K_hat(tuning_freq);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(X_t_X.rows(), X_t_X.rows());
  
  double prev_lambda = curr_lambda;
  
  for(int i=1; i < (n_samp+1); i++){
    
    for(int t=0; t < n_thin; t++){
      
      temp = i*n_thin;
      
      new_alpha_beta = update_alpha_beta(X_t_X, X_t_y, curr_sigma, curr_psi, I).array();
      curr_alpha = new_alpha_beta(0);
      curr_beta = new_alpha_beta.tail(p).matrix();
      curr_psi = update_psi(p, curr_beta, curr_gamma, curr_lambda);
      curr_sigma = update_sigma(n, y, x, curr_alpha, curr_beta);
      curr_gamma = update_gamma(p, M, curr_lambda, curr_psi);
      prev_lambda = curr_lambda;
      curr_lambda = update_lambda(p, curr_psi, curr_gamma, curr_lambda, M, proposal_sd);
      
      if((temp+t) <= tuning_time){
        // store previous values to be used in adaptive metropolis step
        log_lambdas(temp+t) = std::log(curr_lambda);
        tuning_count += 1;
      }
      
      if((tuning_count == tuning_freq) & ((temp+t) <= tuning_time)){
        // Metropolis proposal variance tuning step (Haario)
        K = log_lambdas.segment(temp+t-tuning_freq, tuning_freq);
        K_hat = K - K.mean();
        proposal_sd = 2.4*std::sqrt((K_hat*K_hat).sum()/(tuning_freq-1));
        tuning_count = 0;
      }
      
    }
    
    // store updated parameters
    chains(i, 0) = curr_alpha;
    chains.block(i, 1, 1, p) = curr_beta.transpose();
    chains.block(i, p+1, 1, p ) = curr_psi.transpose();
    chains(i, sigma_ind) = curr_sigma;
    chains(i, lambda_ind) = curr_lambda;
    chains(i, gamma_ind) = curr_gamma;
    
    if((verbose == 1) & (verb_freq_count == 100)){
      // show percentage completed
      Rcpp::Rcout << std::round(100 * double(i+1)/n_samp) << " %\r";
      flush_console();
      verb_freq_count=0;
    }
    
    if(interupt_count == 500){
      // check if user wants to terminate after every 500 iterations
      Rcpp::checkUserInterrupt();
      interupt_count = 0;
    }
    interupt_count += 1;
    verb_freq_count += 1;
    
  }
  
  // return random number state
  PutRNGstate();
  
  flush_console();
  return chains;
}




