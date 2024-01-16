\encoding{UTF-8}
\name{NG}
\alias{NGlm-package}
\alias{NGlm}
\alias{NG}
\title{
Sparse Bayesian Linear Models with Normal-Gamma Priors
}
\description{
  This package exposes a single method, \code{NG()}, which performs sparse Bayesian linear regression using the Normal-Gamma priors described by Griffin & Brown (2010). This is particularly suited to problems where the number of covariates is greater than the number of observations and the number of covariates with non-zero coefficients is small. The MCMC sampler is implemented in C++ for performance, utilizing the template library 'Eigen' for fast linear algebra operations (interfaced to R using Rcpp/RcppEigen). In addition, OpenMP is enabled to allow Eigen to exploit multiple cores (on machines that have them) when performing dense matrix multiplication.
}

\usage{
NG(y, X, n_samples=10000, n_thin=1, init_alpha=0.0, init_beta=NULL, 
   init_psi=NULL, init_sigma=1.0, init_lambda=1.0, init_gamma=1.0,    
   init_prop_sd=1, tuning_time=2000, tuning_freq=200, verbose=TRUE)
}

\arguments{
  \item{y}{a numeric vector representing the response.}
  \item{X}{a numeric data.frame or matrix, representing the matrix of covariates. Should have the same number of rows as y. This should not contain an intercept column of ones as this will be added internally.}
  \item{n_samples}{the number of MCMC samples to be collected.}
  \item{n_thin}{the number of MCMC samples to discard before collecting the next sample. For example, if n_thin=4, then every 4th sample will be stored. By default, n_thin=1, corresponding to no thinning (all samples are collected).}
  \item{init_alpha}{the inital state of the markov chain for the intercept.}
  \item{init_beta}{the inital states of the markov chains for the regression coefficients. If supplied, it should be a numeric vector with length equal to the number of columns of X. If not supplied, then each chain is intialized at 0.}
  \item{init_psi}{the inital states of the markov chains for the latent psi variables. If supplied, it should be a numeric vector with length equal to the number of columns of X. If not supplied, then each chain is intialized to 1.}
  \item{init_sigma}{the inital state of the markov chain for the errors standard deviation. If not supplied, is set to 1.}
  \item{init_lambda}{the inital state of the markov chain for the lambda parameter. If not supplied, is set to 1.}
  \item{init_gamma}{the inital state of the markov chain for the gamma parameter. If not supplied, is set to 1.}
  \item{init_prop_sd}{the inital standard deviation for the Gaussian proposal distribution used in the Metropolis sampling step for lambda.}
  \item{tuning_time}{the number of initial iterations to be used for tuning the proposal standard deviation for the Metropolis step.}
  \item{tuning_freq}{the frequency with which to update the proposal standard deviation during the tuning_time.}
  \item{verbose}{if verbose=TRUE (default) then the percentage of MCMC samples collected will be displayed as sampling takes place, along with the time taken at the end. If FALSE, then nothing will be displayed.}
}

\value{
  \code{NG} returns an object of \link{class} "NG", which is essentially a holder containing the MCMC samples of the model parameters. 

An object of class "NG" contains the following components:
  \item{y}{the response used.}
  \item{X}{the model matrix used.}
  \item{call}{the matched call.}
  \item{n_samples}{the number of samples obtained.}
  \item{samples}{A matrix holding the samples for each chain as columns.}
  \item{n_thin}{the level of thinning.}
  \item{tuning_time}{the number of iterations spent tuning the proposal variance for lambda.}
  \item{tuning_freq}{the tuning frequency used.}
}

\details{
  A summary of the model specification and sampling procedure is given below, however, full details can be found in Griffin & Brown (2010). The response vector, \eqn{\boldsymbol{Y}}, is modeled as: 
  
  \deqn{%
  \boldsymbol{Y} = \alpha + \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon} \\
  }{}

  Where \eqn{\boldsymbol{X}} is an \eqn{n \times p} matrix of covariates, \eqn{\boldsymbol{\beta} = (\beta_1, ..., \beta_p)^T} is a vector of coefficients, \eqn{\boldsymbol{\epsilon} = (\epsilon_1, ..., \epsilon_n)^T}, with \eqn{\epsilon_i \overset{iid}{\sim} N(0, \sigma^2)} and \eqn{\alpha} is an intercept parameter. The hierarchical priors for these parameters are:
  
  \deqn{%
  \beta_j | \Psi_j \sim N(0, \Psi_j) \\
  \Psi_j | \lambda, \gamma \sim Ga(\lambda, \frac{1}{2\gamma^2}) \\
  \lambda \sim Exp(1) \\
  \pi(\alpha) \propto 1 \\
  \pi(\sigma^{-2}) \propto 1
  }{}
  The prior for \eqn{\gamma} given \eqn{\lambda} is given by
  \deqn{%
  \nu_{\beta} = \lambda \gamma^2 \sim IG(2, M)
  }{}
  where \eqn{M = \frac{1}{p}\sum_{j=1}^{p} \hat{\beta_j ^ 2}}. If X is non-singular (when p < n), then \eqn{\hat{\beta}} is the least squares estimate, otherwise it is the Minimum Length Least Squares estimate.
  
  The MCMC sampler implements the Gibbs sampler as described in Griffin & Brown (2010), differing only in the adaptive Metropolis step for lambda. The adaptive step here uses the approach of Haario, Saksman & Tamminen (1999), where the proposal variance at the current iteration is based on the covariance of the previous \eqn{H} sampled values. This \eqn{H}, is denoted by \code{tuning_freq} here.
  
  The Gibbs sampling step for \eqn{\Psi_j} requires sampling from the Generalized Inverse Gaussian distribution. For this, the algorithm developed by Hormann & Leydold (2014) is used.
  
  The function \code{summary} can be used to obtain a summary table containing posterior summaries, see \link{summary.NG} on how to use this.
}

\author{
Carson McKee \email{carson.mckee@kcl.ac.uk}
}

\keyword{Bayesian, Regression, Shrinkage, Sparse, Normal-Gamma}

\examples{
# generate dummy data
set.seed(1234)
p <- 100
n <- 75
X <- matrix(rnorm(p*n, 0, 1), nrow=n, ncol=p)
beta <- rep(0, p)
beta[1:5] <- 5
y <- X\%*\%beta + rnorm(n, 0, 1)

# fit the model
mod <- NG(y, X, verbose=FALSE, n_samples=10000, n_thin=1)

# inspect the beta 1 chain
plot(mod$samples[2000:10000, 'beta_1'], type='l')

# geneate posterior summary estimates
summary(mod, burnin=2000)
}

\references{
Bates, D., & Eddelbuettel, D. (2013). Fast and Elegant Numerical Linear Algebra Using the RcppEigen Package. Journal of Statistical Software, 52(5), p.1–24. https://doi.org/10.18637/jss.v052.i05

Eddelbuettel, D., & Francois, R. (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), p.1–18. https://doi.org/10.18637/jss.v040.i08

Griffin, J. E., & Brown, P. J. (2010). Inference with normal-gamma prior distributions in regression problems. Bayesian Analysis, 5(1), p.171–188. https://doi.org/10.1214/10-BA507

Guennebaud, G., Jacob, B., & others. (2010). Eigen v3. http://eigen.tuxfamily.org. http://eigen.tuxfamily.org

Haario, H., Saksman, E., and Tamminen, J. 1999. Adaptive proposal distribution for random walk Metropolis algorithm. Computational Statistics, 14(3), p.375-395.

Hormann, W., and Leydold, J. 2014. Generating generalized inverse Gaussian random variates. Statistics and Computing, 24(4), p.547-557.

}

\seealso{
\code{\link{summary.NG}}, \code{\link{lm}}, \code{\link{glm}}
}
