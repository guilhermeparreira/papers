// COM-Poisson full
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y1);
  DATA_VECTOR(Y2);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_MATRIX(U);
  PARAMETER(rho);
  PARAMETER_VECTOR(sigma);
  PARAMETER_VECTOR(nu);
  int n = Y1.size();
  int u = U.cols();
  Type nll = 0;
  // Response Variables
  vector<Type> mu1 = X*beta1 + U.col(0).array();
  for(int i = 0; i < n; i++)
    nll -= dcompois2(Y1(i), exp(mu1(i)), exp(nu(0)), true);
  vector<Type> mu2 = X*beta2 + U.col(1).array();
  for(int i = 0; i < n; i++)
    nll -= dcompois2(Y2(i), exp(mu2(i)), exp(nu(1)), true);
  // Reparametrization of input parameters
  sigma(0) = pow(exp(sigma(0)),2);
  sigma(1) = pow(exp(sigma(1)),2);
  rho = -1 + 2*(exp(rho)/(exp(rho) + 1)); // // transf.rho <- function(rho){-1+2*(exp(rho)/(exp(rho)+1))}
  // Random effect
  Type cov_rho = rho*sqrt(sigma(0))*sqrt(sigma(1));
  matrix<Type> Sigma(u,u);
  Sigma.row(0) << sigma(0), cov_rho;
  Sigma.row(1) << cov_rho, sigma(1);
  MVNORM_t<Type> neg_log_density(Sigma); // Create object from covariance matrix Sigma; (neg_log_density could be any name)
  for (int i=0; i<n; i++){
    nll += neg_log_density(U.row(i)); // Process likelihood
    SIMULATE {
      U.row(i) = neg_log_density.simulate();
    }
  }
  SIMULATE {
  vector<Type> mu1 = X*beta1 + U.col(0).array();
  Y1 = rcompois2(exp(mu1), exp(nu(0)));  // Simulate response
  REPORT(Y1);          // Report the simulation
  
  vector<Type> mu2 = X*beta2 + U.col(1).array();
  Y2 = rcompois2(exp(mu2), exp(nu(1)));  // Simulate response
  REPORT(Y2);          // Report the simulation
  }

  ADREPORT(exp(nu(0)));
  ADREPORT(exp(nu(1)));
  ADREPORT(exp(sigma(0))); // Standard Deviation
  ADREPORT(exp(sigma(1))); // Standard Deviation
  ADREPORT(-1 + 2*(exp(rho)/(exp(rho) + 1)));

  // for(int i = 0; i < n; i++){
  //   nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i)); // Likelihood
  // }
  // 
  // matrix<Type> Cor(2,2);
  // Cor = UNSTRUCTURED_CORR(rho).cov();
  // ADREPORT(Cor);
  
  
  return nll;
}
