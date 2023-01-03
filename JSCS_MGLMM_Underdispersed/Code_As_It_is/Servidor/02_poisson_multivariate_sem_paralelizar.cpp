// Poisson Multivariate (n responses) same fixed Effects.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(X);          // n x p
  PARAMETER_MATRIX(beta);  // p x r
  DATA_MATRIX(Y);          // n x r
  PARAMETER_MATRIX(U);     // n x r
  PARAMETER_VECTOR(rho);   // r(r-1)/2
  PARAMETER_VECTOR(sigma); // r
  
  Type nll = 0;
  // parallel_accumulator<Type> nll(this);
  
  int n = Y.rows(); // Number of rows
  int c = Y.cols();  // Number of cols
  matrix<Type> Xbeta(n, c);
  Xbeta = X*beta;
  
  vector<Type> Yj(n);
  vector<Type> mu(n);
  for (int j = 0; j<c; j++){
    Yj = Y.col(j);
    mu = exp(Xbeta.col(j).array() + U.col(j).array());
    nll -= sum(dpois(Yj, mu, true)); // Talvez cada response tenha que ter uma var_name diff
  }
  
  // Sigma as standard deviation only
  int k = sigma.size();
  for (int j = 0; j<k; j++){
    sigma(j) = exp(sigma(j));
  }
  // Rho Rep
  int rr = rho.size();
  for (int l = 0; l<rr; l++){
    rho(l) = -1 + 2*(exp(rho(l))/(exp(rho(l)) + 1));
  }
  // Type cov_rho = rho*sqrt(sigma(0))*sqrt(sigma(1));
  // vector<Type> rho_temp(1);
  // rho_temp(0) = cov_rho;
  
  
  // NORMAL MULTIVARIATE
  matrix<Type> Sigma(c,c); // Sigma has to be a covariance matrix
  int counter = 0;
  for(int i = 0; i < c; i++){
    Sigma(i,i) = pow(sigma(i),2);
    for(int j = 0; j < i; j++){
        Sigma(i,j) = rho(counter)*sigma(i)*sigma(j);
        Sigma(j,i) = Sigma(i,j);
        counter += 1;
      }
    }
  // Sigma(0,0) = pow(exp(sigma(0)),2);
  // Sigma(1,1) = pow(exp(sigma(1)),2);
  // Sigma(2,2) = pow(exp(sigma(2)),2);
  // Sigma(3,3) = pow(exp(sigma(3)),2);
  // Sigma(4,4) = pow(exp(sigma(4)),2);
  
  // Sigma(1,0) = rho_temp(0)*sqrt(exp(sigma(0)))*sqrt(exp(sigma(1)));
  // Sigma(0,1) = rho_temp(0)*sqrt(exp(sigma(0)))*sqrt(exp(sigma(1)));
  
  for (int i = 0; i < n; i++){
    nll += MVNORM(Sigma)(U.row(i));
  }
  
  // ADREPORT((exp(2*rho) - 1)/(exp(2*rho) + 1)); //Correlation
  // ADREPORT(exp(sigma)); // Standard deviation
  return nll;
}
