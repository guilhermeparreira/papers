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
  
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);
  
  int n = Y.rows(); // Number of rows
  int c = Y.cols();  // Number of cols
  matrix<Type> Xbeta(n, c);
  Xbeta = X*beta;
  
  vector<Type> Yj(n);
  vector<Type> mu(n);
  for (int j = 0; j<c; j++){ //Density for response
    Yj = Y.col(j);
    mu = exp(Xbeta.col(j).array() + U.col(j).array());
    nll -= sum(dpois(Yj, mu, true));
  }

  for(int i = 0; i < n; i++) // Density for R.E.
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));

  matrix<Type> Cor(c,c);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);

  return nll;
}
