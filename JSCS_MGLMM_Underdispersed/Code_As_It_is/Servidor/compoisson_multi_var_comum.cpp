// COM-Poisson multivariate
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(X);             // n x p
  DATA_MATRIX(Y);             // n x r
  PARAMETER_MATRIX(beta);     // p x r
  PARAMETER_MATRIX(U);        // n x r
  PARAMETER_VECTOR(rho);      // r(r-1)/2
  PARAMETER_VECTOR(nu);      //r -> F.E.
  PARAMETER(sigma);    //r -> R.E.
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);
  
  int n = Y.rows(); // Number of rows
  int c = Y.cols();  // Number of cols
  matrix<Type> Xbeta(n, c);
  Xbeta = X*beta;
  
  // It makes all standard deviation the same
  vector<Type> sigma_aux(c); 
  for (int k = 0; k<c; k++){
    sigma_aux(k) = sigma;
  }
  
  vector<Type> Yj(n);
  vector<Type> mu(n);
  for (int j = 0; j<c; j++){
    Yj = Y.col(j);
    mu = exp(Xbeta.col(j).array() + U.col(j).array());
    nll -= sum(dcompois2(Yj, mu, exp(nu(j)), true));
  }
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma_aux)(U.row(i));
  matrix<Type> Cor(c,c);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(nu));
  return nll;
}
