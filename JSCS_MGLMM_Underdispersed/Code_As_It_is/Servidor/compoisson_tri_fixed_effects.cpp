// COM-Poisson multivariate
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(X1);             // n x p
  DATA_MATRIX(X2);             // n x p
  DATA_MATRIX(X3);             // n x p
  DATA_MATRIX(Y);             // n x r
  PARAMETER_VECTOR(beta1);     // p x 1
  PARAMETER_VECTOR(beta2);     // p x 2
  PARAMETER_VECTOR(beta3);     // p x 3
  PARAMETER_MATRIX(U);        // n x r
  PARAMETER_VECTOR(rho);      // r(r-1)/2
  PARAMETER_VECTOR(sigma);    //r -> R.E.
  PARAMETER_VECTOR(nu);      //r -> F.E.
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);
  
  int n = Y.rows(); // Number of rows
  int c = Y.cols();  // Number of cols
  matrix<Type> Xbeta(n, c);
  vector<Type> Xbeta1 = X1*beta1;
  vector<Type> Xbeta2 = X2*beta2;
  vector<Type> Xbeta3 = X3*beta3;
  
  Xbeta.col(0) = Xbeta1;
  Xbeta.col(1) = Xbeta2;
  Xbeta.col(2) = Xbeta3;
  
  vector<Type> Yj(n);
  vector<Type> mu(n);
  for (int j = 0; j<c; j++){
    Yj = Y.col(j);
    mu = exp(Xbeta.col(j).array() + U.col(j).array());
    nll -= sum(dcompois2(Yj, mu, exp(nu(j)), true));
  }
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));
  matrix<Type> Cor(c,c);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(nu));
  return nll;
}
