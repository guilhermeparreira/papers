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
  PARAMETER_VECTOR(rho);
  PARAMETER_VECTOR(sigma);
  PARAMETER_VECTOR(nu);
  int n = Y1.size();
  vector<Type> mu1 = exp(X*beta1 + U.col(0).array());
  vector<Type> mu2 = exp(X*beta2 + U.col(1).array());
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);
  
  nll -= sum(dcompois2(Y1, mu1, exp(nu(0)), true));
  nll -= sum(dcompois2(Y2, mu2, exp(nu(1)), true));
  
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));

  matrix<Type> Cor(2,2);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(nu));
  return nll;
}
