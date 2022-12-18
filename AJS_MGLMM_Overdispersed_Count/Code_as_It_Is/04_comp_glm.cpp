// Bivariate Negative Binomial full
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y1);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(nu);
  int n = Y1.size();
  
  parallel_accumulator<Type> nll(this);
  
  vector<Type> mu1 = exp(X*beta1);

  nll -= sum(dcompois2(Y1, mu1, exp(nu(0)), true));

  ADREPORT(exp(nu));

  return nll;
}
