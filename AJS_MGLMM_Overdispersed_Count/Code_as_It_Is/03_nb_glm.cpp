// Bivariate Negative Binomial full
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y1);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(phi);
  int n = Y1.size();
  vector<Type> mu1 = exp(X*beta1);
  vector<Type> var1 = mu1 + (mu1*mu1)/exp(phi(0));

  Type nll1 = -sum(dnbinom2(Y1, mu1, var1, true));

  ADREPORT(exp(phi));

  return nll1;
}
