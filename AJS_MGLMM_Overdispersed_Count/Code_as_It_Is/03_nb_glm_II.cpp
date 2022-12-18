// Bivariate Negative Binomial full
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y);
  PARAMETER(phi);
  PARAMETER(mu);
  Type var = mu + (mu*mu)/phi;

  Type nll1 = sum(dnbinom2(Y, mu, var, false));

  return nll1;
}
