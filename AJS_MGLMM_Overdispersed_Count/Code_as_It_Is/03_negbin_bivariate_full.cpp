// Bivariate Negative Binomial full
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
  PARAMETER_VECTOR(phi);
  int n = Y1.size();
  vector<Type> mu1 = exp(X*beta1 + U.col(0).array());
  vector<Type> mu2 = exp(X*beta2 + U.col(1).array());
  vector<Type> var1 = mu1 + (mu1*mu1)/exp(phi(0));
  vector<Type> var2 = mu2 + (mu2*mu2)/exp(phi(1));
  
  Type nll1 = -sum(dnbinom2(Y1, mu1, var1, true));
  Type nll2 = -sum(dnbinom2(Y2, mu2, var2, true));
  
  Type nll = 0;
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));

  matrix<Type> Cor(2,2);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(phi));

  return nll1 + nll2 + nll;
}
