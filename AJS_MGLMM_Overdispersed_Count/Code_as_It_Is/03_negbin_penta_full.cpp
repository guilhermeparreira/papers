// Bivariate Negative Binomial full
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y1);
  DATA_VECTOR(Y2);
  DATA_VECTOR(Y3);
  DATA_VECTOR(Y4);
  DATA_VECTOR(Y5);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(beta3);
  PARAMETER_VECTOR(beta4);
  PARAMETER_VECTOR(beta5);
  PARAMETER_MATRIX(U);
  PARAMETER_VECTOR(rho);
  PARAMETER_VECTOR(sigma);
  PARAMETER_VECTOR(phi);
  int n = Y1.size();
  
  vector<Type> mu1 = exp(X*beta1 + U.col(0).array());
  vector<Type> mu2 = exp(X*beta2 + U.col(1).array());
  vector<Type> mu3 = exp(X*beta3 + U.col(2).array());
  vector<Type> mu4 = exp(X*beta4 + U.col(3).array());
  vector<Type> mu5 = exp(X*beta5 + U.col(4).array());
  
  vector<Type> var1 = mu1 + (mu1*mu1)/exp(phi(0));
  vector<Type> var2 = mu2 + (mu2*mu2)/exp(phi(1));
  vector<Type> var3 = mu3 + (mu3*mu3)/exp(phi(2));
  vector<Type> var4 = mu4 + (mu4*mu4)/exp(phi(3));
  vector<Type> var5 = mu5 + (mu5*mu5)/exp(phi(4));
  
  Type nll = 0;
  // parallel_accumulator<Type> nll(this);
    
  nll -= sum(dnbinom2(Y1, mu1, var1, true));
  nll -= sum(dnbinom2(Y2, mu2, var2, true));
  nll -= sum(dnbinom2(Y3, mu3, var3, true));
  nll -= sum(dnbinom2(Y4, mu4, var4, true));
  nll -= sum(dnbinom2(Y5, mu5, var5, true));
  
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));

  matrix<Type> Cor(5,5);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(phi));

  return nll;
}
