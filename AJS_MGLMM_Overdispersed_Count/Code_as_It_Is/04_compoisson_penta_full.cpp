// COM-Poisson full
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
  PARAMETER_VECTOR(nu);
  int n = Y1.size();
  
  vector<Type> mu1 = exp(X*beta1 + U.col(0).array());
  vector<Type> mu2 = exp(X*beta2 + U.col(1).array());
  vector<Type> mu3 = exp(X*beta3 + U.col(2).array());
  vector<Type> mu4 = exp(X*beta4 + U.col(3).array());
  vector<Type> mu5 = exp(X*beta5 + U.col(4).array());
  
  Type nll = 0; 
  // parallel_accumulator<Type> nll(this);
  nll -= sum(dcompois2(Y1, mu1, exp(nu(0)), true));
  nll -= sum(dcompois2(Y2, mu2, exp(nu(1)), true));
  nll -= sum(dcompois2(Y3, mu3, exp(nu(2)), true));
  nll -= sum(dcompois2(Y4, mu4, exp(nu(3)), true));
  nll -= sum(dcompois2(Y5, mu5, exp(nu(4)), true));
  
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));

  matrix<Type> Cor(5,5);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(nu));
  return nll;
}
