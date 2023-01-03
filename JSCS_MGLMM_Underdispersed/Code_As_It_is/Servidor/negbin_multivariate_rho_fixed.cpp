// Binomial Negative Multivariate (n responses) same fixed Effects.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(X);             // n x p
  DATA_MATRIX(Y);             // n x r
  PARAMETER_MATRIX(beta);     // p x r
  PARAMETER_MATRIX(U);        // n x r
  DATA_VECTOR(rho);      // r(r-1)/2
  PARAMETER_VECTOR(sigma);    //r -> R.E.
  PARAMETER_VECTOR(phi);      //r -> F.E.
  
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);
  
  int n = Y.rows(); // Number of rows
  int c = Y.cols();  // Number of cols
  matrix<Type> Xbeta(n, c);
  Xbeta = X*beta;
  
  vector<Type> Yj(n);
  vector<Type> mu(n);
  vector<Type> var(n);
  for (int j = 0; j<c; j++){
    Yj = Y.col(j);
    mu = exp(Xbeta.col(j).array() + U.col(j).array());
    var = mu + (mu.array()*mu.array())/exp(phi(j));
    nll -= sum(dnbinom2(Yj, mu, var, true));
  }
  for(int i = 0; i < n; i++)
    nll += VECSCALE(UNSTRUCTURED_CORR(rho), sigma)(U.row(i));
  matrix<Type> Cor(c,c);
  Cor = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Cor);
  ADREPORT(exp(phi));
  return nll;
}
