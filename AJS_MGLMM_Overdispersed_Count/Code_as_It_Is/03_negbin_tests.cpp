// Binomial Negative Multivariate (n responses) same fixed Effects.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_IVECTOR(Y1);
  PARAMETER(mu);
  PARAMETER(var);
  Type nll = 0;
  int n = Y1.size();
  Type var2 = exp(var);
  // Binomial Negative model
  for (int i=0; i<n; i++){
    nll -= dnbinom2(Type(Y1(i)), mu, var2);
  }
  return nll;
  // // TEST 1 -> Understanding 
  // // dnbinom2 is defined as Mu and Var
  // Type mu_v0 = dnbinom2(Type(2), Type(1.285714), Type(1.836735), true); // NAN. SO, phi is size, because it does
  // REPORT(mu_v0);
  // 
  // // Size and prob
  // Type size_prob = dnbinom(Type(2), Type(3.0), Type(.7), true);
  // REPORT(size_prob);
  
  // ------------- TRASH, TO BE DELETED
  
  // Trying as mu and phi
  // mu = 1, r = 1
  // Type mu_var1 = dnbinom2(Type(2), Type(1.0), Type(1.0), true); // NAN. SO, phi is size, because it does
  // REPORT(mu_var1);
  // 
  // mu = 1, r = .5
  // Type mu_var_05 = dnbinom2(Type(2), Type(1.0), Type(.5), true); // NAN
  // REPORT(mu_var_05);
  
  // Type mu_var3 = dnbinom2(Type(2), Type(1.0), Type(3.0), true);
  // REPORT(mu_var3);
  
  // Type mu_var5 = dnbinom2(Type(2), Type(1.0), Type(5.0), true);
  // REPORT(mu_var5);
  
  // mu = 1.285714, r = .5
  // Type mu_128_r_05 = dnbinom2(Type(2), Type(285714), Type(3.0), true); // NAN
  // REPORT(mu_128_r_05);
}

