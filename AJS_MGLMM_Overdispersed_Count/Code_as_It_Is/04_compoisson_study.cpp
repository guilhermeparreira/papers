// COM-Poisson full
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(mu);
  
  Type d0_log = -sum(dcompois2(Y, mu, Type(1.0), true));
  // .13539
  
  // Type d0_prob = dcompois2(Type(1), Type(9.5212), Type(.5), false);
  // REPORT(d0_prob);
  
  // Type d1_exp = dcompois(Type(1), exp(Type(3)), Type(.5), true);
  // REPORT(d1_exp);
  // 
  // 
  // 
  // Type d2 = dcompois2(Type(1), Type(3), Type(.5), true);
  // REPORT(d2); 
  // // .17167
  // 
  // Type d2_exp = dcompois2(Type(1), exp(Type(3)), Type(.5), true);
  // REPORT(d2_exp); 
  // 
  // Type d1 = dcompois(Type(1), Type(3), Type(.5), true);
  // REPORT(d1);
  // // .13539
  // 
  // Type d1_exp = dcompois(Type(1), exp(Type(3)), Type(.5), true);
  // REPORT(d1_exp);
  
  return d0_log;
}
