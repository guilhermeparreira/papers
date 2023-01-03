template<class Type>
Type c_mu_phi_ctt(Type mu, Type sigma){
  // # Input of the function call()
   int lmu = 1;
   int maxV = 500;
   // CppAD::vector<Type> ans(ly);
   Type ans = 0;
   int ly2 = maxV + 1;

   // // Variables declared inside of function
   CppAD::vector<Type> ylogofy(ly2);
   CppAD::vector<Type> lga(ly2);
   CppAD::vector<Type> ym(ly2);

  // Function as it was written in C++
    for (int j=0 ; j <ly2 ; j++) {
      ylogofy[j] = j * ((j==0)? 1 : log(j));
      Type j1 = j + 1.0; // Converting to double (Type)
      lga[j] = lgamma(j1);
      ym[j] = (j - ylogofy[j]);
    }
    for (int i=0; i<lmu; i++){
      Type sumC = 0;
      Type mus = mu / sigma;
      Type lsig2 = -0.5 * log(sigma);
      Type lmus = log(mu) / sigma - 1;
      Type invs = 1 / sigma;
      Type ls = lsig2 - mus;
      for (int j=0 ; j <ly2 ; j++){
        sumC += exp(ls - lga[j] + ylogofy[j] + j * lmus + invs * ym[j]);
      }
      ans = pow(sumC,-1);
      ans = log(ans);
    }
    
   return ans;
}
VECTORIZE2_tt(c_mu_phi_ctt)

template <class Type>
Type doublepoisson(Type y, Type mu, Type sigma)
{
  Type c_mu_phi = c_mu_phi_ctt(mu, sigma); //Já é em log()
  Type logy = 1;
  if (y != 0)
  {
    logy = log(y);
  }
  Type logvero = -.5 * log(sigma) - (mu / sigma) - lgamma(y + 1) + y * logy - y + (y * log(mu)) / sigma + y / sigma - (y * logy) / sigma +
    c_mu_phi;
  return logvero;
}
VECTORIZE3_ttt(doublepoisson)
// VECTORIZE1_t(doublepoisson)

    // template<class Type>
    // Type doublepoisson(Type y, Type mu, Type sigma){
    //   Type c_mu_phi = ((1-(1/sigma))/(12*mu*(1/sigma)))*( 1 + (1/ (mu*(1/sigma)) ));
    //
    //   Type logy = 1;
    //   if (y != 0) {
    //     logy = log(y);
    //   }
    //
    //   Type logvero = -.5 * log(sigma) - (mu/sigma) - lgamma(y + 1) + y *
    //     logy - y + (y * log(mu))/sigma + y/sigma - (y * logy)/sigma +
    //     c_mu_phi;
    //   return logvero;
    // }
    // VECTORIZE3_ttt(doublepoisson)
    // VECTORIZE1_t(doublepoisson)