// Binomial Negative Multivariate (n responses) same fixed Effects.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(X1);             // n x p
  DATA_MATRIX(X2);             // n x p
  DATA_MATRIX(X3);             // n x p
  DATA_MATRIX(X4);             // n x p
  DATA_MATRIX(X5);             // n x p
  DATA_MATRIX(X6);             // n x p
  DATA_MATRIX(X7);             // n x p
  DATA_MATRIX(X8);             // n x p
  DATA_MATRIX(X9);             // n x p
  DATA_MATRIX(X10);             // n x p
  DATA_MATRIX(X11);             // n x p
  DATA_MATRIX(X12);             // n x p
  DATA_MATRIX(X13);             // n x p
  DATA_MATRIX(X14);             // n x p
  DATA_MATRIX(X15);             // n x p
  DATA_MATRIX(X16);             // n x p
  DATA_MATRIX(X17);             // n x p
  DATA_MATRIX(X18);             // n x p
  DATA_MATRIX(X19);             // n x p
  DATA_MATRIX(X20);             // n x p
  DATA_MATRIX(X21);             // n x p
  DATA_MATRIX(X22);             // n x p
  DATA_MATRIX(X23);             // n x p
  DATA_MATRIX(X24);             // n x p
  DATA_MATRIX(X25);             // n x p
  DATA_MATRIX(X26);             // n x p
  DATA_MATRIX(X27);             // n x p
  DATA_MATRIX(X28);             // n x p
  DATA_MATRIX(X29);             // n x p
  DATA_MATRIX(X30);             // n x p
  DATA_MATRIX(X31);             // n x p
  DATA_MATRIX(X32);             // n x p
  DATA_MATRIX(X33);             // n x p
  DATA_MATRIX(X34);             // n x p
  DATA_MATRIX(X35);             // n x p
  DATA_MATRIX(X36);             // n x p
  DATA_MATRIX(X37);             // n x p
  DATA_MATRIX(X38);             // n x p
  DATA_MATRIX(X39);             // n x p
  DATA_MATRIX(X40);             // n x p
  DATA_MATRIX(X41);             // n x p
  
  DATA_MATRIX(Y);             // n x r
  
  
  PARAMETER_VECTOR(beta1);     // p x r
  PARAMETER_VECTOR(beta2);     // p x r
  PARAMETER_VECTOR(beta3);     // p x r
  PARAMETER_VECTOR(beta4);     // p x r
  PARAMETER_VECTOR(beta5);     // p x r
  PARAMETER_VECTOR(beta6);     // p x r
  PARAMETER_VECTOR(beta7);     // p x r
  PARAMETER_VECTOR(beta8);     // p x r
  PARAMETER_VECTOR(beta9);     // p x r
  PARAMETER_VECTOR(beta10);     // p x r
  PARAMETER_VECTOR(beta11);     // p x r
  PARAMETER_VECTOR(beta12);     // p x r
  PARAMETER_VECTOR(beta13);     // p x r
  PARAMETER_VECTOR(beta14);     // p x r
  PARAMETER_VECTOR(beta15);     // p x r
  PARAMETER_VECTOR(beta16);     // p x r
  PARAMETER_VECTOR(beta17);     // p x r
  PARAMETER_VECTOR(beta18);     // p x r
  PARAMETER_VECTOR(beta19);     // p x r
  PARAMETER_VECTOR(beta20);     // p x r
  PARAMETER_VECTOR(beta21);     // p x r
  PARAMETER_VECTOR(beta22);     // p x r
  PARAMETER_VECTOR(beta23);     // p x r
  PARAMETER_VECTOR(beta24);     // p x r
  PARAMETER_VECTOR(beta25);     // p x r
  PARAMETER_VECTOR(beta26);     // p x r
  PARAMETER_VECTOR(beta27);     // p x r
  PARAMETER_VECTOR(beta28);     // p x r
  PARAMETER_VECTOR(beta29);     // p x r
  PARAMETER_VECTOR(beta30);     // p x r
  PARAMETER_VECTOR(beta31);     // p x r
  PARAMETER_VECTOR(beta32);     // p x r
  PARAMETER_VECTOR(beta33);     // p x r
  PARAMETER_VECTOR(beta34);     // p x r
  PARAMETER_VECTOR(beta35);     // p x r
  PARAMETER_VECTOR(beta36);     // p x r
  PARAMETER_VECTOR(beta37);     // p x r
  PARAMETER_VECTOR(beta38);     // p x r
  PARAMETER_VECTOR(beta39);     // p x r
  PARAMETER_VECTOR(beta40);     // p x r
  PARAMETER_VECTOR(beta41);     // p x r

  PARAMETER_MATRIX(U);        // n x r
  PARAMETER_VECTOR(rho);      // r(r-1)/2
  PARAMETER_VECTOR(sigma);    //r -> R.E.
  PARAMETER_VECTOR(phi);      //r -> F.E.
  
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);
  
  int n = Y.rows(); // Number of rows
  int c = Y.cols();  // Number of cols

  // Building the linear predictor
  vector<Type> Xbeta1 = X1*beta1;
  vector<Type> Xbeta2 = X2*beta2;
  vector<Type> Xbeta3 = X3*beta3;
  vector<Type> Xbeta4 = X4*beta4;
  vector<Type> Xbeta5 = X5*beta5;
  vector<Type> Xbeta6 = X6*beta6;
  vector<Type> Xbeta7 = X7*beta7;
  vector<Type> Xbeta8 = X8*beta8;
  vector<Type> Xbeta9 = X9*beta9;
  vector<Type> Xbeta10 = X10*beta10;
  vector<Type> Xbeta11 = X11*beta11;
  vector<Type> Xbeta12 = X12*beta12;
  vector<Type> Xbeta13 = X13*beta13;
  vector<Type> Xbeta14 = X14*beta14;
  vector<Type> Xbeta15 = X15*beta15;
  vector<Type> Xbeta16 = X16*beta16;
  vector<Type> Xbeta17 = X17*beta17;
  vector<Type> Xbeta18 = X18*beta18;
  vector<Type> Xbeta19 = X19*beta19;
  vector<Type> Xbeta20 = X20*beta20;
  vector<Type> Xbeta21 = X21*beta21;
  vector<Type> Xbeta22 = X22*beta22;
  vector<Type> Xbeta23 = X23*beta23;
  vector<Type> Xbeta24 = X24*beta24;
  vector<Type> Xbeta25 = X25*beta25;
  vector<Type> Xbeta26 = X26*beta26;
  vector<Type> Xbeta27 = X27*beta27;
  vector<Type> Xbeta28 = X28*beta28;
  vector<Type> Xbeta29 = X29*beta29;
  vector<Type> Xbeta30 = X30*beta30;
  vector<Type> Xbeta31 = X31*beta31;
  vector<Type> Xbeta32 = X32*beta32;
  vector<Type> Xbeta33 = X33*beta33;
  vector<Type> Xbeta34 = X34*beta34;
  vector<Type> Xbeta35 = X35*beta35;
  vector<Type> Xbeta36 = X36*beta36;
  vector<Type> Xbeta37 = X37*beta37;
  vector<Type> Xbeta38 = X38*beta38;
  vector<Type> Xbeta39 = X39*beta39;
  vector<Type> Xbeta40 = X40*beta40;
  vector<Type> Xbeta41 = X41*beta41;
  
  matrix<Type> Xbeta(n, c);
  Xbeta.col(0) = Xbeta1;
  Xbeta.col(1) = Xbeta2;
  Xbeta.col(2) = Xbeta3;
  Xbeta.col(3) = Xbeta4;
  Xbeta.col(4) = Xbeta5;
  Xbeta.col(5) = Xbeta6;
  Xbeta.col(6) = Xbeta7;
  Xbeta.col(7) = Xbeta8;
  Xbeta.col(8) = Xbeta9;
  Xbeta.col(9) = Xbeta10;
  Xbeta.col(10) = Xbeta11;
  Xbeta.col(11) = Xbeta12;
  Xbeta.col(12) = Xbeta13;
  Xbeta.col(13) = Xbeta14;
  Xbeta.col(14) = Xbeta15;
  Xbeta.col(15) = Xbeta16;
  Xbeta.col(16) = Xbeta17;
  Xbeta.col(17) = Xbeta18;
  Xbeta.col(18) = Xbeta19;
  Xbeta.col(19) = Xbeta20;
  Xbeta.col(20) = Xbeta21;
  Xbeta.col(21) = Xbeta22;
  Xbeta.col(22) = Xbeta23;
  Xbeta.col(23) = Xbeta24;
  Xbeta.col(24) = Xbeta25;
  Xbeta.col(25) = Xbeta26;
  Xbeta.col(26) = Xbeta27;
  Xbeta.col(27) = Xbeta28;
  Xbeta.col(28) = Xbeta29;
  Xbeta.col(29) = Xbeta30;
  Xbeta.col(30) = Xbeta31;
  Xbeta.col(31) = Xbeta32;
  Xbeta.col(32) = Xbeta33;
  Xbeta.col(33) = Xbeta34;
  Xbeta.col(34) = Xbeta35;
  Xbeta.col(35) = Xbeta36;
  Xbeta.col(36) = Xbeta37;
  Xbeta.col(37) = Xbeta38;
  Xbeta.col(38) = Xbeta39;
  Xbeta.col(39) = Xbeta40;
  Xbeta.col(40) = Xbeta41;
    
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
