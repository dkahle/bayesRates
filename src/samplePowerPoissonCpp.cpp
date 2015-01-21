#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double samplePowerPoissonCpp(int t, double a1, double b1, double a2, double b2, double a, double b, double pi0, double pi1, double c){

  double power = 0;
  
  // determine essential support of the poisson (/gamma mixture)
  int y1min = Rf_qnbinom(.00001, a1, b1/(t+b1), 1, 0);
  int y2min = Rf_qnbinom(.00001, a2, b2/(t+b2), 1, 0);
  
  int y1max = Rf_qnbinom(.99999, a1, b1/(t+b1), 1, 0);
  int y2max = Rf_qnbinom(.99999, a2, b2/(t+b2), 1, 0);

  // run for loops
  for(int y1 = y1min; y1 < y1max+1; ++y1){
    for(int y2 = y2min; y2 < y2max+1; ++y2){
      
      if(
        a*log(b) + lgamma(y1+y2+a) - lgamma(a) - (y1+y2+a)*log(2*t+b) + 
        lgamma(a1) + lgamma(a2) + (y1+a1)*log(t+b1) + (y2+a2)*log(t+b2) - 
        a1*log(b1) - a2*log(b2) - lgamma(y1+a1) - lgamma(y2+a2) <= log(c * pi1/pi0)    
      ){
        power += exp( 
          (y1+y2)*log(t) + a1*log(b1) + a2*log(b2) + lgamma(y1+a1) +
          lgamma(y2+a2) - lgamma(y1+1) - lgamma(y2+1) - 
          lgamma(a1) - lgamma(a2) - (y1+a1)*log(t+b1) - (y2+a2)*log(t+b2) 
        );
      }
    
    }
  }

  return power;
}
