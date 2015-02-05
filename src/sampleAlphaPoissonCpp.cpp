#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sampleAlphaPoissonCpp(int t, double a1, double b1, double a2, double b2, double a, double b, double pi0, double pi1, double c){

  double alpha = 0;
  
  // determine essential support of the poisson (/gamma mixture)
  int y1min = Rf_qnbinom(.00001, a, b/(t+b), 1, 0);
  int y2min = Rf_qnbinom(.00001, a, b/(t+b), 1, 0);
  
  int y1max = Rf_qnbinom(.99999, a, b/(t+b), 1, 0);
  int y2max = Rf_qnbinom(.99999, a, b/(t+b), 1, 0);

  // run for loops
  for(int y1 = y1min; y1 < y1max+1; ++y1){
    for(int y2 = y2min; y2 < y2max+1; ++y2){
      
      if(
        a*log(b) + lgamma(y1+y2+a) - lgamma(a) - (y1+y2+a)*log(2*t+b) + 
        lgamma(a1) + lgamma(a2) + (y1+a1)*log(t+b1) + (y2+a2)*log(t+b2) - 
        a1*log(b1) - a2*log(b2) - lgamma(y1+a1) - lgamma(y2+a2) <= log(c * pi1/pi0)    
      ){
        alpha += exp( 
          (y1+y2)*log(t) + a*log(b) + lgamma(y1+y2+a) -
            lgamma(y1+1) - lgamma(y2+1) - lgamma(a) - 
            (y1+y2+a)*log(2*t+b) 
        );
      }
    
    }
  }

  return alpha;
}

