#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double samplePowerBinomialCpp(int n, double a1, double b1, double a2, double b2, double a, double b, double pi0, double pi1, double c){

  double power = 0;

  for(int y1 = 0; y1 < n+1; ++y1){
    for(int y2 = 0; y2 < n+1; ++y2){
      
      if(
        Rf_lbeta(y1 + y2 + a, 2 * n - y1 - y2 + b) - Rf_lbeta(a, b) -
          Rf_lbeta(y1 + a1, n - y1 + b1) + Rf_lbeta(a1, b1) -
          Rf_lbeta(y2 + a2, n - y2 + b2) + Rf_lbeta(a2, b2) <= log(c * pi1/pi0)
      ){
        power +=  exp(
          Rf_lchoose(n, y1) + Rf_lbeta(y1 + a1, n - y1 + b1) -
          Rf_lbeta(a1, b1) + Rf_lchoose(n, y2) + 
          Rf_lbeta(y2 + a2, n - y2 + b2) - Rf_lbeta(a2, b2)
        );
      }
    
    }
  }

  return power;
}
