#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param length_of_chains the number of samples
//' @param from_point how many number ahead are to abandon
//' @param a,b the shape of distribution Y = beta(x+a,n-x+b)
//' @param x_range = n the X values uses in distribution X = Binomial(n, y)
//' @param mu_x,mu_y the expectation of X and Y
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' ZC <- cp_gibbsC()
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix cp_gibbsC(int length_of_chains=1e4,int from_point=1000,int a =1,int b = 1,int x_range=10,float mu_x = -1,float mu_y= -1) {
  // x,y can not be negative, and the NA is not accepted, so I choose -1 to avoid no initial values.
  int to_point = length_of_chains + from_point;
  NumericMatrix Z(to_point, 2);
  // Initial Value
  if (mu_x == -1 || mu_y == -1){
    Z(0,1) = a/(a+b);
    Z(0,0) = x_range * Z(0,1);
  } else {
    Z(0,0) = mu_x;
    Z(0,1) = mu_y;
  }
  for(int i = 1; i < to_point; i++) {
      // Update x
      float Zy = Z(i-1,1);
      Z(i,0) = rbinom(1,x_range,Zy)[0];
      // Update y
      float Zx = Z(i,0);
      Z(i,1) = rbeta(1,Zx+a,x_range-Zx+b)[0];
  }
  return(Z(Range(from_point,to_point-1),_));
  // 1000 - 10999: 10000 in total.
}
