#include <Rcpp.h>
#include <stan/math.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector normLogLikAD(double y_, double mu_, double sigma_) {
  NumericVector ret(4);

  stan::math::var y = y_;
  stan::math::var mu = mu_;
  stan::math::var sigma = sigma_;
  stan::math::var lp = 0;

  lp += -0.5* ((y-mu)/sigma) * ((y-mu)/sigma) - log(sigma) - 0.5 * log(2*PI);

  // Compute Function Value
  ret[0] = lp.val();

  // Compute Gradient
  lp.grad();
  ret[1] = y.adj();
  ret[2] = mu.adj();
  ret[3] = sigma.adj();

  // Memory is allocated on a global stack
  stan::math::recover_memory();
  stan::math::ChainableStack::memalloc_.free_all();

  return ret;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
NumericVector(1,1,1)
*/
