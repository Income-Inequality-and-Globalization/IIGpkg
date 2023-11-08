#include "00_ghk_sampler_helper.h"
// [[Rcpp::export]]
double norm_rej(const double a, const double b) {
  double x;
  do {
    x = R::rnorm(0.0, 1.0);  
  } while (x < a || x > b);
  return x;
}

// [[Rcpp::export]]
double unif_rej(const double a, const double b) {
  double x, u, r;
  do {
    x = R::runif(a, b);
    u = R::runif(0.0, 1.0);
    if (a <= 0.0 && b >= 0.0) 
      r = exp(-x*x/2.0);
    else if (a > 0.0)
      r = exp(-(x*x-a*a)/2.0);
    else 
      r = exp(-(x*x-b*b)/2.0);
  } while (u > r);
  return x;
}

// [[Rcpp::export]]
double halfnorm_rej(const double a, const double b) {
  double x;
  do {
    x = fabs(R::rnorm(0.0, 1.0));  
  } while (x < a || x > b);
  return x;
}

// [[Rcpp::export]]
double exp_rej(const double a, const double b) {
  double lambda = (a+sqrt(a*a+4.0))/2.0;
  double x, u, r;
  do {
    x = a+R::rweibull(1, 1.0/lambda);
    u = R::runif(0.0, 1.0);
    r = exp(-pow(x-lambda,2)/2.0);
  } while (u > r || x > b);
  return x;
}


// [[Rcpp::export]]
arma::vec rtnormcpp(const arma::vec& mean, 
                    const double sd, 
                    const arma::vec& lower, 
                    const arma::vec& upper) {
  const unsigned int n=mean.n_elem, n1=lower.n_elem, n2=upper.n_elem;
  
  // check dimensions
  if (!(n==n1 || n1==1)) {
    Rcpp::stop("lower must be a row vector or have the same number of rows as mean");
  }
  
  if (!(n1==n2)) {
    Rcpp::stop("lower and upper must have the same number of rows");
  }
  
  // check boundary conditions
  if (any(lower >= upper)) {
    Rcpp::stop("lower bound must be smaller than upper bound");
  }
  
  static const double pi = 3.141592653589793238462643383279;
  
  auto imp_case1 = [](double a, double b) { 
    double w;
    if (a < 0.0) w = norm_rej(a, b);
    else if (a < 0.25696) w = halfnorm_rej(a, b);
    else w = exp_rej(a, b);
    return w;
  };
  
  auto imp_case2 = [](double a, double b) { 
    double w;
    if (b <= a + sqrt(2*pi)) w = unif_rej(a, b);
    else w = norm_rej(a, b);
    return w;
  };
  
  auto imp_case3 = [](double a, double b) { 
    double w, lambda;
    if (a <= 0.25696) {
      if (b <= a+sqrt(pi/2)*exp(a*a/2)) w = unif_rej(a,b);
      else w = halfnorm_rej(a,b);
    } else {
      lambda = (a+sqrt(a*a+4.0))/2.0;      
      if (b <= a+1/lambda*exp((a*a-a*sqrt(a*a+4))/4+0.5)) w = unif_rej(a,b);
      else w = exp_rej(a,b);
    }
    return w;
  };
  
  double a, b, w;
  arma::vec x(n);
  for (arma::uword i=0; i<n; i++) {
    if (n1==n) {
      a = (lower(i) - mean(i))/sd;
      b = (upper(i) - mean(i))/sd;
    } else {
      a = (lower(0) - mean(i))/sd;
      b = (upper(0) - mean(i))/sd;
    }
    
    if (std::isinf(a) || std::isinf(b)) {
      if (std::isinf(b)) w = imp_case1(a,b);
      else w = -imp_case1(-b,-a); // case 4
    } else {
      if (a<0.0 && b>0.0) w = imp_case2(a,b);
      else if (a>=0.0) w = imp_case3(a,b);
      else w = -imp_case3(-b,-a); // case 5
    }
    
    x(i) = mean(i)+sd*w;
  }
  return x;
}