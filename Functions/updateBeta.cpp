#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
mat updateBeta(arma::mat xxt, arma::mat x, arma::vec winv, arma::mat z, 
	 arma::mat tildebrho, arma::mat xtw, const double rhoinv){
	mat temp = diagmat(winv);
	mat t1 = xtw*z + tildebrho; 
	mat out = t1 - rhoinv*x.t()*solve(temp + xxt, x*t1, solve_opts::allow_ugly);
	return out*rhoinv; 
}


