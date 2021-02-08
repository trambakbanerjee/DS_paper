// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
using namespace std;
using namespace Rcpp;
using namespace arma;

// Gradient Estimation using Mahalanobis Distance
// [[Rcpp::export]]
mat ksd_md_cpp(const mat & x, const mat & K,double Omega_11,
               double Omega_12,double Omega_22,
               double eta)
{
  int n = x.n_rows;
  mat w_hat = zeros<mat>(n,2);
  mat d_K = zeros<mat>(n,2);
  mat x_i = zeros<mat>(n,2);
  mat I(n,n,fill::eye);
  for (int i=0; i<n; i++) {
    
    x_i.each_row()=x(i,span::all);
    mat x_x = x_i-x;
    double temp_1 = sum(K(span::all,i)%(x_x(span::all,0)*Omega_11+x_x(span::all,1)*Omega_12));
    double temp_2 = sum(K(span::all,i)%(x_x(span::all,1)*Omega_22+x_x(span::all,0)*Omega_12));
    d_K(i,0) = temp_1;
    d_K(i,1) = temp_2;
  }
  //mat A = K+eta*I;//inv_sympd(K+eta*I);
  w_hat = -solve(K+eta*I,d_K,solve_opts::likely_sympd);//A*d_K;
  return w_hat; 
}

// Gradient Estimation using Mahalanobis Distance (3d)
// [[Rcpp::export]]
mat ksd_3d_md_cpp(const mat & x, const mat & K,double Omega_11,
               double Omega_12,double Omega_13,double Omega_22,
               double Omega_23,double eta)
{
  int n = x.n_rows;
  mat w_hat = zeros<mat>(n,2);
  mat d_K = zeros<mat>(n,2);
  mat x_i = zeros<mat>(n,3);
  mat I(n,n,fill::eye);
  for (int i=0; i<n; i++) {
    
    x_i.each_row()=x(i,span::all);
    mat x_x = x_i-x;
    double temp_1 = sum(K(span::all,i)%(x_x(span::all,0)*Omega_11+x_x(span::all,1)*Omega_12+x_x(span::all,2)*Omega_13));
    double temp_2 = sum(K(span::all,i)%(x_x(span::all,1)*Omega_22+x_x(span::all,0)*Omega_12+x_x(span::all,2)*Omega_23));
    d_K(i,0) = temp_1;
    d_K(i,1) = temp_2;
  }
  //mat A = K+eta*I;//inv_sympd(K+eta*I);
  w_hat = -solve(K+eta*I,d_K,solve_opts::likely_sympd);//A*d_K;
  return w_hat; 
}

// Gradient Estimation using Mahalanobis Distance
// [[Rcpp::export]]
mat ksd_mdprod_cpp(const mat & x, const mat & K,double Omega_11,
               double Omega_12,double Omega_22,
               double eta)
{
  int n = x.n_rows;
  mat w_hat = zeros<mat>(n,2);
  mat d_K = zeros<mat>(n,2);
  mat x_i = zeros<mat>(n,2);
  mat I(n,n,fill::eye);
  for (int i=0; i<n; i++) {
    
    x_i.each_row()=x(i,span::all);
    mat x_x = x_i-x;
    double temp_1 = sum(K(span::all,i)%(x_x(span::all,0)*Omega_11+x_x(span::all,1)*Omega_12));
    double temp_2 = sum(K(span::all,i)%(x_x(span::all,1)*Omega_22+x_x(span::all,0)*Omega_12));
    d_K(i,0) = temp_1;
    d_K(i,1) = temp_2;
  }
  //mat A = K+eta*I;//inv_sympd(K+eta*I);
  w_hat = -solve(K+eta*I,d_K,solve_opts::likely_sympd);//A*d_K;
  return w_hat; 
}