#ifndef switchSelection_mvoprobit_H
#define switchSelection_mvoprobit_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

double lnL_mvoprobit(const arma::vec par,
                     const List control_lnL,
                     const String out_type,
                     const int n_sim,
                     const int n_cores,
                     const Nullable<List> regularization);

NumericVector grad_mvoprobit(const arma::vec par,
                             const List control_lnL,
                             const String out_type,
                             const int n_sim,
                             const int n_cores,
                             const Nullable<List> regularization);

#endif
