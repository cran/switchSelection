#ifndef switchSelection_msel_H
#define switchSelection_msel_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

double lnL_msel(const arma::vec par,
                const List control_lnL,
                const String out_type,
                const int n_sim,
                const int n_cores,
                const Nullable<List> regularization);

NumericVector grad_msel(const arma::vec par,
                        const List control_lnL,
                        const String out_type,
                        const int n_sim,
                        const int n_cores,
                        const Nullable<List> regularization);

#endif
