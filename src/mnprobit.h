#ifndef switchSelection_mnprobit_H
#define switchSelection_mnprobit_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

double lnL_mnprobit(const arma::vec par,
                    const List control_lnL,
                    const String out_type,
                    const int n_sim,
                    const int n_cores,
                    const Nullable<List> regularization);

NumericVector grad_mnprobit(const arma::vec par,
                            const List control_lnL,
                            const String out_type,
                            const int n_sim,
                            const int n_cores,
                            const Nullable<List> regularization);

#endif
