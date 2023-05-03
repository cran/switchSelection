#ifndef switchSelection_helpFunctions_H
#define switchSelection_helpFunctions_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

arma::field<arma::vec> findGroup(NumericMatrix z, 
                                 NumericMatrix groups,
                                 const int adj);

LogicalVector matrixInMatrix(NumericMatrix x, 
                             NumericMatrix y);

#endif
