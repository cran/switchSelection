// [[Rcpp::depends(mnorm)]]
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <mnorm.h>
#include "helpFunctions.h"
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export(rng = false)]]
arma::field<arma::vec> findGroup(NumericMatrix z, 
                                 NumericMatrix groups, 
                                 const int adj = 1)
{
  const int n_groups = groups.nrow();
  const int n_obs = z.nrow();
  const int n_eq = z.ncol();
  
  arma::field<arma::colvec> ind_g(n_groups);
  arma::vec n_g(n_groups, arma::fill::zeros);
  for (int i = 0; i < n_groups; i++)
  {
    ind_g(i) = arma::vec(n_obs);
  }

  for (int i = 0; i < n_obs; i++)
  {
    for (int g = 0; g < n_groups; g++)
    {
      int correct = 0;
      for (int j = 0; j < n_eq; j++)
      {
        if (z(i, j) == groups(g, j))
        {
          correct++;
        }
        else
        {
          break;
        }
      }
      if (correct == n_eq)
      {
        ind_g(g).at(n_g.at(g)) = i + adj;
        n_g.at(g)++;
        break;
      }
    }
  }
  
  for (int i = 0; i < n_groups; i++)
  {
    ind_g(i).shed_rows(n_g.at(i), n_obs - 1);
  }
  
  return(ind_g);
}

// [[Rcpp::export(rng = false)]]
LogicalVector matrixInMatrix(NumericMatrix x, 
                             NumericMatrix y)
{
  // Get dimensional data
  const int x_row = x.nrow();
  const int x_col = x.ncol();
  
  const int y_row = y.nrow();
  
  // Create vector to store indexes
  LogicalVector ind(x_row);
  
  // Perform main routine
  for (int i = 0; i < x_row; i++)
  {
    for (int j = 0; j < y_row; j++)
    {
      int counter = 0;
      for (int t = 0; t < x_col; t++)
      {
        if (x(i, t) == y(j, t))
        {
          counter++;
        }
      }
      if (counter == x_col)
      {
        ind[i] = true;
        break;
      }
    }
  }
  
  return(ind);
}
