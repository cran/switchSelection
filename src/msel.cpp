// [[Rcpp::depends(mnorm)]]
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <mnorm.h>
#include "helpFunctions.h"
using namespace Rcpp;

#ifdef _OPENMP
// [[Rcpp::plugins(openmp)]]
#endif

//' Log-likelihood Function of Multivariate Ordered Probit Model
//' @description Calculates log-likelihood function of multivariate ordered
//' probit model.
//' @param par vector of parameters.
//' @param control_lnL list with some additional parameters.
//' @param out_type string represeint the output type of the function.
//' @param n_sim the number of random draws for multivariate 
//' normal probabilities.
//' @param n_cores the number of cores to be used. 
//' @param regularization list of regularization parameters.
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix lnL_msel(const arma::vec par,
                       const List control_lnL,
                       const String out_type = "val",
                       const int n_sim = 1000,
                       const int n_cores = 1,
                       const Nullable<List> regularization = R_NilValue)
{
  // Determine whether gradient and Jacobian should be estimated
  const bool is_grad = (out_type == "grad");
  const bool is_jac  = (out_type == "jac");
  const bool is_diff = is_grad | is_jac;

  // Get the regressors (covariates) for the ordered equations
  const arma::field<arma::mat> W     = control_lnL["W"];
  const arma::field<arma::mat> W_var = control_lnL["W_var"];

  // Get the regressors (covariates) and the dependent variable
  // for the continuous equations
  const arma::mat y              = control_lnL["y"];
  const arma::field<arma::mat> X = control_lnL["X"];

  // Get the regressors (covariates) for the multinomial equations
  const arma::mat W_mn = control_lnL["W_mn"];

  // Get information on the dimensions of the data
  const int n_par                = control_lnL["n_par"];
  const int n_obs                = control_lnL["n_obs"];
  const arma::vec n_regimes      = control_lnL["n_regimes"];
  const arma::vec n_coef         = control_lnL["n_coef"];
  const arma::vec n_coef2        = control_lnL["n_coef2"];
  const int n_coef3              = control_lnL["n_coef3"];
  const int n_eq                 = control_lnL["n_eq"];
  const int n_eq2                = control_lnL["n_eq2"];
  const int n_eq3                = control_lnL["n_eq3"];
  const int n_groups             = control_lnL["n_groups"];
  const arma::mat groups         = control_lnL["groups"];
  const arma::mat groups2        = control_lnL["groups2"];
  const arma::vec groups3        = control_lnL["groups3"];
  const arma::vec n_eq_g         = control_lnL["n_eq_g"];
  const arma::vec n_eq2_g        = control_lnL["n_eq2_g"];
  const arma::vec n_eq_all_g     = control_lnL["n_eq_all_g"];
  const arma::vec n_cuts_eq      = control_lnL["n_cuts_eq"];
  const bool is1                 = control_lnL["is1"];
  const bool is2                 = control_lnL["is2"];
  const bool is3                 = control_lnL["is3"];

  // Create matrices to store the Jacobian if need
  arma::mat jac;
  if (is_jac | is_grad)
  {
    jac = arma::mat(n_obs, n_par);
  }

  // Get indexes of the observations for different groups
  const arma::field<arma::uvec> ind_g = control_lnL["ind_g"];
  const arma::vec n_obs_g             = control_lnL["n_obs_g"];

  // Get indexes of the observed equations for each group
  const arma::field<arma::uvec> ind_eq     = control_lnL["ind_eq"];
  const arma::field<arma::uvec> ind_eq2    = control_lnL["ind_eq2"];
  const arma::field<arma::uvec> ind_eq_all = control_lnL["ind_eq_all"];

  // Get the parameters of the regularization
  List regularization1(regularization);
  arma::uvec ridge_ind;
  arma::uvec lasso_ind;
  arma::vec ridge_scale;
  arma::vec ridge_location;
  arma::vec lasso_scale;
  arma::vec lasso_location;
  bool is_ridge = false;
  bool is_lasso = false;
  const bool is_regularization = regularization1.size() > 0;
  if (is_regularization)
  {
    // Ridge parameters
    is_ridge = regularization1.containsElementNamed("ridge_ind");

    if (is_ridge)
    {
      // Indexes
      arma::uvec ridge_ind_tmp = regularization1["ridge_ind"];
      ridge_ind                = ridge_ind_tmp;

      // Scales
      arma::vec ridge_scale_tmp = regularization1["ridge_scale"];
      ridge_scale               = ridge_scale_tmp;

      // Locations
      arma::vec ridge_location_tmp = regularization1["ridge_location"];
      ridge_location               = ridge_location_tmp;
    }
    
    // Lasso parameters
    is_lasso = regularization1.containsElementNamed("lasso_ind");

    if (is_lasso)
    {
      // Indexes
      arma::uvec lasso_ind_tmp = regularization1["lasso_ind"];
      lasso_ind                = lasso_ind_tmp;
      
      // Scales
      arma::vec lasso_scale_tmp = regularization1["lasso_scale"];
      lasso_scale               = lasso_scale_tmp;
      
      // Locations
      arma::vec lasso_location_tmp = regularization1["lasso_location"];
      lasso_location               = lasso_location_tmp;
    }
  }

  // Determine the indexes of the coefficients and the
  // covariances in the vector of parameters
  const arma::field<arma::uvec> coef_ind       = control_lnL["coef_ind"];
  const arma::field<arma::umat> coef2_ind      = control_lnL["coef2_ind"];
  const arma::umat coef3_ind                   = control_lnL["coef3_ind"];
  const arma::uvec sigma_ind                   = control_lnL["sigma_ind"];
  const arma::field<arma::uvec> sigma2_ind     = control_lnL["sigma2_ind"];
  const arma::uvec sigma3_ind                  = control_lnL["sigma3_ind"];
  const arma::mat sigma_omit                   = control_lnL["sigma_omit"];
  const arma::umat sigma_ind_mat               = control_lnL["sigma_ind_mat"];
  const arma::field<arma::umat> sigma2_ind_mat = control_lnL["sigma2_ind_mat"];
  const arma::umat sigma3_ind_mat              = control_lnL["sigma3_ind_mat"];
  const arma::field<arma::uvec> var2_ind       = control_lnL["var2_ind"];
  const arma::field<arma::umat> cov2_ind       = control_lnL["cov2_ind"];
  const arma::field<arma::uvec> cuts_ind       = control_lnL["cuts_ind"];
  // Store the coefficients for each equation
  arma::field<arma::vec> coef(n_eq);
  if (is1)
  {
    for (int i = 0; i < n_eq; i++)
    {
      coef.at(i) = par.elem(coef_ind(i));
    }
  }

  // Construct the covariance matrix for the ordered equations
  arma::mat sigma(n_eq + n_eq2, n_eq + n_eq2);
  sigma.diag().ones();
  if (is1)
  {
    arma::vec sigma_vec = par.elem(sigma_ind);
    if (n_eq > 1)
    {
      int counter = 0;
      for (int i = 1; i < n_eq; i++)
      {
        for (int j = 0; j < i; j++)
        {
          if (sigma_omit.at(i, j) != 1)
          {
            sigma.at(i, j) = sigma_vec.at(counter);
            sigma.at(j, i) = sigma.at(i, j);
            counter++;
          }
        }
      }
    }
  }
  NumericMatrix sigma_R = wrap(sigma);

  // Assign cuts and check that they are sorted
  arma::field<arma::vec> cuts(n_eq);
  if (is1)
  {
    for (int i = 0; i < n_eq; i++)
    {
      cuts.at(i) = par.elem(cuts_ind(i));
      if (!cuts.at(i).is_sorted("strictascend"))
      {
        if (is_diff)
        {
          NumericMatrix out_tmp(n_par,1);
          std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
          return(out_tmp);
        }
        NumericMatrix out_tmp(1,1);
        out_tmp(0, 0) = -(1e+100);
        return(out_tmp);
      }
    }
  }

  // Get variables related to the heteroscedasticity if need
  const LogicalVector is_het                 = control_lnL["is_het"];
  const arma::field<arma::uvec> coef_var_ind = control_lnL["coef_var_ind"];
  arma::field<arma::vec> coef_var(n_eq);
  if (is1)
  {
    for (int i = 0; i < n_eq; i++)
    {
      if (is_het[i])
      {
        coef_var.at(i) = par.elem(coef_var_ind(i));
      }
    }
  }

  // Assign coefficients associated with the continuous equations
  arma::field<arma::mat> coef2(n_eq2);
  if (is2)
  {
    for (int i = 0; i < n_eq2; i++)
    {
      coef2.at(i) = arma::mat(n_regimes[i], n_coef2[i]);
      for (int j = 0; j < n_regimes[i]; j++)
      {
        coef2.at(i).row(j) = par.elem(coef2_ind(i).row(j)).t();
      }
    }
  }

  // Marginal distributions
  const CharacterVector marginal_names           = control_lnL["marginal_names"];
  const arma::ivec marginal_par_n                = control_lnL["marginal_par_n"];
  const arma::field<arma::uvec> marginal_par_ind = control_lnL["marginal_par_ind"];
  const bool is_marginal                         = marginal_names.size() > 0;
  arma::field<arma::vec> marginal_par(n_eq + n_eq2);
  List marginal_par_list;
  if (is_marginal)
  {
    bool is_marginal_par_valid = true;
    for (int i = 0; i < n_eq; i++)
    {
      if (marginal_par_n.at(i) > 0)
      {
        // Assign marginal distribution parameters
        marginal_par.at(i) = par.elem(marginal_par_ind.at(i));
        
        // Validate marginal distribution parameters
        if (marginal_names[i] == "student")
        {
          if (marginal_par.at(i).at(0) <= 2)
          {
            is_marginal_par_valid = false;
          }
        }
      }
    }
    // Stop if there are incorrect marginal
    // distribution parameters
    if (!is_marginal_par_valid)
    {
      if (is_diff)
      {
        NumericMatrix out_tmp(n_par, 1);
        std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
        return(out_tmp);
      }
      NumericMatrix out_tmp(1,1);
      out_tmp(0, 0) = -(1e+100);
      return(out_tmp);
    }
    // store the results into Rcpp list
    marginal_par_list = wrap(marginal_par);
  }

  // Store the coefficients of the multinomial equations
  arma::mat coef3;
  if (is3)
  {
    coef3 = arma::mat(n_eq3 - 1, n_coef3);
    for (int i = 0; i < (n_eq3 - 1); i++)
    {
      coef3.row(i) = par.elem(coef3_ind.row(i)).t();
    }
  }

  // Estimate the linear indexes of the ordered and
  // continuous equations
  arma::mat li_mean(n_obs, n_eq);
  arma::mat li_var(n_obs, n_eq);
  arma::mat li_lower(n_obs, n_eq);
  arma::mat li_upper(n_obs, n_eq);
  arma::mat li_y(n_obs, n_eq2);
  for (int i = 0; i < n_groups; i++)
  {
    // Ordered equations
    if (n_eq_g(i) > 0)
    {
      for (int j = 0; j < n_eq_g.at(i); j++)
      {
        //Prepare the indexes
        int j_o             = ind_eq(i).at(j);
        int group_j         = groups(i, j_o);
        arma::uvec j_o_uvec = {(unsigned int)j_o};
        // Calculate linear indexes for mean and variance parts
        li_mean.submat(ind_g(i), j_o_uvec) = W(j_o).rows(ind_g(i)) * 
                                             coef.at(j_o);
        if (is_het[j_o])
        {
          li_var.submat(ind_g(i), j_o_uvec) = arma::exp(W_var(j_o).rows(ind_g(i)) * 
                                              coef_var.at(j_o));
        }
        // Adjust lower limits for cuts and heteroscedasticity
        if (group_j == 0)
        {
          li_lower.submat(ind_g(i), j_o_uvec).fill(-arma::datum::inf);
        }
        else
        {
          li_lower.submat(ind_g(i), j_o_uvec) = cuts(j_o).at(group_j - 1) - 
                                                li_mean.submat(ind_g(i), 
                                                               j_o_uvec);
          if (is_het[j_o])
          {
            li_lower.submat(ind_g(i), j_o_uvec) /= li_var.submat(ind_g(i), 
                                                                 j_o_uvec);
          }
        }
        // Adjust upper limits for cuts and heteroscedasticity
        if (group_j == n_cuts_eq.at(j_o))
        {
          li_upper.submat(ind_g(i), j_o_uvec).fill(arma::datum::inf);
        }
        else
        {
          li_upper.submat(ind_g(i), j_o_uvec) = cuts(j_o).at(group_j) -
                                                li_mean.submat(ind_g(i), 
                                                               j_o_uvec);
          if (is_het[j_o])
          {
            li_upper.submat(ind_g(i), j_o_uvec) /= li_var.submat(ind_g(i), 
                                                                 j_o_uvec);
          }
        }
      }
    }
    // Continuous equations
    if (is2)
    {
      for (int j = 0; j < n_eq2_g.at(i); j++)
      {
        //Prepare indexes
        int j_o                         = ind_eq2(i).at(j);
        arma::uvec j_o_uvec             = {(unsigned int)j_o};
        li_y.submat(ind_g(i), j_o_uvec) = y.submat(ind_g(i), j_o_uvec) - 
                                          X(j_o).rows(ind_g(i)) * 
                                          coef2.at(j_o).row(groups2.at(i, j_o)).t();
      }
    }
  }

  // Vector to store the log-likelihoods
  arma::vec lnL;
  if ((out_type == "val") | is3)
  {
    lnL = arma::vec(n_obs);
  }

  // Get the type of the multinomial regression
  const String type3 = control_lnL["type3"];
  
  // Special routine for the multinomial logit
  if (is3 & (type3 == "logit"))
  {
    // Prepare the data
    arma::mat li_mn(n_obs, n_eq3);
    arma::mat li_exp(n_obs, n_eq3);
    
    // Calculate the linear indexes
    for (int i = 0; i < (n_eq3 - 1); i++)
    {
      li_mn.col(i) = W_mn * coef3.row(i).t();
    }
    li_mn.col(n_eq3 - 1).fill(0);
    
    // Get the exponents of the linear indexes
    for (int i = 0; i < (n_eq3 - 1); i++)
    {
      li_exp.col(i) = exp(li_mn.col(i));
    }
    li_exp.col(n_eq3 - 1).fill(1);
    
    // Estimate the probabilities
    arma::mat prob = li_exp.each_col() / sum(li_exp, 1);

    // Calculate the log-likelihood depending on the group
    for (int i = 0; i < n_groups; i++)
    {
      // Check the observability
      int groups3_int = groups3.at(i);
      if (groups3_int != -1)
      {
        // Calculate the log-likelihoods
        arma::uvec group3  = {(unsigned int)groups3_int};
        lnL.elem(ind_g(i)) = log(prob.submat(ind_g(i), group3));

        // Differentiate if need
        if (is_diff)
        {
          arma::mat W_i   = W_mn.rows(ind_g(i));
          arma::mat prob2 = arma::pow(prob, 2);
          for (int j = 0; j < (n_eq3 - 1); j++)
          {
            arma::uvec j_adj = {(unsigned int)j};
            if (groups3_int == j)
            {
              jac.submat(ind_g(i), coef3_ind.row(j)) = 
                W_i.each_col() % 
                ((prob.submat(ind_g(i), j_adj) - 
                prob2.submat(ind_g(i), j_adj)) / 
                prob.submat(ind_g(i), group3));
            }
            else
            {
              jac.submat(ind_g(i), coef3_ind.row(j)) = 
                W_i.each_col() % (-prob.submat(ind_g(i), j_adj));
            }
          }
        }
      }
    }
  }

  // Special routine for the multinomial probit
  if (is3 & (type3 == "probit"))
  {
    // Store the indexes of the worse alternatives for each best alternative
    arma::umat worse_ind(n_eq3 - 1, n_eq3);
    for (int i = 0; i < n_eq3; i++) 
    {
      int counter = 0;
      for (int j = 0; j < n_eq3; j++) 
      {
        if (i != j)
        {
          worse_ind.at(counter, i) = j;
          counter++;
        }
      }
    }
    
    // Construct the covariance matrix
    arma::vec sigma3_vec = par.elem(sigma3_ind);
    arma::mat sigma3(n_eq3 - 1, n_eq3 - 1);
    sigma3(0, 0) = 1;
    if (n_eq3 > 2)
    {
      int counter = 0;
      for (int i = 0; i < (n_eq3 - 1); i++)
      {
        for (int j = 0; j <= i; j++)
        {
          if (!((i == 0) & (j == 0)))
          {
            sigma3.at(i, j) = sigma3_vec.at(counter);
            sigma3.at(j, i) = sigma3.at(i, j);
            counter++;
          }
        }
      }
    }

    // Estimate the linear index for each alternative
    arma::mat li_mn(n_obs, n_eq3 - 1);
    for (int i = 0; i < (n_eq3 - 1); i++)
    {
      li_mn.col(i) = W_mn * coef3.row(i).t();
    }

    // Calculate the log-likelihood depending on the alternative
    for (int i = 0; i < n_groups; i++)
    {
      // Get the alternative
      int alt = groups3.at(i);
      if (alt == -1)
      {
        continue;
      }
      
      // Transformation matrix
      arma::mat transform_mat(n_eq3, n_eq3);
      transform_mat.diag().ones();
      
      if (alt != (n_eq3 - 1))
      {
        transform_mat.col(alt).ones();
        transform_mat.col(alt) = -transform_mat.col(alt);
      }
      transform_mat.shed_row(alt);
      transform_mat.shed_col(n_eq3 - 1);

      // Construct covariance matrix for the alternative
      arma::mat sigma3_transform = transform_mat * sigma3 * transform_mat.t();
      NumericMatrix sigma3_R     = wrap(sigma3_transform);
      
      // Check that covariance matrix is positive defined
      if(!sigma3.is_sympd())
      {
        if (is_diff)
        {
          NumericMatrix out_tmp(n_par, 1);
          std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
          return(out_tmp);
        }
        NumericMatrix out_tmp(1,1);
        out_tmp(0, 0) = -(1e+100);
        return(out_tmp);
      }
      
      // Vector of zero means
      NumericVector mean_zero_R(sigma3_R.ncol());
      
      // Estimate all necessary differences in utilities
      // related to linear indexes
      arma::mat li_diff = arma::mat(n_obs_g.at(i), n_eq3 - 1);
      if (alt == (n_eq3 - 1))
      {
        li_diff = -li_mn.rows(ind_g(i));
      }
      else
      {
        arma::uvec alt_uvec = {(unsigned int)alt};
        li_diff.col(n_eq3 - 2) = li_mn.submat(ind_g(i), alt_uvec);        
        int counter = 0;
        for (int j = 0; j < (n_eq3 - 1); j++)
        {
          arma::uvec j_uvec = {(unsigned int)j};
          if (alt != j)
          {
            li_diff.col(counter) = li_mn.submat(ind_g(i), alt_uvec) - 
              li_mn.submat(ind_g(i), j_uvec);
            counter++;
          }
        }
      }
      NumericMatrix li_diff_R = wrap(li_diff);
      
      // Matrix of negative infinite values
      NumericMatrix lower_neg_inf(n_obs_g.at(i), n_eq3 - 1);
      std::fill(lower_neg_inf.begin(), lower_neg_inf.end(), R_NegInf);
      
      // Calculate the probabilities
      List prob_list     = mnorm::pmnorm(lower_neg_inf, li_diff_R,
                                         NumericVector(), 
                                         mean_zero_R, sigma3_R,
                                         NumericVector(),
                                         n_sim, "default", "NO", true,
                                         is_diff, is_diff, is_diff, false,
                                         false, R_NilValue, n_cores,  
                                         R_NilValue, false, false);
      arma::vec prob_tmp = prob_list["prob"];
      
      // Calculate log-likelihoods
      lnL.elem(ind_g(i)) = prob_tmp;
      
      // Check the validity of the likelihood
      if(lnL.elem(ind_g(i)).has_inf() || lnL.elem(ind_g(i)).has_nan())
      {
        if (is_diff)
        {
          NumericMatrix out_tmp(n_par, 1);
          std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
          return(out_tmp);
        }
        NumericMatrix out_tmp(1,1);
        out_tmp(0, 0) = -(1e+100);
        return(out_tmp);
      }
      
      // Differentiate if need
      if (is_diff)
      {
        // Get some gradients
        arma::mat grad_upper  = prob_list["grad_upper"];
        arma::cube grad_sigma = prob_list["grad_sigma"];
        // Сoefficients
        // First (n_eq3 - 1) alternatives
        if (alt < (n_eq3 - 1))
        {
          for (int i1 = 0; i1 < (n_eq3 - 1); i1++)
          {
            for (int j1 = 0; j1 < n_coef3; j1++)
            {
              arma::uvec j1_uvec = {(unsigned int)j1};
              arma::uvec coef3_ind_uvec = {coef3_ind.at(alt, j1)};
              jac.submat(ind_g(i), coef3_ind_uvec) = 
                jac.submat(ind_g(i), coef3_ind_uvec) + 
                grad_upper.col(i1) % W_mn.submat(ind_g(i), j1_uvec);
            }
          }
          for (int i1 = 0; i1 < (n_eq3 - 2); i1++)
          {
            int worse_ind_i1 = worse_ind.at(i1, alt);
            for (int j1 = 0; j1 < n_coef3; j1++)
            {
              arma::uvec j1_uvec = {(unsigned int)j1};
              arma::uvec  coef3_ind_uvec = {coef3_ind.at(worse_ind_i1, j1)};
              jac.submat(ind_g(i), coef3_ind_uvec) = 
                jac.submat(ind_g(i), coef3_ind_uvec) -
                grad_upper.col(i1) % W_mn.submat(ind_g(i), j1_uvec);
            }
          }
        }
        // The last alternative
        else
        {
          for (int i1 = 0; i1 < (n_eq3 - 1); i1++)
          {
            for (int j1 = 0; j1 < n_coef3; j1++)
            {
              arma::uvec j1_uvec = {(unsigned int)j1};
              arma::uvec  coef3_ind_uvec = {coef3_ind.at(i1, j1)};
              jac.submat(ind_g(i), coef3_ind_uvec) = 
                jac.submat(ind_g(i), coef3_ind_uvec) -
                grad_upper.col(i1) % W_mn.submat(ind_g(i), j1_uvec);
            }
          }
        }
        // Сovariances
        int n_sigma3 = sigma3_vec.size();
        for (int i1 = 0; i1 < n_sigma3; i1++)
        {
          arma::uvec sigma3_ind_uvec_i1 = {sigma3_ind(i1)};
          arma::mat sigma3_alt_diff(n_eq3 - 1, n_eq3 - 1);
          sigma3_alt_diff.at(sigma3_ind_mat.at(i1, 0), 
                             sigma3_ind_mat.at(i1, 1)) = 1;
          sigma3_alt_diff.at(sigma3_ind_mat.at(i1, 1), 
                             sigma3_ind_mat.at(i1, 0)) = 1;
          sigma3_alt_diff = transform_mat * sigma3_alt_diff * transform_mat.t();
          for (int j1 = 0; j1 < (n_eq3 - 1); j1++)
          {
            for (int j2 = 0; j2 <= j1; j2++)
            {
              if (sigma3_alt_diff.at(j1, j2) != 0)
              {
                arma::vec mat_tmp = grad_sigma.tube(j1, j2);
                jac.submat(ind_g(i), sigma3_ind_uvec_i1) =
                  jac.submat(ind_g(i), sigma3_ind_uvec_i1) +
                  mat_tmp * sigma3_alt_diff.at(j1, j2);
              }
            }
          }
        }
      }
    }
  }

  // Routine for the ordered equations
  for (int i = 0; i < n_groups; i++)
  {
    // Leave the cycle if there is only multinomial equation
    if (ind_eq_all(i).size() == 0)
    {
      continue;
    }

    // Transform index
    arma::uvec i_uvec = {(unsigned int)i};
    
    // Modify covariance matrix addressing for the
    // continuous equations and their regimes
    if (is2)
    {
      for (int j1 = 0; j1 < n_eq2_g.at(i); j1++)
      {
        int j1_o = ind_eq2(i).at(j1);
        sigma.at(n_eq + j1_o, n_eq + j1_o) = par.at(var2_ind(j1_o).at(groups2.at(i, j1_o)));
        for (int j2 = 0; j2 < n_eq; j2++)
        {
          sigma.at(j1_o + n_eq, j2) = par.at(cov2_ind(j1_o).at(groups2.at(i, j1_o), j2));
          sigma.at(j2, j1_o + n_eq) = par.at(cov2_ind(j1_o).at(groups2.at(i, j1_o), j2));
        }
      }
      if (n_eq2 >= 2)
      {
        for (int j1 = 1; j1 < n_eq2_g.at(i); j1++)
        {
          for (int j2 = 0; j2 < j1; j2++)
          {
            int j1_o                           = ind_eq2(i).at(j1);
            int j2_o                           = ind_eq2(i).at(j2);
            arma::uvec uvec_tmp                = {sigma2_ind_mat(i).at(j1_o, 
                                                                       j2_o)};
            arma::vec vec_tmp                  =  par.elem(uvec_tmp);
            sigma.at(j1_o + n_eq, j2_o + n_eq) = vec_tmp.at(0);
            sigma.at(j2_o + n_eq, j1_o + n_eq) = vec_tmp.at(0);
          }
        }
      }
    }

    // Construct covariance matrix for the observable equations
    arma::mat sigma_g       = sigma.submat(ind_eq_all(i), ind_eq_all(i));
    NumericMatrix sigma_g_R = wrap(sigma_g);

    // Check that matrix is positive defined after accounting for 
    // the elements associated with the continuous equations
    if(!sigma_g.is_sympd())
    {
      if (is_diff)
      {
        NumericMatrix out_tmp(n_par, 1);
        std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
        return(out_tmp);
      }
      NumericMatrix out_tmp(1,1);
      out_tmp(0, 0) = -(1e+100);
      return(out_tmp);
    }
    

    // Vector of zero means
    NumericVector mean_zero_R(n_eq_all_g.at(i));

    // Get lower and upper integration limits
    arma::mat li_lower_tmp   = li_lower.submat(ind_g(i), ind_eq(i));
    arma::mat li_upper_tmp   = li_upper.submat(ind_g(i), ind_eq(i));
    NumericMatrix li_lower_R = wrap(li_lower_tmp);
    NumericMatrix li_upper_R = wrap(li_upper_tmp);
    arma::mat li_y_tmp;
    NumericMatrix li_y_R;
    if (n_eq2_g(i) > 0)
    {
      li_y_tmp = li_y.submat(ind_g(i), ind_eq2(i));
      li_y_R   = wrap(li_y_tmp);
    }

    // Create a list of marginal distributions parameters
    // accounting for the unobservable equations
    Nullable<List> marginal_par_list_i = R_NilValue;
    if (is_marginal)
    {
      int n_list = n_eq_g.at(i) + n_eq2_g.at(i);
      List marginal_par_list_tmp(n_list);
      CharacterVector marginal_names_tmp(n_list);
      for (int j = 0; j < n_eq_g.at(i); j++)
      {
        int j_o = ind_eq(i).at(j);
        marginal_par_list_tmp[j] = marginal_par_list[j_o];
        marginal_names_tmp[j]    = marginal_names[j_o];
      }
      if (is2)
      {
        for (int j = 0; j < n_eq2_g.at(i); j++)
        {
          int j_adj                    = j + n_eq_g.at(i);
          marginal_par_list_tmp[j_adj] = R_NilValue;
          marginal_names_tmp[j_adj]    = "normal";
        }
      }
      marginal_par_list_tmp.names() = marginal_names_tmp;
      Nullable<List> marginal_par_list_tmp2(marginal_par_list_tmp);
      marginal_par_list_i           = marginal_par_list_tmp2;
    }

    // Calculate the probabilities and differentiate respect to them
    List prob_list;
    if (n_eq_g(i) > 0)
    {
      if (is2 & (n_eq2_g(i) > 0))
      {
        IntegerVector given_x_tmp = Rcpp::seq(n_eq_g(i) + 1, n_eq_all_g(i));
        NumericVector given_x     = as<NumericVector>(given_x_tmp);
        prob_list                 = mnorm::pmnorm(li_lower_R, li_upper_R,
                                                  li_y_R,
                                                  mean_zero_R, sigma_g_R,
                                                  given_x,
                                                  n_sim, "default", "NO", true,
                                                  is_diff, is_diff, 
                                                  is_diff, is_diff,
                                                  false, R_NilValue, n_cores,
                                                  marginal_par_list_i, 
                                                  is_diff, is_marginal);
      }
      else
      {
        prob_list = mnorm::pmnorm(li_lower_R, li_upper_R,
                                  NumericVector(),
                                  mean_zero_R, sigma_g_R,
                                  NumericVector(),
                                  n_sim, "default", "NO", true,
                                  is_diff, is_diff, 
                                  is_diff, false,
                                  false, R_NilValue, n_cores,
                                  marginal_par_list_i, 
                                  is_diff, is_marginal);
      }
    }
    else
    {
      arma::cube grad_sigma_tmp = arma::cube();
      arma::mat grad_given_tmp  = arma::mat();
      prob_list["grad_sigma"]   = grad_sigma_tmp;
      prob_list["grad_given"]   = grad_given_tmp;
      
    }

    // Calculate the densities
    List den_list;
    if (is2)
    {
      if (n_eq2_g(i) > 0)
      {
        NumericVector mean2_zero_R(n_eq2_g.at(i));
        arma::mat sigma2_g       = sigma.submat(ind_eq2(i) + n_eq,
                                                ind_eq2(i) + n_eq);
        NumericMatrix sigma2_g_R = wrap(sigma2_g);
        den_list                 = mnorm::dmnorm(li_y_R,
                                                 mean2_zero_R, sigma2_g_R,
                                                 NumericVector(), 
                                                 true, is_diff, is_diff,
                                                 false,
                                                 R_NilValue, n_cores);
      }
    }

    // Store the log-likelihood contributions
    if (out_type == "val")
    {
      if (is1)
      {
        if (n_eq_g(i) > 0)
        {
          NumericVector prob = prob_list["prob"];
          arma::vec prob_arma(prob.begin(), n_obs_g.at(i), false);
          lnL.elem(ind_g(i)) = lnL.elem(ind_g(i)) + prob_arma;
        }
      }
      if (is2)
      {
        if (n_eq2_g(i) > 0)
        {
          NumericVector den  = den_list["den"];
          arma::vec den_arma(den.begin(), n_obs_g.at(i), false);
          lnL.elem(ind_g(i)) = lnL.elem(ind_g(i)) + den_arma;
        }
      }

      // Check the validity of the likelihood
      if(lnL.elem(ind_g(i)).has_inf() || lnL.elem(ind_g(i)).has_nan())
      {
        if (is_diff)
        {
          NumericMatrix out_tmp(n_par, 1);
          std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
          return(out_tmp);
        }
        NumericMatrix out_tmp(1,1);
        out_tmp(0, 0) = -(1e+100);
        return(out_tmp);
      }
    }

    // Store the Jacobian if need
    if (is_diff)
    {
      arma::cube grad_sigma   = prob_list["grad_sigma"];
      if (is1 & (n_eq_g.at(i) > 0))
      {
        arma::mat grad_lower = prob_list["grad_lower"];
        arma::mat grad_upper = prob_list["grad_upper"];
        arma::mat grad_sum   = grad_upper + grad_lower;
        // for the mean equation coefficients
        for (int j = 0; j < n_eq_g.at(i); j++)
        {
          int j_o             = ind_eq(i).at(j);
          arma::uvec j_o_uvec = {(unsigned int)j_o};
          arma::mat W_tmp     = W(j_o).rows(ind_g(i));
          arma::mat jac_tmp   = W_tmp.each_col() % (-grad_sum.col(j));
          if (is_het[j_o])
          {
            jac.submat(ind_g(i), coef_ind(j_o)) = jac_tmp.each_col() /
                                                  li_var.submat(ind_g(i), j_o_uvec);
          }
          else
          {
            jac.submat(ind_g(i), coef_ind(j_o)) = jac_tmp;
          }
        }
        // for the covariances
        arma::cube grad_sigma      = prob_list["grad_sigma"];
        arma::umat sigma_ind_mat_o = sigma_ind_mat.submat(ind_eq(i), ind_eq(i));
        for (int j1 = 1; j1 < n_eq_g.at(i); j1++)
        {
          for (int j2 = 0; j2 < j1; j2++)
          {
            arma::vec grad_sigma_tmp             = grad_sigma.tube(j1, j2);
            arma::uvec sigma_ind_uvec            = {sigma_ind_mat_o.at(j1, j2)};
            jac.submat(ind_g(i), sigma_ind_uvec) = grad_sigma_tmp;
          }
        }
        // for the cuts
        for (int j = 0; j < n_eq_g.at(i); j++)
        {
          int j_o             = ind_eq(i).at(j);
          arma::uvec j_o_uvec = {(unsigned int)j_o};
          int group_j         = groups(i, j_o);
          // lower cuts
          if (group_j != 0)
          {
            arma::uvec j_tmp = {cuts_ind(j_o).at(group_j - 1)};
            if (is_het[j_o])
            {
              jac.submat(ind_g(i), j_tmp) = grad_lower.col(j) / 
                                            li_var.submat(ind_g(i), j_o_uvec);
            }
            else
            {
              jac.submat(ind_g(i), j_tmp) = grad_lower.col(j);
            }
          }
          // upper cuts
          if (group_j != n_cuts_eq.at(j_o))
          {
            arma::uvec j_tmp = {cuts_ind(j_o).at(group_j)};
            if (is_het[j_o])
            {
              jac.submat(ind_g(i), j_tmp) = grad_upper.col(j) / 
                                            li_var.submat(ind_g(i), j_o_uvec);
            }
            else
            {
              jac.submat(ind_g(i), j_tmp) = grad_upper.col(j);
            }
          }
        }
        // for the variance equation coefficients
        for (int j = 0; j < n_eq_g.at(i); j++)
        {
          int j_o             = ind_eq(i).at(j);
          arma::uvec j_o_uvec = {(unsigned int)j_o};
          if (is_het[j_o])
          {
            arma::vec li_lower_j = li_lower.submat(ind_g(i), j_o_uvec);
            arma::vec li_upper_j = li_upper.submat(ind_g(i), j_o_uvec);
            li_lower_j.replace(-arma::datum::inf, 0);
            li_upper_j.replace(arma::datum::inf, 0);
            arma::mat W_var_tmp = W_var(j_o).rows(ind_g(i));
            jac.submat(ind_g(i), coef_var_ind(j_o)) = W_var_tmp.each_col() % 
                                                      (-((grad_lower.col(j) % li_lower_j) +
                                                      (grad_upper.col(j) % li_upper_j)));
          }
        }
      }
      // Jacobian for the parameters related to the continuous equations
      if (is2 & (n_eq2_g.at(i) > 0))
      {
        arma::mat grad_given      = prob_list["grad_given"];
        arma::mat grad_x          = den_list["grad_x"];
        arma::cube grad_sigma_den = den_list["grad_sigma"];
        // for the coefficients, variances and covariances of the continuous
        // equations with the ordered equations
        for (int j = 0; j < n_eq2_g.at(i); j++)
        {
          int j_o             = ind_eq2(i).at(j);
          arma::uvec j_o_uvec = {(unsigned int)j_o};
          arma::mat X_tmp     = X(j_o).rows(ind_g(i));
          // part related to the probability
          arma::vec jac_tmp_1_var = arma::vec(n_obs_g.at(i));
          if (n_eq_g.at(i) > 0)
          {
            jac_tmp_1_var = grad_sigma.tube(n_eq_g(i) + j, n_eq_g(i) + j);
            for (int j1 = 0; j1 < n_eq_g.at(i); j1++)
            {
              int j1_o                 = ind_eq(i).at(j1);
              arma::uvec cov2_ind_uvec = {cov2_ind(j_o).at(groups2.at(i, j_o), 
                                                                      j1_o)};
              arma::vec tube_tmp       = grad_sigma.tube(n_eq_g(i) + j, j1);
              jac.submat(ind_g(i), cov2_ind_uvec) = tube_tmp;
            }
          }
          // part related to density
          arma::mat jac_tmp_2     = X_tmp.each_col() % (-grad_x.col(j));
          arma::vec jac_tmp_2_var = grad_sigma_den.tube(j, j);
          // combine both parts
          arma::uvec var2_ind_uvec = {var2_ind.at(j_o).at(groups2(i, j_o))};
          if (n_eq_g.at(i) > 0)
          {
            jac_tmp_2 = jac_tmp_2 + X_tmp.each_col() % (-grad_given.col(j));
          }
          jac.submat(ind_g(i), coef2_ind(j_o).row(groups2(i, j_o))) = jac_tmp_2;
          jac.submat(ind_g(i), var2_ind_uvec) = jac_tmp_1_var + jac_tmp_2_var;
        }
        // for the covariances between the continuous equations
        if (n_eq2_g.at(i) > 1)
        {
          arma::umat sigma2_ind_mat_o = sigma2_ind_mat(i).submat(ind_eq2(i), 
                                                                 ind_eq2(i));
          for (int j1 = 1; j1 < n_eq2_g.at(i); j1++)
          {
            for (int j2 = 0; j2 < j1; j2++)
            {
              arma::vec grad_sigma2_den = grad_sigma_den.tube(j1, j2);
              if (n_eq_g.at(i) > 0)
              {
                arma::vec grad_sigma2_prob = grad_sigma.tube(n_eq_g(i) + j1,
                                                             n_eq_g(i) + j2);
                grad_sigma2_den            = grad_sigma2_den + grad_sigma2_prob;
                                                                     
              }
              arma::uvec sigma2_ind_uvec            = {sigma2_ind_mat_o.at(j1, j2)};
              jac.submat(ind_g(i), sigma2_ind_uvec) = grad_sigma2_den;
            }
          }
        }
      }
      // Jacobian for the parameters related to the marginal 
      // distributions parameters
      if (is_marginal)
      {
        arma::field<arma::mat> grad_marginal = prob_list["grad_marginal"];
        for (int j = 0; j < n_eq_g.at(i); j++)
        {
          int j_o = ind_eq(i).at(j);
          if (marginal_par_n.at(j_o) > 0)
          {
            jac.submat(ind_g(i), marginal_par_ind.at(j_o)) = grad_marginal.at(j);
          }
        }
      }
    }
  }
  
  // Return the Jacobian if need
  if (is_jac)
  {
    return(wrap(jac));
  }

  // Return the gradient if need
  if (is_grad)
  {
    arma::rowvec grad = sum(jac, 0);
    // Account for the regularization if need
    if (is_regularization)
    {
      // Ridge
      if (is_ridge)
      {
        arma::vec par_ridge  = par.elem(ridge_ind) - ridge_location;
        grad.elem(ridge_ind) = grad.elem(ridge_ind) - 
                               2 * par_ridge % ridge_scale;
      }
      // Lasso
      if (is_lasso)
      {
      arma::vec par_lasso  = par.elem(lasso_ind) - lasso_location;
      grad.elem(lasso_ind) = grad.elem(lasso_ind) - 
                             arma::sign(par_lasso) % lasso_scale;
      }
    }
    return(wrap(grad));
  }

  // Calculate the log-likelihood
  double lnL_val = sum(lnL);
  
  // Perform the regularization if need
  if (is_regularization)
  {
    if (is_ridge)
    {
      arma::vec par_ridge = par.elem(ridge_ind) - ridge_location;
      lnL_val            -= sum(pow(par_ridge, 2) % ridge_scale);
    }
    if (is_lasso)
    {
      arma::vec par_lasso = par.elem(lasso_ind) - lasso_location;
      lnL_val            -= sum(arma::abs(par_lasso) % lasso_scale);
    }
  }
  
  // Check the validity of the likelihood
  NumericMatrix lnL_mat(1, 1);
  lnL_mat(0, 0) = lnL_val;

  return(lnL_mat);
}

//' Gradient of the Log-likelihood Function of Multivariate Ordered
//' Probit Model
//' @description Calculates gradient of the log-likelihood function of 
//' multivariate ordered probit model.
//' @param par vector of parameters.
//' @param control_lnL list with some additional parameters.
//' @param out_type string representing the output type of the function.
//' @param n_sim the number of random draws for multivariate 
//' normal probabilities.
//' @param n_cores the number of cores to be used. 
//' @param regularization list of regularization parameters.
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix grad_msel(const arma::vec par,
                             const List control_lnL,
                             const String out_type = "grad",
                             const int n_sim = 1000,
                             const int n_cores = 1,
                             const Nullable<List> regularization = R_NilValue)
{
  return(lnL_msel(par, control_lnL, "grad", 
                       n_sim, n_cores, regularization));
}
