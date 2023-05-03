// [[Rcpp::depends(mnorm)]]
#include <RcppArmadillo.h>
#include <mnorm.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

//' Log-likelihood Function of Multinomial Probit Model
//' @description Calculates log-likelihood function of multinomial probit model.
//' @param par vector of parameters.
//' @param control_lnL list with some additional parameters.
//' @param out_type string represeint the output type of the function.
//' @param n_sim the number of random draws for multivariate 
//' normal probabilities.
//' @param n_cores the number of cores to be used. 
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix lnL_mnprobit(const arma::vec par,
                           const List control_lnL,
                           const String out_type = "val",
                           const int n_sim = 1000,
                           const int n_cores = 1)
{
  // Determine whether gradient and Jacobian should be estimated
  const bool is_grad = (out_type == "grad");
  const bool is_jac = (out_type == "jac");
  const bool is_diff = is_grad | is_jac;
  
  // Get some dimensional data
  const int n_par = control_lnL["n_par"];
  const int n_alt = control_lnL["n_alt"];
  const int n_obs = control_lnL["n_obs"];
  const int n_coef = control_lnL["n_coef"];
  const int n_coef2 = control_lnL["n_coef2"];
  const int n_coef_total = control_lnL["n_coef_total"];
  const int n_sigma = control_lnL["n_sigma"];
  
  // Get data
  const arma::vec z = control_lnL["z"];
  const arma::mat W = control_lnL["W"];
  const arma::vec y = control_lnL["y"];
  const arma::mat X = control_lnL["X"];

  // Check whether there are continuous equations
  const int n_regimes = control_lnL["n_regimes"];
  const bool is2 = n_regimes > 0;
  const arma::vec regimes = control_lnL["regimes"];

  // Get indexes of observations associated with different
  // alternatives and regimes
  const arma::field<arma::uvec> ind_alt = control_lnL["ind_alt"];
  const arma::field<arma::uvec> ind_regime = control_lnL["ind_regime"];

  // Determine indexes of coefficients and 
  // covariances in the vector of parameters
  const arma::umat coef_ind_alt = control_lnL["coef_ind_alt"];
  const arma::uvec sigma_ind = control_lnL["sigma_ind"];
  const arma::umat sigma_ind_mat = control_lnL("sigma_ind_mat");
  const arma::umat coef2_ind_regime = control_lnL["coef2_ind_regime"];
  const arma::umat cov2_ind_regime = control_lnL["cov2_ind_regime"];
  const arma::uvec var2_ind_regime = control_lnL["var2_ind_regime"];
  
  // Initialize Jacobian matrix
  arma::mat jac;
  if (is_diff)
  {
    jac = arma::mat(n_obs, n_par);
  }

  // Store indexes of worse alternatives for each best alternative
  arma::umat worse_ind(n_alt - 1, n_alt);
  for (int i = 0; i < n_alt; i++) 
  {
    int counter = 0;
    for (int j = 0; j < n_alt; j++) 
    {
      if (i != j)
      {
        worse_ind.at(counter, i) = j;
        counter++;
      }
    }
  }

  // Store the coefficients for each alternative
  arma::mat coef(n_coef, n_alt - 1);
  for (int i = 0; i < (n_alt - 1); i++) 
  {
    coef.col(i) = par.elem(coef_ind_alt.col(i));
  }

  // Store the coefficients for each regime
  arma::mat coef2;
  if (is2)
  {
    coef2 = arma::mat(n_coef2, n_regimes);
    for (int i = 0; i < n_regimes; i++)
    {
      coef2.col(i) = par.elem(coef2_ind_regime.col(i));
    }
  }
  
  // Construct the covariance matrix
  arma::vec sigma_vec = par.elem(sigma_ind);
  arma::mat sigma(n_alt - 1, n_alt - 1);
  sigma(0, 0) = 1;
  if (n_alt > 2)
  {
    int counter = 0;
    for (int i = 0; i < (n_alt - 1); i++)
    {
      for (int j = 0; j <= i; j++)
      {
        if (!((i == 0) & (j == 0)))
        {
          sigma.at(i, j) = sigma_vec.at(counter);
          sigma.at(j, i) = sigma.at(i, j);
          counter++;
        }
      }
    }
  }

  // Check that matrix is positive defined
  if(!sigma.is_sympd())
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
  NumericVector mean_zero_R(n_alt - 1);
  
  // Estimate the linear index for each alternative
  arma::mat li(n_obs, n_alt - 1);
  for (int i = 0; i < (n_alt - 1); i++)
  {
    li.col(i) = W * coef.col(i);
  }
  
  // Estimate the linear index for each regime
  // and corresponding densities
  arma::vec li_y;
  arma::vec den;
  arma::vec var2;
  if (is2)
  {
    li_y = arma::vec(n_obs);
    den = arma::vec(n_obs);
    var2 = arma::vec(n_regimes);
    for (int i = 0; i < n_regimes; i++)
    {
      arma::vec li_y_tmp = y.elem(ind_regime(i)) - 
                           X.rows(ind_regime(i)) * coef2.col(i);
      li_y.elem(ind_regime(i)) = li_y_tmp;
      var2.at(i) = par.at(var2_ind_regime.at(i));
      den.elem(ind_regime(i)) = arma::log_normpdf(li_y_tmp, 
                                                  0.0, sqrt(var2.at(i)));
    }
  }

  // Vector to store log-likelihood information
  arma::vec lnL(n_obs);

  // Calculate log-likelihood depending on the alternative
  for (int i = 0; i < n_alt; i++)
  {
    // Get indexes of observations related
    // to the alternative
    int n_alt_obs = ind_alt(i).size();
    
    // Transformation matrix
    arma::mat transform_mat(n_alt, n_alt);
    transform_mat.diag().ones();
    
    if (i != (n_alt - 1))
    {
      transform_mat.col(i).ones();
      transform_mat.col(i) = -transform_mat.col(i);
    }
    transform_mat.shed_row(i);
    transform_mat.shed_col(n_alt - 1);
    
    // Account for continuous equations
    const int regime = regimes.at(i);
    arma::mat sigma2 = arma::mat(n_alt, n_alt);
    arma::mat transform_mat2;
    if (regime != -1)
    {
      // extend covariance matrix
      sigma2.submat(0, 0, n_alt - 2, n_alt - 2) = sigma;
      sigma2.at(n_alt - 1, n_alt - 1) = par.at(var2_ind_regime.at(regime));
      sigma2.submat(n_alt - 1, 0, n_alt - 1, n_alt - 2) = par.elem(
        cov2_ind_regime.col(regime)).t();
      sigma2.submat(0, n_alt - 1, n_alt - 2, n_alt - 1) = par.elem(
        cov2_ind_regime.col(regime));
      // extend transformation matrix
      transform_mat2 = arma::mat(n_alt, n_alt);
      transform_mat2.submat(0, 0, n_alt - 2, n_alt - 2) = transform_mat;
      transform_mat2.at(n_alt - 1, n_alt - 1) = 1;
    }
    else
    {
      sigma2 = sigma;
      transform_mat2 = transform_mat;
    }

    // Check that covariance matrix is positive
    if (regime != -1)
    {
      if(!sigma2.is_sympd())
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

    // Construct covariance matrix for alternative
    arma::mat sigma_alt = transform_mat2 * sigma2 * transform_mat2.t();
    NumericMatrix sigma_alt_R = wrap(sigma_alt);

    // Estimate all necessary differences in utilities
    // related to linear indexes
    arma::mat li_diff = arma::mat(n_alt_obs, n_alt - 1);
    if (i == (n_alt - 1))
    {
      li_diff = -li.rows(ind_alt(i));
    }
    else
    {
      arma::uvec i_uvec = {(unsigned int)i};
      li_diff.col(n_alt - 2) = li.submat(ind_alt(i), i_uvec);        
      int counter = 0;
      for (int j = 0; j < (n_alt - 1); j++)
      {
        arma::uvec j_uvec = {(unsigned int)j};
        if (i != j)
        {
          li_diff.col(counter) = li.submat(ind_alt(i), i_uvec) - 
            li.submat(ind_alt(i), j_uvec);
          counter++;
        }
      }
    }
    NumericMatrix li_diff_R = wrap(li_diff);
    
    // Matrix of negative infinite values
    NumericMatrix lower_neg_inf(n_alt_obs, n_alt - 1);
    std::fill(lower_neg_inf.begin(), lower_neg_inf.end(), R_NegInf);
    
    // Values of continuous variable
    arma::mat li_y_tmp;
    NumericVector li_y_R;
    if (regime != -1)
    {
      li_y_tmp = li_y.elem(ind_alt(i));
      li_y_R = wrap(li_y_tmp);
    } 

    // Calculate probabilities
    NumericVector given_x;
    if (regime != -1)
    {
      given_x = NumericVector::create(n_alt);
    }
    List prob_list = mnorm::pmnorm(lower_neg_inf, li_diff_R,
                                   li_y_R, 
                                   mean_zero_R, sigma_alt_R,
                                   given_x,
                                   n_sim, "default", "NO", true,
                                   is_diff, is_diff, is_diff, is_diff & is2,
                                   false, R_NilValue, n_cores,
                                   R_NilValue);
    arma::vec prob_tmp = prob_list["prob"];

    // Calculate log-likelihoods
    lnL.elem(ind_alt(i)) = prob_tmp;
    if (regime != -1)
    {
      lnL.elem(ind_alt(i)) = lnL.elem(ind_alt(i)) + den(ind_alt(i));
    }

    // Differentiate if need
    if (is_diff)
    {
      // Get some gradients
      arma::mat grad_upper = prob_list["grad_upper"];
      arma::cube grad_sigma = prob_list["grad_sigma"];

      // Сoefficients
      // First (n_alt - 1) alternatives
      if (i < (n_alt - 1))
      {
        for (int i1 = 0; i1 < (n_alt - 1); i1++)
        {
          for (int j1 = 0; j1 < n_coef; j1++)
          {
            arma::uvec j1_uvec = {(unsigned int)j1};
            arma::uvec  coef_ind_alt_uvec = {coef_ind_alt.at(j1, i)};
            jac.submat(ind_alt(i), coef_ind_alt_uvec) = 
              jac.submat(ind_alt(i), coef_ind_alt_uvec) + 
              grad_upper.col(i1) % W.submat(ind_alt(i), j1_uvec);
          }
        }
        for (int i1 = 0; i1 < (n_alt - 2); i1++)
        {
          int worse_ind_i1 = worse_ind.at(i1, i);
          for (int j1 = 0; j1 < n_coef; j1++)
          {
            arma::uvec j1_uvec = {(unsigned int)j1};
            arma::uvec  coef_ind_alt_uvec = {coef_ind_alt.at(j1, worse_ind_i1)};
            jac.submat(ind_alt(i), coef_ind_alt_uvec) = 
              jac.submat(ind_alt(i), coef_ind_alt_uvec) -
              grad_upper.col(i1) % W.submat(ind_alt(i), j1_uvec);
          }
        }
      }
      // Last alternative
      else
      {
        for (int i1 = 0; i1 < (n_alt - 1); i1++)
        {
          for (int j1 = 0; j1 < n_coef; j1++)
          {
            arma::uvec j1_uvec = {(unsigned int)j1};
            arma::uvec  coef_ind_alt_uvec = {coef_ind_alt.at(j1, i1)};
            jac.submat(ind_alt(i), coef_ind_alt_uvec) = 
              jac.submat(ind_alt(i), coef_ind_alt_uvec) -
              grad_upper.col(i1) % W.submat(ind_alt(i), j1_uvec);
          }
        }
      }

      // Сovariances
      for (int i1 = 0; i1 < n_sigma; i1++)
      {
        arma::uvec sigma_ind_uvec_i1 = {sigma_ind(i1)};
        arma::mat sigma_alt_diff(n_alt - 1, n_alt - 1);
        sigma_alt_diff.at(sigma_ind_mat.at(i1, 0), sigma_ind_mat.at(i1, 1)) = 1;
        sigma_alt_diff.at(sigma_ind_mat.at(i1, 1), sigma_ind_mat.at(i1, 0)) = 1;
        sigma_alt_diff = transform_mat * sigma_alt_diff * transform_mat.t();
        for (int j1 = 0; j1 < (n_alt - 1); j1++)
        {
          arma::uvec j1_uvec = {(unsigned int)j1};
          for (int j2 = 0; j2 <= j1; j2++)
          {
            arma::uvec j2_uvec = {(unsigned int)j2};
            if (sigma_alt_diff.at(j1, j2) != 0)
            {
              arma::vec mat_tmp = grad_sigma.tube(j1, j2);
              jac.submat(ind_alt(i), sigma_ind_uvec_i1) =
                jac.submat(ind_alt(i), sigma_ind_uvec_i1) +
                mat_tmp * sigma_alt_diff.at(j1, j2);
            }
          }
        }
      }
      
      // Part related to continuous equation
      if (regime != -1)
      {
        // Coefficients
        arma::colvec grad_given = prob_list["grad_given"];
        arma::mat X_tmp = X.rows(ind_alt(i));  
        arma::uvec coef2_ind_regime_uvec = {coef2_ind_regime.col(regime)};
        arma::vec li_y_adj = li_y.elem(ind_alt(i)) / var2.at(regime);
        jac.submat(ind_alt(i), coef2_ind_regime_uvec) = 
          X_tmp.each_col() % (li_y_adj - grad_given);
        
        // Covariances
        arma::cube grad_sigma = prob_list["grad_sigma"];
        for (int j = 0; j < (n_alt - 1); j++)
        {
          arma::uvec uvec_tmp = {cov2_ind_regime.at(j, regime)};
          for (int j1 = 0; j1 < (n_alt - 1); j1++)
          {
            if (transform_mat.at(j1, j) != 0)
            {
              arma::vec grad_sigma_tmp = grad_sigma.tube(n_alt - 1, j1);
              jac.submat(ind_alt(i), uvec_tmp) = 
                jac.submat(ind_alt(i), uvec_tmp) + 
                transform_mat.at(j1, j) * 
                grad_sigma_tmp;
            }
          }
        }
        
        // Variance
        arma::vec grad_sigma_tmp = grad_sigma.tube(n_alt - 1, n_alt - 1);
        arma::uvec uvec_tmp = {var2_ind_regime.at(regime)};
        jac.submat(ind_alt(i), uvec_tmp) =  grad_sigma_tmp +
                                            (arma::pow(li_y_adj, 2) - 
                                             1 / var2.at(regime)) / 2;
                                            //(arma::pow(li_y.elem(ind_alt(i)), 2) -
                                            //var2.at(regime)) / 
                                            //(2 * pow(var2.at(regime) , 2));
      }
    }
  }
  
  // Return jacobian if need
  if (is_jac)
  {
    return(wrap(jac));
  }

  // Return gradient if need
  if (is_grad)
  {
    arma::rowvec grad = sum(jac, 0);
    return(wrap(grad));
  }

  // By default return value of the log-likelihood function
  NumericMatrix out(1,1);
  out(0, 0) = sum(lnL);
  return(out);
}

//' Gradient of the Log-likelihood Function of Multinomial Probit Model
//' @description Calculates gradient of the log-likelihood function of 
//' multinomial probit model.
//' @param par vector of parameters.
//' @param control_lnL list with some additional parameters.
//' @param out_type string represeint the output type of the function.
//' @param n_sim the number of random draws for multivariate 
//' normal probabilities.
//' @param n_cores the number of cores to be used. 
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix grad_mnprobit(const arma::vec par,
                            const List control_lnL,
                            const String out_type = "grad",
                            const int n_sim = 1000,
                            const int n_cores = 1)
{
  return(lnL_mnprobit(par, control_lnL, "grad", n_sim, n_cores));
}
