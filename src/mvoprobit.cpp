// [[Rcpp::depends(mnorm)]]
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <mnorm.h>
#include "helpFunctions.h"
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

//' Log-likelihood Function of Multivariate Ordered Probit Model
//' @description Calculates log-likelihood function of multivariate ordered
//' probit model.
//' @param par vector of parameters.
//' @param control_lnL list with some additional parameters.
//' @param out_type string represeint the output type of the function.
//' @param n_sim the number of random draws for multivariate 
//' normal probabilities.
//' @param n_cores the number of cores to be used. 
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix lnL_mvoprobit(const arma::vec par,
                            const List control_lnL,
                            const String out_type = "val",
                            const int n_sim = 1000,
                            const int n_cores = 1)
{
  // Determine whether gradient and Jacobian should be estimated
  const bool is_grad = (out_type == "grad");
  const bool is_jac = (out_type == "jac");
  const bool is_diff = is_grad | is_jac;
  
  // Get regressors for ordered equations
  const arma::field<arma::mat> W = control_lnL["W"];
  const arma::field<arma::mat> W_var = control_lnL["W_var"];

  // Get some dimensional data
  const int n_par = control_lnL["n_par"];
  const int n_obs = control_lnL["n_obs"];
  const arma::vec n_regimes = control_lnL["n_regimes"];
  const arma::vec n_coef = control_lnL["n_coef"];
  const arma::vec n_coef2 = control_lnL["n_coef2"];
  const int n_eq = control_lnL["n_eq"];
  const int n_eq2 = control_lnL["n_eq2"];
  const int n_groups = control_lnL["n_groups"];
  const arma::mat groups = control_lnL["groups"];
  const arma::mat groups2 = control_lnL["groups2"];
  const arma::vec n_eq_g = control_lnL["n_eq_g"];
  const arma::vec n_eq2_g = control_lnL["n_eq2_g"];
  const arma::vec n_eq_all_g = control_lnL["n_eq_all_g"];
  const arma::vec n_cuts_eq = control_lnL["n_cuts_eq"];
  const bool is2 = control_lnL["is2"];
  const arma::mat y = control_lnL["y"];
  const arma::field<arma::mat> X = control_lnL["X"];
  
  // Create matrices to store Jacobian if need
  arma::mat jac;
  if (is_jac | is_grad)
  {
    jac = arma::mat(n_obs, n_par);
  }

  // Get indexes of observations for different groups
  const arma::field<arma::uvec> ind_g = control_lnL["ind_g"];
  const arma::vec n_obs_g = control_lnL["n_obs_g"];
  
  // Get indexes of observed equations for each group
  const arma::field<arma::uvec> ind_eq = control_lnL["ind_eq"];
  const arma::field<arma::uvec> ind_eq2 = control_lnL["ind_eq2"];
  const arma::field<arma::uvec> ind_eq_all = control_lnL["ind_eq_all"];

  // Determine indexes of coefficients and 
  // covariances in the vector of parameters
  const arma::field<arma::uvec> coef_ind = control_lnL["coef_ind"];
  const arma::field<arma::umat> coef2_ind = control_lnL["coef2_ind"];
  const arma::uvec sigma_ind = control_lnL["sigma_ind"];
  const arma::field<arma::uvec> sigma2_ind = control_lnL["sigma2_ind"];
  const arma::umat sigma_ind_mat = control_lnL["sigma_ind_mat"];
  const arma::field<arma::umat> sigma2_ind_mat = control_lnL["sigma2_ind_mat"];
  const arma::field<arma::uvec> var_y_ind = control_lnL["var_y_ind"];
  const arma::field<arma::umat> cov_y_ind = control_lnL["cov_y_ind"];
  const arma::field<arma::uvec> cuts_ind = control_lnL["cuts_ind"];

  // Store the coefficients for each equation
  arma::field<arma::vec> coef(n_eq);
  for (int i = 0; i < n_eq; i++)
  {
    coef.at(i) = par.elem(coef_ind(i));
  }

  // Construct the covariance matrix for ordered equations
  arma::vec sigma_vec = par.elem(sigma_ind);
  arma::mat sigma(n_eq + n_eq2, n_eq + n_eq2);
  sigma.diag().ones();
  if (n_eq > 1)
  {
    int counter = 0;
    for (int i = 1; i < n_eq; i++)
    {
      for (int j = 0; j < i; j++)
      {
        sigma.at(i, j) = sigma_vec.at(counter);
        sigma.at(j, i) = sigma.at(i, j);
        counter++;
      }
    }
  }
  NumericMatrix sigma_R = wrap(sigma);

  // Check that matrix is positive defined
  if(!sigma.is_sympd() & !is2)
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
  
  // Assign cuts and check that they are sorted
  arma::field<arma::vec> cuts(n_eq);
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

  // Get variables related to heteroscedasticity if need
  const LogicalVector is_het = control_lnL["is_het"];
  const arma::field<arma::uvec> coef_var_ind = control_lnL["coef_var_ind"];
  arma::field<arma::vec> coef_var(n_eq);
  for (int i = 0; i < n_eq; i++)
  {
    if (is_het[i])
    {
      coef_var.at(i) = par.elem(coef_var_ind(i));
    }
  }
  
  // Assign coefficients associated with continious equations
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
  const CharacterVector marginal_names = control_lnL["marginal_names"];
  const arma::ivec marginal_par_n = control_lnL["marginal_par_n"];
  const arma::field<arma::uvec> marginal_par_ind = control_lnL["marginal_par_ind"];
  const bool is_marginal = marginal_names.size() > 0;
  arma::field<arma::vec> marginal_par(n_eq + n_eq2);
  List marginal_par_list;
  if (is_marginal)
  {
    for (int i = 0; i < n_eq; i++)
    {
      if (marginal_par_n.at(i) > 0)
      {
        marginal_par.at(i) = par.elem(marginal_par_ind.at(i));
        // return if coefficients are too large
        // if (any(abs(marginal_par.at(i)) > 10))
        // {
        //   if (is_diff)
        //   {
        //     NumericMatrix out_tmp(n_par, 1);
        //     std::fill(out_tmp.begin(), out_tmp.end(), -(1e+100));
        //     return(out_tmp);
        //   }
        //   NumericMatrix out_tmp(1,1);
        //   out_tmp(0, 0) = -(1e+100);
        //   return(out_tmp);
        // }
      }
    }
    // store the results into Rcpp list
    marginal_par_list = wrap(marginal_par);
  }

  // Estimate the linear index for each alternative
  arma::mat li_mean(n_obs, n_eq);
  arma::mat li_var(n_obs, n_eq);
  arma::mat li_lower(n_obs, n_eq);
  arma::mat li_upper(n_obs, n_eq);
  arma::mat li_y(n_obs, n_eq2);
  for (int i = 0; i < n_groups; i++)
  {
    // Ordered equations
    for (int j = 0; j < n_eq_g.at(i); j++)
    {
      //Prepare indexes
      int j_o = ind_eq(i).at(j);
      int group_j = groups(i, j_o);
      arma::uvec j_o_uvec = {(unsigned int)j_o};
      // Calculate linear indexes for mean and variance parts
      li_mean.submat(ind_g(i), j_o_uvec) = W(j_o).rows(ind_g(i)) * coef.at(j_o);
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
                                              li_mean.submat(ind_g(i), j_o_uvec);
        if (is_het[j_o])
        {
          li_lower.submat(ind_g(i), j_o_uvec) /= li_var.submat(ind_g(i), j_o_uvec);
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
                                              li_mean.submat(ind_g(i), j_o_uvec);
        if (is_het[j_o])
        {
          li_upper.submat(ind_g(i), j_o_uvec) /= li_var.submat(ind_g(i), j_o_uvec);
        }
      }
    }
    // Continuous equations
    if (is2)
    {
      for (int j = 0; j < n_eq2_g.at(i); j++)
      {
        //Prepare indexes
        int j_o = ind_eq2(i).at(j);
        arma::uvec j_o_uvec = {(unsigned int)j_o};
        li_y.submat(ind_g(i), j_o_uvec) = y.submat(ind_g(i), j_o_uvec) - 
                                          X(j_o).rows(ind_g(i)) * 
                                          coef2.at(j_o).row(groups2.at(i, j_o)).t();
      }
    }
  }

  // Vector to store log-likelihoods
  arma::vec lnL;
  if (out_type == "val")
  {
    lnL = arma::vec(n_obs);
  }

  // Calculate probabilities for each group
  for (int i = 0; i < n_groups; i++)
  {
    // Transform index
    arma::uvec i_uvec = {(unsigned int)i};
    
    // Modify covariance matrix accounting for 
    // continuous equations and their regimes
    if (is2)
    {
      for (int j1 = 0; j1 < n_eq2_g.at(i); j1++)
      {
        int j1_o = ind_eq2(i).at(j1);
        // arma::uvec j1_o_uvec = {j1_o};
        sigma.at(n_eq + j1_o, n_eq + j1_o) = par.at(var_y_ind(j1_o).at(groups2.at(i, j1_o)));
        for (int j2 = 0; j2 < n_eq; j2++)
        {
          sigma.at(j1_o + n_eq, j2) = par.at(cov_y_ind(j1_o).at(groups2.at(i, j1_o), j2));
          sigma.at(j2, j1_o + n_eq) = par.at(cov_y_ind(j1_o).at(groups2.at(i, j1_o), j2));
        }
      }
      if (n_eq2 >= 2)
      {
        for (int j1 = 1; j1 < n_eq2_g.at(i); j1++)
        {
          for (int j2 = 0; j2 < j1; j2++)
          {
            int j1_o = ind_eq2(i).at(j1);
            int j2_o = ind_eq2(i).at(j2);
            arma::uvec uvec_tmp = {sigma2_ind_mat(i).at(j1_o, j2_o)};
            arma::vec vec_tmp =  par.elem(uvec_tmp);
            sigma.at(j1_o + n_eq, j2_o + n_eq) = vec_tmp.at(0);
            sigma.at(j2_o + n_eq, j1_o + n_eq) = vec_tmp.at(0);
          }
        }
      }
    }

    // Construct covariance matrix for observable equations
    arma::mat sigma_g = sigma.submat(ind_eq_all(i), ind_eq_all(i));
    NumericMatrix sigma_g_R = wrap(sigma_g);

    // Check that matrix is positive defined after
    // accounting for the elements associated with
    // continuous equations
    if (is2)
    {
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
    }

    // Vector of zero means
    NumericVector mean_zero_R(n_eq_all_g.at(i));

    // Get lower and upper integration limits
    arma::mat li_lower_tmp = li_lower.submat(ind_g(i), ind_eq(i));
    arma::mat li_upper_tmp = li_upper.submat(ind_g(i), ind_eq(i));
    NumericMatrix li_lower_R = wrap(li_lower_tmp);
    NumericMatrix li_upper_R = wrap(li_upper_tmp);
    arma::mat li_y_tmp;
    NumericMatrix li_y_R;
    if (n_eq2_g(i) > 0)
    {
      li_y_tmp = li_y.submat(ind_g(i), ind_eq2(i));
      li_y_R = wrap(li_y_tmp);
    }

    // Create a list of marginal distributions parameters
    // accounting for unobservable equations
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
        marginal_names_tmp[j] = marginal_names[j_o];
      }
      if (is2)
      {
        for (int j = 0; j < n_eq2_g.at(i); j++)
        {
          int j_adj = j + n_eq_g.at(i);
          marginal_par_list_tmp[j_adj] = R_NilValue;
          marginal_names_tmp[j_adj] = "normal";
        }
      }
      marginal_par_list_tmp.names() = marginal_names_tmp;
      Nullable<List> marginal_par_list_tmp2(marginal_par_list_tmp);
      marginal_par_list_i = marginal_par_list_tmp2;
    }

    // Calculate probabilities and differentiate (if need) them
    List prob_list;
    if (is2 & (n_eq2_g(i) > 0))
    {
      IntegerVector given_x_tmp = Rcpp::seq(n_eq_g(i) + 1, n_eq_all_g(i));
      NumericVector given_x = as<NumericVector>(given_x_tmp);
      prob_list = mnorm::pmnorm(li_lower_R, li_upper_R,
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

    // Calculate densities if need
    List den_list;
    if (n_eq2_g(i) > 0)
    {
      NumericVector mean2_zero_R(n_eq2_g.at(i));
      arma::mat sigma2_g = sigma.submat(ind_eq2(i) + n_eq,
                                        ind_eq2(i) + n_eq);
      NumericMatrix sigma2_g_R = wrap(sigma2_g);
      den_list = mnorm::dmnorm(li_y_R,
                               mean2_zero_R, sigma2_g_R,
                               NumericVector(), 
                               true, is_diff, is_diff,
                               false,
                               R_NilValue, n_cores);
    }
    
    // // Store gradient information
    // if (is_jac)
    // {
    //   arma::mat jac_g = prob_list["grad_upper"];
    //   if (n_eq2_g(i) > 0)
    //   {
    //     arma::mat jac2_g = prob_list["grad_x"];
    //     jac.submat(ind_g(i), ind_eq(i)) = jac_g + jac2_g;
    //   }
    //   else
    //   {
    //     jac.submat(ind_g(i), ind_eq(i)) = jac_g;
    //   }
    // }
    
    // Store log-likelihood contributions
    if (out_type == "val")
    {
      NumericVector prob = prob_list["prob"];
      arma::vec prob_arma(prob.begin(), n_obs_g.at(i), false);
      if (n_eq2_g(i) > 0)
      {
        NumericVector den = den_list["den"];
        arma::vec den_arma(den.begin(), n_obs_g.at(i), false);
        lnL.elem(ind_g(i)) = prob_arma + den_arma;
      }
      else
      {
        lnL.elem(ind_g(i)) = prob_arma;
      }
    }

    // Store Jacobian if need
    if (is_diff)
    {
      arma::mat grad_lower = prob_list["grad_lower"];
      arma::mat grad_upper = prob_list["grad_upper"];
      arma::mat grad_sum = grad_upper + grad_lower;
      // for mean equation coefficients
      for (int j = 0; j < n_eq_g.at(i); j++)
      {
        int j_o = ind_eq(i).at(j);
        arma::uvec j_o_uvec = {(unsigned int)j_o};
        arma::mat W_tmp = W(j_o).rows(ind_g(i));
        arma::mat jac_tmp = W_tmp.each_col() % (-grad_sum.col(j));
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
      // for covariances
      arma::cube grad_sigma = prob_list["grad_sigma"];
      arma::umat sigma_ind_mat_o = sigma_ind_mat.submat(ind_eq(i), ind_eq(i));
      for (int j1 = 1; j1 < n_eq_g.at(i); j1++)
      {
        for (int j2 = 0; j2 < j1; j2++)
        {
          arma::vec grad_sigma_tmp = grad_sigma.tube(j1, j2);
          arma::uvec sigma_ind_uvec = {sigma_ind_mat_o.at(j1, j2)};
          jac.submat(ind_g(i), sigma_ind_uvec) = grad_sigma_tmp;
        }
      }
      // for cuts
      for (int j = 0; j < n_eq_g.at(i); j++)
      {
        int j_o = ind_eq(i).at(j);
        arma::uvec j_o_uvec = {(unsigned int)j_o};
        int group_j = groups(i, j_o);
        // lower cut
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
        // upper cut
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
      // for variance equation coefficients
      for (int j = 0; j < n_eq_g.at(i); j++)
      {
        int j_o = ind_eq(i).at(j);
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
      // Jacobian for parameters related to continuous equations
      if (is2 & (n_eq2_g.at(i) > 0))
      {
        arma::mat grad_given = prob_list["grad_given"];
        arma::mat grad_x = den_list["grad_x"];
        arma::cube grad_sigma_den = den_list["grad_sigma"];
        // for coefficients, variances and covariances of continuous
        // equations with ordered equations
        for (int j = 0; j < n_eq2_g.at(i); j++)
        {
          int j_o = ind_eq2(i).at(j);
          arma::uvec j_o_uvec = {(unsigned int)j_o};
          arma::mat X_tmp = X(j_o).rows(ind_g(i));
          // part related to probability
          arma::mat jac_tmp_1 = X_tmp.each_col() % (-grad_given.col(j));
          arma::vec jac_tmp_1_var = grad_sigma.tube(n_eq_g(i) + j,
                                                    n_eq_g(i) + j);
          for (int j1 = 0; j1 < n_eq_g.at(i); j1++)
          {
            int j1_o = ind_eq(i).at(j1);
            arma::uvec j1_o_uvec = {(unsigned int)j1_o};
            arma::uvec cov_y_ind_uvec = {cov_y_ind(j_o).at(groups2.at(i, j_o), j1_o)};
            arma::vec tube_tmp = grad_sigma.tube(n_eq_g(i) + j, j1);
            jac.submat(ind_g(i), cov_y_ind_uvec) = tube_tmp;
          }
          // part related to density
          arma::mat jac_tmp_2 = X_tmp.each_col() % (-grad_x.col(j));
          arma::vec jac_tmp_2_var = grad_sigma_den.tube(j, j);
          // combine both parts
          arma::uvec var_y_ind_uvec = {var_y_ind.at(j_o).at(groups2(i, j_o))};
          jac.submat(ind_g(i), coef2_ind(j_o).row(groups2(i, j_o))) = jac_tmp_1 + jac_tmp_2;
          jac.submat(ind_g(i), var_y_ind_uvec) = jac_tmp_1_var + jac_tmp_2_var;
        }
        // for covariances between continuous equations
        if (n_eq2_g.at(i) > 1)
        {
          arma::umat sigma2_ind_mat_o = sigma2_ind_mat(i).submat(ind_eq2(i), ind_eq2(i));
          for (int j1 = 1; j1 < n_eq2_g.at(i); j1++)
          {
            for (int j2 = 0; j2 < j1; j2++)
            {
              arma::vec grad_sigma2_prob = grad_sigma.tube(n_eq_g(i) + j1,
                                                           n_eq_g(i) + j2);
              arma::vec grad_sigma2_den = grad_sigma_den.tube(j1, j2);
              arma::uvec sigma2_ind_uvec = {sigma2_ind_mat_o.at(j1, j2)};
              jac.submat(ind_g(i), sigma2_ind_uvec) = grad_sigma2_prob + 
                                                      grad_sigma2_den;
            }
          }
        }
      }
      // Jacobian for parameters related to marginal distributions parameters
      if (is_marginal)
      {
        arma::field<arma::mat> grad_marginal = prob_list["grad_marginal"];
        for (int j = 0; j < n_eq_g.at(i); j++)
        {
          int j_o = ind_eq(i).at(j);
          if (marginal_par_n.at(j_o) > 0)
          {
            jac.submat(ind_g(i),  marginal_par_ind.at(j_o)) = grad_marginal.at(j);
          }
        }
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

  // Calculate log-likelihood
  double lnL_val = sum(lnL);
  if(std::isnan(lnL_val))
  {
    NumericMatrix out_tmp(1, 1);
    out_tmp(0, 0) = -(1e+100);
    return(out_tmp);
  }
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
//' @param out_type string represeint the output type of the function.
//' @param n_sim the number of random draws for multivariate 
//' normal probabilities.
//' @param n_cores the number of cores to be used. 
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix grad_mvoprobit(const arma::vec par,
                             const List control_lnL,
                             const String out_type = "grad",
                             const int n_sim = 1000,
                             const int n_cores = 1)
{
  return(lnL_mvoprobit(par, control_lnL, "grad", n_sim, n_cores));
}
