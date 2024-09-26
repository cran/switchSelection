#' Multivariate and multinomial sample selection and endogenous switching models
#' with multiple outcomes.
#' @description This function allows to estimate parameters of the 
#' multivariate and multinomial sample selection and endogenous switching models
#' with multiple outcomes. Both maximum-likelihood and two-step estimators are
#' implemented.
#' @template msel_param_Template
#' @template msel_details_Template
#' @template msel_return_Template
#' @template msel_references_Template
#' @template msel_examples_cps_Template
#' @template msel_examples1_Template
#' @template msel_examples2_Template
#' @template msel_examples3_Template
#' @template msel_examples4_Template
#' @template msel_examples5_Template
#' @template msel_examples6_Template
#' @template msel_examples7_Template
#' @template msel_examples8_Template
#' @template msel_examples9_Template
#'
msel <- function(formula        = NA,
                 formula2       = NA,
                 formula3       = NA,
                 data           = NULL,
                 groups         = NA,
                 groups2        = NA,
                 groups3        = NA,
                 marginal       = list(),
                 opt_type       = "optim",
                 opt_args       = NA,
                 start          = NULL,
                 estimator      = "ml",
                 cov_type       = "mm",
                 degrees        = NA,
                 degrees3       = NA,
                 n_sim          = 1000, 
                 n_cores        = 1,
                 control        = list(),
                 regularization = list(),
                 type3          = "logit")
{
  # List to store the output
  out        <- list()
  class(out) <- "msel"
  out$other  <- list(n_cores = n_cores, n_sim = n_sim)
  
  # Get the first step model if it has been provided
  model1 <- NULL
  if (is(object = formula, class2 = "msel"))
  {
    model1  <- formula
    formula <- model1$formula
  }
  if (is(object = formula3, class2 = "msel"))
  {
    model1   <- formula3
    formula3 <- model1$formula3
  }
  if (!is.null(model1))
  {
    formula  <- model1$formula
    formula3 <- model1$formula3
    marginal <- model1$marginal
    if (!any(is.na(groups2)))
    {
      groups  <- model1$groups
      groups3 <- model1$groups3
    }
    estimator <- "2step"
  }
  
  # Find the groups which have not been provided by the users
  is_na_group           <- c(any(is.na(groups)), 
                             any(is.na(groups2)), 
                             any(is.na(groups3)))
  out$other$is_na_group <- is_na_group
  
  # Deal with the control variables
  out_type <- "default"
  if (hasName(x = control, name = "out_type"))
  {
    out_type <- control$out_type
  }
  
  # -------------------------------------------------------
  # Validation
  # -------------------------------------------------------
  
  # Validate data
  if (!is.data.frame(data))
  {
    stop(paste0("Argument 'data' is missing or wrong. ",
                "Please, insure that 'data' is a dataframe ",
                "containing variables described in 'formula', 'formula2', ",
                "and 'formula3'.\n"))
  }
  
  # Validate estimator
  estimator <- tolower(estimator)
  if (estimator %in% c("two-step", "twostep","2-step", "2st", "2", 'mm', 'gmm'))
  {
    warning(paste0("It is assumed that the 'estimator' '", estimator,
                   "' is '2step'.\n"))
    estimator <- "2step"
  }
  if (estimator %in% c("mle", "likelihood","maximum-likelihood", 
                       "1step", "maxlik", "fiml"))
  {
    warning(paste0("It is assumed that the 'estimator' '", estimator,
                   "' is 'ml'.\n"))
    estimator <- "ml"
  }
  if (!(estimator %in% c("2step", "ml")))
  {
    stop(paste0("Wrong 'estimator' argument. ",
                "Please, insure that it is either '2step' or 'ml'.\n"))
  }
  out$estimator <- estimator
  
  # Validate formula
  is1 <- is.list(formula) | is.language(formula)
  if (is.null(formula))
  {
    formula <- NA
  }
  if (is1)
  {
    if (!is.list(formula))
    {
      tryCatch(formula <- as.formula(formula),
               error = function(e) {stop("Invalid 'formula' argument.\n")})
      formula <- list(formula)
    }
  }
  out$other$is1 <- is1
  
  # Validate formula2
  is2 <- is.list(formula2) | is.language(formula2)
  if (is.null(formula2))
  {
    formula2 <- NA
  }
  if (is2)
  {
    if (!is.list(formula2))
    {
      tryCatch(formula2 <- as.formula(formula2),
               error = function(e) {stop("Invalid 'formula2' argument.\n")})
      formula2 <- list(formula2)
    }
  }
  if (!is2 & (estimator != "ml"))
  {
    estimator <- "ml"
    warning(paste0("Since 'formula2' have not been provided it is assumed ",
                   "that 'estimator' is  'ml'.\n"))
  }
  out$other$is2 <- is2
  
  # Validate formula3
  is3 <- is.list(formula3) | is.language(formula3)
  if (is.null(formula3))
  {
    formula3 <- NA
  }
  if (is3)
  {
    tryCatch(formula3 <- as.formula(formula3),
             error = function(e) {stop("Invalid 'formula3' argument.\n")})
  }
  out$other$is3 <- is3
  
  # Validate the estimator one more time
  if (is2 & is3 & (estimator == "ml"))
  {
    estimator     <- "2step"
    out$estimator <- estimator
    warning("It is assumed that the 'estimator' is '2step'.")
  }
    
  # Check that at least some formulas have been provided
  if (!is1 & !is2 & !is3)
  {
    stop("No 'formula', 'formula2', or 'formula3' have been provided.\n")
  }
  
  # Get the number of equations
  n_eq  <- 0
  n_eq2 <- 0
  if (is1)
  {
    n_eq <- length(formula)
  }
  if (is2)
  {
    n_eq2 <- length(formula2)
  }
  out$other$n_eq  <- n_eq
  out$other$n_eq2 <- n_eq2
  
  # Validate groups
  if (is.null(groups))
  {
    groups <- NA
  }
  if (!is_na_group[1])
  {
    if (!is.matrix(groups))
    {
      if (is.vector(groups))
      {
        groups <- matrix(groups, ncol = 1)
      }
      else
      {
        groups <- as.matrix(groups)
      }
    }
    if (ncol(groups) != length(formula))
    {
      stop(paste0("Argument 'groups' is wrong. ",
                  "Please, insure that 'groups' is a matrix ",
                  "which number of columns equals to the length ",
                  "of 'formula'.",
                  "\n"))
    }
  }
  
  # Validate groups2
  if (is.null(groups2))
  {
    groups2 <- NA
  }
  if (!is_na_group[2])
  {
    if (!is.matrix(groups2))
    {
      if (is.vector(groups2))
      {
        groups2 <- matrix(groups2, ncol = 1)
      }
      else
      {
        groups2 <- as.matrix(groups2)
      }
    }
    if (ncol(groups2) != length(formula2))
    {
      stop(paste0("Argument 'groups2' is wrong. ",
                  "Please, insure that 'groups' is a matrix ",
                  "which number of columns equals to the length ",
                  "of 'formula2'.",
                  "\n"))
    }
  }
  
  # Validate opt_type
  opt_type_vec <- c("optim", "gena", "pso")
  if (!(opt_type %in% opt_type_vec))
  {
    stop(paste0("Argument 'opt_type' is wrong. ",
                "Please, insure that it is one of: ",
                paste0(opt_type_vec, collapse = ", "),
                ".\n"))
  }
  
  # Validate cov_type
  cov_type     <- tolower(cov_type)
  cov_type_vec <- NULL
  if (estimator == "ml")
  {
    cov_type_vec <- c("sandwich", "hessian", "gop", "no", "mm")
  }
  else
  {
    cov_type_vec <- c("mm", "gmm", "no")
  }
  if (!is.matrix(cov_type))
  {
    if (!(cov_type %in% cov_type_vec))
    {
      stop(paste0("Argument 'cov_type' is wrong. ",
                  "Please, insure that it is one of: ",
                  paste0(cov_type_vec, collapse = ", "),
                  ".\n"))
    }
  }
  if (cov_type == "gmm")
  {
    cov_type <- "mm"
  }
  if (cov_type == "opg")
  {
    cov_type <- "gop"
    warning("It is assumed that cov_type = 'gop'.\n")
  }
  
  # Validate marginal
  n_marginal <- length(marginal)
  if (n_marginal > 0)
  {
    # transform vector of characters into the list of NULL values
    if (!is.list(marginal))
    {
      if (is.character(marginal))
      {
        marginal_list        <- vector(mode = "list", length = n_marginal)
        names(marginal_list) <- marginal
        marginal             <- marginal_list
      }
      else
      {
        stop("Argument 'marginal' should be either character vector or a list.")
      }
    }
    if (n_marginal != length(formula))
    {
      stop(paste0("The number of elements in 'marginal' should be the same ",
                  "as the number of ordered equations. ",
                  "Please, insure that 'length(marginal) == length(formula)' ",
                  "or 'marginal' is an empty character vector."))
    }
  }
  
  # Validate degrees
  if ((estimator == "2step") & is1)
  {
    if (is.null(degrees))
    {
      degrees <- NA
    }
    if (any(is.na(degrees)))
    {
      degrees <- matrix(data = 1, nrow = n_eq2, ncol = n_eq)
    }
    if (is.vector(degrees))
    {
      if (length(degrees) != n_eq)
      {
        stop(paste0("Invalid 'degrees' argument. Please, insure that ", 
                    "if 'degrees' is a vector then its length equals to ",
                    "the number of selection equations."))
      }
      degrees <- matrix(rep(degrees, n_eq2), nrow = n_eq2, byrow = TRUE)
    }
    if (!is.matrix(degrees))
    {
      stop(paste0("Invalid 'degrees' argument. Please, insure that ", 
                  "it is a matrix."))
    }
    if (nrow(degrees) > n_eq2)
    {
      stop(paste0("Invalid 'degrees' argument since it has ", nrow(degrees),
                  " rows. Please, insure that it has ", n_eq2, "rows."))
    }
    if (ncol(degrees) > n_eq)
    {
      stop(paste0("Invalid 'degrees' argument since it has ", ncol(degrees),
                  " columns. Please, insure that it has ", n_eq, "columns."))
    }
    if (any((degrees %% 1) != 0) | any(degrees < 0))
    {
      stop(paste0("Invalid 'degrees' argument. Please, insure that ", 
                  "all elements of degress are non-negative integers."))
    }
    for (v in seq_len(n_eq2))
    {
      formula2_terms <- attr(terms(formula2[[v]]), "term.labels")
      if (any(grepl("lambda", formula2_terms, fixed = TRUE)))
      {
        degrees[v, ] <- 0
      }
    }
  }
  out$other$degrees   <- degrees
  
  # Validate type3
  type3_vec <- c("logit", "probit")
  if (type3 %in% c("logistic", "mlogit", "mnlogit", "multinomial logit"))
  {
    type3 <- "logit"
    warning("It is assumed that 'type3' is 'logit'.")
  }
  if (type3 %in% c("normal", "mprobit", "mnprobit", "normal", "gaussian",
                       "norm", "multinomial probit"))
  {
    type3 <- "probit"
    warning("It is assumed that 'type3' is 'probit'.")
  }
  if (!(type3 %in% type3_vec))
  {
    stop(paste0("Wrong 'type3' argument. ",
                "Please, insure that it is either 'logit' or 'probit'.\n"))
  }
  out$type3 <- type3
  
  # Validate degrees3
  n_degrees3 <- 0
  if ((estimator == "2step") & is3)
  {
    z_mn_name  <- all.vars(formula3)[1]
    n_degrees3 <- sum(unique(na.omit(data[, z_mn_name])) >= 0)
    if (type3 == "probit")
    {
      n_degrees3 <- n_degrees3 - 1
    }
    if (is.null(degrees3))
    {
      degrees3 <- NA
    }
    if (any(is.na(degrees3)))
    {
      degrees3 <- matrix(data = 1, nrow = n_eq2, ncol = n_degrees3)
    }
    if (is.vector(degrees3))
    {
      if (length(degrees3) != n_degrees3)
      {
        stop(paste0("Invalid 'degrees3' argument. Please, insure that ", 
                    "if 'degrees3' is a vector then its length equals to ",
                    "the number of multinomial equations."))
      }
      degrees3 <- matrix(rep(degrees3, n_eq2), nrow = n_eq2, byrow = TRUE)
    }
    if (!is.matrix(degrees3))
    {
      stop(paste0("Invalid 'degrees3' argument. Please, insure that ", 
                  "it is a matrix."))
    }
    if (nrow(degrees3) > n_eq2)
    {
      stop(paste0("Invalid 'degrees3' argument since it has ", nrow(degrees3),
                  " rows. Please, insure that it has ", n_eq2, "rows."))
    }
    if (ncol(degrees3) > n_degrees3)
    {
      stop(paste0("Invalid 'degrees3' argument since it has ", ncol(degrees3),
                  " columns. Please, insure that it has ", 
                  model1$other$n_eq3, "columns."))
    }
    if (any((degrees3 %% 1) != 0) | any(degrees3 < 0))
    {
      stop(paste0("Invalid 'degrees3' argument. Please, insure that ", 
                  "all elements of degress3 are non-negative integers."))
    }
    for (v in seq_len(n_eq2))
    {
      formula2_terms <- attr(terms(formula2[[v]]), "term.labels")
      if (any(grepl("lambda_mn", formula2_terms, fixed = TRUE)))
      {
        degrees3[v, ] <- 0
      }
    }
  }
  out$other$degrees3   <- degrees3
  
  # Validate groups3
  if (is.matrix(groups3))
  {
    groups3 <- as.vector(groups3)
  }

  # -------------------------------------------------------
  # Formulas
  # -------------------------------------------------------
  
  # Adjust the formulas
  formula_list <- formula_msel(object   = out,      formula  = formula, 
                               formula2 = formula2, formula3 = formula3,
                               degrees  = degrees,  degrees3 = degrees3)
  
  # Get the results
  is_het          <- formula_list$is_het
  formula         <- formula_list$formula
  formula2        <- formula_list$formula2
  formula3        <- formula_list$formula3
  formula_mean    <- formula_list$formula_mean
  formula_var     <- formula_list$formula_var
  coef_lambda_ind <- formula_list$coef_lambda_ind
  
  # Clear the memory
  formula_list <- NULL
  
  # Save the results to the output
  out$other$is_het          <- is_het
  out$formula               <- formula
  out$formula2              <- formula2
  out$formula3              <- formula3
  out$other$formula_mean    <- formula_mean
  out$other$formula_var     <- formula_var
  out$other$coef_lambda_ind <- coef_lambda_ind
  
  # -------------------------------------------------------
  # The first step
  # -------------------------------------------------------
  
  # Prepare the variables to store the selectivity terms
  lambda    <- NULL
  lambda_mn <- NULL
  
  # Get estimates of the first step
  if ((estimator == "2step") & (is1 | is3))
  {
    if (is.null(model1))
    {
      # Multivariate and multinomial models
      model1 <- msel(formula   = formula,  formula3  = formula3,
                     groups    = groups,   groups3   = groups3,
                     data      = data,     estimator = "ml",
                     marginal  = marginal, opt_type  = opt_type,
                     opt_args  = opt_args, n_sim     = n_sim,
                     n_cores   = n_cores,  cov_type  = "mm",
                     type3     = type3)
      if (!is(object = model1, class2 = "msel"))
      {
        return(model1)
      }
    }
    out$model1 <- model1
    data       <- model1$data
    # Selectivity terms of the ordinal equations
    if (is1)
    {
      lambda <- model1$lambda
      for (i in seq_len(model1$other$n_eq))
      {
        data[[paste0("lambda", i)]] <- lambda[, i]
      }
    }
    # Selectivity term of the multinomial equation
    if (is3)
    {
      lambda_mn <- model1$lambda_mn
      for (i in seq_len(model1$other$n_eq3 - (type3 == "probit")))
      {
        data[[paste0("lambda", i, "_mn")]] <- lambda_mn[, i]
      }
    }
  }

  # -------------------------------------------------------
  # Missings
  # -------------------------------------------------------
  
  # Determine the groups before removing the missing values
  if (all(is_na_group))
  {
    groups_list           <- groups_msel(object  = out,    data    = data, 
                                         groups  = groups, groups2 = groups2, 
                                         groups3 = groups3)
    groups                <- groups_list$groups
    groups2               <- groups_list$groups2
    groups3               <- groups_list$groups3
    is_na_group           <- c(!is1, !is2, !is3)
    out$other$is_na_group <- is_na_group
    groups_list           <- NULL
  }

  # Get data including only the complete observations 
  data <- complete_msel(object = out, data = data)
  
  # Need run second time because of the endogenous regressors
  data <- complete_msel(object = out, data = data)

  # -------------------------------------------------------
  # Groups
  # -------------------------------------------------------

  # Get the groups associated with the data
  groups_list <- groups_msel(object  = out,    data    = data, 
                             groups  = groups, groups2 = groups2, 
                             groups3 = groups3)

  # Store to the variables
  groups     <- groups_list$groups
  groups2    <- groups_list$groups2
  groups3    <- groups_list$groups3
  ind_g      <- groups_list$ind_g
  ind_eq     <- groups_list$ind_eq
  ind_eq2    <- groups_list$ind_eq2
  ind_eq3    <- groups_list$ind_eq3
  ind_eq_all <- groups_list$ind_eq_all
  n_cuts_eq  <- groups_list$n_cuts_eq
  n_eq3      <- groups_list$n_eq3
  n_groups   <- groups_list$n_groups
  data       <- groups_list$data
  n_obs      <- groups_list$n_obs
  n_obs_g    <- groups_list$n_obs_g
  n_eq_g     <- groups_list$n_eq_g
  n_eq2_g    <- groups_list$n_eq2_g
  n_eq_all_g <- groups_list$n_eq_all_g
  ind_g_all  <- groups_list$ind_g_all
  
  # Clear the memory
  groups_list <- NULL
  
  # Store to the outup
  out$groups           <- groups
  out$groups2          <- groups2
  out$groups3          <- groups3
  out$other$ind_g      <- ind_g
  out$other$ind_eq     <- ind_eq
  out$other$ind_eq2    <- ind_eq2
  out$other$ind_eq3    <- ind_eq3
  out$other$ind_eq_all <- ind_eq_all
  out$other$n_cuts_eq  <- n_cuts_eq
  out$other$n_eq3      <- n_eq3
  out$other$n_groups   <- n_groups
  out$data             <- data
  out$other$n_obs      <- n_obs
  out$other$n_obs_g    <- n_obs_g
  out$other$n_eq_g     <- n_eq_g
  out$other$n_eq2_g    <- n_eq2_g
  out$other$n_eq_all_g <- n_eq_all_g
  out$other$ind_g_all  <- ind_g_all
  
  # -------------------------------------------------------
  # Variables
  # -------------------------------------------------------
  
  # Get dependent and independent variables
  data_list <- data_msel(object = out, data = data)
  
  # Store to the variables
  W_mean <- data_list$W_mean
  W_var  <- data_list$W_var
  X      <- data_list$X
  W_mn   <- data_list$W_mn
  z      <- data_list$z
  y      <- data_list$y
  z_mn   <- data_list$z_mn
  
  # Clear the memory
  data_list <- NULL
  
  # Store to the output
  out$W_mean <- W_mean
  out$W_var  <- W_var
  out$X      <- X
  out$W_mn   <- W_mn
  out$z      <- z
  out$y      <- y
  out$z_mn   <- z_mn
  
  # -------------------------------------------------------
  # Names of the variables
  # -------------------------------------------------------
  
  names_list <- names_msel(object = out)
  
  # Store to the variables
  z_names         <- names_list$z_names
  y_names         <- names_list$y_names
  z_mn_names      <- names_list$z_mn_names
  colnames_X      <- names_list$colnames_X 
  colnames_W_mean <- names_list$colnames_W_mean
  colnames_W_var  <- names_list$colnames_W_var
  colnames_W_mn   <- names_list$colnames_W_mn
  
  # Clear the memory
  names_list <- NULL
  
  # Store to the output
  out$other$z_names         <- z_names
  out$other$y_names         <- y_names
  out$other$z_mn_names      <- z_mn_names
  out$other$colnames_X      <- colnames_X 
  out$other$colnames_W_mean <- colnames_W_mean
  out$other$colnames_W_var  <- colnames_W_var
  out$other$colnames_W_mn   <- colnames_W_mn
  
  # All the names
  all_names           <- c(z_names,        y_names,        z_mn_names, 
                           unlist(colnames_X),     unlist(colnames_W_mean), 
                           unlist(colnames_W_var), colnames_W_mn)
  all_names           <- na.omit(unique(all_names))
  out$other$all_names <- all_names
  
  # Names of the selectivity correction terms
  lambda_names <- numeric()
  if (estimator == "2step")
  {
    lambda_names <- all_names[grepl("lambda", all_names, fixed = TRUE)]
  }
  out$other$lambda_names <- lambda_names
  
  # Name the groups
  if (is1)
  {
    colnames(groups) <- z_names
    out$groups       <- groups
  }
  if (is2)
  {
    colnames(groups2) <- y_names
    out$groups2       <- groups2
  }
  
  # -------------------------------------------------------
  # Indexes of the parameters
  # -------------------------------------------------------
  
  # Get the number of coefficients (regressors) for each equation and 
  # determine their indexes in the parameters vector
  n_coef   <- vector("mode" = "numeric", length = n_eq)
  coef_ind <- vector("mode" = "list",    length = n_eq)
  if (is1)
  {
    for(i in seq_len(n_eq))
    {
      n_coef[i] <- ncol(W_mean[[i]])
      if (i != 1)
      {
        coef_ind[[i]] <- (coef_ind[[i - 1]][n_coef[i - 1]] + 1):
                         ((coef_ind[[i - 1]][n_coef[i - 1]] + n_coef[i]))
      }
      else
      {
        coef_ind[[i]] <- seq_len(n_coef[i])
      }
    }
  }
  
  # Get sigma elements which are not identified
  # and store them into the matrix
  sigma_omit <- matrix(0, nrow = n_eq, ncol = n_eq)
  if (is1)
  {
    if (n_eq > 1)
    {
      for (i in seq_len(n_eq - 1))
      {
        for(j in (i + 1):n_eq)
        {
          sigma_omit[i, j] <- as.numeric(all((groups[, i, drop = FALSE] == -1) | 
                                             (groups[, j, drop = FALSE] == -1)))
          sigma_omit[j, i] <- sigma_omit[i, j]
        }
      }
    }
  }
  out$other$sigma_omit <- sigma_omit
  
  # The number of the unidentified sigma elements
  n_sigma_omit <- sum(sigma_omit) / 2
  
  # Calculate the total number of the estimated coefficients
  n_coef_total <- sum(n_coef)
  
  # Get the number of the covariance matrix elements to be estimated
  n_sigma           <- (n_eq ^ 2 - n_eq) / 2 - n_sigma_omit
  out$other$n_sigma <- n_sigma
  
  # Get the number of the cut points to be considered
  n_cuts <- sum(n_cuts_eq)
  
  # Get indexes of the sigma elements in the parameters vector
  sigma_ind <- is1
  if (n_sigma > 0)
  {
    sigma_ind <- (n_coef_total + 1):(n_coef_total + n_sigma)
  }

  # Store the indexes of the sigma elements in a matrix form
  sigma_ind_mat <- matrix(1, nrow = n_eq, ncol = n_eq)
  if (is1)
  {
    if (n_eq >= 2)
    {
      counter <- sigma_ind[1]
      for (i in 2:n_eq)
      {
        for (j in seq_len(i - 1))
        {
          if (sigma_omit[i, j] == 0)
          {
            sigma_ind_mat[i, j] <- counter
            sigma_ind_mat[j, i] <- sigma_ind_mat[i, j]
            counter             <- counter + 1
          }
          else
          {
            # special number for omitted covariances
            sigma_ind_mat[i, j] <- 0
          }
        }
      }
    }
  }

  # Store the indexes of the cuts in a list form
  cuts_ind <- vector(mode = "list", length = n_eq)
  n_par    <- n_coef_total + n_sigma
  if (is1)
  {
    for(i in seq_len(n_eq))
    {
      for (j in seq_len(n_cuts_eq[i]))
      {
        n_par            <- n_par + 1
        cuts_ind[[i]][j] <- n_par
      }
    }
  }
  
  # Deal with the coefficients of the variance equation
  n_coef_var   <- vector("mode" = "numeric", length = n_eq)
  coef_var_ind <- vector("mode" = "list",    length = n_eq)
  if (is1)
  {
    for(i in seq_len(n_eq))
    {
      if(is_het[i])
      {
        n_coef_var[i]     <- ncol(W_var[[i]])
        coef_var_ind[[i]] <- (n_par + 1):(n_par + n_coef_var[i])
        n_par             <- n_par + n_coef_var[i]
      }
      else
      {
        coef_var_ind[[i]] <- as.vector(1)
      }
    }
  }
  
  # Calculate the number of the regimes for each continuous variable
  n_regimes <- vector(mode = "numeric", length = 0)
  if (is2)
  {
    n_regimes <- vector(mode = "numeric", length = n_eq2)
    for (i in seq_len(n_eq2))
    {
      n_regimes[i] <- max(groups2[, i] + 1)
    }
  }
  out$other$n_regimes <- n_regimes

  # Deal with the coefficients of the continuous equations
  n_coef2   <- vector(mode = "numeric", length = 0)
  coef2_ind <- list(matrix(1))
  if (is2)
  {
    n_coef2   <- vector(mode = "numeric", length = n_eq2)
    coef2_ind <- vector(mode = "list",    length = n_eq2)
    for(i in seq_len(n_eq2))
    {
      n_coef2[i] <- ncol(X[[i]])
      # coefficients for each regime are located in the different rows
      coef2_ind[[i]] <- matrix((n_par + 1):(n_par + n_coef2[i] * n_regimes[i]), 
                               nrow = n_regimes[i], ncol = n_coef2[i], 
                               byrow = TRUE)
      n_par          <- n_par + length(coef2_ind[[i]])
    }
  }
  
  # Determine the structure of the covariances between the continuous
  # equations for different regimes
  regimes         <- as.matrix(t(hpa::polynomialIndex(n_regimes - 1)))
  n_regimes_total <- nrow(regimes)
  
  # If need add the covariances related to the continuous equations
  var2_ind       <- list(as.vector(1))
  cov2_ind       <- list(matrix(1))
  sigma2_ind     <- list(as.vector(1))
  sigma2_ind_mat <- list(matrix(1))
  if (is2 & (estimator == "ml"))
  {
    # variances
    for (i in seq_len(n_eq2))
    {
      var2_ind[[i]] <- (n_par + 1):(n_par + n_regimes[i])
      n_par         <- n_par + n_regimes[i]
    }
      # ordinal equations
    if (is1)
    {
      cov2_ind <- vector(mode = "list", length = n_eq2)
      for (i in seq_len(n_eq2))
      {
        cov2_ind[[i]] <- matrix(ncol = n_eq, nrow = n_regimes[i])
        for (j in seq_len(n_regimes[i]))
        {
          cov2_ind[[i]][j, ] <- (n_par + 1):(n_par + n_eq)
          n_par              <- n_par + n_eq
        }
      }
    }
      # multinomial equations
    if (is3 & (type3 == "probit"))
    {
      for (i in seq_len(n_eq2))
      {
        cov2_ind[[i]] <- matrix(ncol = n_eq3 - 1, nrow = n_regimes[i])
        for (j in seq_len(n_regimes[i]))
        {
          cov2_ind[[i]][j, ] <- (n_par + 1):(n_par + n_eq3 - 1)
          n_par              <- n_par + n_eq3 - 1
        }
      }
    }
    # covariances between the continuous equations
    if (n_eq2 >= 2)
    {
      n_sigma2           <- (n_eq2 ^ 2 - n_eq2) / 2
      out$other$n_sigma2 <- n_sigma2
      sigma2_ind         <- vector(mode = "list", length = n_sigma2)
      regimes_pair       <- vector(mode = "list", length = n_sigma2)
      sigma2_ind_mat     <- vector(mode = "list", length = n_groups)
      for (g in seq_len(n_groups))
      {
        sigma2_ind_mat[[g]] <- matrix(1, nrow = n_eq2, ncol = n_eq2)
      }
      # indexes for each covariance depending on the regime
      counter <- 1
      for (i in 2:n_eq2)
      {
        for (j in seq_len(i - 1))
        {
          sigma2_ind[[counter]]   <- numeric()
          regimes_pair[[counter]] <- numeric()
          # note that i > j so in regimes_pair equation with
          # the greater index goes first
          pairs          <- unique(groups2[, c(i, j), drop = FALSE])
          pairs          <- pairs[(pairs[, 1, drop = FALSE] != -1) & 
                                  (pairs[, 2, drop = FALSE] != -1), , 
                                  drop = FALSE]
          regimes_pair_n <- ifelse(is.matrix(pairs), nrow(pairs), 0)
          if (regimes_pair_n > 0)
          {
            regimes_pair[[counter]] <- pairs
            sigma2_ind[[counter]]   <- (n_par + 1):(n_par + regimes_pair_n)
            n_par                   <- n_par + regimes_pair_n
            # indexes of covariances depending on group
            for (g in seq_len(n_groups))
            {
              cond <- which((regimes_pair[[counter]][, 1] == groups2[g, i]) &
                            (regimes_pair[[counter]][, 2] == groups2[g, j]))
              if (length(cond) > 0)
              {
                sigma2_ind_mat[[g]][i, j] <- sigma2_ind[[counter]][cond]
                sigma2_ind_mat[[g]][j, i] <- sigma2_ind[[counter]][cond]
              }
            }
          }
          counter <- counter + 1
        }
      }
      out$other$regimes_pair <- regimes_pair
    }
  }

  # Parameters associated with the marginal distributions
  is_marginal              <- length(marginal) > 0
  out$other$is_marginal    <- is_marginal
  marginal_par_ind         <- list(as.vector(1))
  marginal_par_n           <- rep(0, n_eq)
  out$other$marginal_par_n <- marginal_par_n
  marginal_names           <- character()
  if (is_marginal)
  {
    marginal_names <- names(marginal)
    for (i in seq_len(n_eq))
    {
      marginal_par_n[i] <- ifelse((length(marginal[[i]]) != 0) &
                                  !is.null(marginal[[i]]), 
                                  as.numeric(marginal[[i]]), 0)
      if (marginal_par_n[i] > 0)
      {
        marginal_par_ind[[i]] <- (n_par + 1):(n_par + marginal_par_n[i])
        n_par                 <- n_par + length(marginal_par_ind[[i]])
      }
    }
  }
  out$marginal <- marginal
  
  # Round the number of parameters to prevent the bugs
  marginal_par_n           <- round(marginal_par_n)
  out$other$marginal_par_n <- marginal_par_n
  
  # Coefficients of the multinomial equation
  n_coef3   <- 0
  coef3_ind <- matrix(1)
  if (is3)
  {
    n_coef3   <- ncol(W_mn)
    coef3_ind <- matrix(nrow = n_eq3 - 1, ncol = n_coef3)
    for(i in seq_len(n_eq3 - 1))
    {
      coef3_ind[i, ] <- n_par + (((i - 1) * n_coef3 + 1):(i * n_coef3))
    }
    n_par <- n_par + length(coef3_ind)
  }
  out$other$n_coef3 <- n_coef3
  
  # Covariances of the multinomial equations
  # Note: if sigma3_ind_mat[k, ] = c(i, j) it means
  # that sigma3_ind[k] refers to the sigma3[i, j]
  n_sigma3       <- ((n_eq3 - 1) ^ 2 - n_eq3 + 1) / 2 + (n_eq3 - 2)
  sigma3_ind     <- 1
  sigma3_ind_mat <- matrix(1)
  if (is3 & (type3 == "probit"))
  {
    if (n_sigma3 > 0)
    {
      sigma3_ind     <- (n_par + 1):(n_par + n_sigma3)
      sigma3_ind_mat <- matrix(1, nrow = n_sigma3, ncol = 2)
    }
    if (n_eq3 > 2)
    {
      counter <- 1
      for (i in seq_len(n_eq3 - 1))
      {
        for (j in 1:i)
        {
          if (!((i == 1) & (j == 1)))
          {
            sigma3_ind_mat[counter, 1] <- i
            sigma3_ind_mat[counter, 2] <- j
            counter                    <- counter + 1
          }
        }
      }
    }
    n_par <- n_par + length(sigma3_ind)
  }
  
  # Save the number of the estimated parameters
  out$other$n_par <- n_par
  
  # Save the indexes
  out$ind <- list(coef         = coef_ind,      coef_var     = coef_var_ind,
                  cuts         = cuts_ind,      sigma        = sigma_ind,
                  g            = ind_g,         eq           = ind_eq,
                  var2         = var2_ind,      sigma2       = sigma2_ind,
                  coef2        = coef2_ind,     cov2         = cov2_ind,
                  sigma2       = sigma2_ind,    marginal_par = marginal_par_ind,
                  sigma_mat    = sigma_ind_mat, coef3        = coef3_ind,
                  sigma3       = sigma3_ind,    sigma3_mat   = sigma3_ind_mat)
  
  # -------------------------------------------------------
  # Maximum-likelihood estimator
  # -------------------------------------------------------

  # Store the information need for the maximum-likelihood 
  # estimator into the list
  groups_tmp <- groups
  if (any(is.na(groups)))
  {
    groups_tmp <- matrix()
  }
  groups2_tmp <- groups2
  if (any(is.na(groups2)))
  {
    groups2_tmp <- matrix()
  }
  groups3_tmp <- groups3
  if (any(is.na(groups3)))
  {
    groups3_tmp <- vector(mode = "numeric")
  }
  control_lnL <- list(n_par            = n_par,
                      n_obs            = n_obs,
                      n_eq             = n_eq,
                      n_eq2            = n_eq2,
                      n_eq3            = n_eq3,
                      n_coef           = n_coef,
                      n_coef2          = n_coef2,
                      n_coef3          = n_coef3,
                      n_regimes        = n_regimes,
                      coef_ind         = lapply(coef_ind, function(x){x - 1}),
                      coef_var_ind     = lapply(coef_var_ind, 
                                                function(x){x - 1}),
                      coef2_ind        = lapply(coef2_ind, function(x){x - 1}),
                      sigma_ind        = sigma_ind - 1,
                      sigma2_ind       = lapply(sigma2_ind, function(x){x - 1}),
                      sigma_ind_mat    = sigma_ind_mat - 1,
                      sigma_omit       = sigma_omit,
                      sigma2_ind_mat   = lapply(sigma2_ind_mat, 
                                                function(x){x - 1}),
                      sigma3_ind       = sigma3_ind - 1,
                      sigma3_ind_mat   = sigma3_ind_mat - 1,
                      var2_ind         = lapply(var2_ind, function(x){x - 1}),
                      cov2_ind         = lapply(cov2_ind, function(x){x - 1}),
                      cuts_ind         = lapply(cuts_ind, function(x){x - 1}),
                      marginal_par_ind = lapply(marginal_par_ind, 
                                                function(x){x - 1}),
                      coef3_ind        = coef3_ind - 1,
                      n_groups         = n_groups,
                      ind_g            = lapply(ind_g, function(x){x - 1}),
                      ind_eq           = lapply(ind_eq, function(x){x - 1}),
                      ind_eq2          = lapply(ind_eq2, function(x){x - 1}),
                      ind_eq_all       = lapply(ind_eq_all, function(x){x - 1}),
                      n_cuts_eq        = n_cuts_eq,
                      groups           = groups_tmp,
                      groups2          = groups2_tmp,
                      groups3          = groups3_tmp,
                      n_obs_g          = n_obs_g,
                      n_eq_g           = n_eq_g,
                      n_eq2_g          = n_eq2_g,
                      n_eq_all_g       = n_eq_all_g,
                      is_het           = is_het,
                      is1              = is1,
                      is2              = is2,
                      is3              = is3,
                      y                = y,
                      X                = X,
                      W                = W_mean,
                      W_var            = W_var,
                      W_mn             = W_mn,
                      marginal_names   = marginal_names,
                      marginal_par_n   = marginal_par_n,
                      type3            = type3)
  out$control_lnL <- control_lnL
  
  # Validate the regularization
  regularization <- regularization_validate(regularization = regularization, 
                                            n_par          = n_par, 
                                            estimator      = estimator)
  out$other$regularization <- regularization
  
  # Vector of the estimates of the estimated parameters
  par <- rep(NA, n_par)
  
  # Create a starting point for the numeric optimization
  if (is.null(start) & (estimator == "ml"))
  {
    start <- rep(0, n_par)
    if (is1)
    {
      for (i in seq_len(n_eq))
      {
        # adjusted ordinal outcome
        z_tmp                <- z[, i]
        z_tmp[z_tmp == -1]   <- NA
        is_na_z              <- is.na(z[, i])
        z_tmp                <- na.omit(z[, i])
        # cuts
        for (j in seq_len(n_cuts_eq[i]))
        {
          start[cuts_ind[[i]]][j] <- qnorm(mean(na.omit(z_tmp) <= (j - 1)))
        }
        # coefficients
        z_tmp_table          <- table(z_tmp)
        z_tmp_ind            <- which.min(abs((cumsum(z_tmp_table) / 
                                               sum(z_tmp_table)) - 0.5)) - 1
        z_tmp_cond           <- z_tmp <= z_tmp_ind
        z_tmp[z_tmp_cond]    <- 0
        z_tmp[!z_tmp_cond]   <- 1
        data_tmp             <- data[!is_na_z, ]
        data_tmp[z_names[i]] <- z_tmp
        model_tmp            <- glm(formula_mean[[i]], data = data_tmp,
                                    family = binomial(link = "logit"))
        start[coef_ind[[i]]] <- coef(model_tmp)[-1] * (sqrt(3) / pi)
        # heteroscedasticity
        if (is_het[i])
        {
          n_coef_var_ind_i <- length(coef_var_ind[[i]])
          for (j in 1:n_coef_var_ind_i)
          {
            W_var_j_mean <- mean(W_var[[i]][, j], na.rm = TRUE)
            val_tmp <- 1
            if (abs(W_var_j_mean) > 1)
            {
              val_tmp <- val_tmp / W_var_j_mean
            }
            val_tmp <- val_tmp / sd(W_var[[i]][, j], na.rm = TRUE)
            start[coef_var_ind[[i]][j]] <- val_tmp
          }
          start[coef_var_ind[[i]]] <- start[coef_var_ind[[i]]] / 
                                      (5 * n_coef_var_ind_i)
        }
      }
    }
    for (i in seq_len(n_eq2))
    {
      data_tmp <- data.frame(y = y[, i], X[[i]])
      data_tmp <- data_tmp[rowSums(is.na(X[[i]])) == 0, ]
      for (j in seq_len(n_regimes[i]))
      {
        ind_g_include_tmp <- which(groups2[, i] == (j - 1))
        ing_g_tmp         <- NULL
        for (t in ind_g_include_tmp)
        {
          ing_g_tmp <- c(ing_g_tmp, ind_g[[t]])
        }
        data_tmp_g <- data_tmp[ing_g_tmp, ]
        model_tmp  <- lm(y ~ . + 0, 
                         data = data_tmp_g[!is.na(data_tmp_g$y), ])
        start[coef2_ind[[i]][j, ]] <- coef(model_tmp)
        start[var2_ind[[i]][j]]    <- sigma(model_tmp) ^ 2
      }
    }
    if (is_marginal)
    {
      for (i in seq_len(n_eq))
      {
        if ((marginal_names[i] == "PGN") | (marginal_names[i] == "hpa"))
        {
          start[marginal_par_ind[[i]]] <- 1e-8
        }
        if (marginal_names[i] == "student")
        {
          start[marginal_par_ind[[i]]] <- 30
        }
      }
    }
    if (is3 & (type3 == "probit"))
    {
      if (n_eq3 > 2)
      {
        for (i in seq_len(n_sigma3))
        {
          if (sigma3_ind_mat[i, 1] == sigma3_ind_mat[i, 2])
          {
            start[sigma3_ind[i]] <- 1
          }
        }
      }
    }
  }
  out$start <- start
  
  # Code to test the derivatives
  if (hasName(control, "test_grad"))
  {
    f0 <- lnL_msel(par = start,
                   n_sim = n_sim, n_cores = n_cores,
                   control_lnL = control_lnL)
    grad.a <- lnL_msel(par = start,
                       n_sim = n_sim, n_cores = n_cores,
                       control_lnL = control_lnL,
                       out_type = "grad")
    grad.n <- rep(NA, n_par)
    for (i in seq_len(n_par))
    {
      delta <- 1e-6
      start.delta <- start
      start.delta[i] <- start[i] + delta
      f1 <- lnL_msel(par = start.delta,
                          n_sim = n_sim, n_cores = n_cores,
                          control_lnL = control_lnL)
      grad.n[i] <- (f1 - f0) / delta
    }
    return(cbind(analytical = as.vector(grad.a), numeric = as.vector(grad.n)))
  }

  # Optimization
  opt <- NULL
  if (estimator == "ml")
  {
    opt <- opt_switchSelection(opt_args       = opt_args, 
                               control_lnL    = control_lnL,
                               n_sim          = n_sim, 
                               n_cores        = n_cores,
                               type           = "msel", 
                               start          = start, 
                               opt_type       = opt_type, 
                               regularization = regularization)
    par <- opt$par
  }

  # -------------------------------------------------------
  # Two-step estimator
  # -------------------------------------------------------

  # List to store the least squares estimators
  model2_list <- NULL
  
  # List to store the indexes of the observations for different regimes
  ind_regime  <- NULL
  
  # Matrix to store the predictions of the continuous equations
  # used by V-stage least squares estimator
  y_pred  <- NULL

  # The main routine
  if (estimator == "2step")
  {
    # Initialize the matrix to store the predictions
    y_pred           <- matrix(NA, nrow = n_obs, ncol = n_eq2)
    colnames(y_pred) <- y_names
    
    # Store the estimates of the first step
    if (is1)
    {
      for (i in seq_len(n_eq))
      {
        par[cuts_ind[[i]]] <- model1$cuts[[i]]
        par[coef_ind[[i]]] <- model1$coef[[i]]
        if (is_het[i])
        {
          par[coef_var_ind[[i]]] <- model1$coef_var[[i]]
        }
        if (is_marginal)
        {
          if (model1$other$marginal_par_n[i] > 0)
          {
            par[marginal_par_ind[[i]]] <- model1$par[
              model1$ind$marginal_par[[i]]]
          }
        }
      }
      if (n_sigma > 0)
      {
        par[sigma_ind] <- model1$par[model1$ind$sigma]
      }
    }
    if (is3)
    {
      for (i in seq_len(n_eq3 - 1))
      {
        par[coef3_ind[i, ]] <- model1$par[model1$ind$coef3[i, ]]
      }
      par[sigma3_ind] <- model1$par[model1$ind$sigma3]
    }
    
    # List to store the least squares models
    model2_list <- vector(mode = "list", length = n_eq2)
    
    # List to store the indexes of the observations in each regime
    ind_regime <- vector(mode = "list", length = n_eq2)
    
    # Two-step procedure for each equation
    for (v in seq_len(n_eq2))
    {
      # Conduct two-step procedure for each regime
      ind_regime[[v]] <- vector(mode = "list", length = n_regimes[v])
      for (i in seq_len(n_regimes[v]))
      {
        # Indexes of the observations corresponding to the regime
        for (j in seq_len(n_groups))
        {
          if (groups2[j, v] == (i - 1))
          {
            ind_regime[[v]][[i]] <- c(ind_regime[[v]][[i]], ind_g[[j]])
          }
        }
        
        # Regime specific data
        data_regime <- data[ind_regime[[v]][[i]], ]
        if (v >= 2)
        {
          # predictions of the previous equation
          for (v0 in seq_len(v - 1))
          {
            if (y_names[v0] %in% colnames(data_regime))
            {
              data_regime[, y_names[v0]] <- y_pred[ind_regime[[v]][[i]], v0]
            }
          }
        }
        
        # Estimate least squares regression
        model2_list[[v]][[i]] <- lm(formula2[[v]], data = data_regime, 
                                    na.action = na.exclude)
        
        # Save predictions of the regression
        y_pred[ind_regime[[v]][[i]], v] <- predict(model2_list[[v]][[i]])
        
        # Store the coefficients
        par[coef2_ind[[v]][i, ]] <- coef(model2_list[[v]][[i]])
      }
    }
  }
  out$twostep <- model2_list
  out$y_pred  <- y_pred
  
  # Store parameters into the variables
    # get the list of the parameters
  par_list <- par_msel(object = out, par = par)

    # assign them to the variables
  par             <- par_list$par
  coef            <- par_list$coef
  coef_var        <- par_list$coef_var
  cuts            <- par_list$cuts
  coef2           <- par_list$coef2
  sigma           <- par_list$sigma
  var2            <- par_list$var2
  cov2            <- par_list$cov2
  sigma2          <- par_list$sigma2
  marginal_par    <- par_list$marginal_par
  sigma_vec_ind   <- par_list$sigma_vec_ind
  coef3           <- par_list$coef3
  sigma3          <- par_list$sigma3
  
  # Store the variables in the output list
  out$par                 <- par
  out$coef                <- coef
  out$coef_var            <- coef_var
  out$cuts                <- cuts
  out$coef2               <- coef2
  out$sigma               <- sigma
  out$var2                <- var2
  out$cov2                <- cov2
  out$sigma2              <- sigma2
  out$marginal_par        <- marginal_par
  out$other$sigma_vec_ind <- sigma_vec_ind
  out$coef3               <- coef3
  out$sigma3              <- sigma3
  
  # Clear the memory
  par_list <- NULL

  # Estimate lambda for the sample selection models
  if (is1)
  {
    out$lambda <- predict(out, type = "lambda")
  }
  if (is3)
  {
    out$lambda_mn <- predict(out, type = "lambda_mn")
  }

  # Estimate the asymptotic covariance matrix
    # Prepare some values
  H     <- NULL                # Hessian
  H_inv <- NULL                # Inverse Hessian
  J     <- NULL                # Jacobian (or scores for MM)
  cov   <- diag(rep(1, n_par)) # Asymptotic covariance matrix
    # Manual covariance matrix
  if (is.matrix(cov_type))
  {
    cov <- cov_type
  }
    # Maximum-likelihood asymptotic covariance matrix estimator
  if (!is.matrix(cov_type) & (estimator == "ml"))
  {
    vcov_object <- vcov_ml(object = out, type = cov_type, 
                           n_cores = n_cores, n_sim = n_sim)
    cov <- vcov_object$vcov
    if (hasName(vcov_object, "H"))
    {
      H     <- vcov_object$H
      out$H <- H
    }
    if (hasName(vcov_object, "J"))
    {
      J     <- vcov_object$J
      out$J <- J
    }
  }
    # Two-step asymptotic covariance matrix estimator
  if (estimator == "2step" & (cov_type != "no"))
  {
    vcov_object <- vcov_2step(object  = out,
                              n_cores = n_cores, 
                              n_sim   = n_sim)
    cov         <- vcov_object$vcov
    out$J       <- vcov_object$scores
    out$H       <- vcov_object$scores_jac
  }
  out$cov_type <- cov_type
  out$cov      <- cov

  # Create tables to store the results in a format used by summary function
  tbl_list    <- tbl_msel(object = out)
  out$tbl     <- tbl_list$tbl
  out$se      <- tbl_list$se
  out$p_value <- tbl_list$p_value
  
  # Calculate the log-likelihood value
  logLik_val <- NA
  if (estimator == "ml")
  {
    if (length(regularization) == 0)
    {
      logLik_val <- opt$value
    }
    else
    {
      # because of the regularization
      logLik_val <- lnL_msel(par = par, control_lnL = control_lnL)
    }
  }
  out$logLik <- logLik_val

  return(out)
}