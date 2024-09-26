# Optimization routine
opt_switchSelection <- function(opt_args = NULL, control_lnL,
                                n_sim = 1000, n_cores = 1,
                                type = "msel", start, 
                                opt_type = "optim", regularization = NULL)
{
  # Adjustment for Rcpp code
  if (!is.null(regularization))
  {
    if (hasName(x = regularization, name = "ridge_ind"))
    {
      regularization$ridge_ind <- regularization$ridge_ind - 1
    }
    if (hasName(x = regularization, name = "lasso_ind"))
    {
      regularization$lasso_ind <- regularization$lasso_ind - 1
    }
  }
  
  # Initialize variable to store optimization results
  opt <- NULL
  
  # Initialize list to store arguments
  optim_args <- list()
  if (hasName(opt_args, "optim"))
  {
    optim_args <- opt_args$optim
  }
    
  # obligatory arguments
  optim_args$par <- start
  optim_args$fn  <- lnL_msel
  optim_args$gr  <- grad_msel
    
  # technical arguments
  if (!hasName(optim_args, "method"))
  {
    optim_args$method <- "BFGS"
  }
  if (!hasName(optim_args, "hessian"))
  {
    optim_args$hessian <- FALSE
  }
  if (!hasName(optim_args, "control"))
  {
    optim_args$control <- list(maxit = 10000000,
                               fnscale = -1,
                               reltol = 1e-10,
                               abstol = 1e-10)
  }
    
  # likelihood arguments
  optim_args$n_sim <- n_sim
  optim_args$n_cores <- n_cores
  optim_args$control_lnL <- control_lnL
  optim_args$regularization <- regularization
  # Perform the optimization routine
  opt <- do.call(what = optim, args = optim_args)
    
  # Add genetic optimization if need
  if (opt_type == "gena")
  {
    # Initialize list to store arguments
    gena_args <- list()
    if (hasName(opt_args, "gena"))
    {
      gena_args <- opt_args$gena
    }
      
    # obligatory arguments
    gena_args$fn <- lnL_msel
    gena_args$gr <- grad_msel
      
    # technical arguments
    if (!hasName(gena_args, "pop.initial"))
    {
      gena_args$pop.initial <- opt$par
    }
    if (!hasName(gena_args, "mutation.method"))
    {
      gena_args$mutation.method <- "percent"
    }
    if (!hasName(gena_args, "maxiter"))
    {
      gena_args$maxiter <- 100
    }
    if (!hasName(gena_args, "info"))
    {
      gena_args$info <- TRUE
    }
    if (!hasName(gena_args, "lower"))
    {
      gena_args$lower <- -2 * abs(opt$par)
    }
    if (!hasName(gena_args, "upper"))
    {
      gena_args$upper <- -2 * abs(opt$par)
    }
    if (!hasName(gena_args, "hybrid.prob"))
    {
      gena_args$hybrid.prob <- 0.2
    }
      
    # likelihood arguments
    gena_args$n_sim <- n_sim
    gena_args$n_cores <- n_cores
    gena_args$control_lnL <- control_lnL
    gena_args$regularization <- regularization
      
    # Perform the optimization routine
    opt <- do.call(what = gena::gena, args = gena_args)
  }
    
  # Add particle swarm optimization if need
  if (opt_type == "pso")
  {
    # Initialize list to store arguments
    pso_args <- list()
    if (hasName(opt_args, "pso"))
    {
      pso_args <- opt_args$pso
    }
      
    # obligatory arguments
    if (type == "msel")
    pso_args$fn <- lnL_msel
    pso_args$gr <- grad_msel
      
    # technical arguments
    if (!hasName(pso_args, "pop.initial"))
    {
      pso_args$pop.initial <- opt$par
    }
    if (!hasName(pso_args, "maxiter"))
    {
      pso_args$maxiter <- 100
    }
    if (!hasName(pso_args, "info"))
    {
      pso_args$info <- TRUE
    }
    if (!hasName(pso_args, "lower"))
    {
      pso_args$lower <- -2 * abs(opt$par)
    }
    if (!hasName(pso_args, "upper"))
    {
      pso_args$upper <- 2 * abs(opt$par)
    }
    if (!hasName(pso_args, "hybrid.prob"))
    {
      pso_args$hybrid.prob <- 0
    }
      
    # likelihood arguments
    pso_args$n_sim <- n_sim
    pso_args$n_cores <- n_cores
    pso_args$control_lnL <- control_lnL
    pso_args$regularization <- regularization
      
    # Perform the optimization routine
    opt <- do.call(what = gena::pso, args = pso_args)
  }
    
  return (opt)
}