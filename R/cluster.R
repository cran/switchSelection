# Create the clusters
cluster_msel <- function(object, cluster = NULL)
{
  # Take the default cluster
  if (is.null(cluster))
  {
    cluster <- object$cluster
  }
  
  # Check for character
  if (is.character(cluster))
  {
    cluster_str <- "~"
    for (i in 1:length(cluster))
    {
      cluster_str <- paste0(cluster_str, "+", cluster[i])
    }
    cluster <- as.formula(cluster_str)
  }
  
  # Check whether clustering is required
  is_cluster <- is(object = cluster, class2 = "formula")
  
  # Set the default values
  n_cluster   <- object$other$n_obs
  ind_cluster <- 1:n_cluster
  
  # Transform formula to variable
  if (is_cluster)
  {
    data_cluster        <- model.frame(as.formula(cluster), object$data, 
                                       na.action = na.pass)
    ind_cluster         <- as.factor(do.call(paste0, data_cluster))
    n_cluster           <- length(levels(ind_cluster))
    levels(ind_cluster) <- 1:n_cluster
    ind_cluster         <- as.numeric(ind_cluster)
  }
  
  # Validate clustering
  if (is_cluster & any(is.na(ind_cluster)))
  {
    stop (paste0("There are undefined clusters. ",
                 "Please, check that variables in 'cluster' have ",
                 "no missing values."))
  }
  
  # Store to the output
  out <- list(ind_cluster = ind_cluster, 
              n_cluster   = n_cluster, 
              is_cluster  = is_cluster)
  
  # Return the results
  return (out)
}