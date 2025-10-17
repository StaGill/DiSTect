#' Title
#'
#' @param list_y a list of binary variables to indicate the spot status
#' @param coords a matrix with two columns to indicate the coordinate information for each spot/cell
#'
#' @return imputed outcome dataset
#' @export
#' @importFrom rjags jags.model
#' @importFrom rjags coda.samples
#' @importFrom sp spDistsN1
#' @examples \donttest{}


missing_imputation <- function(list_y, coords) {
  library(sp)      # for spDistsN1()
  library(rjags)   # for running JAGS
  library(coda)    # for processing MCMC samples
  missing_index <- which(is.na(list_y))
  obs_idx <- c(1:length(list_y))

  if (length(missing_index) == 0) {
    stop("There are no missing values to impute.")
  }



  neighbors_list <- lapply(missing_index, function(miss_i) {

    dists <- spDistsN1(coords[obs_idx, ], coords[miss_i, ], longlat = FALSE)

    neighbor_idx <- which(dists <= 1)
    neighbor_idx <- neighbor_idx[neighbor_idx != miss_i]
    return(neighbor_idx)
  }) # find neighbors for each missing spot


  split_neighbors <- function(neighbor_index,missing_index) {
    new_list1 <- list()
    new_list2 <- list()

    for (i in seq_along(neighbor_index)) {
      neighbors <- neighbor_index[[i]]

      if (any(missing_index %in% neighbors)) {
        new_list1 <- append(new_list1, list(setdiff(neighbors,missing_index)))
        new_list2 <- append(new_list2, intersect(neighbors,missing_index))
      } else {
        new_list1 <- append(new_list1, list(neighbors))
        new_list2 <- append(new_list2, 0)
      }
    }

    return(list("Without_NA" = new_list1, "With_NA" = new_list2))
  }

  list_split<-split_neighbors(neighbors_list,missing_index)

  replace_with_positions <- function(vec, ref_vec) {
    match(vec, ref_vec)
  }
  position_NA <- lapply(list_split$With_NA, replace_with_positions, ref_vec=missing_index)



  N_miss <- length(missing_index)
  n_neighbors <- sapply(neighbors_list, length)
  n_neighbors_withoutNA<-sapply(list_split$Without_NA,length)
  n_neighbors_withNA<-sapply(position_NA,length)

  if(any(n_neighbors == 0)){
    stop("One or more missing observations have no neighbors within the specified distance.")
  }

  max_neighbors_withoutNA <- max(sapply(list_split$Without_NA,length))
  neighbor_matrix1 <- matrix(NA, nrow = N_miss, ncol = max_neighbors_withoutNA )

  max_neighbors_withNA <- max(sapply(list_split$With_NA,length))
  neighbor_matrix2<- matrix(NA, nrow = N_miss, ncol = max_neighbors_withNA )







  model_string <- "
  model {
    for(j in 1:N_miss) {
      y_missing[j] ~ dbern((mean(y_obs)))

    }
      }


  "


  data_list <- list(
      N_miss=N_miss,
      y_obs=list_y[-missing_index]
  )

  init_values <- function() {
    list(y_missing = rep(0, N_miss))
  }

  jags_mod <- jags.model(file = textConnection(model_string),
                         data = data_list,
                         inits = init_values,
                         n.chains = 3,
                         n.adapt = 1000)

  update(jags_mod, 1000)

  samples <- coda.samples(jags_mod, variable.names = c("y_missing"), n.iter = 5000)
  list_y[missing_index]<-colMeans(as.matrix(samples))
  return(list_y)
}



