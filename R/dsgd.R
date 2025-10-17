#' Disease-Specific Gene Detection (Wrapper Function)
#' Fitting the model with the option of single or multiple correlated tissues and the option of adding interaction terms.
#' 
#' @param list_y a list of binary variables to indicate the spot status
#' @param matrix_x a matrix of covariates, with the coordinates placed in the final two columns of the matrix
#' @param interaction optional, character vector specifying the interaction terms
#' 
#' @return a fitted stan mode
#' @export
#' @importFrom rstan stan_model
#' @importFrom rstan vb
#'
#' @examples \donttest{}

#' @export
dsgd_single<-function(list_y, matrix_x, method="VI", niter=2000, nwarmup=1000, nchain=3){
  
  tmp_program<-"
data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N,P+2] x;
  int y[N];
}
parameters {
  vector[P] beta;
  real<lower=0,upper=8> eta;
  vector <lower=0>[P] beta_gamma;
  real <lower=0,upper=1> w;
}

model{
vector[N] mu;
matrix[N,N] prob_neigh;

for(i in 1:P){
  w ~ uniform(0,1);
  beta_gamma[i] ~ inv_gamma(5,50);
  target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]), log(w)+normal_lpdf(beta[i]|0,beta_gamma[i]));
}
for(i in 1:N) {
  for(j in 1:N){
    eta ~ uniform(0,8);
    if (j != i && sqrt(square(x[i,P+1]-x[j,P+1])+square(x[i,P+2]-x[j,P+2]))<=1){
      prob_neigh[i,j]=eta*y[j];
    }else{
      prob_neigh[i,j]=0;
    }}
  mu[i] = dot_product(x[i,1:P],beta)+sum(prob_neigh[i,]);

}
y ~ bernoulli_logit(mu);
}"



stan_program <- tmp_program

stan_data <- list(
  N      = nrow(matrix_x),
  P      = ncol(matrix_x)-2,
  x      = matrix_x,
  y      = list_y
)
library(rstan)
autologistic_model<-stan_model(model_code = stan_program)

if (method == "VI"){
  fit <- vb(autologistic_model,data = stan_data, algorithm = "fullrank", seed = 128, tol_rel_obj = 0.00001)
  } else if (method == "NUTS") {
    fit <- sampling(
      object = autologistic_model,
      data = stan_data,
      iter = niter,
      warmup = nwarmup,
      chains = nchain,
      seed = seed)
}

return(fit)
}


#' @export
dsgd_multiple<-function(list_y, matrix_x, label_list, method="VI", niter=2000, nwarmup=1000, nchain=3){
  tmp_program<-"
data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> S;
  matrix[N,P+2] x;
  int y[N];
  vector[S] u_mu;
  int label[N];
}
parameters {
  vector[P] beta;
  real<lower=0,upper=8> eta;
  vector<lower=0>[P] beta_gamma;
  real<lower=0,upper=1> w;
  vector<lower=-1,upper=1>[S] U;
  real<lower=0,upper=1> sigma_square;
  real<lower=0,upper=1> rho;
}
model{

vector[N] mu;
matrix[N,N] prob_neigh;
matrix[S,S] covar_matrix;
rho ~ uniform(0,1);
sigma_square ~ inv_gamma(5,50);
for (i in 1:S){
  for (j in 1:S){
    if (i == j){
      covar_matrix[i,j]=sigma_square;
    }else{
      covar_matrix[i,j]=sigma_square*rho;
    }
  }
}

U ~ multi_normal(u_mu, covar_matrix);


for(i in 1:P){
  w ~ uniform(0,1);
  beta_gamma[i] ~ inv_gamma(5,50);
  target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]), log(w)+normal_lpdf(beta[i]|0,beta_gamma[i]));
}
for(i in 1:N) {
  for(j in 1:N){
    eta ~ uniform(0,8);
    if (label[i]==label[j] && j != i && sqrt(square(x[i,P+1]-x[j,P+1])+square(x[i,P+2]-x[j,P+2]))<=1){
      prob_neigh[i,j]=eta*y[j];
    }else{
      prob_neigh[i,j]=0;
    }}
  mu[i] = dot_product(x[i,1:P],beta)+sum(prob_neigh[i,])+U[label[i]];

}
y ~ bernoulli_logit(mu);
}"

stan_program <- tmp_program

stan_data <- list(
  N      = nrow(matrix_x),
  P      = ncol(matrix_x)-2,
  x      = matrix_x,
  y      = list_y,
  label  = label_list,
  S      = length(unique(label_list)),
  u_mu=rep(0,length(unique(label_list)))
)
autologistic_model<-stan_model(model_code = stan_program)

if (method=="VI"){
  fit <- vb(autologistic_model,data = stan_data, algorithm = "fullrank", seed = 128, tol_rel_obj = 0.00001)
  } else if (method=="NUTS") {
    fit <- sampling(
      object = autologistic_model,
      data = stan_data,
      iter = niter,
      warmup = nwarmup,
      chains = nchain,
      seed = seed)
    }

return(fit)
}

#' @export
add_interaction <- function(matrix_x, interaction) {
  x_df <- as.data.frame(matrix_x)
  coord_cols <- tail(names(x_df), 2)
  x_df_no_coords <- x_df[, 1:(ncol(x_df) - 2)]
  
  if (!all(interaction %in% names(x_df_no_coords))) {
    stop("Error: Some variables do not exist.")
  }
  
  # Create all unique two-way combinations
  interaction_pairs <- combn(interaction, 2, simplify = FALSE)
  
  # For each interaction pair, create a new column
  for (pair in interaction_pairs) {
    new_col <- x_df_no_coords[[pair[1]]] * x_df_no_coords[[pair[2]]]
    col_name <- paste(pair, collapse = "*")
    x_df_no_coords[[col_name]] <- new_col
  }
  
  x_df_final <- cbind(x_df_no_coords, x_df[, coord_cols])
  return(as.matrix(x_df_final))
}


#' @export
dsgd <- function(list_y, matrix_x, label_list = NULL, interaction = NULL, method="VI", niter=2000, nwarmup=1000, nchain=3) {
  
  # Add interaction terms if specified
  if (!is.null(interaction)) {
    matrix_x <- add_interaction(matrix_x, interaction)
  }
  
  # Check the validity of method argument
  if (method !="VI" & method != "NUTS"){
    stop("Error: method must be either 'VI' or 'NUTS'")
  }
  
  # Choose model based on label_list presence
  if (is.null(label_list)) {
    return(dsgd_single(list_y, matrix_x, method="VI", niter=2000, nwarmup=1000, nchain=3))
  } else {
    return(dsgd_multiple(list_y, matrix_x, label_list, method="VI", niter=2000, nwarmup=1000, nchain=3))
  }
}


