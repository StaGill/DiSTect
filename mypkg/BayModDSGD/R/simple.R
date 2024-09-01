
#' @title Disease-Specific Gene Detection in Simple Scenario
#'
#' @param list_y a list of binary variables to indicate the spot status
#' @param matrix_x a martix of covaraites (Important: coordinates must be placed at the last columns of matrix)
#' @param num_samples number of spots
#' @param num_covariates number of genes
#'
#' @return a fitted stan model
#' @export
#'
#' @examples \donttest{}
dsgd_simple<-function(list_y,matrix_x,num_samples,num_covariates){

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
  N      = num_samples,
  P      = num_covariates,
  x      = matrix_x,
  y      = list_y
)
library(rstan)
autologistic_model<-stan_model(model_code = stan_program)
fit <- vb(autologistic_model,data = stan_data, algorithm = "fullrank",seed=128,tol_rel_obj=0.00001)
return(fit)
}

