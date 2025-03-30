#' Negative Binomial Log-Likelihood
#'
#' Computes the negative log-likelihood for the Negative Binomial distribution.
#'
#' @param x Vector of observed counts.
#' @param mu Vector of expected means.
#' @param disp Dispersion parameter.
#'
#' @return Negative log-likelihood value.
#' @keywords internal
nb_lik = function (x, mu, disp)
{
  return(-sum(dnbinom(x = x, size = disp, mu = mu, log = TRUE)))
}


#' Gaussian Log-Likelihood
#'
#' Computes the negative log-likelihood for the Gaussian distribution.
#'
#' @param x Vector of observed values.
#' @param mu Vector of predicted means.
#' @param sigma Standard deviation.
#'
#' @return Negative log-likelihood value.
#' @keywords internal
gaussian_lik = function (x, mu, sigma)
{
  return(-sum(dnorm(x = x, sd = sigma, mean = mu, log = TRUE)))
}



#' Poisson Log-Likelihood
#'
#' Computes the negative log-likelihood for the Poisson distribution.
#'
#' @param x Vector of observed counts.
#' @param mu Vector of predicted means (Poisson rates).
#'
#' @return Negative log-likelihood value.
#' @keywords internal
poisson_lik = function (x, mu)
{
  return(-sum(dpois(x = x, lambda = mu, log = TRUE)))
}


#' Binomial Log-Likelihood
#'
#' Computes the negative log-likelihood for the Binomial distribution.
#'
#' @param x Vector of observed successes.
#' @param mu Vector of predicted probabilities.
#' @param num_obs Number of trials for each observation.
#'
#' @return Negative log-likelihood value.
#' @keywords internal
binomial_lik = function (x, mu, num_obs)
{
  return(-sum(dbinom(x = x, size = num_obs,prob = mu,log = TRUE)))
}




#' Custom GLM Object: Spot Poisson
#'
#' Defines a Poisson-based GLM family object tailored for modeling spatial transcriptomics
#' count data using deconvolution. Includes link functions, marginal mean, likelihood,
#' gradient computation, and prediction utilities.
#'
#' @format A list of class-specific functions:
#' \describe{
#'   \item{\code{link}}{Link and inverse link functions with derivatives.}
#'   \item{\code{marginal_mu}}{Marginal mean function.}
#'   \item{\code{original_class}}{Base GLM family used (Poisson).}
#'   \item{\code{predict}}{Predict means from \code{X}, \code{beta}, and \code{lambda}.}
#'   \item{\code{likelihood}}{Negative log-likelihood function.}
#'   \item{\code{grad}}{Gradient function with respect to beta.}
#' }
#'
#' @export
#spot poisson class
spot_poisson = vector("list",3)
names(spot_poisson) = c("link","marginal_mu","original_class")

#link function
spot_poisson[["link"]] <- function(C,lambda,offset = rep(0,length(C))) {
  ## link y is eta = linkfun(E(Y|X))
  linkfun <- function(mu){
    A = mu-C*exp(offset)
    A = pmax(A,1e-3)
    B = lambda*exp(offset)
    if(mean(is.na(log(A/B))) > 0){
      #print("errors")
    }
    return(log(A/B))
  }

  ## inverse link
  linkinv <- function(eta){
    eta <- pmax(eta, -20 - offset)
    eta = pmin(eta,-1)
    #check = linkfun((lambda*exp(eta)+C)*exp(offset))
    return ((lambda*exp(eta)+C)*exp(offset))
  }

  ## derivative of dmu/deta
  mu.eta <- function(eta) {
    eta <- pmax(eta, -10)
    mu = exp(eta)
    return ((lambda*mu)*exp(offset))
  }
  valideta <- function(eta) TRUE
  link <- "spot_poisson"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),
            class = "link-glm")
}

#marginal mean function
spot_poisson[['marginal_mu']] = function(eta){
  return(exp(eta))
}

#original class
spot_poisson[["original_class"]] = stats:::poisson


#predict from X and beta, predict y
spot_poisson[["predict"]] = function(X,beta,lambda,offset = rep(0,nrow(X))){
  #Calculate eta for all other cell types
  eta_rest = X%*%beta + offset
  #Calculate mu for all other cell types
  mu_rest = exp(eta_rest)
  #Calculate C for all spots
  C = rowSums(lambda*mu_rest)
  C[is.na(C)] = 1e-20
  C[C == 0] = 1e-20
  #return expectation
  return(list(total = C, individual = mu_rest))
}

spot_poisson[["likelihood"]] = poisson_lik


spot_poisson[["grad"]] = function(X,y,beta,lambda,offset = rep(0,length(y))){
  means = spot_poisson[["predict"]](X,beta,lambda,offset)
  #get cell type specifici expression means
  CT_means = means$individual
  #get total mean spot expression
  total_means = means$total
  #get number of cell types
  nCT = ncol(lambda)
  #initialize gradient
  grad_beta = beta * 0
  # Compute weights matrix for all observations and components
  weights = -lambda * CT_means + y * (1 / total_means) * lambda * CT_means  # Vectorized computation for weights
  # Batch computation of gradients
  grad_beta = t(X) %*% weights  # Computes all gradients at once
  return(list(grad = grad_beta,weights = weights))
}







#' Custom GLM Object: Spot Gaussian
#'
#' Defines a Gaussian-based GLM family object for spatial transcriptomics data,
#' incorporating link functions and gradient-based optimization for continuous-valued
#' expression deconvolution.
#'
#' @format A list of class-specific functions:
#' \describe{
#'   \item{\code{link}}{Link and inverse link functions with derivatives.}
#'   \item{\code{marginal_mu}}{Marginal mean function.}
#'   \item{\code{original_class}}{Base GLM family used (Gaussian).}
#'   \item{\code{predict}}{Predict means from \code{X}, \code{beta}, and \code{lambda}.}
#'   \item{\code{likelihood}}{Negative log-likelihood function.}
#' }
#'
#' @export
#spotgaussian  class
spot_gaussian = vector("list",3)
names(spot_gaussian) = c("link","marginal_mu","original_class")

#link function
spot_gaussian[["link"]] <- function(C,lambda) {
  ## link y is E(Y|X)
  linkfun <- function(mu){
    A = mu-C
    B = lambda
    return(A/B)
  }
  ## inverse link
  linkinv <- function(eta){
    return (lambda*eta + C)
  }
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) {
    return ((lambda))
  }
  valideta <- function(eta) TRUE
  link <- "spot_gaussian"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),
            class = "link-glm")
}

#marginal mean
spot_gaussian[['marginal_mu']] = function(eta){
  return(eta)
}

#link function
spot_gaussian[["original_class"]] = stats:::gaussian


#predict mean from X beta and lambda
spot_gaussian[["predict"]] = function(X,beta,lambda){
  #Calculate eta for all other cell types
  eta_rest = X%*%beta
  #Calculate mu for all other cell types
  mu_rest = eta_rest
  #Calculate C for all spots
  C = rowSums(lambda*mu_rest)
  #return expectation
  return(list(total = C, individual = mu_rest))
}

#likelihood function
spot_gaussian[["likelihood"]] = gaussian_lik





#' Custom GLM Object: Spot Binomial
#'
#' Defines a Binomial-based GLM family object for modeling bounded expression data
#' or proportions in spatial transcriptomics. Supports deconvolution and optimization
#' with logistic link functions.
#'
#' @format A list of class-specific functions:
#' \describe{
#'   \item{\code{link}}{Logistic-based link and inverse link functions with derivatives.}
#'   \item{\code{marginal_mu}}{Marginal mean function (inverse logit).}
#'   \item{\code{original_class}}{Base GLM family used (Binomial).}
#'   \item{\code{predict}}{Predict probabilities from \code{X}, \code{beta}, and \code{lambda}.}
#'   \item{\code{likelihood}}{Negative log-likelihood function.}
#'   \item{\code{grad}}{Gradient function with respect to beta.}
#' }
#'
#' @export
#spot binomial class
spot_binomial = vector("list",3)
names(spot_binomial) = c("link","marginal_mu","original_class")

#link function
spot_binomial[["link"]] <- function(C,lambda) {
  ## link y is eta = linkfun(E(Y|X))
  linkfun <- function(mu){
    #print(mu)
    A = mu-C
    B = lambda
    return(LaplacesDemon::logit(A/B))
  }

  ## inverse link
  linkinv <- function(eta){
    eta = pmax(eta,-10)
    return (lambda*LaplacesDemon::invlogit(eta)+C)
  }

  ## derivative of dmu/deta
  mu.eta <- function(eta) {
    A = 1/(exp(eta) + 1) - 1/(exp(eta) + 1)^2
    return ((lambda*A))
  }

  valideta <- function(eta) TRUE
  link <- "spot_binomial"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),
            class = "link-glm")
}

#marginal mean
spot_binomial[['marginal_mu']] = function(eta){
  return(LaplacesDemon::invlogit(eta))
}

#original class
spot_binomial[["original_class"]] = stats:::binomial

#get predicted p from X beta and lambda
spot_binomial[["predict"]] = function(X,beta,lambda){
  #Calculate eta for all other cell types
  eta_rest = X%*%beta
  #Calculate mu for all other cell types
  mu_rest = LaplacesDemon::invlogit(eta_rest)
  #Calculate C for all spots
  C = rowSums(lambda*mu_rest)
  return(list(total = C, individual = mu_rest))
}

#likelihood function
spot_binomial[["likelihood"]] = binomial_lik




#gradient
spot_binomial[["grad"]] = function(X,y,beta,lambda,weights){
  means = spot_binomial[["predict"]](X,beta,lambda)
  #get cell type specifici expression means
  CT_means = means$individual
  CT_means = pmax(CT_means,0.000001)
  CT_means = pmin(CT_means,0.999999)
  #get total mean spot expression
  total_means = means$total
  total_means = pmax(total_means,0.000001)
  total_means = pmin(total_means,0.999999)
  #get number of cell types
  nCT = ncol(lambda)
  #initialize gradient
  grad_beta = beta * 0
  # Precompute invariant terms
  common_term = (y - weights * total_means) * (1 / total_means) * (1 / (1 - total_means))
  # Compute weights matrix for all observations and components
  weights = common_term * lambda * CT_means * (1 - CT_means)  # Produces a matrix of weights
  # Batch computation of gradients
  grad_beta = t(X) %*% weights  # This computes gradients for all j in one operation

  return(list(grad = grad_beta,weights = weights))
}






#' Custom GLM Object: Spot Negative Binomial
#'
#' Defines a Negative Binomial-based GLM family object suitable for overdispersed
#' spatial transcriptomics count data. Incorporates flexible link functions and
#' gradient computations for both regression coefficients and dispersion parameters.
#'
#' @format A list of class-specific functions:
#' \describe{
#'   \item{\code{link}}{Negative binomial link and inverse link functions with derivatives.}
#'   \item{\code{marginal_mu}}{Marginal mean function.}
#'   \item{\code{original_class}}{Base GLM family used (Poisson or Negative Binomial).}
#'   \item{\code{predict}}{Predict means from \code{X}, \code{beta}, and \code{lambda}.}
#'   \item{\code{likelihood}}{Negative log-likelihood function.}
#'   \item{\code{grad}}{Gradient function with respect to beta and dispersion.}
#' }
#'
#' @export

#spot negative binomial class
spot_negative_binomial = vector("list",3)
names(spot_negative_binomial) = c("link","marginal_mu","original_class")

#link function
spot_negative_binomial[["link"]] <- function(C,lambda,gamma,offset = rep(0,length(C))) {
  base = MASS::negative.binomial(theta = gamma)
  #edit link parameters
  base$linkfun = function(mu){
    A = mu-C*exp(offset)
    B = lambda*exp(offset)
    if(mean(is.na(log(A/B))) > 0){
    }
    return(log(A/B))
  }
  #edit inverse link
  base$linkinv = function(eta){
    return ((lambda*exp(eta)+C)*exp(offset))
  }
  #edit dmu/deta
  base$mu.eta = function(eta){
    eta <- pmax(eta, -10)
    mu = exp(eta)
    return ((lambda*mu)*exp(offset))
  }
  return (base)
}


#marginal mean
spot_negative_binomial[['marginal_mu']] = function(eta){
  return(exp(eta))
}

#original class
#spot_negative_binomial[["original_class"]] = MASS::negative.binomial
spot_negative_binomial[["original_class"]] = stats:::poisson


#predict y from X beta and lambda
spot_negative_binomial[["predict"]] = function(X,beta,lambda,offset = rep(0,nrow(X))){
  #Calculate eta for all other cell types
  eta_rest = X%*%beta + offset
  #Calculate mu for all other cell types
  mu_rest = exp(eta_rest)
  #Calculate C for all spots
  C = rowSums(lambda*mu_rest)
  #return expectation
  return(list(total = C, individual = mu_rest))
}

#likelihood
spot_negative_binomial[["likelihood"]] = nb_lik



#gradient
spot_negative_binomial[["grad"]] =
  function(X,y,beta,lambda,disp,offset = rep(0,length(y))){
  means = spot_negative_binomial[["predict"]](X,beta,lambda,offset)
  #get cell type specifici expression means
  CT_means = means$individual
  #get total mean spot expression
  total_means = means$total
  #get number of cell types
  nCT = ncol(lambda)
  #initialize gradient
  grad_beta = beta * 0
  #get offset value
  A = exp(offset)
  # Precompute invariant terms
  common_term = (y - total_means) * (1 / total_means) * disp / (total_means + disp)
  # Element-wise product across lambda and CT_means
  weights = common_term * lambda * CT_means  # Produces a matrix of weights
  # Matrix multiplication for gradient computation
  grad_beta = t(X) %*% weights  # This computes gradients for all j in a single operation
  #compute dispersion gradient
  disp_grad = sum(digamma(total_means+disp) - digamma(disp) + (total_means-y)/(total_means + disp) + log(disp/(disp + total_means)))
  return(list(grad = grad_beta,disp_grad = disp_grad, weights = weights))
}








