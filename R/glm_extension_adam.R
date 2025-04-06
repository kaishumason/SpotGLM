#' Fit a Spatial GLM with Deconvolution
#'
#' Fits a generalized linear model (GLM) with deconvolution for spatial transcriptomics data,
#' supporting Poisson, Gaussian, Negative Binomial, and Binomial families. Optimization is
#' performed via gradient descent (or closed-form for Gaussian).
#'
#' @param X Design matrix of covariates (spots × covariates).
#' @param y Response vector (e.g., gene expression).
#' @param lambda Deconvolution matrix (spots × cell types).
#' @param family GLM family: one of `"spot poisson"`, `"spot gaussian"`, `"spot binomial"`, `"spot negative binomial"`.
#' @param beta_0 Initial coefficients (matrix: covariates × cell types).
#' @param offset Optional offset term (numeric vector).
#' @param weights Optional observation-level weights.
#' @param fix_coef Logical matrix indicating which coefficients to fix.
#' @param learning_rate Initial learning rate for gradient descent.
#' @param n_epochs Number of epochs (iterations over the full data).
#' @param batch_size Size of mini-batches for gradient descent.
#' @param max_diff Convergence threshold on likelihood improvement ratio.
#' @param improvement_threshold Minimum change in ratio between epochs to signal progress.
#' @param max_conv Number of consecutive small-improvement epochs to trigger convergence.
#'
#' @return A list with:
#' \describe{
#'   \item{beta}{Estimated coefficient matrix.}
#'   \item{vcov}{Variance-covariance matrix of estimates.}
#'   \item{dispersion}{Estimated dispersion (for NB models).}
#'   \item{likelihood}{Final negative log-likelihood.}
#'   \item{converged}{Logical indicating if convergence was reached.}
#'   \item{num_epoch}{Number of epochs used.}
#'   \item{fixed_coef}{Final fix coefficient matrix.}
#' }
#'
#' @keywords internal
spot_glm = function(
    X, y, lambda, family,
    beta_0 = matrix(0, ncol(X), ncol(lambda)),
    offset = rep(0, nrow(X)),
    weights = NULL,
    fix_coef = matrix(FALSE, ncol(X), ncol(lambda)),
    learning_rate = 100,
    n_epochs = 50,            # <-- NEW parameter: # of epochs
    batch_size = 128,         # <-- mini-batch size
    max_diff = 1 - 1e-6,       # threshold for relative improvement
    improvement_threshold = 1e-6,
    max_conv = 10
) {
  # ----------------------------------
  # 0. Checks and setup
  # ----------------------------------
  if(ncol(X) != nrow(beta_0)){
    stop("#Rows of initial beta does not match #cols of covariate matrix X")
  }
  if(ncol(lambda) != ncol(beta_0)){
    stop("#Cols of initial beta does not match #cols of deconvolution matrix lambda")
  }
  if(nrow(lambda) != nrow(X)){
    stop("#Rows of deconvolution matrix lambda and covariate matrix X must match")
  }
  
  if(is.null(fix_coef)){
    print("Fix coefficients matrix is not supplied. No coefficients will be fixed")
    fix_coef = matrix(FALSE,ncol(X),ncol(lambda))
  } else if(dim(fix_coef)[1] != dim(beta_0)[1] & dim(fix_coef)[2] != dim(beta_0)[2]){
    print("Fix coefficients matrix is of incorrect dimension. No coefficients will be fixed")
    fix_coef = matrix(FALSE,ncol(X),ncol(lambda))
  }
  
  if(mean(is.na(y)) != 0){
    stop("Y has na entries. Remove prior to running spotGLM")
  }
  if(mean(is.na(X)) != 0){
    stop("X has na entries. Remove prior to running spotGLM")
  }
  if(mean(is.na(lambda)) != 0){
    stop("Lambda has na entries. Remove prior to running spotGLM")
  }
  if(mean(is.na(beta_0)) != 0){
    stop("Initial beta has na entries. Remove prior to running spotGLM")
  }
  
  if(is.null(weights)){
    weights = rep(1, length(y))
  } else if(length(weights) != length(y)){
    stop("Weight vector must be NULL or same length as response vector.")
  }
  
  # Choose family model
  if(family == "spot poisson"){
    family_model = spot_poisson
  } else if (family == "spot gaussian"){
    family_model = spot_gaussian
  } else if (family == "spot binomial"){
    family_model = spot_binomial
  } else if (family == "spot negative binomial"){
    family_model = spot_negative_binomial
  } else if (family == "spot gaussian power"){
    family_model = spot_gaussian
  } else {
    stop("Family not found. Family should be one of gaussian,gaussian_power,binomial,poisson,or negative binomial")
  }
  
  if((sum(offset != 0) > 0) & !(family %in% c("spot poisson","spot negative binomial"))){
    stop("Offsets can only be used with spot poisson and spot negative binomial models")
  }
  
  # ----------------------------------
  # 1. Basic initialization
  # ----------------------------------
  nCT = ncol(lambda)        # # of cell types
  p   = ncol(X)             # # of covariates
  beta_curr = beta_0
  
  # Add col/rownames for clarity
  if(!is.null(colnames(lambda))){
    colnames(beta_curr) = colnames(lambda)
  }
  if(!is.null(colnames(X))){
    rownames(beta_curr) = colnames(X)
  }
  
  # "Good" betas or spots
  good_beta = 0 * beta_0
  good_spots_ct = vector("list", nCT)
  for(j in seq_len(nCT)){
    if(is.null(weights)){
      good_spots_ct[[j]] = which(lambda[, j] * exp(offset) > 1e-3)
    } else {
      good_spots_ct[[j]] = which(lambda[, j] * exp(offset) > 1e-3 & (weights != 0))
    }
  }
  # If fix_coef is TRUE, that means don't update. So good_beta = 1 - fix_coef
  good_beta = matrix(as.numeric(1 - fix_coef), nrow(fix_coef), ncol(fix_coef))
  
  # For Gaussian, we just solve closed form and return
  if(family == "spot gaussian"){
    model_outcome = spot_glm_gaussian(X, y, lambda, fix_coef, beta_0)
    beta_curr     = model_outcome$beta
    V             = model_outcome$vcov
    sigma.sq      = model_outcome$sigma.sq
    fitted_vals   = family_model[["predict"]](X, beta_curr, lambda)$total
    lik           = -gaussian_lik(x = y, mu = fitted_vals, sigma = sqrt(sigma.sq))
    return(list(
      beta       = beta_curr,
      vcov       = V,
      dispersion = sigma.sq,
      likelihood = lik,
      R          = 0,
      converged  = TRUE
    ))
  }
  
  # Negative binomial => estimate initial dispersion
  if(family == "spot negative binomial"){
    fitted_vals = family_model[["predict"]](X, beta_0, lambda, offset)$total
    A = optimize(nb_lik, x = y, mu = fitted_vals, lower = 0.05, upper = 100)
    dispersion = A$minimum
  } else {
    dispersion = 1e8
  }
  
  # ----------------------------------
  # 2. Setup for mini-batch gradient descent
  # ----------------------------------
  beta_new = beta_0
  N        = length(y)
  vec      = seq_len(N)     # all observation indices
  
  # Adam parameters
  m_beta = 0
  v_beta = 0
  m_disp = 0
  v_disp = 0
  beta_1 = 0.9
  beta_2 = 0.999
  epsilon = 1e-8
  
  t = 0  # Adam timestep
  
  # Compute initial full-data likelihood
  if(family == "spot poisson"){
    pred = family_model[["predict"]](X, beta_0, lambda, offset)$total
    lik_old = -poisson_lik(y, pred)
  } else if(family == "spot negative binomial"){
    pred = family_model[["predict"]](X, beta_0, lambda, offset)$total
    lik_old = -nb_lik(y, pred, dispersion)
  } else if(family == "spot binomial"){
    pred = family_model[["predict"]](X, beta_0, lambda)$total
    lik_old = -binomial_lik(y, pred, weights)
  }
  
  converged = FALSE
  conv_counter = 0
  lik_ratio = 0
  
  # --------------------------------------------------
  # 3. Outer loop: run up to n_epochs
  #    Each epoch: shuffle data, do mini-batch updates,
  #    then compute full-data likelihood to check
  # --------------------------------------------------
  for(epoch in seq_len(n_epochs)) {
    #print(epoch)
    # Shuffle the data indices for this epoch
    shuffled_indices = sample(vec)
    
    # Split into contiguous mini-batches of size batch_size
    batch_starts = seq(1, N, by = batch_size)
    
    for(start_idx in batch_starts) {
      end_idx = min(start_idx + batch_size - 1, N)
      spots   = shuffled_indices[start_idx:end_idx]
      
      # -----------------------------------------
      # (A) Compute gradient on these 'spots'
      # -----------------------------------------
      recompute_gradients = TRUE
      G = gradient_descent_update(
        family,
        X[spots,,drop = F],
        y[spots],
        beta_new,
        lambda[spots,,drop = F],
        offset[spots],
        disp = dispersion,
        weights = weights[spots]
      )
      while(recompute_gradients) {
        # Scale by batch size to keep consistent step
        grad      = (1 / length(spots)) * G$beta_grad
        disp_grad = (1 / length(spots)) * G$disp_grad
        
        # Fix coefficients that should not update
        grad[is.na(grad)] = 0
        grad[fix_coef == TRUE] = 0
        
        # Adam moment updates
        t_temp        = t + 1
        m_beta_temp   = beta_1 * m_beta + (1 - beta_1) * grad
        v_beta_temp   = beta_2 * v_beta + (1 - beta_2) * (grad^2)
        m_beta_hat    = m_beta_temp / (1 - beta_1^t_temp)
        v_beta_hat    = v_beta_temp / (1 - beta_2^t_temp)
        
        grad_update   = learning_rate * m_beta_hat / (sqrt(v_beta_hat) + epsilon)
        beta_new_temp = beta_new + grad_update
        
        #if(family == "spot negative binomial"){
        #m_disp_temp  = beta_1 * m_disp + (1 - beta_1) * disp_grad
        #v_disp_temp  = beta_2 * v_disp + (1 - beta_2) * (disp_grad^2)
        #m_disp_hat   = m_disp_temp / (1 - beta_1^t_temp)
        #v_disp_hat   = v_disp_temp / (1 - beta_2^t_temp)
        #disp_update  = learning_rate * m_disp_hat / (sqrt(v_disp_hat) + epsilon)
        #dispersion_temp = dispersion + disp_update
        #print(dispersion_temp)
        #dispersion_temp = max(dispersion_temp, 0.05)
        #} else {
        #dispersion_temp = dispersion
        #}
        
        # Here we do a line-search style step check using the FULL data
        # If you'd rather skip re-checking on full data every mini-batch
        # (for performance), you can just accept the step or limit check
        # with the mini-batch. But we'll do as your code does:
        
        if (family == "spot poisson") {
          # Predictions for mini-batch, old params
          pred_spots_old = family_model[["predict"]](
            X[spots,,drop = F],
            beta_new,
            lambda[spots,,drop = F],
            offset[spots]
          )$total
          
          # "Old" mini-batch likelihood
          lik_old_batch = -poisson_lik(y[spots], pred_spots_old)
          
        } else if (family == "spot negative binomial") {
          pred_spots_old = family_model[["predict"]](
            X[spots, ,drop = F],
            beta_new,
            lambda[spots, ,drop = F],
            offset[spots]
          )$total
          
          lik_old_batch = -nb_lik(y[spots], pred_spots_old, dispersion)
          
        } else if (family == "spot binomial") {
          pred_spots_old = family_model[["predict"]](
            X[spots, ,drop = F],
            beta_new,
            lambda[spots, ,drop = F]
          )$total
          
          lik_old_batch = -binomial_lik(y[spots], pred_spots_old, weights[spots])
        }
        
        # ... Then you compute gradient, do Adam step, get beta_new_temp, etc. ...
        # e.g., beta_new_temp = beta_new + grad_update
        
        # Now compute new mini-batch likelihood with the *tentative* parameters
        if (family == "spot poisson") {
          pred_spots_new = family_model[["predict"]](
            X[spots, ,drop = F],
            beta_new_temp,
            lambda[spots, ,drop = F],
            offset[spots]
          )$total
          
          lik_new_batch = -poisson_lik(y[spots], pred_spots_new)
          
        } else if (family == "spot negative binomial") {
          pred_spots_new = family_model[["predict"]](
            X[spots, ,drop = F],
            beta_new_temp,
            lambda[spots, ,drop = F],
            offset[spots]
          )$total
          
          lik_new_batch = -nb_lik(y[spots], pred_spots_new, dispersion)
          
        } else if (family == "spot binomial") {
          pred_spots_new = family_model[["predict"]](
            X[spots, ,drop = F],
            beta_new_temp,
            lambda[spots, ,drop = F]
          )$total
          
          lik_new_batch = -binomial_lik(y[spots], pred_spots_new, weights[spots])
        }
        
        # Now you can compare the new vs. old batch likelihood:
        lik_diff = lik_new_batch / lik_old_batch
        if(is.na(lik_diff)){
          learning_rate = learning_rate * 0.5
          break
        }
        
        
        # If overshoot or small improvement, reduce step and try again
        if((lik_diff > 1 + 1e-6) | (lik_diff > 1 + 1e-6 & epoch > 1)) {
          learning_rate = learning_rate * 0.5

          recompute_gradients = FALSE
        } else {
          # Accept update
          lik_diff           = min(lik_diff, 1 / lik_diff)
          beta_new           = beta_new_temp
          #dispersion         = dispersion_temp
          m_beta             = m_beta_temp
          v_beta             = v_beta_temp
          
          #if(family == "spot negative binomial"){
          #m_disp          = m_disp_temp
          #v_disp          = v_disp_temp
          #}
          
          t = t_temp
          recompute_gradients = FALSE
        }
      } # end while(recompute_gradients)
      
    } # end for each mini-batch
    
    # --------------------------------------------------
    # (B) After completing all mini-batches in this epoch,
    #     compute full-data likelihood again
    # --------------------------------------------------
    if(family == "spot poisson"){
      pred = family_model[["predict"]](X, beta_new, lambda, offset)$total
      lik_epoch = -poisson_lik(y, pred)
    } else if(family == "spot negative binomial"){
      pred = family_model[["predict"]](X, beta_new, lambda, offset)$total
      lik_epoch = -nb_lik(y, pred, dispersion)
    } else if(family == "spot binomial"){
      pred = family_model[["predict"]](X, beta_new, lambda)$total
      lik_epoch = -binomial_lik(y, pred, weights)
    } else {
      # fallback if needed
      lik_epoch = lik_old
    }
    
    # Evaluate improvement in "full-data" likelihood
    # e.g., check ratio or difference
    ratio = lik_epoch / lik_old
    lik_ratio_test = min(ratio,1/ratio)
    
    #check if convergence occurs from slowing down
    if((abs(lik_ratio_test - lik_ratio) < improvement_threshold)&lik_ratio_test > 0.99){
      conv_counter = conv_counter + 1
    }else{
      conv_counter = 0
    }
    
    lik_ratio = lik_ratio_test
    
    #compute new dispersion
    if(family == "spot negative binomial"){
      fitted_vals = family_model[["predict"]](X, beta_new, lambda, offset)$total
      A = optimize(nb_lik, x = y, mu = fitted_vals, lower = 0.05, upper = 100)
      dispersion = A$minimum
    } else {
      dispersion = 1e8
    }
    
    
    # If ratio is near 1, improvement is small => possible convergence
    #  e.g. if ratio >= max_diff => improvement < some threshold
    if(lik_ratio >= max_diff | conv_counter > max_conv) {
      # Converged
      converged = TRUE
      # cat("Converged at epoch:", epoch, " ratio:", lik_ratio, "\n")
      break
    }
    
    # Update lik_old to the new value for next epoch
    lik_old = lik_epoch
  } # end for epoch
  
  # Final standard errors
  if(sum(good_beta) != 0){
    V = compute_standard_errors(X, y, lambda, family, beta_new,
                                good_beta, dispersion, weights, offset)
  }else{
    V = matrix(NA,ncol(X)*ncol(lambda),ncol(X)*ncol(lambda))
  }
  
  # Final likelihood
  if(family == "spot poisson"){
    pred = family_model[["predict"]](X, beta_new, lambda, offset)$total
    lik = -poisson_lik(y, pred)
  } else if(family == "spot negative binomial"){
    pred = family_model[["predict"]](X, beta_new, lambda, offset)$total
    lik = -nb_lik(y, pred, dispersion)
  } else if(family == "spot binomial"){
    pred = family_model[["predict"]](X, beta_new, lambda)$total
    lik = -binomial_lik(y, pred, weights)
  } else {
    lik = NA
  }
  
  colnames(beta_new) = colnames(lambda)
  rownames(beta_new) = colnames(X)
  
  return(list(
    beta       = beta_new,
    vcov       = V,
    dispersion = dispersion,
    likelihood = lik,
    converged  = converged,
    num_epoch = epoch,
    fixed_coef = fix_coef
    # Could return # of epochs used, final learning rate, etc.
  ))
}













#' Internal Gaussian Spot-GLM Solver
#'
#' Fits a Gaussian GLM for spatial transcriptomics data with optional fixed coefficients.
#' Used internally by \code{\link{spot_glm}} when the family is "spot gaussian".
#'
#' @param X Design matrix (spots × covariates).
#' @param y Response vector.
#' @param lambda Deconvolution matrix (spots × cell types).
#' @param fix_coef Logical matrix indicating which coefficients to fix (TRUE = do not update).
#' @param beta_0 Initial coefficient matrix (covariates × cell types).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{Estimated coefficient matrix.}
#'   \item{vcov}{Variance-covariance matrix of estimated coefficients.}
#'   \item{sigma.sq}{Residual variance estimate.}
#' }
#'
#' @keywords internal
spot_glm_gaussian = function(X,y,lambda,fix_coef,beta_0){
  #make big X
  nct = ncol(lambda)
  #initialize our big covaraite matrix
  big_X = X*lambda[,1]
  #add other covariates
  if(nct > 1){
    for(j in c(2:nct)){
      big_X = cbind(big_X,X*lambda[,j])
    }
  }
  #remove some columns based on fix_coef
  counter = 1
  rem_ind = c()
  keep_ind = c()
  beta_fixed = c()
  for (j in c(1:ncol(fix_coef))){
    for(k in c(1:nrow(fix_coef))){
      #if we need to fix coef, add current column to list of those to remove
      if (fix_coef[k,j] == TRUE){
        rem_ind = c(rem_ind,counter)
        beta_fixed = c(beta_fixed,beta_0[k,j])
      }else{
        keep_ind = c(keep_ind,counter)
      }
      #add to counter
      counter = counter + 1
    }
  }
  #remove fix columns
  fixed_effects = rep(0,nrow(big_X))
  if(length(rem_ind) > 0){
    #get covariates taht we are removing
    X_to_remove = big_X[,rem_ind]
    #update big X matrix
    big_X = big_X[,-rem_ind]
    #account for effect from fixed coefficients
    fixed_effects = X_to_remove%*%matrix(beta_fixed)
  }
  #run regression
  lm1 = lm(y~offset(fixed_effects) + big_X - 1)
  B = coef(lm1)
  #make B into a matrix
  counter = 1
  for (j in c(1:ncol(fix_coef))){
    for(k in c(1:nrow(fix_coef))){
      #if we need to fix coef, add current column to list of those to remove
      if (fix_coef[k,j] != TRUE){
        beta_0[k,j] = B[counter]
        #add to counter
        counter = counter + 1
      }
    }
  }
  #get vcvov matrix
  N = ncol(fix_coef)*nrow(fix_coef)
  V = matrix(NA,N,N)
  V[keep_ind,keep_ind] = vcov(lm1)

  colnames(beta_0) = colnames(lambda)
  rownames(beta_0) = colnames(X)

  #return matrix
  return(list(beta = beta_0,vcov = V,sigma.sq = summary(lm1)$sigma^2))
}


#' Fit Spatial Linear Model (Spot-GLM Gaussian)
#'
#' Fits a linear model for spatial transcriptomics data using deconvolution and fixed coefficients.
#'
#' @param y Response vector.
#' @param X Covariate matrix (spots × covariates).
#' @param lambda Deconvolution matrix (spots × cell types).
#' @param big_X Optional full design matrix (covariates × cell types).
#' @param fix_coef Logical matrix specifying which coefficients to fix.
#' @param beta_0 Initial coefficients.
#'
#' @return A list with:
#' \describe{
#'   \item{beta}{Estimated coefficients.}
#'   \item{vcov}{Variance-covariance matrix.}
#'   \item{std_err_mat}{Standard error matrix.}
#'   \item{sigma.sq}{Residual variance estimate.}
#' }
#'
#' @export
spot_lm = function(y,X,lambda,big_X = NULL,fix_coef = NULL,beta_0 = NULL){
  if(is.null(fix_coef)){
    fix_coef = matrix(FALSE,ncol(X),ncol(lambda))
  }else if(nrow(fix_coef)!= ncol(X)){
    stop("Fix coef must have same #rows as #cols of X")
  }else if(ncol(fix_coef)!= ncol(lambda)){
    stop("Fix coef must have same #cols as #cols of lambda")
  }
  
  
  if(is.null(beta_0)){
    beta_0 = matrix(0,ncol(X),ncol(lambda))
  }else if(nrow(beta_0)!= ncol(X)){
    stop("Initial beta must have same #rows as #cols of X")
  }else if(ncol(beta_0)!= ncol(lambda)){
    stop("Initial beta must have same #cols as #cols of lambda")
  }
  
  #make big X
  nct = ncol(lambda)
  if(is.null(big_X) == T){
    #initialize our big covaraite matrix
    big_X = X*lambda[,1]
    #add other covariates
    if(nct > 1){
      for(j in c(2:nct)){
        big_X = cbind(big_X,X*lambda[,j])
      }
    }
  }
  
  if(ncol(big_X)!= (ncol(X) * nct)){
    stop("Big X must contain same number of rows as #CT times #columns of X")
  }
  if(nrow(big_X)!= (nrow(X) )){
    stop("Big X must contain same numebr of rows as X")
  }
  
  #remove some columns based on fix_coef
  counter = 1
  rem_ind = c()
  keep_ind = c()
  beta_fixed = c()
  for (j in c(1:ncol(fix_coef))){
    for(k in c(1:nrow(fix_coef))){
      #if we need to fix coef, add current column to list of those to remove
      if (fix_coef[k,j] == TRUE){
        rem_ind = c(rem_ind,counter)
        beta_fixed = c(beta_fixed,beta_0[k,j])
      }else{
        keep_ind = c(keep_ind,counter)
      }
      #add to counter
      counter = counter + 1
    }
  }
  #remove fix columns
  fixed_effects = rep(0,nrow(big_X))
  if(length(rem_ind) > 0){
    #get covariates taht we are removing
    X_to_remove = big_X[,rem_ind]
    #update big X matrix
    big_X = big_X[,-rem_ind]
    #account for effect from fixed coefficients
    fixed_effects = X_to_remove%*%matrix(beta_fixed)
  }
  #run regression
  y_shift = y - fixed_effects
  lm1 = lm(y_shift~big_X - 1)
  B = coef(lm1)
  #make B into a matrix
  counter = 1
  for (j in c(1:ncol(fix_coef))){
    for(k in c(1:nrow(fix_coef))){
      #if we need to fix coef, add current column to list of those to remove
      if (fix_coef[k,j] != TRUE){
        beta_0[k,j] = B[counter]
        #add to counter
        counter = counter + 1
      }
    }
  }
  #get vcvov matrix
  N = ncol(fix_coef)*nrow(fix_coef)
  V = matrix(NA,N,N)
  V[keep_ind,keep_ind] = vcov(lm1)
  #save standard error matrix 
  standard_error_mat = matrix(sqrt(diag(V)),ncol(X),ncol(lambda))
  

  colnames(beta_0) = colnames(lambda)
  rownames(beta_0) = colnames(X)
  
  colnames(standard_error_mat) = colnames(lambda)
  rownames(standard_error_mat) = colnames(X)

  #return matrix
  return(list(beta = beta_0,vcov = V,std_err_mat = standard_error_mat,sigma.sq = summary(lm1)$sigma^2))
}



#' Compute Gradient Descent Updates for Spot-GLM
#'
#' Computes gradients of the negative log-likelihood for supported GLM families.
#'
#' @keywords internal
gradient_descent_update = function(family,X,y,beta,lambda,offset,disp = 1e8,weights = rep(1,length(y))){
  if(family == "spot poisson"){
    beta_grad = spot_poisson[["grad"]](X = X,y = y,beta = beta,lambda = lambda,offset = offset)
    return(list(beta_grad = beta_grad$grad,disp_grad = 0,weights = beta_grad$weights))
  }else if(family == "spot negative binomial"){
    beta_grad = spot_negative_binomial[["grad"]](X = X,y = y,beta = beta,lambda = lambda,disp = disp,offset = offset)
    return(list(beta_grad = beta_grad$grad,disp_grad = beta_grad$disp_grad,weights = beta_grad$weights))
  }else if(family == "spot binomial"){
    beta_grad = spot_binomial[["grad"]](X = X,y = y,beta = beta,lambda = lambda,weights = weights)
    return(list(beta_grad = beta_grad$grad,disp_grad = 0,weights = beta_grad$weights))
  }
}


#' Compute Standard Errors via Fisher Information
#'
#' Computes the variance-covariance matrix of estimates from the Fisher Information Matrix.
#'
#' @keywords internal
compute_standard_errors = function(X,y,lambda,family,beta,good_beta,dispersion,weights,offset){
  #get vcov matrix
  fisher_info = glm_fisher(X,y,lambda,family,beta,good_beta,dispersion = dispersion,weights = weights,offset = offset)
  #get degenerate columns
  degen_ind = which(diag(fisher_info) == 0)
  #update good_beta
  counter = 0
  for(j_ind in c(1:ncol(good_beta))){
    for(k_ind in c(1:nrow(good_beta))){
      #skip to next index if it is 0
      if(good_beta[k_ind,j_ind] == 1){
        counter = counter + 1
      }else{
        next
      }
      #check if counter is the value of a degenerate
      if(counter %in% degen_ind){
        good_beta[k_ind,j_ind] = 0
      }
    }
  }
  #remove degen columns
  if(length(degen_ind) > 0){
    fisher_info = fisher_info[-degen_ind,-degen_ind]
  }
  invert = FALSE
  #invert fisher information
  tryCatch({
    A = Matrix::chol(fisher_info,LDL = FALSE,perm = FALSE)
    #get covariance matrix
    vcov = Matrix::solve(A,tol = 1e-500)%*%Matrix::t(Matrix::solve(A,tol = 1e-500))
    invert = TRUE
  }, error = function(e) {
    #print(paste("An error occurred:", e))
    #print(c(loop_iter,j))
    #nothing
  })
  if(invert == FALSE){
    tryCatch({
      #get covariance matrix
      vcov = solve(fisher_info,tol = 1e-500)
      invert = TRUE
    }, error = function(e) {
      #print(paste("An error occurred:", e))
      #print(c(loop_iter,j))
      #nothing
    })
  }
  if(invert == FALSE){
    print("Failed to generate standard errors")
    return (matrix(NA,length(beta),length(beta)))
  }


  #get number of covariates
  p = nrow(beta)
  #get number of cell types
  nCT = ncol(beta)
  #get full covariance matrix
  V = matrix(NA,p*nCT,p*nCT)
  #get indices per cell type
  fill_ind = c()
  for(j in c(1:nCT)){
    #get indices for cell type j
    true_start = p*(j-1) + 1
    true_end = p*j
    total_ind = (true_start:true_end)
    #get which indices are actually filled
    cell_ind = which(good_beta[,j] == 1)
    if(length(cell_ind) > 0){
      cell_ind = total_ind[cell_ind]
    }
    fill_ind = append(fill_ind,cell_ind)
  }

  #fill it with values
  V[fill_ind,fill_ind] = vcov

  #turn to sparce matrix
  V = Matrix::Matrix(V, sparse = TRUE)
}


#' Fisher Information Matrix for Spot-GLM
#'
#' Computes the Fisher Information matrix for a given GLM family and model.
#'
#' @keywords internal
glm_fisher = function(X,y,lambda,family,beta,good_beta,dispersion,weights = NULL,offset = rep(0,length(y))){
  #Step 0:making sure family is correct
  if(family == "spot poisson"){
    family_model = spot_poisson
    predictions = family_model[['predict']](X,beta,lambda,offset)
  }else if (family == "spot gaussian"){
    family_model = spot_gaussian
    predictions = family_model[['predict']](X,beta,lambda)
  }else if (family == "spot binomial"){
    family_model = spot_binomial
    predictions = family_model[['predict']](X,beta,lambda)
  }else if (family == "spot negative binomial"){
    family_model = spot_negative_binomial
    predictions = family_model[['predict']](X,beta,lambda,offset)
  }else{
    stop("Family not found. Family should be one of gaussian,binomial,poisson, or negative binomial")
  }

  #get number of cell types
  nCT = ncol(lambda)
  #get number of covariates
  p = ncol(X)
  #get number of true covariates
  total_cov = sum(good_beta)
  #initialize covariance matrix
  V = matrix(NA,total_cov,total_cov)
  #get indices per cell type
  CT_ind = list()
  ind_start = 1
  for(j in c(1:ncol(good_beta))){
    #get indices for cell type j
    ind_end = ind_start + sum(good_beta[,j]) - 1
    if(ind_end >= ind_start){
      CT_ind[[j]] = ind_start:ind_end
    }else{
      CT_ind[[j]] = c(NA)
    }
    #reset ind_start
    ind_start = ind_end + 1
  }



  if(family == "spot gaussian"){
    #iterate over each cell type combo
    for(j in c(1:nCT)){
      if(is.na(CT_ind[[j]][1])){
        next
      }
      #get number of cell type j for spots
      N_t = lambda[,j]
      #get indicies for this cell type
      t_ind = which(good_beta[,j] == 1)
      X_t = X[,t_ind]
      for(k in c(j:nCT)){
        if(is.na(CT_ind[[k]][1])){
          next
        }
        #get number of cell type j for spots
        N_t2 = lambda[,k]
        #get indices for this cell type
        t_ind2 = which(good_beta[,k] == 1)
        X_t2 = X[,t_ind2]
        #get A matrix
        A = N_t * N_t2
        #get indices for each cell type
        j_ind = CT_ind[[j]]

        k_ind = CT_ind[[k]]
        #compute matrix
        fill = (1/dispersion) * t(X_t)%*%(A*X_t2)
        V[(j_ind),(k_ind)] = fill
        V[(k_ind),(j_ind)] = t(fill)
      }
    }
    return(V)
  }


  if(family == "spot poisson"){
    #iterate over each cell type combo
    for(j in c(1:nCT)){
      if(is.na(CT_ind[[j]][1])){
        next
      }
      #get number of cell type j for spots
      N_t = lambda[,j]
      mu_t = predictions$individual[,j]
      #get indicies for this cell type
      t_ind = which(good_beta[,j] == 1)
      X_t = X[,t_ind]
      for(k in c(j:nCT)){
        if(is.na(CT_ind[[k]][1])){
          next
        }
        #get number of cell type j for spots
        N_t2 = lambda[,k]
        mu_t2 = predictions$individual[,k]
        #get indices for this cell type
        t_ind2 = which(good_beta[,k] == 1)
        X_t2 = X[,t_ind2]
        #get A matrix
        A = (N_t * N_t2 * mu_t * mu_t2)/(predictions$total)
        #make NA 0
        A[is.na(A)] = 0
        #get indices for each cell type
        j_ind = CT_ind[[j]]
        k_ind = CT_ind[[k]]
        #compute matrix
        #compute matrix
        fill = t(X_t)%*%(A*X_t2)
        V[(j_ind),(k_ind)] = fill
        V[(k_ind),(j_ind)] = t(fill)
      }
    }
    return(V)
  }

  if(family == "spot negative binomial"){
    #iterate over each cell type combo
    for(j in c(1:nCT)){
      if(is.na(CT_ind[[j]][1])){
        next
      }
      #get number of cell type j for spots
      N_t = lambda[,j]
      mu_t = predictions$individual[,j]
      #get indicies for this cell type
      t_ind = which(good_beta[,j] == 1)
      X_t = X[,t_ind]
      for(k in c(j:nCT)){
        if(is.na(CT_ind[[k]][1])){
          next
        }
        #get number of cell type j for spots
        N_t2 = lambda[,k]
        mu_t2 = predictions$individual[,k]
        #get indices for this cell type
        t_ind2 = which(good_beta[,k] == 1)
        X_t2 = X[,t_ind2]
        #get A matrix
        A = (dispersion * N_t * N_t2 * mu_t * mu_t2)/(predictions$total * (predictions$total + dispersion))
        #make NA 0
        A[is.na(A)] = 0
        #get indices for each cell type
        j_ind = CT_ind[[j]]
        k_ind = CT_ind[[k]]
        #compute matrix
        #compute matrix
        fill = t(X_t)%*%(A*X_t2)
        V[(j_ind),(k_ind)] = fill
        V[(k_ind),(j_ind)] = t(fill)
      }
    }
    return(V)
  }


  if(family == "spot binomial"){
    #iterate over each cell type combo
    for(j in c(1:nCT)){
      if(is.na(CT_ind[[j]][1])){
        next
      }
      #get number of cell type j for spots, add offset to avoid A being 0
      N_t = lambda[,j]
      mu_t = 0.9999*predictions$individual[,j] + 0.0001*0.5
      #get indicies for this cell type
      t_ind = which(good_beta[,j] == 1)
      X_t = X[,t_ind]
      for(k in c(j:nCT)){
        if(is.na(CT_ind[[k]][1])){
          next
        }
        #get number of cell type j for spots
        N_t2 = lambda[,k]
        mu_t2 = 0.9999*predictions$individual[,k] + 0.0001*0.5
        #get indices for this cell type
        t_ind2 = which(good_beta[,k] == 1)
        X_t2 = X[,t_ind2]
        #get A matrix
        A1 = (weights * N_t * N_t2)
        A2 = (mu_t/(1+mu_t)^2) * (mu_t2/(1+mu_t2)^2)
        A3 = predictions$total * (1-predictions$total)
        A = A1*A2/A3

        #make NA 0
        A[is.na(A)] = 0
        #make infinites 0 (basically where we predict 0 or 1)
        A[is.infinite(A)] = 0
        #get indices for each cell type
        j_ind = CT_ind[[j]]

        k_ind = CT_ind[[k]]
        #compute matrix
        fill = t(X_t)%*%(A*X_t2)
        V[(j_ind),(k_ind)] = fill
        V[(k_ind),(j_ind)] = t(fill)

      }
    }
    return(V)
  }
}





#' Run Single-Cell GLM per Cell Type
#'
#' Fits individual GLMs for each cell type using pseudo-bulked or single-cell data.
#'
#' @keywords internal
run_single_cell = function(y,X,lambda,sc_family = "gaussian",offset = rep(0,length(y)),weights =rep(1,length(y)),
                           fix_coef = matrix(FALSE,ncol(X),ncol(lambda)),CT = NULL){
  if(is.null(weights)){
    weights = rep(1,length(y))
  }
  
  #get main cell types
  if(is.null(CT)){
    CT = apply(lambda,1,function(x){which.max(x)})
  }
  #initialize beta and std_err matrices
  beta = matrix(0,nrow = ncol(X),ncol = ncol(lambda))
  std_err_mat = matrix(NA,nrow = ncol(X),ncol = ncol(lambda))
  t1 = Sys.time()
  #iterate model over each cell type 
  for(r in c(1:ncol(lambda))){
    #get spots classified as this cell type 
    cells = which(CT == r)
    if(length(cells) < 2500){
      cells = order(lambda[,r],decreasing = T)[1:2500]
      cells = cells[lambda[cells,r] >= 0.1]
    }
    cells = cells[is.na(cells) == F]
    #get how often a covaraite appears
    freq = apply(X[cells,,drop = F],2,function(x){sum(x!=0)})
    #update fix_coef
    bad_cov = which(freq < 5)
    if(length(bad_cov) > 0){
      fix_coef[bad_cov,r] = TRUE
    }
    if(sum(fix_coef[,r,drop = F]) == nrow(fix_coef)){
      next
    }
    #get covariates for these spots 
    features = X[cells,which(fix_coef[,r,drop = F] == FALSE),drop = FALSE]
    if(ncol(features) == 0){
      next
    }
    #get observations 
    obs = y[cells]
    #get offset
    offset_CT = offset[cells]
    #get weights
    weights_CT = weights[cells]
    #run regular regression 
    if(sc_family == "negative binomial"){
      lm1 = suppressWarnings(glm(obs/weights_CT~offset(offset_CT) + features-1, family = "poisson",weights = weights_CT*lambda[cells,r]))
    }else if(sc_family == "binomial"){
      lm1 = suppressWarnings(glm(obs/weights_CT~ features-1, family = sc_family,weights = weights_CT*lambda[cells,r]))
    }else{
      lm1 = suppressWarnings(glm(obs/weights_CT~offset(offset_CT) + features-1, family = sc_family,weights = weights_CT*lambda[cells,r]))
    }
    #get vcov matrix
    V = diag(vcov(lm1))
    #update beta and std err matrices
    beta[which(fix_coef[,r] == FALSE),r] = coef(lm1)
    std_err_mat[which(fix_coef[,r] == FALSE),r] = sqrt(V)
  }
  t2 = Sys.time()
  #remove NAs
  beta[is.na(beta)] = 0
  #print(t2-t1)
  return (list(beta_est = beta, stand_err_mat = std_err_mat,fix_coef = fix_coef))
}













#' Initialize Fixed Coefficient Matrix
#'
#' Automatically sets which coefficients to fix based on low signal or coverage.
#'
#' @keywords internal
initialize_fix_coef = function(X,lambda,fix_coef = matrix(FALSE,ncol(X),ncol(lambda)),
                               weights = NULL,min_deconv = 0.1,min_freq = 100){
  #number of cell types
  nCT = ncol(lambda)
  p = ncol(X)
  #get valid betas 
  good_beta = 0*fix_coef
  #get spots that are good for each cell type
  good_spots_ct = vector("list",nCT)
  for (j in c(1:nCT)){
    if(is.null(weights)){
      good_spots_ct[[j]] = which(lambda[,j] > min_deconv)
    }else{
      good_spots_ct[[j]] = which(lambda[,j] > min_deconv & (weights != 0))
    }
  }
  
  good_cov_ct = vector("list",nCT)
  for (j in c(1:nCT)){
    if(length(good_spots_ct[[j]]) == 0){
      fix_coef[, j]
      next
    }
    # Determine covariates with sufficient variation
    good_cov <- which(apply(X[good_spots_ct[[j]],,drop = FALSE],2,function(x){sum(x!=0)}) > min_freq)
    bad_cov <- setdiff(1:p, good_cov)
    
    # Update fixed coefficients
    if (length(bad_cov) > 0) {
      fix_coef[bad_cov, j] = TRUE
    }
  }
  
  #return fix_coef and good_beta
  return(list(fix_coef = fix_coef))
}
















#' Fit a Spatial GLM with Initialization and Optimization
#'
#' Fits a generalized linear model (GLM) with spatial deconvolution for a single response variable 
#' (e.g., gene expression), supporting Poisson, Gaussian, Binomial, and Negative Binomial families.
#' This function handles coefficient initialization, model fitting via mini-batch gradient descent, 
#' and automatic coefficient filtering for weak covariates or poorly represented cell types.
#'
#' @param y Numeric response vector (e.g., gene expression for one gene across spots).
#' @param X Covariate matrix (spots × covariates).
#' @param lambda Deconvolution matrix (spots × cell types).
#' @param family GLM family: `"spot gaussian"`, `"spot poisson"`, `"spot negative binomial"`, or `"spot binomial"`.
#' @param beta_0 Optional initial coefficient matrix (covariates × cell types).
#' @param fix_coef Logical matrix (covariates × cell types) indicating coefficients to fix during optimization.
#' @param offset Optional numeric vector (same length as `y`), used for Poisson or NB normalization.
#' @param initialization Boolean if initialization via single cell approximation should be performed. Default TRUE.
#' @param CT Optional vector of dominant cell type labels per spot.
#' @param weights Observation weights (same length as `y`).
#' @param ct_cov_weights Optional vector of cell-type–specific weights (length = number of cell types).
#' @param n_epochs Number of training epochs for gradient descent.
#' @param batch_size Size of mini-batches used during gradient descent.
#' @param learning_rate Initial learning rate for optimization.
#' @param max_diff Convergence threshold based on likelihood ratio.
#' @param improvement_threshold Minimum required improvement in likelihood ratio between epochs.
#' @param max_conv Number of consecutive low-improvement epochs before convergence is assumed.
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_est}{Estimated coefficient matrix (covariates × cell types).}
#'   \item{stand_err_mat}{Standard error matrix for each coefficient.}
#'   \item{time}{Elapsed fitting time (in seconds).}
#'   \item{disp}{Estimated dispersion (for NB models).}
#'   \item{converged}{Logical indicating if convergence was reached.}
#'   \item{likelihood}{Final negative log-likelihood.}
#'   \item{vcov}{Variance-covariance matrix.}
#'   \item{niter}{Number of optimization epochs completed.}
#'   \item{fixed_coef}{Final matrix indicating fixed coefficients.}
#' }
#'
#' @export

run_model = function(y,X,lambda,family = "spot gaussian",beta_0 = NULL,fix_coef = NULL,
                     offset = rep(0,length(y)),initialization = T,
                     CT = NULL, weights = rep(1,length(y)),ct_cov_weights = rep(1,ncol(lambda)),
                     n_epochs = 100,batch_size = 500,learning_rate = 1,max_diff = 1-1e-6, improvement_threshold = 1e-6,
                     max_conv = 10){
  #Step 0: Pre-processing 
  if(is.null(weights)){
    weights = rep(1,length(y))
  }else if(length(weights)!= length(y)){
    stop("Weights must be same length as observations")
  }
  
  if(is.null(offset)){
    offset = rep(0,length(y))
  }else if(length(offset)!= length(y)){
    stop("Offsets must be same length as observations")
  }
  
  #remove spots with no weight
  bad_spots = which(weights == 0)
  if(length(bad_spots) > 0){
    y = y[-bad_spots]
    X = X[-bad_spots,,drop = F]
    lambda = lambda[-bad_spots,,drop = F]
    offset = offset[-bad_spots]
    weights = weights[-bad_spots]
  }
  
  
  #weight lambda by cov weights and normalize
  if(is.null(ct_cov_weights) == F){
    if(length(ct_cov_weights) != ncol(lambda)){
      stop("Cell type covariate weights must be the same length as #col lambda")
    }
    lambda = sweep(lambda,2,ct_cov_weights,"*")
    lambda = sweep(lambda,1,rowSums(lambda),"/")
    lambda[is.na(lambda)] = 1/ncol(lambda)
  }
  
  #get family
  if(family == "spot gaussian"){
    model_family = "spot_gaussian"
    sc_family = "gaussian"
  }else if(family == "spot poisson"){
    model_family = "spot_poisson"
    sc_family = "poisson"
  }else if(family == "spot negative binomial"){
    model_family = "spot_negative_binomial"
    sc_family = "poisson"
  }else if(family == "spot binomial"){
    model_family = "spot_binomial"
    sc_family = "binomial"
  }else{
    stop("Family must be one of spot gaussian, spot poisson, spot negative binomial, or spot binomial")
  }
  
  
  if(is.null(fix_coef)){
    fix_coef = matrix(FALSE,ncol(X),ncol(lambda))
  }else if( (nrow(fix_coef)!= ncol(X)) |(ncol(fix_coef)!= ncol(lambda))){
    stop("Fixed coefficients matrix must be of dimension ncol(X) by ncol(lambda)")
  }
  if(family != "spot gaussian"){
    #Step 1: get initial beta
    if(is.null(beta_0)){
      beta_0 = matrix(0,ncol(X),ncol(lambda))
    }else if( (nrow(beta_0)!= ncol(X)) |(ncol(beta_0)!= ncol(lambda))){
      stop("Initial beta matrix must be of dimension ncol(X) by ncol(lambda)")
    }
    
    
    
    
    if(initialization == T){
      initial= spotglm:::run_single_cell(y = y, X = X, lambda = lambda,sc_family = sc_family,offset = offset,
                                                                   weights = weights,fix_coef = fix_coef,CT = CT)
      beta_0 = initial$beta
      fix_coef = initial$fix_coef
    }
    #Step 2: Get fix coef and good beta
    #get fixed coefficients 
    #coef_info = spotglm:::initialize_fix_coef(X = X, lambda = lambda,fix_coef = fix_coef,
                                              #weights = weights,min_deconv = min_deconv, min_freq = min_freq)
    #fix_coef = coef_info$fix_coef
    
  }
 
  #Step 3:Run model 
  
  t1 = Sys.time()
  LR = learning_rate
  conv = FALSE
  result = spot_glm(X = X,y = y,lambda = lambda,beta_0 = beta_0,family = family,offset = offset,fix_coef = fix_coef,
          weights = weights,learning_rate = LR,
          n_epochs = n_epochs,batch_size = batch_size, max_diff = max_diff,improvement_threshold = improvement_threshold,max_conv = max_conv)
  
  conv = result$converged
  t2 = Sys.time()
  #get standard errors for coefs
  std_err = rep(NA,nrow(result$vcov))
  for (j in c(1:nrow(result$vcov))){
    std_err[j] = suppressWarnings(sqrt(result$vcov[j,j]))
  }
  
  stand_err_mat = matrix(std_err,nrow = ncol(X),byrow = F)
  colnames(stand_err_mat) = colnames(lambda)
  rownames(stand_err_mat) = colnames(X)
                         
  return (list(beta_est = result$beta, stand_err_mat = stand_err_mat,time = t2 - t1,
               disp = result$dispersion,converged = result$converged,likelihood = result$likelihood,vcov = result$vcov,
               niter = result$num_epoch,fixed_coef = result$fixed_coef))
}







#' Parallelized Spot-GLM Model Fitting (Windows)
#'
#' Fits a Spot-GLM model for multiple response variables (e.g., genes) in parallel using `foreach` and
#' `doParallel`. Optimized for Windows systems where `mclapply` is not available.
#'
#' @param Y Response matrix (spots × responses).
#' @param X Covariate matrix (spots × covariates).
#' @param lambda Deconvolution matrix (spots × cell types).
#' @param family The GLM family to use. One of: `"spot gaussian"`, `"spot poisson"`, 
#'   `"spot negative binomial"`, or `"spot binomial"`.
#' @param beta_0 Optional initial coefficient matrix (covariates × cell types).
#' @param fix_coef Optional logical matrix indicating which coefficients to fix (same dimensions as `beta_0`).
#' @param offset Optional numeric vector (length equal to number of spots).
#' @param initialization Boolean if initialization via single cell approximation should be performed. Default TRUE.
#' @param G Maximum chunk size (in GB) to control memory usage during parallelization.
#' @param num_cores Number of CPU cores to use in parallel.
#' @param CT Optional vector of dominant cell types per spot.
#' @param weights Optional observation-level weight matrix (spots × genes).
#' @param ct_cov_weights Optional cell-type-specific weight matrix (cell types × genes).
#' @param n_epochs Number of training epochs.
#' @param batch_size Size of each mini-batch.
#' @param learning_rate Initial learning rate.
#' @param max_diff Convergence threshold based on likelihood improvement ratio.
#' @param improvement_threshold Minimum improvement ratio between epochs.
#' @param max_conv Number of low-improvement epochs before stopping.
#'
#' @return A named list of model results (one per gene), each containing:
#' \describe{
#'   \item{beta_est}{Estimated coefficients.}
#'   \item{stand_err_mat}{Standard error matrix.}
#'   \item{disp}{Dispersion estimate (if applicable).}
#'   \item{likelihood}{Final log-likelihood.}
#'   \item{converged}{Convergence status.}
#'   \item{niter}{Number of epochs run.}
#'   \item{vcov}{Variance-covariance matrix.}
#'   \item{fixed_coef}{Final fixed coefficients matrix.}
#' }
#'
#' @details
#' This function splits the gene expression matrix into memory-safe chunks, 
#' then evaluates each chunk in parallel using `foreach` and `doParallel`. 
#' For Mac/Linux, use \code{\link{run_spot_glm_mac}}.
#'
#' External dependencies include \code{Matrix}, \code{MASS}, \code{LaplacesDemon}.
#'
#' @importFrom foreach %dopar%
#' @import parallel
#' @import doParallel
#' @import Matrix
#' @import MASS
#' @importFrom LaplacesDemon invlogit logit
#' @export
run_model_parallel_windows = function(Y,X,lambda,family = "spot",beta_0 = NULL,fix_coef = NULL,offset = NULL,
                        initialization = T,G = 0.1,num_cores = 1,
                        CT = NULL, weights = NULL,ct_cov_weights = NULL,
                        n_epochs = 100,batch_size = 500,learning_rate = 1,max_diff = 1-1e-6, improvement_threshold = 1e-6,
                        max_conv = 10){
  
  #get offset values
  if(is.null(offset) == TRUE & family %in% c("spot poisson","spot negative binomial")){
    offset = log(Matrix::rowSums(Y))
    if(sum(is.infinite(offset)) > 0){
      stop("Remove spots with 0 counts")
    }
  }else if(length(offset) < nrow(Y)){
    stop("Length of offsets must be the same as the number of observations")
  }
  #check if ct_cov_weights and weights are proper matrices
  if(is.null(ct_cov_weights) == F){
    if( (nrow(ct_cov_weights) != ncol(lambda)) |(ncol(ct_cov_weights) != ncol(y)) ){
      stop("Cell type covariate weights should be a vector of dimension #cell types by #responses(e.g. genes)")
    }
  }
  
  if(is.null(weights) == F){
    if( (nrow(weights) != nrow(y)) |(ncol(weights) != ncol(y)) ){
      stop("Weights should be a vector of dimension #observations (e.g. spots) by #responses(e.g. genes)")
    }
  }
  
  
  #estimate size of full data matrix 
  data_size = 8 * prod(dim(Y))/1e+09
  #number of chunks needed
  nchunks = ceiling(data_size/G)
  
  print(paste0("Splitting Data into ", nchunks, " chunks in order to avoid memory overload. Each chunk is less than ", 
               G, " gigabytes."))
  #get number of chunks 
  chunk_size = ceiling(ncol(Y)/nchunks)
  #group data into chunks 
  grouping <- rep(1:nchunks, each = chunk_size, length.out = ncol(Y))
  index_chunks = split(1:ncol(Y), grouping)
  #start chunk counter
  chunk_counter = 1
  #start spotGLM
  T_1 = Sys.time()
  for (I in c(1:length(index_chunks))) {
    t1 = Sys.time()
    print(paste0("Evaluating chunk ", I, " out of ", 
                 nchunks))
    counts_chunk = as.matrix(Y[, index_chunks[[I]]])
    
    if(is.null(ct_cov_weights)){
      ct_cov_weights_chunk = NULL
    }else{
      ct_cov_weights_chunk = ct_cov_weights[,index_chunks[[I]]]
    }
    
    if(is.null(weights)){
      weights_chunk = NULL
    }else{
      weights_chunk = weights[,index_chunks[[I]]]
    }
    
    #make cluster 
    print("Initializing cluster")
    cluster <- parallel::makeCluster(num_cores, outfile = "")
    doParallel::registerDoParallel(cluster)
    # Ensure the cluster stops even if an error or interrupt occurs
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    #iterate over chunk
    NC = ncol(counts_chunk)
    results_chunk = foreach::foreach(i = 1:NC,.export = c("X", "lambda", "offset", "CT", "initialization", 
                                                          "n_epochs","batch_size", "learning_rate", "max_diff", "family",
                                                          "counts_chunk", "weights_chunk", "ct_cov_weights_chunk","improvement_threshold","max_conv"), .packages = c("Matrix", "MASS", "LaplacesDemon","spotglm")) %dopar% {
      tryCatch({
        print(paste0("On iteration ", i))
        t1 = Sys.time()
        y <- counts_chunk[, i]
        
        CT_COVARIATE_WEIGHTS <- if (is.null(ct_cov_weights_chunk)) {
          NULL
        } else {
          ct_cov_weights_chunk[, i]
        }
        
        WEIGHTS <- if (is.null(weights_chunk)) {
          NULL
        } else {
          weights_chunk[, i]
        }
        
        run_model(
          y = y, X = X, lambda = lambda, family = family,
          beta_0 = NULL, fix_coef = NULL, offset = offset,
          initialization = initialization, 
          CT = CT, weights = WEIGHTS, ct_cov_weights = CT_COVARIATE_WEIGHTS,
          n_epochs = n_epochs,batch_size = batch_size, learning_rate = learning_rate,
          max_diff = max_diff, improvement_threshold = improvement_threshold,
          max_conv = max_conv
        )
        print(Sys.time() - t1)
      }, error = function(e) {
        print(e)
        warning(paste("Gene", i, "failed:", conditionMessage(e)))
        return(NULL)
      })
    }   
    
    if (chunk_counter == 1) {
      results = results_chunk
    }
    else {
      results = c(results, results_chunk)
    }
    rm(counts_chunk)
    gc()
    chunk_counter = chunk_counter + 1
    print("Closing cluster")
    parallel::stopCluster(cluster)
    print(paste0("Chunk took ",Sys.time() - t1))
  }
  names(results) = colnames(Y)
  T_2 = Sys.time()
  
  return(results)
  
}

#' Parallelized Spot-GLM Model Fitting (macOS / Linux)
#'
#' Fits a Spot-GLM model for multiple responses (e.g., genes) in parallel using memory-safe chunking
#' and `pbmclapply`, which relies on `mclapply`. Only available on Unix-based systems.
#'
#' @inheritParams run_model_parallel_windows
#'
#' @details
#' This version uses `pbmcapply::pbmclapply` for parallelism. On Windows systems,
#' please use \code{\link{run_spot_glm_windows}}.
#'
#' @import parallel
#' @import Matrix
#' @import MASS
#' @importFrom LaplacesDemon invlogit logit
#' @importFrom pbmcapply pbmclapply
#' @export
run_model_parallel_mac = function(Y, X, lambda, family = "spot gaussian", beta_0 = NULL, fix_coef = NULL,
                                initialization = T, G = 0.1, num_cores = 1,offset = NULL,
                                 CT = NULL, weights = NULL, ct_cov_weights = NULL,
                                 n_epochs = 100,batch_size = 500, learning_rate = 1, max_diff = 1 - 1e-6, improvement_threshold = 1e-6,
                                max_conv = 10) {
  #get offset values
  if(is.null(offset) == TRUE & family %in% c("spot poisson","spot negative binomial")){
    offset = log(Matrix::rowSums(Y))
    if(sum(is.infinite(offset)) > 0){
      stop("Remove spots with 0 counts")
    }
  }else if(length(offset) < nrow(Y)){
    stop("Length of offsets must be the same as the number of observations")
  }
  
  
  if (!is.null(ct_cov_weights)) {
    if ((nrow(ct_cov_weights) != ncol(lambda)) || (ncol(ct_cov_weights) != ncol(Y))) {
      stop("Cell type covariate weights should be a matrix of dimension #cell types by #responses (e.g., genes)")
    }
  }
  
  if (!is.null(weights)) {
    if ((nrow(weights) != nrow(Y)) || (ncol(weights) != ncol(Y))) {
      stop("Weights should be a matrix of dimension #observations (e.g., spots) by #responses (e.g., genes)")
    }
  }
  
  data_size = 8 * prod(dim(Y)) / 1e+09
  nchunks = max(ceiling(data_size / G), num_cores)
  
  message("Splitting Data into ", nchunks, " chunks to avoid memory overload. Each chunk < ", G, " GB.")
  
  chunk_size = ceiling(ncol(Y) / nchunks)
  grouping <- rep(1:nchunks, each = chunk_size, length.out = ncol(Y))
  index_chunks = split(1:ncol(Y), grouping)
  
  process_chunk <- function(indices, chunk_id) {
    message("Processing chunk ", chunk_id, " with ", length(indices), " genes")
    counts_chunk = as.matrix(Y[, indices])
    if (!is.null(ct_cov_weights)) {
      ct_cov_weights_chunk = ct_cov_weights[, indices]
    } else {
      ct_cov_weights_chunk = NULL
    }
    
    if (!is.null(weights)) {
      weights_chunk = weights[, indices]
    } else {
      weights_chunk = NULL
    }
    
    NC = length(indices)
    results_chunk = vector("list", NC)
    for (i in seq_len(NC)) {
      if (i %% 50 == 0 || i == 1 || i == NC) {
        message("Chunk ", chunk_id, ": Processing gene ", i, "/", NC)
      }
      
      y = counts_chunk[, i]
      CT_COVARIATE_WEIGHTS <- if (is.null(ct_cov_weights_chunk)) NULL else ct_cov_weights_chunk[, i]
      WEIGHTS <- if (is.null(weights_chunk)) NULL else weights_chunk[, i]
      
      result <- tryCatch({
        spotglm::run_model(
          y = y, X = X, lambda = lambda, family = family,
          beta_0 = beta_0, fix_coef = fix_coef, offset = offset,
          initialization = initialization,  
          CT = CT, weights = WEIGHTS, ct_cov_weights = CT_COVARIATE_WEIGHTS,
          n_epochs = n_epochs,batch_size = batch_size, learning_rate = learning_rate,
          max_diff = max_diff,improvement_threshold = improvement_threshold,max_conv = max_conv
        )
      }, error = function(e) {
        warning(paste("Gene", indices[i], "failed:", conditionMessage(e)))
        return(NULL)
      })
      results_chunk[[i]] <- result
    }
    names(results_chunk) <- colnames(Y)[indices]
    results_chunk
  }
  
  results_list <- pbmcapply::pbmclapply(seq_along(index_chunks), function(i) {
    tryCatch({
      process_chunk(index_chunks[[i]], chunk_id = i)
    }, error = function(e) {
      message("Chunk ", i, " failed: ", conditionMessage(e))
      return(NULL)
    })
  }, mc.cores = num_cores,mc.preschedule = FALSE, ignore.interactive = TRUE)
  
  results <- do.call(c, results_list)
  
  return(results)
}



#' @keywords internal
compute_effective_niche = function(coord,data,sigma,batch_size = 1000, cutoff = 0.05)
{
  #initialize EN matrix
  EN = matrix(NA, nrow(coord), ncol(data))
  #initialize KNN list
  rownames(EN) = rownames(coord)
  colnames(EN) = colnames(data)
  #get number of batches
  C = nrow(coord)
  #split cells into rectangles. #rectangles = C/batch_size
  num_iter = ceiling(C/batch_size)
  #get number of rows and number of columns
  rows = ceiling(sqrt(num_iter))
  row_blocks = seq(min(coord[,2]),max(coord[,2]),length.out = rows+1)
  cols = ceiling(num_iter/rows)
  col_blocks = seq(min(coord[,1]),max(coord[,1]),length.out = cols+1)
  print(paste0("Number of rows: ",rows))
  print(paste0("Number of cols: ",cols))
  possible_inds = c(1:nrow(coord))
  for(row in c(1:rows)){
    for(col in c(1:cols)){
      if(length(possible_inds) == 0){
        next
      }
      print(paste0("Row: ",row ," Col ",col))
      #subset coord into rows and columns
      row_block = c(row_blocks[row],row_blocks[row+1])
      col_block = c(col_blocks[col],col_blocks[col+1])
      #get indicies that are encompassed
      inds_ = which(spatstat.utils::inside.range(coord[possible_inds,1],col_block) & spatstat.utils::inside.range(coord[possible_inds,2],row_block))
      if(length(inds_) == 0){
        next
      }
      inds = possible_inds[inds_]
      possible_inds = possible_inds[-inds_]
      #get range of coordiantes
      range_x = range(coord[inds,1])
      range_y = range(coord[inds,2])
      #add sigma to each side (i.e. we are guarnateed that matched points are within the new interval)
      range_x[1] = range_x[1] - sigma
      range_x[2] = range_x[2] + sigma
      
      range_y[1] = range_y[1] - sigma
      range_y[2] = range_y[2] + sigma
      #get candidate cells; cells which may be near a cell in set
      candidate = which(spatstat.utils::inside.range(coord[,1],range_x) & spatstat.utils::inside.range(coord[,2],range_y))
      #compute distance between cells in set and candidate cells
      D = Rfast::dista(xnew = as.matrix(coord[inds,], ncol = 2, byrow = T), x = coord[candidate,], type = "euclidean", trans = T)
      #if we want to use nearest neighbors definitely use this one instead
      #FNN::get.knnx(deconv_quad$reference_embed,deconv_quad$query_embed, 100)
      #set those that are closer than sigma equal to 1
      D[D < sigma] = 1
      #else equal to 0
      D[D>sigma] = 0
      #D[D < cutoff] = 0
      #call this matrix neighbors
      neighbors = D
      #compute sum of expression in neighborhood
      A = (neighbors%*%data[candidate,])
      for(k in c(1:length(inds))){
        EN[inds[k],] = A[k,]
      }
    }
  }
  return(EN)
}










