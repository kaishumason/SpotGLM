
#glm spot model
#' @export
spot_glm = function(X,y,lambda,family,beta_0 = matrix(0,ncol(X),ncol(lambda)),offset = rep(0,nrow(X)),weights = NULL,M = 50,
                    fix_coef = matrix(FALSE,ncol(X),ncol(lambda)),learning_rate = 100, max_gd_steps = 500,max_diff = 1-1e-6,intercept = TRUE){

  #checks for dimensionality
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
  }else if(dim(fix_coef)[1] != dim(beta_0)[1] & dim(fix_coef)[2] != dim(beta_0)[2]){
    print("Fix coefficients matrix is of incorrect dimension. No coefficients will be fixed")
    fix_coef = matrix(FALSE,ncol(X),ncol(lambda))
  }

  #remove NA rows
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

  #Step 0:making sure family is correct
  if(family == "spot poisson"){
    family_model = spot_poisson
  }else if (family == "spot gaussian"){
    family_model = spot_gaussian
  }else if (family == "spot binomial"){
    family_model = spot_binomial
  }else if (family == "spot negative binomial"){
    family_model = spot_negative_binomial
  }else if (family == "spot gaussian power"){
    family_model = spot_gaussian
  }else{
    stop("Family not found. Family should be one of gaussian,gaussian_power,binomial,poisson,or negative binomial")
  }
  #add names to y
  if(is.null(names(y))){
    names(y) = c(1:length(y))
  }
  if((sum(offset != 0) > 0) & (family %in% c("spot poisson","spot negative binomial") == F)){
    stop("Offsets can only be used with spot poisson and spot negative binomial models")
  }

  #Step 1: Get number of cell types
  nCT = ncol(lambda)
  #get number of covariates
  p = ncol(X)
  #get current beta
  beta_curr = beta_0
  #add colnames and rownames to beta
  if (is.null(colnames(lambda)) == F){
    colnames(beta_curr) = colnames(lambda)
  }
  if (is.null(colnames(X)) == F){
    rownames(beta_curr) = colnames(X)
  }
  #get sigma
  sigma.sq = 0
  #While loop to get betas
  loop_iter = 1
  #initialize whihc betas are "good", i.e will be fitted
  good_beta = 0*beta_0

  #get spots that are good for each cell type
  good_spots_ct = vector("list",nCT)
  for (j in c(1:nCT)){
    if(is.null(weights)){
      good_spots_ct[[j]] = which(lambda[,j]*exp(offset)  > 1e-3)
    }else{
      good_spots_ct[[j]] = which(lambda[,j]*exp(offset)  > 1e-3 & (weights != 0))
    }
  }

  good_cov_ct = vector("list",nCT)
  for (j in c(1:nCT)){
    # Determine covariates with sufficient variation
    good_cov <- which(apply(X[good_spots_ct[[j]],,drop = FALSE],2,function(x){length(x) - max(table(x))}) > M)
    bad_cov <- setdiff(1:p, good_cov)

    # Add intercept if needed
    if (intercept) {
      good_cov <- unique(c(1, good_cov))
      bad_cov <- setdiff(1:p, good_cov)
    }

    # Update fixed coefficients
    if (length(bad_cov) > 0) {
      fix_coef[bad_cov, j] = TRUE
    }
    #filter the covariate matrix
    fixed_betas = which(fix_coef[,j] == TRUE)
    #remove those that align with
    if(length(fixed_betas) > 0){
      good_cov_copy = c()
      for(l in c(1:length(good_cov))){
        val = good_cov[l]
        if(val %in% fixed_betas == F){
          good_cov_copy = c(good_cov_copy,val)
        }
      }
      #update good covs
      good_cov = good_cov_copy
    }
    #save good covariates for this ct
    if(is.null(good_cov)){
      good_cov_ct[[j]] = c(-1)
    }else{
      good_cov_ct[[j]] = good_cov
      good_beta[good_cov,j] = 1
    }
  }

  if(family == "spot gaussian"){
    #get model
    model_outcome = spot_glm_gaussian(X,y,lambda,fix_coef,beta_0)
    #save beta
    beta_curr = model_outcome$beta
    #save vcvov matrix
    V = model_outcome$vcov
    #save sigma square
    sigma.sq = model_outcome$sigma.sq
    #get fitted values
    fitted_vals = family_model[['predict']](X,beta_curr,lambda)$total
    #save likelihood and sigma parameter
    lik = -gaussian_lik(x = y,mu = fitted_vals,sigma = sqrt(sigma.sq)) #get variance of observations
    #return list
    return (list(beta = beta_curr,vcov = V,dispersion = sigma.sq,likelihood = lik,R = 0,converged = TRUE))
  }

  if(family == "spot negative binomial"){
    fitted_vals = family_model[['predict']](X,beta_0,lambda,offset)$total
    #maximize likelihood using gamma
    #get initial overdispersion parameter
    A = optimize(nb_lik,x = y,mu = fitted_vals, lower = 0.05, upper = 10)
    dispersion = A$minimum
  }else{
    dispersion = 1e8
  }
  #gradient descent step
  beta_new = beta_0

  #get number of data points
  N = length(y)
  #get weights for spots
  batch_size = N
  #split data into batches
  vec <- 1:N
  #initialize counter
  counter = 1
  #get difference in likelihood across successive steps of gd
  lik_diff = -1
  # Adam parameters
  m_beta = 0  # First moment for beta
  v_beta = 0  # Second moment for beta
  m_disp = 0  # First moment for dispersion (if needed)
  v_disp = 0  # Second moment for dispersion (if needed)
  beta_1 = 0.9
  beta_2 = 0.999
  epsilon = 1e-8  # Small constant to prevent division by zero

  # Initialize timestep
  t = 0
  if (family == "spot poisson") {
    pred = family_model[["predict"]](X, beta_0, lambda, offset)$total
    lik_old = -poisson_lik(y, pred)
  } else if (family == "spot negative binomial") {
    pred = family_model[["predict"]](X, beta_0, lambda, offset)$total
    lik_old = -nb_lik(y, pred, dispersion)
  } else if (family == "spot binomial") {
    pred = family_model[["predict"]](X, beta_0, lambda)$total
    lik_old = -binomial_lik(y, pred, weights)
  }
  while (lik_diff < max_diff & counter < max_gd_steps) {
    # Select spots
    spots = vec

    # Compute gradients (initially or after learning rate adjustment)
    recompute_gradients = TRUE
    while (recompute_gradients) {
      G = gradient_descent_update(family, X[spots, ], y[spots], beta_new, lambda[spots, ],
                                  offset[spots], disp = dispersion, weights = weights[spots])
      grad = (1 / length(spots)) * G$beta_grad
      disp_grad = (1 / length(spots)) * G$disp_grad


      # Handle fixed coefficients
      grad[is.na(grad)] = 0
      grad[fix_coef == TRUE] = 0


      # Tentative moment updates
      t_temp = t + 1
      m_beta_temp = beta_1 * m_beta + (1 - beta_1) * grad
      v_beta_temp = beta_2 * v_beta + (1 - beta_2) * (grad^2)

      # Bias correction
      m_beta_hat = m_beta_temp / (1 - beta_1^t_temp)
      v_beta_hat = v_beta_temp / (1 - beta_2^t_temp)

      # Compute adaptive learning rate
      grad_update = learning_rate * m_beta_hat / (sqrt(v_beta_hat) + epsilon)
      # Tentative updates for beta and loop_dev
      beta_new_temp = beta_new + grad_update

      # Dispersion gradient update (if applicable)
      if (family == "spot negative binomial") {
        m_disp_temp = beta_1 * m_disp + (1 - beta_1) * disp_grad
        v_disp_temp = beta_2 * v_disp + (1 - beta_2) * (disp_grad^2)

        m_disp_hat = m_disp_temp / (1 - beta_1^t_temp)
        v_disp_hat = v_disp_temp / (1 - beta_2^t_temp)

        disp_update = learning_rate * m_disp_hat / (sqrt(v_disp_hat) + epsilon)
        dispersion_temp = dispersion + disp_update
        dispersion_temp = max(dispersion_temp, 0.05)
      } else {
        dispersion_temp = dispersion
      }

      # Compute likelihood
      if (family == "spot poisson") {
        pred = family_model[["predict"]](X, beta_new_temp, lambda, offset)$total
        lik_new = -poisson_lik(y, pred)
      } else if (family == "spot negative binomial") {
        pred = family_model[["predict"]](X, beta_new_temp, lambda, offset)$total
        lik_new = -nb_lik(y, pred, dispersion_temp)
      } else if (family == "spot binomial") {
        pred = family_model[["predict"]](X, beta_new_temp, lambda)$total
        lik_new = -binomial_lik(y, pred, weights)
      }

      lik_diff = lik_new / lik_old
      #if we are greater than some cutoff (i.e bad step) or we are already starting to converge and take a small bad step
      if ((lik_diff > 1+1e-3) | (lik_diff > 1+1e-6 & counter > 20)) {
        # If likelihood difference is too large, reduce step size and recompute gradients
        learning_rate = learning_rate * 0.5
        recompute_gradients = TRUE  # Recompute gradients with updated learning rate
      } else {
        lik_diff = min(lik_diff,1/lik_diff)
        # Accept updates if likelihood difference is acceptable
        beta_new = beta_new_temp
        dispersion = dispersion_temp
        m_beta = m_beta_temp
        v_beta = v_beta_temp
        if (family == "spot negative binomial") {
          m_disp = m_disp_temp
          v_disp = v_disp_temp
        }
        t = t_temp
        lik_old = lik_new
        recompute_gradients = FALSE  # Exit gradient recomputation loop
      }
    }

    # Update counter
    counter = counter + 1
  }

  #print if not converged
  conv = TRUE
  if(counter == max_gd_steps){
    #print("Model did not converge. Try increasing max iterations, max step size, or decreasing max difference.")
    conv = FALSE
  }
  #get standard error matrix
  V = compute_standard_errors(X,y,lambda,family,beta_new,good_beta,dispersion,weights,offset)

  #compute last likelihood
  lik = NA
  if(family == "spot poisson"){
    pred = family_model[["predict"]](X,beta_new,lambda,offset)$total
    lik = -poisson_lik(y,pred)
  }else if(family == "spot negative binomial"){
    pred = family_model[["predict"]](X,beta_new,lambda,offset)$total
    lik = -nb_lik(y,pred,dispersion)
  }else if(family == "spot binomial"){
    pred = family_model[["predict"]](X,beta_new,lambda)$total
    lik = -binomial_lik(y,pred,weights)
  }

  #rename dimnames of beta
  colnames(beta_new) = colnames(lambda)
  rownames(beta_new) = colnames(X)
  return (list(beta = beta_new,vcov = V,dispersion = dispersion,likelihood = lik,
               converged = conv,counter = counter,learning_rate = learning_rate))
}




#Gaussian spotglm
#' @export
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


#Gaussian spotglm
#' @export
spot_glm_gaussian_new = function(X,y,big_X,lambda,fix_coef,beta_0){
  #make big X
  nct = ncol(lambda)
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

  colnames(beta_0) = colnames(lambda)
  rownames(beta_0) = colnames(X)

  #return matrix
  return(list(beta = beta_0,vcov = V,sigma.sq = summary(lm1)$sigma^2))
}


#gradient descent
#' @export
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


#compute var-cov matrix
#' @export
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


#get covariance matrix for models
#' @export
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



#run spotGLM
run_spot_glm = function(y,X,lambda,family = "spot gaussian",beta_0 = matrix(0,ncol(X),ncol(lambda)),fix_coef = NULL,
                        offset = rep(0,length(y)),M = 50,weights = rep(1,length(y)),
                        max_gd_steps = 5000,learning_rate = 1,max_diff = 1-1e-6){
  t1 = Sys.time()
  LR = learning_rate
  conv = FALSE
  while(LR > 0.0001 & conv == FALSE){
    result = spot_glm(X,y,lambda,beta_0 = beta_0,family = family,fix_coef = fix_coef,offset = offset,M = M,
                      weights = weights,learning_rate = LR,
                      max_gd_steps = max_gd_steps,max_diff = max_diff)
    conv = result$converged
    LR = LR/2
  }
  t2 = Sys.time()
  #get standard errors for coefs
  std_err = rep(NA,nrow(result$vcov))
  for (j in c(1:nrow(result$vcov))){
    std_err[j] = sqrt(result$vcov[j,j])
  }
  return (list(beta_est = result$beta, stand_err_mat = matrix(std_err,nrow = ncol(X) ,byrow = F),time = t2 - t1,
               disp = result$dispersion,converged = result$converged,likelihood = result$likelihood,vcov = result$vcov,
               niter = result$counter,LR = result$learning_rate))
}





