# Complex sampling analysis of SEM models
# Daniel Oberski, 2013-09-25

lavaan.survey <- 
  function(lavaan.fit, survey.design, 
           estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
           estimator.gamma=c("default","Yuan-Bentler")) {
  
  # Not all estimators in lavaan make sense to use here, therefore matching args
  estimator <- match.arg(estimator) 
  if(estimator=="ML") warning("Estimator 'ML' will not correct standard errors and chi-square statistic.")
  estimator.gamma <- match.arg(estimator.gamma) # Smoothing Gamma or not
  
  # Names of the observed variables (same for each group)
  ov.names <- lavaanNames(lavaan.fit, type="ov", group=1)
  
  # The MP-inverse duplication matrix is handy for removing redundancy
  Dplus <- ginv(lavaan::duplicationMatrix(length(ov.names)))
  # Create a formula that includes all observed variables for svymean
  ov.formula <- as.formula(paste("~", paste(ov.names, collapse="+")))
  
  # <no. group>-sized lists that will contain the asy. covariance matrix,
  #  and sample covariance matrix and mean vector
  Gamma <- vector("list", lavaan.fit@Data@ngroups)
  sample.cov <- vector("list", lavaan.fit@Data@ngroups)
  sample.mean <- vector("list", lavaan.fit@Data@ngroups)
  
  for(g in seq(lavaan.fit@Data@ngroups)) {
    if(lavaan.fit@Data@ngroups > 1) {
      # Use survey::subset to create data groups
      survey.design.g <- 
        subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                                                      lavaan.fit@call$group, 
                                                      lavaan.fit@Data@group.label[[g]]))))
    } 
    else { # In case of no groups, just use the original survey design object.
      survey.design.g <- survey.design  
    }
    
    # Function that takes survey design and returns the Gamma & observed moments
    get.stats.design <- function(survey.design.g, sample.nobs) {
      sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
      # survey package returns the variance matrix of the (co)variances as attr:
      Gamma.cov.g <- attr(sample.cov.g, "var")
      # Remove (co)variances wrt redundant elements of covs; not used by lavaan. 
      Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
      
      # Same for mean vector
      sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
      Gamma.mean.g <- attr(sample.mean.g, "var")
      
      # Join asy. variance matrices for means and covariances
      # TODO add offdiag
      Gamma.g <- as.matrix(Matrix::bdiag(Gamma.mean.g, Gamma.cov.g))
      
      Gamma.g <- Gamma.g * sample.nobs[g] # lavaan wants nobs * Gamma.
      
      # Since the above nonparametric estimate of Gamma can be unstable, Yuan
      # and Bentler suggested a model-smoothed estimate of it, optional here:
      if(estimator.gamma == "Yuan-Bentler") {
        r <- get.residuals(lavaan.fit) # Iff these asy = 0, all will be well...
        Gamma.g <- Gamma.g + (sample.nobs[g]/(sample.nobs[g] - 1)) * (r %*% t(r))
      }
      # Get rid of attribute, preventing errors caused by lazy evaluation
      # (This has to be at the end or lazy evaluation mayhem will ensue)
      attr(sample.cov.g, "var") <- NULL
      tmp  <- as.vector(sample.mean.g)
      names(tmp) <- names(sample.mean.g)
      sample.mean.g <- tmp	
      
      list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }
    # The data may be a list of multiply imputed datasets
    if(!any(class(survey.design.g) == "svyimputationList")) {
      # If no imputations, just use usual no. observations and asy variance
      sample.nobs <- unlist(lavaan.fit@Data@nobs)
      stats <- get.stats.design(survey.design.g, sample.nobs)
    } 
    else { # In case of multiply imputed data
      # Not only can nobs differ from lavaan.fit, but also per imputation
      sample.nobs <- get.sample.nobs(survey.design.g, lavaan.fit@call$group)
      # Retrieve point and variance estimates per imputation
      stats.list <- lapply(survey.design.g[[1]], get.stats.design, sample.nobs=sample.nobs)
      m  <- length(stats.list) # no. imputation
      
      # Point estimates are average over imputations
      sample.cov.list <- lapply(stats.list, `[[`, 'sample.cov.g')
      sample.cov.g <- Reduce(`+`, sample.cov.list) / m
      cov.df <- Reduce(`rbind`, lapply(sample.cov.list, vech))
      sample.mean.list <- lapply(stats.list, `[[`, 'sample.mean.g')
      sample.mean.g <- Reduce(`+`, sample.mean.list) / m
      mean.df <- Reduce(`rbind`, sample.mean.list)
      
      # Variance estimates depend on within- and between-imputation variance:
      Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
      Gamma.between <- cov(cbind(mean.df, cov.df))
      Gamma.g <- Gamma.within + ((m + 1)/m) * Gamma.between
      
      # set stats with multiple imputation point and variance estimates
      stats <- list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }

    # Augment the list for this group
    Gamma[[g]] <- stats$Gamma.g
    sample.cov[[g]] <- stats$sample.cov.g
    sample.mean[[g]] <- stats$sample.mean.g
  } # End of loop over groups

  new.call <- lavaan.fit@call
  new.call$data <- NULL                # Remove any data argument
  new.call$sample.cov <- sample.cov    # Set survey covariances
  new.call$sample.mean <- sample.mean  # Set survey means
  new.call$sample.nobs <- sample.nobs  
  new.call$estimator <- estimator  # Always use Satorra-Bentler or WLS estimator

  if(substr(estimator, 1, 2) == "ML") { # ML, robust options
    # Set asymptotic covariance matrix of sample means and covariances
    new.call$NACOV <- Gamma      
  }
  if(estimator %in% c("WLS", "DWLS")) {
    # Weighted Least Squares, adjust the weight matrix: MP inverse of Gamma
    # Note that Gamma may be singular.
    new.call$WLS.V <- lapply(Gamma, ginv)
  }
  new.fit <- eval(new.call) # Run lavaan with the new arguments
  
  if(estimator %in% c("WLS", "DWLS")) return(new.fit) # We are done for WLS

  # For ML with robust se's, check that a possibly singular Gamma has not
  # created dependencies in the parameter estimates.
  # (Code below should really be implemented in lavaan...)
  evs.too.small <- sapply(Gamma, function(Gamma.g) {
    any(eigen(Gamma.g, only.values=TRUE)$values < .Machine$double.eps*10)
  })
  if(any(evs.too.small)) {
    V.est <- vcov(new.fit)
    if(any(eigen(V.est, only.values=TRUE)$values < (.Machine$double.eps*10))) {
      long.string  <- sprintf("Some of the standard errors may not be trustworthy.
        Some of the observed covariances or means are
        collinear, and this has generated collinearity in your
        parameter estimates.  This may be a sample size issue,
        missing data problem, or due to having too few
        clusters relative to the number of parameters. Problem
        encountered in group(s) %s",
        paste(which(evs.too.small), collapse=", "))
  
      warning(strwrap(long.string, width=9999, simplify=TRUE))#gotta love it
    }
  }

  new.fit
}

# Obtain residuals from a lavaan fit object, concatenating means w/ covariances
# (used in Yuan-Bentler correction)
get.residuals <- function(fit) {
    r  <- residuals(fit)
    c(r$mean, vech(r$cov))
}

# Obtain sample size from multiply imputed svydesign object.
# In case sample size differs over imputations, takes median over imputations.
# TODO: Does not work with multiple group yet.
get.sample.nobs  <- function(svy.imp.design, group=NULL) {
  nobs.imp <- lapply(svy.imp.design[[1]], function(des) {nrow(des$variables)})
  return(median(unlist(nobs.imp)))
}

# Use the pFsum function from the survey package to obtain p value for the 
#   overall model fit using an F reference distribution where the 
#   denominator degrees of freedom is the design degrees of freedom.  
# An anonymous reviewer for J Stat Software suggested that 
#  "in surveys with small numbers of primary sampling units this sort of 
#   correction has often improved the 
#   behavior of tests in other contexts."
# The eigenvalues of the U.Gamma matrix will be the coefficients in the 
#   mixture of F's distribution (Skinner, Holt & Smith, pp. 86-87).
pval.pFsum <- function(lavaan.fit, survey.design, method = "saddlepoint") {
  # Check that Satorra-Bentler or Satterthwaite adjustment is present
  if(!lavaan.fit@Options$test %in% 
              c("satorra.bentler", "mean.var.adjusted", "Satterthwaite")) {
    stop("Please refit the model with Satorra-Bentler (MLM) or Satterthwaite (MLMVS) adjustment.") 
  }
  test <- lavaan.fit@Fit@test
  if("UGamma.eigenvalues" %in% names(test[[2]])) {
    real.eigen.values <- test[[2]]$UGamma.eigenvalues
    return(survey::pFsum(x=test[[1]]$stat, df=rep(1, length(real.eigen.values)), 
                  a=real.eigen.values, ddf=degf(survey.design), lower.tail=FALSE,
                  method=method))
  } 
  warning("U Gamma eigenvalues not available from this version of lavaan, defaulting to Satterthwaite method.")
  if(!lavaan.fit@Options$test %in% c("mean.var.adjusted", "Satterthwaite")) {
    # (Needed in case eigenvalues are not available)
    stop("Please refit the model with Satterthwaite (MLMVS) adjustment.") 
  }
  # Eigenvalues of U%*%Gamma only available in dev version of lavaan, 
  #     if n/a use Satterthwaite adjustments (same as method="Satterthwaite"):
  pf(test[[2]]$stat / test[[2]]$df, df1=test[[2]]$df, df2=degf(survey.design), lower.tail=FALSE)
}
