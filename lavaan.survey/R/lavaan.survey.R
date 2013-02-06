get.sample.nobs  <- function(svy.imp.design, group=NULL) {
    if(is.null(group)) {
	nobs.imp <- lapply(svy.imp.design[[1]], function(des) {nrow(des$variables)})
	return(median(unlist(nobs.imp)))
    }
    # Multiple group analysis
    # TODO
}

lavaan.survey <- function(lavaan.fit, survey.design, 
	     estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
	     estimator.gamma=c("default","Yuan-Bentler")) {
  estimator <- match.arg(estimator)
  estimator.gamma <- match.arg(estimator.gamma)
  ov.names <- lavaan.fit@Data@ov.names[[1]]
  Dplus <- ginv(lavaan::duplicationMatrix(length(ov.names)))
  ov.formula <- as.formula(paste("~",paste(ov.names, collapse="+")))
  
  Gamma <- vector("list", lavaan.fit@Data@ngroups)
  sample.cov <- vector("list", lavaan.fit@Data@ngroups)
  sample.mean <- vector("list", lavaan.fit@Data@ngroups)
  
  for(g in seq(lavaan.fit@Data@ngroups)) {
    if(lavaan.fit@Data@ngroups > 1) {
      survey.design.g <- subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                              lavaan.fit@call$group, lavaan.fit@Data@group.label[[g]]))))
    } else { survey.design.g <- survey.design  }
    

    get.stats.design <- function(survey.design.g, sample.nobs) {
	sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
	Gamma.cov.g <- attr(sample.cov.g, "var")
	Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
      
	sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
	Gamma.mean.g <- attr(sample.mean.g, "var")
	Gamma.g <- as.matrix(Matrix::bdiag(Gamma.mean.g, Gamma.cov.g)) # TODO add offdiag
	
	Gamma.g <- Gamma.g * sample.nobs[g]
	
	if(estimator.gamma == "Yuan-Bentler") {
	    r <- get.residuals(lavaan.fit)
	    Gamma.g <- Gamma.g + (sample.nobs[g]/(sample.nobs[g] - 1)) * (r %*% t(r))
	}
	# This has to be at the end or lazy evaluation mayhem will ensue:
	attr(sample.cov.g, "var") <- NULL
	tmp  <- as.vector(sample.mean.g)
	names(tmp) <- names(sample.mean.g)
	sample.mean.g <- tmp	

	list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }
    if(!any(class(survey.design.g) == "svyimputationList")) {
        sample.nobs <- unlist(lavaan.fit@Data@nobs)
	stats <- get.stats.design(survey.design.g, sample.nobs)
    } else {
	# Not only can nobs differ from lavaan.fit, but also per imputation
	sample.nobs <- get.sample.nobs(survey.design.g, lavaan.fit@call$group)
	stats.list <- lapply(survey.design.g[[1]], get.stats.design, sample.nobs=sample.nobs)
	m  <- length(stats.list)
	sample.cov.list <- lapply(stats.list, `[[`, 'sample.cov.g')
	sample.cov.g <- Reduce(`+`, sample.cov.list) / m
	cov.df <- Reduce(`rbind`, lapply(sample.cov.list, vech))
	sample.mean.list <- lapply(stats.list, `[[`, 'sample.mean.g')
	sample.mean.g <- Reduce(`+`, sample.mean.list) / m
	mean.df <- Reduce(`rbind`, sample.mean.list)

	Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
	Gamma.between <- cov(cbind(mean.df, cov.df))
	Gamma.g <- Gamma.within + ((m + 1)/m) * Gamma.between
	
	stats <- list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }


    Gamma[[g]] <- stats$Gamma.g
    sample.cov[[g]] <- stats$sample.cov.g
    sample.mean[[g]] <- stats$sample.mean.g
  }
  
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
     new.call$WLS.V <- lapply(Gamma, ginv)
  }
  new.fit <- eval(new.call) # Run lavaan with the new arguments
  
  if(estimator %in% c("WLS", "DWLS")) return(new.fit) # We are done for WLS

  # Code below should really be implemented in lavaan...

  # For ML with robust se's, check that a possibly singular Gamma has not
  # created dependencies in the parameter estimates
  evs.too.small <- sapply(Gamma, function(Gamma.g) 
			  any(eigen(Gamma.g, only.values=TRUE)$values < .Machine$double.eps*10))
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


get.residuals <- function(fit) {
    r  <- residuals(fit)
    c(r$mean, vech(r$cov))
}
