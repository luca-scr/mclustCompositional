#' Fit Gaussian mixtures to compositional data
#'
#' `MclustCompositional()` performs model-based clustering for compositional
#' data. Observations are transformed with isometric log-ratio (ILR) coordinates,
#' Gaussian mixture models are fitted with `mclust`, and densities are
#' evaluated on the original compositional scale using the transformation
#' Jacobian.
#'
#' `densityMclustCompositional()` exposes the density-estimation interface used
#' by `MclustCompositional()`.
#'
#' @param data A numeric matrix, data frame, or vector of strictly positive
#'   compositions. Observations are rows and compositional parts are columns.
#'   Rows with missing values are omitted before fitting.
#' @param criterion Model-selection criterion. One of `"ICL"` or `"BIC"` for
#'   `MclustCompositional()` and one of `"BIC"` or `"ICL"` for
#'   `densityMclustCompositional()`.
#' @param ... Additional arguments passed to `densityMclustCompositional()` or
#'   to the internal mclust fitting step.
#' @param G Integer vector with the numbers of mixture components to consider.
#'   Defaults to `1:9`, truncated to the available sample size.
#' @param modelNames Character vector of mclust covariance models. If `NULL`,
#'   univariate transformed data use `"E"` and `"V"`; multivariate transformed
#'   data use `mclust.options("emModelNames")`, with high-dimensional safeguards.
#' @param prior Optional prior specification passed to the mclust fitting step.
#' @param nstart Number of random starts used by k-means initialization when
#'   `G > 1`.
#' @param parallel Logical, numeric, character, or cluster object controlling
#'   parallel model fitting. `FALSE` fits sequentially; `TRUE` uses all detected
#'   cores; a number uses that many cores; `"snow"`/`"PSOCK"` or
#'   `"multicore"`/`"FORK"` selects the cluster type.
#' @param seed Optional integer seed used for reproducible starts.
#'
#' @return
#' `MclustCompositional()` returns an object of class
#' `"MclustCompositional"` and `"densityMclustCompositional"` containing the
#' selected model, posterior probabilities, classification, uncertainty,
#' entropy, and fitted parameters.
#'
#' `densityMclustCompositional()` returns an object of class
#' `"densityMclustCompositional"` containing the selected density model and the
#' full BIC/ICL tables.
#'
#' @details
#' Compositional observations must be strictly positive. The closure constant is
#' arbitrary because the ILR transform is invariant to multiplication of a
#' composition by a positive scalar.
#'
#' @examples
#' set.seed(1)
#'
#' y1 <- matrix(rnorm(80, mean = -2, sd = 0.25), ncol = 2)
#' y2 <- matrix(rnorm(80, mean =  2, sd = 0.25), ncol = 2)
#' x <- ilr_inverse(rbind(y1, y2))
#'
#' fit <- MclustCompositional(
#'   x,
#'   G = 1:2,
#'   modelNames = "EEE",
#'   nstart = 1,
#'   seed = 10
#' )
#'
#' fit$G
#' table(fit$classification)
#' head(predict(fit, what = "z"))
#'
#' @seealso [ilr_transform()], [predict.densityMclustCompositional()]
#' @export
MclustCompositional <- function(data, 
                                criterion = c("ICL", "BIC"), 
                                ...)
{
  mc <- match.call()
  #
  data <- na.omit(data.matrix(data))
  n <- nrow(data)
  d <- ncol(data)
  varname <- deparse(mc$data)
  if(is.null(colnames(data)))
    { if(d == 1) colnames(data) <- varname
      else       colnames(data) <- paste0(varname, seq(d)) }
  # check criterion
  criterion <- match.arg(criterion, several.ok = FALSE)
  #
  obj <- densityMclustCompositional(data, criterion = criterion, ...)
  if(is.null(obj)) return(obj)
  obj$call <- mc
  obj$density <- NULL
  obj$entropy <- with(obj, -rowSums(z * ifelse(z > 0, log(z), 0)))
  obj$nce     <- with(obj, ifelse(G == 1, 0, mean(entropy, na.rm = TRUE)/log(G)))
  class(obj) <- c("MclustCompositional", "densityMclustCompositional")
  return(obj)
}

#' @rdname MclustCompositional
#' @export
densityMclustCompositional <- function(data,
                                       G = NULL,
                                       modelNames = NULL,
                                       criterion = c("BIC", "ICL"),
                                       prior = NULL,
                                       nstart = 25,
                                       parallel = FALSE,
                                       seed = NULL,
                                       ...)
{
  mc <- match.call()
  data <- na.omit(data.matrix(data))
  n <- nrow(data)
  d <- ncol(data)
  varname <- deparse(mc$data)
  if(is.null(colnames(data)))
    { if(d == 1) colnames(data) <- varname
      else       colnames(data) <- paste0(varname, seq(d)) }

  # check G
  G <- if(is.null(G)) 1L:9L else sort(as.integer(unique(G)))
  
  # check modelNames
  if(is.null(modelNames)) 
    { if(d == 1) 
        { modelNames <- c("E", "V") }
      else 
       { modelNames <- mclust::mclust.options("emModelNames")
      if(n <= d) 
        { # select only spherical and diagonal models
          m <- match(modelNames, c("EII", "VII", "EEI", 
                                   "VEI", "EVI", "VVI"),
                     nomatch = 0)
          modelNames <- modelNames[m]
        }
    }
  }

  # check criterion
  criterion <- match.arg(criterion, several.ok = FALSE)
  
  nG <- length(G)
  nM <- length(modelNames)
  if(nG*nM < 2) parallel <- FALSE

  # Start parallel computing (if needed)
  if(is.logical(parallel))
    { if(parallel) 
        { parallel <- startParallel(parallel)
          stopCluster <- TRUE }
      else
      { parallel <- stopCluster <- FALSE } 
    }
  else
    { stopCluster <- if(inherits(parallel, "cluster")) FALSE else TRUE
      parallel <- startParallel(parallel) 
    }
  on.exit(if(parallel & stopCluster)
          stopParallel(attr(parallel, "cluster")) )
  # define operator to use depending on parallel being TRUE or FALSE
  `%DO%` <- if(parallel && requireNamespace("doRNG", quietly = TRUE)) 
               doRNG::`%dorng%`
            else if(parallel) foreach::`%dopar%` else foreach::`%do%`
  # set seed for reproducibility  
  if(is.null(seed)) seed <- sample(1e5, size = 1)
  seed <- as.integer(seed)
  set.seed(seed)
  
  # Run models fitting 
  grid <- expand.grid(modelName = modelNames, G = G)
  fit <- foreach(i = 1:nrow(grid)) %DO%
  { # fit model
    fitMixdensCompositional(data,
                            G = grid$G[i],
                            modelName = grid$modelName[i], 
                            prior = prior,
                            nstart = nstart,
                            ...)
  }
  BIC <- sapply(fit, function(mod) if(is.null(mod)) NA else mod$bic)
  ICL <- sapply(fit, function(mod) if(is.null(mod)) NA else mod$icl)
  i <- if(criterion == "BIC") 
  {
    which(BIC == max(BIC, na.rm = TRUE))[1] 
  } else
  {
    which(ICL == max(ICL, na.rm = TRUE))[1] 
  }
  
  mod <- fit[[i]]
  mod <- append(mod, list(call = mc), after = 0)
  BIC <- matrix(BIC, length(G), length(modelNames), byrow = TRUE,
                dimnames = list(G, modelNames))
  class(BIC) <- "mclustBIC"
  attr(BIC, "prior") <- mod$prior
  attr(BIC, "control") <- mod$control
  mod$BIC <- BIC
  ICL <- matrix(ICL, length(G), length(modelNames), byrow = TRUE,
                dimnames = list(G, modelNames))
  class(ICL) <- "mclustICL"
  mod$ICL <- ICL
  mod$seed <- seed
  class(mod) <- "densityMclustCompositional"
  return(mod)
}

fitMixdensCompositional <- function(data,
                                    kappa = 1,
                                    G,
                                    modelName,
                                    prior = NULL,
                                    control = mclust::emControl(),
                                    nstart = 25,
                                    warn = mclust::mclust.options("warn"),
                                    verbose = FALSE,
                                    eps = sqrt(.Machine$double.eps),
                                    ...)
{
  x <- as.matrix(data)
  kappa <- as.numeric(kappa)
  n <- nrow(x)
  D <- ncol(x)
  d <- D - 1
  G <- as.integer(G)
  modelName <- as.character(modelName)
  # set EM iterations parameters
  tol <- control$tol[1]
  itmax <- min(control$itmax[1], 1000)
  # ilr transform
  y <- ilr_transform(x)
  # initialization using k-means with given G on the transformed variables
  km <- stats::kmeans(y, centers = G, 
                      iter.max = itmax,
                      nstart = ifelse(G > 1, nstart, 1))
  z  <- mclust::unmap(km$cluster)
  # fit GMM model
  # start algorithm
  ME_step <- mclust::me(data = y, modelName = modelName, 
                        z = z, prior = prior, warn = warn, ...)
  if(is.na(ME_step$loglik))
  { 
    if(warn) warning("ME init problems...")
    ME_step$bic <- ME_step$icl <- NA
    return(ME_step) 
  }
  #
  ME_step <- c(ME_step, list(data = x))
  loglik <- do.call("tloglik", ME_step)
  if(is.na(loglik))
  { 
    if(warn) warning("EM init problems...")
    ME_step$bic <- ME_step$icl <- NA
    return(ME_step) 
  }
  
  loglik0 <- loglik - 0.5*abs(loglik)
  iter <- 1
  if(verbose) 
  { 
    cat("\nG =", G, "  Model =", modelName)
    cat("\niter =", iter, "  loglik =", loglik)
  }
  
  while((loglik - loglik0)/(1+abs(loglik)) > tol & iter < itmax)
  { 
    loglik0 <- loglik
    iter <- iter + 1
    # compute ME-step
    ME_step <- mclust::me(data = y, modelName = modelName, 
                          z = ME_step$z, 
                          prior = prior, 
                          warn = warn) # , ...)
    ME_step <- c(ME_step, list(data = x))
    loglik <- do.call("tloglik", ME_step)
    #
    if(is.na(loglik)) 
    { 
      if(warn) warning("EM convergence problems...")
      break 
    }
    #
    if(verbose)
      cat("\niter =", iter, "  loglik =", loglik)
  }
      
  # collect info & estimates  
  mod <- ME_step
  mod$data <- x
  mod$ilr <- y
  mod$loglik <- loglik
  mod$iter <- iter
  cl <- seq(G)
  names(mod$parameters$pro) <- cl
  if(d > 1)
  { 
    varnames <- paste0("ilr", seq(ncol(y)))
    dimnames(mod$parameters$mean)[1] <- list(varnames)
    dimnames(mod$parameters$variance$sigma)[1:2] <- list(varnames, varnames)
  } 
  mod$df <- nMclustParams(modelName, d, G)
  mod$bic <- 2*loglik - mod$df*log(n)
  C <- matrix(0, n, ncol(mod$z))
  for(i in 1:n) C[i, which.max(mod$z[i, ])] <- 1
  mod$icl <- mod$bic + 2 * sum(C * ifelse(mod$z > 0, log(mod$z), 0))
  mod$classification <- cl[map(mod$z)]
  mod$uncertainty <- c(1 - mclustAddons:::rowMax(mod$z))
  mod$density <- do.call("tdens", mod)
  orderedNames <- c("data", "n", "d", "modelName", "G",
                    "ilr", "loglik", "iter", "df", "bic", "icl",
                    "parameters", "z", 
                    "classification", "uncertainty", "density")
  return(mod[orderedNames])
}

# loglik for data-transformed mixture
tloglik <- function(data,
                    modelName,
                    G,
                    parameters = NULL,
                    ...)
{
  l <- sum(tdens(data = data, 
                 modelName = modelName, G = G, 
                 parameters = parameters, 
                 logarithm = TRUE, 
                 what = "dens", ...))
  return(l)
}

# density on the transformed data
tdens <- function(data,
                  modelName,
                  G,
                  parameters,
                  logarithm = FALSE,
                  what = c("dens", "cdens", "z"),
                  warn = mclust::mclust.options("warn"),
                  ...)
{
  D <- parameters$variance$d
  x <- as.matrix(data)
  what <- match.arg(what)
  pro <- parameters$pro; pro <- pro/sum(pro)
  cl <- seq(pro)
  
  # transform data
  y <- ilr_transform(data)
  # log-jacobian of transformation
  logJ <- 0.5 * ( log(rowSums(x^2)) - log(D) - rowSums(2*log(x)) )

  # compute mixture components density 
  logcden <- cdens(modelName = modelName, 
                   data = if(D > 1) y else as.vector(y),
                   logarithm = TRUE, 
                   parameters = parameters, 
                   warn = warn)
	if(attr(logcden, "returnCode") != 0) 
	  return(NA) 
  logcden <- sweep(logcden, 1, FUN = "+", STATS = logJ)
  logcden <- cbind(logcden) # drop redundant attributes
  colnames(logcden) <- cl
           
  if(what == "cdens")
  { 
    # return mixture components density
    cden <- if(logarithm) logcden else exp(logcden)
    return(cden)
  }
  
  if(what == "z")
  { 
    # return probability of belong to mixture components
    z <- mclust::softmax(logcden, log(pro))
    colnames(z) <- cl
    if(logarithm) z <- log(z)
    return(z) 
  }
  
  logden <- mclust::logsumexp(logcden, log(pro))
  den <- if(logarithm) logden else exp(logden)
  return(den)
}

#' Predict densities or posterior probabilities
#'
#' @param object A fitted object returned by [densityMclustCompositional()] or
#'   [MclustCompositional()].
#' @param newdata Optional positive compositional observations. If `NULL`, the
#'   fitted training data are used.
#' @param what Quantity to return: mixture density (`"dens"`), component
#'   densities (`"cdens"`), or posterior probabilities (`"z"`).
#' @param kappa Retained for compatibility with the original function
#'   signature. The fitted density is evaluated on the scale supplied in
#'   `newdata`.
#' @param logarithm Logical; if `TRUE`, return log-densities or log-posterior
#'   probabilities.
#' @param ... Further arguments passed to the density evaluator.
#'
#' @return A numeric vector for `"dens"` or a matrix for `"cdens"` and `"z"`.
#' @export
predict.densityMclustCompositional <- function(object,
                                               newdata = NULL,
                                               what = c("dens", "cdens", "z"),
                                               kappa = 1,
                                               logarithm = FALSE,
                                               ...)
{
  stopifnot(inherits(object, "densityMclustCompositional"))
  D <- object$d+1
  what <- match.arg(what)
  pro <- object$parameters$pro; pro <- pro/sum(pro)
  cl <- seq(pro)
  
  # ilr data
  if(is.null(newdata)) 
  {
    x <- object$data
    y <- object$ilr
  } else
  {
    x <- if(is.vector(newdata)) 
         matrix(newdata, nrow = 1) else as.matrix(newdata)
    y <- ilr_transform(x)
  }
  
  # compute mixture components log-density 
  logcden <- cdens(modelName = object$modelName, 
                   data = y, 
                   logarithm = TRUE, 
                   parameters = object$parameters)
	if(attr(logcden, "returnCode") != 0) 
	  return(NA) 

  # log-jacobian of transformation
  # old:
  # logJ <- log((sqrt(D) / kappa) * apply(x, 1, prod)^(-1 / D))
  # logJ <- 0.5*log(D) -log(kappa) - rowSums(log(x))/D
  # new:
  # logJ <- log( sqrt( rowSums(x^2) / (D * apply(x^2, 1, prod)) ) )
  logJ <- 0.5 * ( log(rowSums(x^2)) - log(D) - rowSums(2*log(x)) )

  logcden <- sweep(logcden, 1, FUN = "+", STATS = logJ)
  logcden <- cbind(logcden) # drop redundant attributes
  colnames(logcden) <- cl
  
  if(what == "cdens")
  { 
    # return mixture components density
    cden <- if(logarithm) logcden else exp(logcden)
    return(cden)
  }
  
  if(what == "z")
  { 
    # return probability of belong to mixture components
    z <- mclust::softmax(logcden, log(pro))
    colnames(z) <- cl
    if(logarithm) z <- log(z)
    return(z) 
  }
  
  logden <- mclust::logsumexp(logcden, log(pro))
  den <- if(logarithm) logden else exp(logden)
  return(den)
}

#-- Print and summary methods -----------------------------------------

catwrap <- function(x,
                    width = getOption("width"),
                    ...)
{
  cat(paste(strwrap(x, width = width, ...), collapse = "\n"), "\n", sep = "")
}

#' Print and summarize fitted compositional mixture models
#'
#' @param x A fitted object returned by [MclustCompositional()] or
#'   [densityMclustCompositional()], or a corresponding summary object.
#' @param object A fitted object returned by [MclustCompositional()] or
#'   [densityMclustCompositional()].
#' @param digits Number of significant digits to use for numeric summaries.
#' @param classification Logical; if `TRUE`, the clustering table is included
#'   in summaries of [MclustCompositional()] objects.
#' @param parameters Logical; if `TRUE`, fitted mixing probabilities, means,
#'   and variances are included in the printed summary.
#' @param ... Further arguments, currently ignored.
#'
#' @return
#' Print methods return the input object invisibly. Summary methods return an
#' object of class `"summary.MclustCompositional"` or
#' `"summary.densityMclustCompositional"`.
#' @export
print.MclustCompositional <- function(x,
                                      digits = getOption("digits"),
                                      ...)
{
  stopifnot(inherits(x, "MclustCompositional"))
  txt <- paste0("'", class(x)[1L], "' model object: (",
                x$modelName, ",", x$G, ")")
  catwrap(txt)
  cat("\n")
  catwrap("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}

#' @rdname print.MclustCompositional
#' @export
print.densityMclustCompositional <- function(x,
                                             digits = getOption("digits"),
                                             ...)
{
  stopifnot(inherits(x, "densityMclustCompositional"))
  txt <- paste0("'", class(x)[1L], "' model object: (",
                x$modelName, ",", x$G, ")")
  catwrap(txt)
  cat("\n")
  catwrap("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}

#' @rdname print.MclustCompositional
#' @export
summary.MclustCompositional <- function(object,
                                        classification = TRUE,
                                        parameters = FALSE,
                                        ...)
{
  stopifnot(inherits(object, "MclustCompositional"))
  summary_compositional_model(
    object,
    classification = classification,
    parameters = parameters,
    density = FALSE,
    summary_class = "summary.MclustCompositional"
  )
}

#' @rdname print.MclustCompositional
#' @export
summary.densityMclustCompositional <- function(object,
                                               parameters = FALSE,
                                               ...)
{
  stopifnot(inherits(object, "densityMclustCompositional"))
  summary_compositional_model(
    object,
    classification = FALSE,
    parameters = parameters,
    density = TRUE,
    summary_class = "summary.densityMclustCompositional"
  )
}

#' @rdname print.MclustCompositional
#' @export
print.summary.MclustCompositional <- function(x,
                                              digits = getOption("digits"),
                                              ...)
{
  stopifnot(inherits(x, "summary.MclustCompositional"))
  print_summary_compositional_model(x, digits = digits)
}

#' @rdname print.MclustCompositional
#' @export
print.summary.densityMclustCompositional <- function(x,
                                                     digits = getOption("digits"),
                                                     ...)
{
  stopifnot(inherits(x, "summary.densityMclustCompositional"))
  print_summary_compositional_model(x, digits = digits)
}

summary_compositional_model <- function(object,
                                        classification,
                                        parameters,
                                        density,
                                        summary_class)
{
  classification <- as.logical(classification)
  parameters <- as.logical(parameters)
  G <- object$G
  pro <- object$parameters$pro
  if(is.null(pro)) {
    pro <- 1
  }
  names(pro) <- seq_len(G)
  mean <- object$parameters$mean
  if(object$d > 1L) {
    sigma <- object$parameters$variance$sigma
  } else {
    sigma <- rep(object$parameters$variance$sigmasq, G)[seq_len(G)]
    names(sigma) <- names(mean)
  }

  printClassification <- isTRUE(classification) && !isTRUE(density)
  classification <- if(printClassification) {
    factor(object$classification, levels = seq_len(G))
  } else {
    NULL
  }
  title <- if(isTRUE(density)) {
    "Density estimation via GMMs on ILR for compositional data"
  } else {
    "Clustering via GMMs on ILR for compositional data"
  }

  obj <- list(title = title,
              objectName = class(object)[1L],
              n = object$n,
              d = object$d,
              G = G,
              modelName = object$modelName,
              loglik = object$loglik,
              df = object$df,
              bic = object$bic,
              icl = object$icl,
              pro = pro,
              mean = mean,
              variance = sigma,
              prior = attr(object$BIC, "prior"),
              printParameters = parameters,
              printClassification = printClassification,
              classification = classification)
  class(obj) <- summary_class
  obj
}

print_summary_compositional_model <- function(x,
                                              digits = getOption("digits"))
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))),
               collapse = "")
  catwrap(txt)
  catwrap(x$title)
  catwrap(txt)

  cat("\n")
  catwrap(paste0(
    x$objectName, " ", x$modelName, " (",
    mclust::mclustModelNames(x$modelName)$type, ") model with ",
    x$G, ifelse(x$G > 1L, " components", " component"), ":"
  ))
  cat("\n")

  if(!is.null(x$prior)) 
  {
    catwrap(paste0(
      "Prior: ", x$prior$functionName, "(",
      paste(names(x$prior[-1L]), x$prior[-1L], sep = " = ",
            collapse = ", "),
      ")"
    ))
    cat("\n")
  }

  tab <- data.frame("log-likelihood" = x$loglik,
                    "n" = x$n,
                    "df" = x$df,
                    "BIC" = x$bic,
                    "ICL" = x$icl,
                    check.names = FALSE)
  print(tab, row.names = FALSE, digits = digits)

  if(x$printClassification) 
  {
    cat("\nClustering table:")
    print(table(x$classification), digits = digits)
  }

  if(x$printParameters) 
  {
    cat("\nMixing probabilities:\n")
    print(x$pro, digits = digits)
    cat("\nMeans:\n")
    print(x$mean, digits = digits)
    cat("\nVariances:\n")
    if(x$d > 1L) 
    {
      for(g in seq_len(x$G)) {
        cat("[,,", g, "]\n", sep = "")
        print(x$variance[, , g], digits = digits)
      }
    } else {
      print(x$variance, digits = digits)
    }
  }

  invisible(x)
}

