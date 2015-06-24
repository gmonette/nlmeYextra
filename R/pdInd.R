#' Factor for pdInd
#'
#' Returns a factor for a 'left' log-Cholesky object for positive-definite variance matrix with zero covariances except in the first row and column. i.e. 
#' $$
#' V = L'L
#' $$
#' with $L$ a lower-triangular matrix.
#' 
#' @param object a pdInd object representing a positive-definite variance matrix with covariances in the first row and column but zero covariances elsewhere. 
#' @return the columns of a lower triangular Choleski factor of the positive-definite matrix as a vector. 
#' @examples
#' mat <- pdInd(diag(1:4))
#' pdFactor(mat)
pdFactor.pdInd <-  
function (object) {
  # Note that pdMatrix.Symm uses pdFactor
  Ncol <- (length(object) + 1)/2
  ret <- matrix(0, Ncol,Ncol)
  diag(ret) <- exp( object[1:Ncol])
  if(Ncol > 1) ret[2:Ncol,1] <- object[(Ncol+1):length(object)]
  ret
}

#' Positive-Definite Matrix With Zero Covariances Between Predictor Random Effects
#' 
#' This function is a constructor for the \code{pdInd} class, representing a positive-definite matrix with zero covariances except possibly in the first row and column. If the matrix associated with \code{object} is of dimension $n$, it is represented by $n + (n-1)$ unrestricted parameters representing a lower-triangular log-Cholesky decomposition. The first $n$ parameters are the logs of the diagonal elements of the matrix and the last $n-1$ components are the $n-1$ remaining elements of the lower-triangular decomposition corresponding the to the possibly non-zero covariances in the first row.
#' @export
#' 
pdInd <-
  function (value = numeric(0), form = NULL, nam = NULL, data = sys.parent()) 
{
 object <- numeric(0)
 class(object) <- c("pdInd", "pdMat")
 pdConstruct(object, value, form, nam, data)
}

#' Construct pdInd object
#' 
#' This function is a constructor for a pdInd object.
#' 
#' @param object an object inheriting from the class \code{pdInd}, representing a positive definite matrix with zero covariances except in the first row and column.
#' @param value and option initialization value, which can be any of the following ...
#' @param form an optional one-sided linear formula specifying the row/column names for the matrix represented by \code{object}.
#' @param nam and optional vector of character strings specifying the row/column names for the matrix represented by \code{object}.
#' @param data and optional data frame i which to evaluate the variables names in \code{value} and \code{form}. ...
#' @export
pdConstruct.pdInd <-
  function (object, value = numeric(0), form = formula(object), 
         nam = Names(object), data = sys.parent(), ...) 
{
    # note that pdConstruct.pdMat return an upper-triangular R factor, i.e. chol(value)
  val <- nlme:::pdConstruct.pdMat(object, value, form, nam, data)
  if (length(val) == 0) {
    class(val) <- c("pdInd", "pdMat")
    return(val)
  }
  isRmat <- function(x) all( x[row(x) > col(x)] == 0)
  if (is.matrix(val)) {
    if( isRmat(val) ){
        # extract the L parameters 
        L <- R2L(val)
        value <- c(log(diag(L)), L[row(L)>1 & col(L)==1])
    } else stop("matrix should be an upper triangular matrix")
    attributes(value) <- 
      attributes(val)[names(attributes(val)) != "dim"]
    class(value) <- c("pdInd", "pdMat")
    attr(value,"invert") <- FALSE
    return(value)
  }
  Ncol <- (length(val) + 1)/2
  if (length(val) != 2*round(Ncol) - 1) {
    stop(gettextf("an object of length %d does not match a pdInd factor (diagonal + covariances with intercept", 
                  length(val)), domain = NA)
  }
  class(val) <- c("pdInd", "pdMat")
  val
}

#' Factor of a pdInd object.
#' 
#' Function to compute the upper triangular factor of a pdInd object.
#' 
#' @param object a 'pdInd' object from which the right-triangular factor of the variance matrix it represents will be extracted
#' @return the full right-triangular factor, including zeros in the lower triangle, is returned as a vector in column order     
#' @export
pdFactor.pdInd <- 
function (object) 
{
  invert <- attr(object,"invert")
  object <- as.vector(object)
  Ncol <- round( (length(object) +1)/2)
  L <- matrix(0,Ncol,Ncol)
  diag(L) <- exp( object[1:Ncol])
  if ( Ncol > 1 ) L[row(L)>1 & col(L)==1] <- 
      object[(Ncol+1):length(object)]
  if(invert) c(t(solve(L))) else c(L2R(L)) 
}
#' pdMatrix method for pdInd objects
#' 
#' This function is the pdMatrix method for pdInd objects
#' 
#' @param object a pdInd object
#' @param factor should an upper-triangular factor be returned of the variance matrix
#' @return a variance matrix or it upper-triangular factor.
pdMatrix.pdInd <-
function (object, factor = FALSE) 
{
  if (!isInitialized(object)) {
    stop("cannot extract matrix from an uninitialized object")
  }
  Ncol <- Dim(object)[2]
  value <- array(pdFactor(object), c(Ncol, Ncol), 
                 attr(object, "Dimnames"))
  ob <- as.vector(object) # subsetting object calls pdMatrix!
  attr(value, "logDet") <- 2*sum(ob[1:Ncol])
  if (factor) value else  crossprod(value)
}
#' solve method for pdInd objects.
#' 
#' This produces a pdInd object corresponding to the inverse of its argument.
#' 
#' @param a the pdInd object to invert.
#' @param b is unused but copied from \code{solve.pdLogChol}.
#' @return a pdInd object corresponding to the matrix inverse of \code{a}.
solve.pdInd <-
function (a, b, ...) 
  {
   if (!isInitialized(a)) {
     stop("cannot get the inverse of an uninitialized object")
   }
   attr(a, 'invert') <- !attr(a, 'invert')
   a
#    Ncol <- (length(a) + 1)/2
#    ob <- as.vector(a)
#    if( Ncol == 1) ret <- -ob[1]
#    else ret <- 
#      c( -ob[1:Ncol] ,
#        - exp(ob[1])*ob[(Ncol+1):length(ob)]/exp(ob[2:Ncol]))
#    attributes(ret) <- attributes(a)
#    ret
}


### Tests
TESTS <- FALSE
if(TESTS) {
  V <- diag(3:6)
  V[1,2:4] <- 1:3
  V[2:4,1] <- 1:3
  V
  eigen(V)
  f <- pdInd(V)
  f
  max( abs(pdMatrix(f)-V)) < 10*.Machine$double.eps
  
  
  
  
  fac <- pdMatrix(f, factor = TRUE)
  fac
  fac2 <- pdFactor(f)
  fac2
  all(c(fac) == fac2)
  
  max(abs(crossprod( pdMatrix(f, factor = TRUE)) - V)) < 10*.Machine$double.eps
  
  finv <- solve(f)
  unclass(f)
  unclass(finv)
  f.mat <- pdMatrix(f)
  finv.mat <- pdMatrix(finv)
  zapsmall(f.mat)
  zapsmall(finv.mat)
  zapsmall(f.mat %*% finv.mat) 
  diag(f.mat)
  diag(finv.mat)
}
