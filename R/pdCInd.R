#' Factor for pdCInd
#' 
#' Originally inteded to produce an independent parametrization it actually produced a conditionally independent parametrization REFLECTING A SERIOUS MISUNDERSTANDING of what it is that is actually parametrized.  It seems to be Sigma inverse but I don't see how that jives with other stuff.
#'
#' All that follows in the way of documentation is wrong.
#' 
#' Returns a factor for a 'left' log-Cholesky object for positive-definite variance matrix with zero covariances except in the first row and column. i.e. 
#' $$
#' V = L'L
#' $$
#' with $L$ a lower-triangular matrix.
#' 
#' @param object a pdCInd object representing a positive-definite variance matrix with covariances in the first row and column but zero covariances elsewhere. 
#' @return the columns of a lower triangular Choleski factor of the positive-definite matrix as a vector. 
#' @examples
#' mat <- pdCInd(diag(1:4))
#' pdFactor(mat)
#' @export
#######################################  SUPERSEDED LATER BUT INTEGRATE DOCUMENTATION
pdFactor.pdCInd <-  
function (object) {
  # Unchanged
  # Note that pdMatrix.Symm uses pdFactor
  Ncol <- (length(object) + 1)/2
  ret <- matrix(0, Ncol,Ncol)
  diag(ret) <- exp( object[1:Ncol])
  if(Ncol > 1) ret[2:Ncol,1] <- object[(Ncol+1):length(object)]
  ret
}

#' Positive-Definite Matrix With Zero Covariances Between Predictor Random Effects
#' 
#' This function is a constructor for the \code{pdCInd} class, representing a positive-definite matrix with zero covariances except possibly in the first row and column. If the matrix associated with \code{object} is of dimension $n$, it is represented by $n + (n-1)$ unrestricted parameters representing a lower-triangular log-Cholesky decomposition. The first $n$ parameters are the logs of the diagonal elements of the matrix and the last $n-1$ components are the $n-1$ remaining elements of the lower-triangular decomposition corresponding the to the possibly non-zero covariances in the first row.
#' @export
pdCInd <-
  function (value = numeric(0), form = NULL, nam = NULL, data = sys.parent()) 
{
  # unchanged  
 object <- numeric(0)
 class(object) <- c("pdCInd", "pdMat")
 pdConstruct(object, value, form, nam, data)
}

#' Construct pdCInd object
#' 
#' This function is a constructor for a pdCInd object.
#' 
#' @param object an object inheriting from the class \code{pdCInd}, representing a positive definite matrix with zero covariances except in the first row and column.
#' @param value and option initialization value, which can be any of the following ...
#' @param form an optional one-sided linear formula specifying the row/column names for the matrix represented by \code{object}.
#' @param nam and optional vector of character strings specifying the row/column names for the matrix represented by \code{object}.
#' @param data and optional data frame i which to evaluate the variables names in \code{value} and \code{form}. ...
#' @export
pdConstruct.pdCInd <-
  function (object, value = numeric(0), form = formula(object), 
         nam = Names(object), data = sys.parent(), ...) 
{
    # START MAKING CHANGES HERE
    # note that pdConstruct.pdMat return an upper-triangular R factor, i.e. chol(value)
  val <- nlme:::pdConstruct.pdMat(object, value = value, form = form, nam = nam, data = data)
  if (length(val) == 0) {
    class(val) <- c("pdCInd", "pdMat")
    return(val)
  }
  isRmat <- function(x) all( x[row(x) > col(x)] == 0)
  if (is.matrix(val)) {
    if( isRmat(val) ){
        # extract the L parameters 
# was:       L <- R2L(val)
      # was:       value <- c(log(diag(L)), L[row(L)>1 & col(L)==1])
      L <- f2L(val)
      value <- c(log(diag(L)), L[col(L)==1 & row(L)>1])
    } else stop("matrix should be an upper triangular matrix")
    attributes(value) <- 
      attributes(val)[names(attributes(val)) != "dim"]
    class(value) <- c("pdCInd", "pdMat")
    attr(value,"invert") <- FALSE
    return(value)
  }
  Ncol <- (length(val) + 1)/2
  if (length(val) != 2*round(Ncol) - 1) {
    stop(gettextf("an object of length %d does not match a pdCInd factor (diagonal + covariances with intercept", 
                  length(val)), domain = NA)
  }
  class(val) <- c("pdCInd", "pdMat")
  val
}

#' Factor of a pdCInd object.
#' 
#' Function to compute the upper triangular factor of a pdCInd object.
#' 
#' @param object a 'pdCInd' object from which the right-triangular factor of the variance matrix it represents will be extracted
#' @return the full right-triangular factor, including zeros in the lower triangle, is returned as a vector in column order     
#' @export
pdFactor.pdCInd <- 
function (object) 
{
  invert <- attr(object,"invert")
  object <- as.vector(object)
  Ncol <- round( (length(object) +1)/2)
# was:
#   L <- matrix(0,Ncol,Ncol)
#   diag(L) <- exp( object[1:Ncol])
#   if ( Ncol > 1 ) L[row(L)>1 & col(L)==1] <- 
#     object[(Ncol+1):length(object)]
#   if(invert) c(t(solve(L))) else c(L2R(L)) 
  
  L <- matrix(0,Ncol,Ncol)
  diag(L) <- exp( object[1:Ncol])
  if ( Ncol > 1 ) L[row(L)>1 & col(L)==1] <- 
    object[(Ncol+1):length(object)]
  if(invert) c(t(solve(L))) else c(f2R(L)) 
}
#' pdMatrix method for pdCInd objects
#' 
#' This function is the pdMatrix method for pdCInd objects
#' 
#' @param object a pdCInd object
#' @param factor should an upper-triangular factor be returned of the variance matrix
#' @return a variance matrix or it upper-triangular factor.
#' @export
pdMatrix.pdCInd <-
function (object, factor = FALSE) 
{
  if (!isInitialized(object)) {
    stop("cannot extract matrix from an uninitialized object")
  }
  Ncol <- (length(object)+1)/2
  value <- array(pdFactor(object), c(Ncol, Ncol), 
                 attr(object, "Dimnames"))
  ob <- as.vector(object) # subsetting object calls pdMatrix!
  attr(value, "logDet") <- 2*sum(ob[1:Ncol])
  if (factor) value else  crossprod(value)
}

#' solve method for pdCInd objects.
#' 
#' This produces a pdCInd object corresponding to the inverse of its argument.
#' 
#' @param a the pdCInd object to invert.
#' @param b is unused but copied from \code{solve.pdLogChol}.
#' @return a pdCInd object corresponding to the matrix inverse of \code{a}.
#' @export
solve.pdCInd <-
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

#' Turns a factor of a PD matrix to a right-triangular factor
#' 
#' This function takes a factor of a positive-definite matrix ($V = A'A$
#' ) and returns the equivalent right-triangular factor ($V = R'R$).
#' @param x a factor of a positive definite matrix $V$
#' @return a right-triangular factor of the same $V$.
#' @seealso \code{f2R}
#' @export 
f2R <- function(A) {
  val <- qr.R(qr(A))
  sgn <- sign(diag(val))
  sgn * val
}
  
#' Turns a factor of a PD matrix to a left-triangular factor
#' 
#' This function takes a factor of a positive-definite matrix ($V = A'A$
#' ) and returns the equivalent left-triangular factor ($V = L'L$).
#' @param A a factor of a positive definite matrix $V$
#' @return a right-triangular factor of the same $V$.
#' @seealso \code{f2L}
#' @export 
f2L <- function(A) {
  A[] <- rev(A)
  A[] <- rev(f2R(A))
  sign(diag(A)) * A
}

### Tests
TESTS <- FALSE
if(TESTS) {
  
  
  
  library(nlme)
  library(magrittr)
  V <- diag(3:6)
  V[1,2:4] <- 1:3
  V[2:4,1] <- 1:3
  V
  eigen(V)  # check that it's pd
  # test f2L and f2R
  
  fc <- chol(V)
  (V - crossprod(fc)) %>% abs %>% max
  
  faceigen <- function(V) {
      ei <- eigen(V)
      sqrt(ei$values) * t(ei$vectors)
  }
  
  fe <- faceigen(V)
  (V - crossprod(fe)) %>% abs %>% max
  
  (V - crossprod(f2L(fe))) %>% abs %>% max
  (V - crossprod(f2R(fe))) %>% abs %>% max
  (V - crossprod(f2L(fc))) %>% abs %>% max
  (V - crossprod(f2R(fc))) %>% abs %>% max
  
  crossprod(fe) - crossprod(fc)
  crossprod(f2L(fe)) - crossprod(f2L(fc))
  f2L(fe) %>% diag %>% sign
  f2L(fe) - f2L(fc)
  f2R(fe) - f2R(fc)
  
 
  f2R(fe)
  f2R(fc)
  f2L(fe)
  f <- pdCInd(V)
  V - pdMatrix(f)
  solve(f)
  pdMatrix(f)
  solve(V)
  V-round(pdMatrix(f),6) 
  pdMatrix(solve(f))
  
  solve(V) - pdMatrix(solve(f))
  V - pdMatrix( solve(solve(f)))
  pf <- pdFactor(f)
  pf <- matrix(pf, 4)
  crossprod(pf) - pdMatrix(f)
  
  # debug(pdMatrix.pdCInd)
  pdMatrix(f)
  max( abs(pdMatrix(f)-V)) < 10*.Machine$double.eps
  pdMatrix(f)-V
  
  summary(f)
  
  
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

# test with lme
  library(spida)
  head(hs)
  hs$sex <- 1*(hs$Sex == 'Female')
fit1 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdDiag( ~ 1 + ses + sex)))
fit2 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdSymm( ~ 1 + ses + sex)),
            control = list(msVerbose=T,returnObject=T))
fit3 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdLogChol( ~ 1 + ses + sex)),
            control = list(msVerbose=T,returnObject=T))
hs$sex <- with(hs, 1*(Sex == "Male"))
fit4 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdCInd( ~ 1 + ses + sex)),
            control = list(msVerbose=T,returnObject=T,msMaxIter=1000))
fit5 <- lme(mathach ~ ses + Sex, hs,
            random = list( school = pdCInd( ~ 1 + ses + Sex)),
            control = list(msVerbose=T,returnObject=T,msMaxIter=1000))
summary(fit4)
summary(fit5)
summary(fit2)
ff4 <- (fit4$modelStruct$reStruct$school)
solve(ff4)
summary(ff4)
unclass(ff4)
summary(solve(ff4))

ff3 <- fit3$modelStruct$reStruct$school
ff2 <- fit2$modelStruct$reStruct$school
summary(ff3)
summary(solve(ff3))
summary(ff3)
summary(ff4)
AIC(fit3)
AIC(fit4)
AIC(fit2)
AIC(fit1)

unclass(ff2)
unclass(ff3)

}
