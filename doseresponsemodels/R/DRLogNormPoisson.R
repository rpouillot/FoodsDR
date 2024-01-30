#' @title Marginal probability of invasive listeriosis over r
#' 
#' @description The function [DRLogNormPoisson()] provides the marginal probability of invasive listeriosis
#' in a given population for a given `Dose` in `CFU`. this function is not vectorized.
#' 
#' @param Dose (`CFU/serving`) Dose (scalar or vector). It should be integers if `Poisson` is `FALSE`. 
#' @param meanlog10 the meanlog10 parameter of the distribution of `r` (parameter of the exponential model).
#' @param sdlog10 the sdlog10 parameter of the distribution of `r` (parameter of the exponential model).
#' @param Poisson if `TRUE`, assume that `Dose` is the mean of a Poisson distribution. (actual LogNormal Poisson).
#'  If `FALSE` (default), assume that `Dose` is the actual number of bacteria.  
#' @param low lower value for the integration.
#' @param up upper value for the integration.
#' @param silent silent the error-try function.
#' @param tol relative tolerance. Note: for `method = "cubature"`, the tolerance will be set to \eqn{1E-05}.
#' @param method either "integrate" (default) or `"cubature"` to specify the integration method.
#' @param ... further arguments to pass to the integrate function.
#' @return Probability of invasive listeriosis integrated over `r`.
#' 
#' @details The function evaluates 
#' \deqn{\int_{low}^{inf} \Phi(x, mulog_{10}, sdlog_{10})\cdot(1-e^{(-Dose \cdot 10^{r})}) dr}
#' using the \code{\link{integrate}} function, with a relative tolerance equals to \code{tol}
#' if \code{Poisson} is \code{TRUE}. If \code{Poisson} is \code{FALSE}, it evaluates 
#' \deqn{\int_{low}^{inf} \Phi(x, mulog_{10}, sdlog_{10})\cdot(1-(1-10^{r})^{Dose}) dr}. 
#' 
#' For `method = "cubature"`, the tolerance will be set to \eqn{1E-5}.
#' `method = "cubature"` will use the `\link[cubature]{hcubature}` function that
#' is much slower but guarantees a tolerance of \eqn{1E-5}.
#' 
#' @note 
#' This function is used by the [DR()] function, a wrapper of [DRLogNormPoisson()].
#' For a quick, vectorized version of it, use [DRQuick()].
#'  
#' @author Régis Pouillot
#'
#' @keywords dose-response model
#'
#' @references
#' \insertRef{Pouillot2015}{doseresponsemodels}
#'
#' @importFrom stats integrate
#' @importFrom cubature hcubature
#'
#' @export
#'
#' @examples
#' DRLogNormPoisson(2, -14.11, 1.62,low=-Inf,up=Inf)
#' DRLogNormPoisson(2, -14.11, 1.62, Poisson=TRUE)
#'
DRLogNormPoisson <- function(Dose, meanlog10, sdlog10, Poisson=FALSE, low=-Inf, up=Inf, silent=TRUE, tol=1E-20,
                             method = "integrate", ...){
  
  if(length(meanlog10) > 1 | length(sdlog10) > 1 ) stop("DRLogNormPoisson is not vectorized for meanlog10 or sdlog10")
  # If length Dose > 1, self call using sapply
  if(length(Dose) > 1) return(sapply(Dose, function(x) DRLogNormPoisson(x, 
                                                                        meanlog10=meanlog10, 
                                                                        sdlog10=sdlog10, 
                                                                        Poisson=Poisson, 
                                                                        low=low, 
                                                                        up=up, 
                                                                        silent=silent, 
                                                                        tol=tol,
                                                                        method=method,...)))
  if(Dose == 0) return(0)
  if(Poisson){
    # dnorm*(1-exp(-r*d)) <=>
    # exp(log(dnorm) + log(1-exp(-r*d))) <=>
    # exp(log(dnorm) + log(-expm1(-exp(log(log10(r)*log(10)+log(d)))))
    
    DRLNDose <- function(log10r, Dose, meanlog, sdlog) {
       exp(stats::dnorm(log10r, meanlog, sdlog, log=TRUE) + 
        log(- expm1(- exp(log(Dose) + pmin(log10r,0) * log(10)))))
    }
  } else {
    # dnorm*(1-(1-r)^d) <=>
    # exp(log(dnorm) + log(1-(1-r)^d)) <=>
    # exp(log(dnorm) + log(1-exp(d*log(1-r)))) <=>
    # exp(log(dnorm) + log(-expm1(d*(1-r))))
    
    DRLNDose <- function(log10r, Dose, meanlog, sdlog) {
      exp(stats::dnorm(log10r, meanlog, sdlog, log=TRUE) + 
        log(-expm1(Dose * log(1-10^pmin(log10r,0)))))
    }
  }
  if(method == "integrate") {
    res <- "try-error"
    printWarn <- TRUE
    while(res=="try-error"){
        Int <- try(stats::integrate(DRLNDose,
                                  lower=low, upper=up,
                                  abs.tol=tol,
                                  Dose=Dose, meanlog=meanlog10, sdlog=sdlog10, ...),silent=silent)

      res <- class(Int)
      tol <- tol * 10
      if(tol > 1E-5 & printWarn) {
        # Provide a warning if precision >1E-5
        warning("DR precision > 1E-5")
        printWarn <- FALSE
      }
      if(tol > 1E-1) {
        # Provide an error if precision >1E-1
        print(Int)
        stop(paste("Can not integrate the DR for dose",Dose,"meanlog:",meanlog10,"sdlog10:",sdlog10))
      }
  }
  } else if(method == "cubature") {
    # One try with cubature
    Int <- try(cubature::hcubature(DRLNDose,
                                   lowerLimit=low, upperLimit=up,
                                   Dose=Dose, meanlog=meanlog10, sdlog=sdlog10,
                                   vectorInterface = TRUE),silent=silent)
    res <- class(Int)
      if(class(res)[1] == "try-error"){
        print(Int)
        stop(paste("Can not integrate the DR for dose",Dose,"meanlog:",meanlog10,"sdlog10:",sdlog10))
      }
    # different name
    Int$value <- Int$integral
  } else stop("wrong specification of method in the DR")
  
    return(Int$value)
}
