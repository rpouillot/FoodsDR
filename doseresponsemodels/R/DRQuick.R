#' @title Dose-response model function for listeriosis - quick version
#' 
#' @description
#' This function provides the marginal probability of invasive listeriosis
#' in a given `population` for a given `Dose` in `CFU` using the 
#' `JEMRA`, the `Pouillot`, the `Fritsch` or the `EFSA` dose-response models 
#' or the model developed within this project (`EFSAMV`,`EFSAV`,`EFSALV`) (see References).
#' 
#' @param Dose (`CFU/serving`) Dose (scalar or vector).
#' @param model either `JEMRA`, `Pouillot`, `Fritsch`, `EFSA`, `EFSAMV`,`EFSAV` or `EFSALV`
#' @param population considered population (scalar or vector).
#' @param Poisson if `TRUE`, assume that `Dose` is the mean of a Poisson distribution.
#' (actual LogNormal Poisson). If `FALSE` (default), assume that `Dose` is the actual number of bacteria.  
#' 
#' @details
#' 
#'   | Model    | Population | Characteristics              | 
#'   |----------|------------|------------------------------|
#'   | JEMRA    | 1          | Healthy population           | 
#'   | JEMRA    | 2          | Increased susceptibility     | 
#'   | Pouillot | 1          | Less than 65 years old       | 
#'   | Pouillot | 2          | More than 65 years old       | 
#'   | Pouillot | 3          | Pregnancy                    |  
#'   | Pouillot | 4          | Nonhematological Cancer      |  
#'   | Pouillot | 5          | Hematological cancer         |  
#'   | Pouillot | 6          | Renal or Liver failure       |  
#'   | Pouillot | 7          | Solid organ transplant       |  
#'   | Pouillot | 8          | Inflammatory diseases        | 
#'   | Pouillot | 9          | HIV/AIDS                     | 
#'   | Pouillot | 10         | Diabetes                     | 
#'   | Pouillot | 11         | Hear diseases                | 
#'   | Fritsch  | 1          | Highly virulent              |
#'   | Fritsch  | 2          | Medium virulent              |
#'   | Fritsch  | 3          | Hypovirulent                 |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 1          | Female 1-4 yo                |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 2          | Male 1-4 yo                  |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 3          | Female 5-14 yo               |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 4          | Male 5-14 yo                 |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 5          | Female 15-24 yo              |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 6          | Male 15-24 yo                |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 7          | Female 25-44 yo              |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 8          | Male 25-44 yo                |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 9          | Female 45-64 yo              |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 10         | Male 45-64 yo                |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 11         | Female 65-74 yo              |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 12         | Male 65-74 yo                |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 13         | Female >75 yo                |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 14         | Male >75 yo                  |
#'
#' See the parameters in the JEMRA report.
#'    
#' @return A vector of size `Dose` (if `population` is a scalar) or a matrix of
#' dimension (length of the `Dose` vector x length of the `population` vector)
#' 
#' @author Regis Pouillot
#'
#' @keywords dose-response model
#'
#' @references
#' 
#' \insertRef{EFSA2018}{doseresponsemodels}
#' 
#' \insertRef{FAO-WHO2004}{doseresponsemodels}
#' 
#' \insertRef{Fritsch2018}{doseresponsemodels}
#' 
#' \insertRef{Pouillot2015}{doseresponsemodels}
#' 
#' @seealso [doseresponsemodels::DR()],  [doseresponsemodels::DRLogNormPoisson()].
#' 
#' @export
#'
#' @note This function uses (for all model but `JEMRA`) a linear approximation (`approxfun`) 
#' from the exact [DR()] model evaluated on \eqn{Dose = c(0,10^{seq(-5,12,length=1701)})} 
#' (if `Poisson=TRUE`) or \eqn{c(0,10^{seq(0,12,length=2000)})} (if `Poisson=FALSE`).
#' Any Dose lower or higher than these ranges will lead to `NA`. 
#'
#' @examples
#' # Compare DR and DRQuick
#' cbind(DR(1:10, model="Pouillot", population = 5), 
#'       DRQuick(1:10,  model="Pouillot", population = 5))
#' DRQuick(1:10,  model="Pouillot", population = 2)
#'
DRQuick <- function(Dose, model="JEMRA", population = 1, Poisson = FALSE){
  if(length(model) > 1) stop("model is not vectorized in DRQuick")
  # if JEMRA: use the original model
  if(model == "JEMRA") {
    risk <- DR(Dose, model=model, population = population, Poisson=Poisson)
  } else {
    # else use the" approx
    risk <- matrix(NA, ncol=length(population), nrow=length(Dose))
    for(i in 1:length(population)){
      f <- listDR[[paste0("DR", model, population[i], Poisson)]]
      if(is.null(f)) stop("Bad specification of the parameters in the DR model")
      risk[,i] <- f(Dose)
    }
  }
  names <- DRParam[DRParam$Model == model,]
  colnames(risk) <- names$Characteristics[population]
  return(risk)
}
