#' @title Dose-response model function for listeriosis
#' 
#' @description The [DR()] function provides the marginal probability of invasive listeriosis
#' in a given `population` for a given `Dose` in CFU using the 
#' `JEMRA`, the `Pouillot`, the `Fritsch`, the `EFSA` dose-response models
#' or the model developed within this project (`EFSAMV`,`EFSAV`,`EFSALV`) (see References).
#' @param Dose average dose in `CFU/serving` assuming a Poisson distribution from serving to serving (scalar or vector)
#' @param model either `JEMRA`, `Pouillot`, `Fritsch`, `EFSA`, `EFSAMV`,`EFSAV` or `EFSALV`
#' @param population considered population (scalar or vector) (see below)
#' @param Poisson if `TRUE`, assume that `Dose` is the mean of a Poisson distribution.
#' (actual LogNormal Poisson). If `FALSE` (default), assume that `Dose` is the actual number of bacteria. 
#' @param method either `integrate` or `cubature` to specify the integration method.
#' @param ... further arguments to pass to the [DRLogNormPoisson()] function
#' @details
#' 
#' For the `Pouillot`, `Fritsch`, `EFSA` `EFSAMV`,`EFSAV` or `EFSALV` model,
#' it integrates a log10Normal Poisson distribution \insertCite{Pouillot2015;textual}{doseresponsemodels}. 
#' or a log10Normal-binomial distribution. For the `JEMRA` model
#' an exponential dose-response model or a binomial model is used. 
#' 
#' This function is slow. Use [DRQuick()] in production.
#' 
#'   | Model    | Population | Characteristics              | 
#'   |----------|------------|------------------------------|
#'   | JEMRA    | 0          | Marginal over subpopulations | 
#'   | JEMRA    | 1          | Healthy population           | 
#'   | JEMRA    | 2          | Increased susceptibility     | 
#'   | Pouillot | 0          | Marginal over subpopulations | 
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
#'   | Fritsch  | 0          | Marginal over subpopulations | 
#'   | Fritsch  | 1          | Highly virulent              |
#'   | Fritsch  | 2          | Medium virulent              |
#'   | Fritsch  | 3          | Hypovirulent                 |
#'   | EFSA-EFSALV-EFSAV-EFSAMV     | 0          | Marginal over subpopulations |
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
#' `method = "cubature"` will use the \link[cubature]{hcubature} function that
#' is much slower but guarantees a tolerance of \eqn{1E-5}.
#' 
#' @return A vector of size Dose (if `population` is a scalar) or a matrix of
#' dimension (length of the `Dose` vector times length of the `population` vector)
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
#'
#' @export
#'
#' @examples
#' DR(5:10, "Pouillot", 1:11)
#' DR(5:10, "Pouillot", 0)
#' 
DR <- function(Dose, model="JEMRA", population = 1, Poisson = FALSE, method="integrate",...){
  
  Res <- matrix(NA, ncol=length(population), nrow=length(Dose))
  for(i in 1:length(population)){
    
    if(population[i] == 0){
      # Deals with marginal over population
      if(model == "JEMRA") {
        # From the report
        pop <- c(1:2) 
        npop <- c(.175, 1-.175)
      } else if(model == "Pouillot") {
        pop <- 1:11
        # Population size from Goulet et al
        npop <- c(48909403, 7038068, 774000, 2065000, 160000, 284000, 25300, 300674, 120000, 2681000, 1400000)
        npop <- npop/sum(npop)
      } else if(model == "Fritsch") {
        # Evaluated from EFSA: sum(RTE * Virulence proportion)
        pop <- 1:3 
        npop <- c(0.205318301,	0.258047568,	0.536634131) 
        } else if(substr(model, 1, 4) == "EFSA"){
          # EFSA
          pop <- 1:14
          npop <- c(2903379572, 3072394715, 6636582087, 7577303434, 6266862469, 8739815828,
                    18101808476, 24976657631, 20243930616, 26798995233, 8978420743, 9748354305, 
                    10090156085, 9056135019)
          npop <- npop/sum(npop)
        } else stop("Bad specification of model in DR")
      # Find the DR for all populations
      DRp <- DR(Dose, model=model, population = pop, Poisson = Poisson, method=method,...) 
      # weighted sum per subpopulations 
      Res[,i] <- DRp %*% matrix(npop, ncol=1)
      
    } else {
      # For subpopulations
      param <- DRParam[DRParam$Model == model & DRParam$Population == population[i],]
      if(nrow(param) !=1) stop("Bad specification of the parameters in the DR model")
      meanLog10 <- param$meanLog10 
      sdLog10 <- param$sdLog10
      if(model == "JEMRA") {
        if(Poisson){
          Res[,i] <- -expm1(-Dose * 10^meanLog10)
        } else {
          Res[,i] <- 1-(1-10^meanLog10)^Dose 
        }
        Res[Dose < 0 ,i] <- NA
      } else {
        Res[,i] <- DRLogNormPoisson(Dose, 
                                    meanlog10 = meanLog10, 
                                    sdlog10 = sdLog10,
                                    Poisson = Poisson,
                                    method=method,...)
      }
    }
  }
  
  names <- DRParam[DRParam$Model == model,]
  colnames(Res) <- names$Characteristics[population]
  return(Res)
}

