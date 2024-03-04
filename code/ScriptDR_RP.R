# define a working directory. The needed files should be stored in this directory
#setwd("path")
rm(list=ls())
library(readxl)
library(tidyverse)
#Option 1 BLS concentration data <- inputs2.xlsx
#Option 2 US concentration data <- inputs.xlsx
#Option 3 mix BLS and US data <- inputs3.xlsx
#Option 4 mix BLS and US data + Specifgic data on virulence <- inputs4ForJEMRA.xlsx
xfile="inputs4ForJEMRA.xlsx"
step=0.1
DoseCont <- seq(0, 12, step)
# Script1 derived DR from the european report and built save_df_DR. 
# Not needed here because we have saved save_df_DR
# source("DR/script1.R") #Stochastic dose response model to be run one time
# # script2_RP.R is a modified version of script2.R from EFSA
# To get the doses as output and to include Virulence prevalence
source("script2_RP.R") #full script from initial conc to number of cases

# IF RERUN EXPO set to TRUE. If not, set to FALSE
reRunExpo <- FALSE

# Load Data
ProbViru <- read_excel(xfile,sheet="ProbVirulence")
DataBase <- read_excel(xfile,sheet="Population")
Fritsch <- read_excel(xfile,sheet="Fritsch")

# Get virulence proportions
# See VirulenceEU.R
print(ProbViru)
# MV is More Virulent
# V is Virulent
# LV is less virulent
# U is unknown (will be set to Virulent)

# Set the various variable
# Don't change names
ProbVirulenceAll   <- c(SF=1,Meat=1,Cheese=1) # ALL
ProbVirulenceMV <- c(SF=ProbViru$pMV[1], Meat = ProbViru$pMV[2], Cheese = ProbViru$pMV[3])
# Unknown are set to Virulent
ProbVirulenceV <- c( SF=ProbViru$pV[1]+ProbViru$pU[1], 
                     Meat = ProbViru$pV[2]+ProbViru$pU[2], 
                     Cheese = ProbViru$pV[3]+ProbViru$pU[3])
ProbVirulenceLV <- c(SF=ProbViru$pLV[1], Meat = ProbViru$pLV[2], Cheese = ProbViru$pLV[3])
ProbHuman <- data.frame(All = 1,
                        MV = ProbViru$pMV[4], 
                        V = ProbViru$pV[4] + ProbViru$pU[4],
                        LV = ProbViru$pLV[4])

############################################################
# EVALUATE THE EXPOSURE FOR ALL virulent groups and for each virulent groups.
# Using EFSA Data and Model
# contamfun is defined in script2_RP.R (slightly modified by RP)
# Additional parameter: prevVirul with the prevalence of each virulence group
# for the three food categories
# Used the default for EFSA
# Save the results for future use
if(reRunExpo){
  for(Virul in c("All","MV","V","LV")){
    cat("running ",Virul,"\n")
    set.seed(666)
    
    rescases=contamfun(runs=1000000, #1000000
                       shift=0,meanTemp=5.9,
                       sdTemp=2.9,
                       Mode_prop_rtime=0.3,
                       Max_prop_rtime=1.1,
                       prevVirul = get(paste0("ProbVirulence",Virul)))
    
    saveRDS(rescases, paste0("rescases",Virul,".rds"))
  }
}


############################################################
# 
# Some Functions from Pouillot et al 2015
#
##################################
NcaseLNDose <- function(r, meanlog, sdlog, i, Meals = MealsCont, Dose=DoseCont) {
  # For a given parameter r | r~ log10normal(meanlog, sdlog), 
  # and a given bin of dose, i, this function 
  # evaluate the number of cases 
  # using N(meals at dose x) * Prob(Cases|r,x) * Prob(x) 
  # protected for r > 1 (log10r > 0)
  # Meals: number of Meals for dose DoseCont
  exp(log(Meals[i]) + dnorm(r, meanlog, sdlog, log=TRUE) + 
    log((r>=0) * (1-exp(-10^Dose[i])) + 
        (r<0) * (-expm1(-10^Dose[i] * 10^r))))
}

NcaseIntLnDose <- function(meanlog,sdlog,Dose,low=-Inf,up=Inf,Print=TRUE){
  #This function integrate NcaseLNDose over r from -Inf to Inf 
  Int <- integrate(NcaseLNDose,
                   lower=low, upper=up, 
                   meanlog=meanlog, sdlog=sdlog, i=Dose,
                   rel.tol = 1E-10
  )
  if(Print) print(Int)
  return(Int$value)
}

f <- function(Mean, Sd = 1, target=0, Prop=1, Print=FALSE){
  # This function evaluate the number of cases for a given mean and sd 
  # and remove a "target". It is used in a root seeking (uniroot) function
  res <- sapply(1:length(DoseCont),
                function(x) NcaseIntLnDose(Mean,Sd,x,low=-Inf,up=Inf,Print=Print))
  return(Prop*sum(res)-target)
}

############################################################
# 
# Specify Standard deviation (see Pouillot et al, 2015)
#
############################################################
# get directly the sd from Fritsch
# Evaluate the intra group standard deviation  (see Pouillot et al, 2015)
SdVirAll <- (5/2)/qnorm(0.95) # Virulence variability
SdVirMV <- sd(Fritsch$r[Fritsch$VirulInd == "MV"])
SdVirV <- sd(Fritsch$r[Fritsch$VirulInd == "V"])
SdVirLV <- sd(Fritsch$r[Fritsch$VirulInd == "LV"])

cat("SdVirAll: ",SdVirAll,"SdVirMV: ",SdVirMV,"SdVirV: ",SdVirV,"SdVirLV: ",SdVirLV, "\n")

# SD intra group for medium (SdSM) or low (SdSL) variability (see Pouillot et al, 2015)
# Note: only SdSM will be used
# SdSH <- (2.9/2)/qnorm(0.95) # High variability
SdSM <- (1.8/2)/qnorm(0.95) # Medium variability
# SdSL <- (0.8/2)/qnorm(0.95) # Low variability (more heterogeneous group)


################################################################################
#
# Solve the function
# REPRODUCE EFSA if Virul = All
#
################################################################################

AllResults <- NULL

for(Virul in c("All","MV","V","LV")){
  
  Data <- DataBase
  
  # Read specific exposure
  rescases <- readRDS(paste0("rescases",Virul,".rds"))
  
  Res <- PerDoseSameSd <- NULL

  Data$Virulence <- Virul
  
  # Extract the prevalence for the subpopulations
  Data$Prev <- rescases$overall_prev$prev 
  # Extract the number of eating occasions for subpopulations
  Data$EO <- rescases$overall_prev$teo
  
  # The number of cases is multiplied by the proportion of this virulence 
  # type
  Data$CasesVirul <- Data$Cases * ProbHuman[[Virul]]
  
  # Choice of the SD
  Data$SdIntraGroup <- SdSM # Middle variability, as in EFSA and Pouillot et al, 2015
  Data$SdVir <- get(paste0("SdVir",Virul))
  Data$RefSdLog <- sqrt(Data$SdVir^2 + Data$SdIntraGroup^2)

  for(i in 1:nrow(Data)){ # For each sup Population
    
    Path <- rescases$overall_prev$Path[i]
    
    PrevDose <- rescases$dose$prob[rescases$dose$Path == Path] 
    # The number of Meals is multiplied by 8 (number of EFSA survey, 2008-2015)
    MealsCont <- Data$Prev[i] * Data$EO[i] * PrevDose * 8 # years of surveys

    cat("\n\n",Path,"\n")
    
    
    # Prop = 1 because the proportion of individuals is already considered in the number of Eating Occasion
    # DataCases as a function of virul
    root <- uniroot(f, interval = c(-40,-5), Sd=Data$RefSdLog[i], target=Data$CasesVirul[i], Prop=1, Print=FALSE)$root
    Res <- c(Res, root)
    RootPerDose <- sapply(1:length(DoseCont),
                          function(x) NcaseIntLnDose(root,Data$RefSdLog[i],x,low=-Inf,up=Inf,Print=FALSE))
    PerDoseSameSd <- rbind(PerDoseSameSd, RootPerDose)
    
    cat("The mean is",root, "the sd is",Data$RefSdLog[i], ", the expected number of cases is ",round(sum(PerDoseSameSd)),
        "\n")
  }
  
  Data$mean <- Res
  Data$CasesVirul <- round(Data$CasesVirul)
  Data$predictedCases <- round(rowSums(PerDoseSameSd))
  print(Data)
  AllResults <- bind_rows(AllResults, Data)
}

write.csv2(AllResults, file = "AllResults.csv")

