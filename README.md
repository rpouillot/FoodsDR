# FoodsDR

This repository contains the R code and data necessary to derive the manuscript
**"Updated parameters for the dose-response model for *Listeria monocytogenes* considering pathogen virulence and age and sex of consumer"**. from 
Régis Pouillot, Andreas Kiermeier, Laurent Guillier, Vasco Cadavez, and Moez Sanaa.
*Foods* **2024**, *13*(5), 751 
(https://doi.org/10.3390/foods13050751), 
as well as the link to install the `doseresponsemodels` package.


## R package to use the Dose Response models

### Package installation

``` r
install.packages("devtools")
devtools::install_github("rpouillot/doseresponsemodels")
```

### Usage

The `doseresponsemodels::DRQuick()` function provides a "quick" version of
the function to derive the marginal probability of invasive listeriosis in a given population for a given dose in CFU (actual dose if the argument `Poisson = FALSE` or average dose
if the argument `Poisson = TRUE`) using the "JEMRA" 2004, the "Pouillot" *et al.*, 2015, the "Fritsch" *et al.* 2018, 
the "EFSA", 2018 dose-response models or the model developed within this project ("EFSAMV" for more virulent strains, "EFSAV" for virulent strains, or "EFSALV" for the less virulent strains).

``` r
library("doseresponsemodels")
help('DRQuick')
DRQuick(1:10,  model="JEMRA", population = 1:2)
DRQuick(1:10,  model="Pouillot", population = 1:11)
DRQuick(1:10,  model="EFSA", population = 1:14)
DRQuick(1:10,  model="EFSAMV", population = 1:14)
```

This function uses (for all model but `JEMRA`) a linear approximation (`approxfun`) 
from the exact `doseresponsemodels::DR()` function. 

The marginal dose-response over all sub-populations (population: 0) is obtained by weighting 
the estimated probabilities according to the proportion of each sub-populations
as provided by JEMRA 2004 (82.5% healthy, 17.5%
with increased susceptibility), Pouillot *et al.* 2015,
(their Table I, based on frequency of each sub-populations in France),
Fritsch *et al.* 2018, and EFSA 2018 (their Table 2, based on the number of eating occasions per sub-population).


   | Model    | Population | Characteristics              |
   |----------|------------|------------------------------|
   | JEMRA    | 0          | Marginal over sub-populations    |  
   | JEMRA    | 1          | Healthy population           |
   | JEMRA    | 2          | Increased susceptibility     |
   | Pouillot | 0          | Marginal over populations    |  
   | Pouillot | 1          | Less than 65 years old       |
   | Pouillot | 2          | More than 65 years old       |
   | Pouillot | 3          | Pregnancy                    |
   | Pouillot | 4          | Nonhematological Cancer      |
   | Pouillot | 5          | Hematological cancer         |
   | Pouillot | 6          | Renal or Liver failure       |
   | Pouillot | 7          | Solid organ transplant       |
   | Pouillot | 8          | Inflammatory diseases        |
   | Pouillot | 9          | HIV/AIDS                     |
   | Pouillot | 10         | Diabetes                     |
   | Pouillot | 11         | Heart diseases               |
   | Fritsch  | 0          | Marginal over virukence      |
   | Fritsch  | 1          | Highly virulent              |
   | Fritsch  | 2          | Medium virulent              |
   | Fritsch  | 3          | Hypovirulent                 |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 0          | Marginal over sub-populations  |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 1          | Female 1-4 yo                |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 2          | Male 1-4 yo                  |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 3          | Female 5-14 yo               |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 4          | Male 5-14 yo                 |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 5          | Female 15-24 yo              |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 6          | Male 15-24 yo                |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 7          | Female 25-44 yo              |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 8          | Male 25-44 yo                |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 9          | Female 45-64 yo              |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 10         | Male 45-64 yo                |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 11         | Female 65-74 yo              |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 12         | Male 65-74 yo                |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 13         | Female >75 yo                |
   | EFSA-EFSALV-EFSAV-EFSAMV     | 14         | Male >75 yo                  |


## Code to derive the dose response

This code provides the necessary material to derive the dose-response model
for *Listeria monocytogenes* dose-response. 

Download the [code](https://github.com/rpouillot/FoodsDR/tree/main/code) repository.

Run `ScriptDR_RP.R`

## References
EFSA (2018). “Scientific opinion on the Listeria monocytogenes contamination of ready-to-eat foods and the risk from human health in the EU.” EFSA Journal, 16(1), 5134 (https://doi.org/10.2903/j.efsa.2018.5134).

FAO-WHO (2004). “Risk assessment of Listeria monocytogenes in ready-to-eat foods: Technical report.” World Health Organization and Food and Agriculture Organization of the United Nations
(https://www.fao.org/3/y5394e/y5394e00.pdf).

Fritsch L, Guillier L, Augustin J (2018). “Next generation quantitative microbiological risk assessment: Refinement of the cold smoked salmon-related listeriosis risk model by integrating genomic data.” Microbial Risk Analysis, 10, 20–27. (https://doi.org/10.1016/j.mran.2018.06.003).

Pouillot R, Hoelzer K, Chen Y, Dennis SB (2015). “*Listeria monocytogenes* dose response revisited–incorporating adjustments for variability in strain virulence and host susceptibility.” Risk Analysis, 35(1), 90–108. (https://doi.org/10.1111/risa.12235).
