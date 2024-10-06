#Stochastic model - Salmonella Dublin in HRO

####PACKAGE DESCRIPTION####
#Packages for modeling
library(deSolve) #Used to solve ordinary differential equations (ODEs). More info at https://cran.r-project.org/web/packages/deSolve/deSolve.pdf
library(mc2d) #Used to apply Monte-Carlo Simulations
#Packages for analysis
library(epiR) #Used to apply Partial rank correlation analysis. More info at https://cran.r-project.org/web/packages/epiR/epiR.pdf
library(EnvStats) #Used to perform statistical analysis. More info at https://cran.r-project.org/web/packages/EnvStats/EnvStats.pdf
library(tidyr) #Used to organize data. More info at https://cran.r-project.org/web/packages/tidyr/tidyr.pdf
library(tidyverse) #Set of packages used to organize data, process and visualize data (also includes ggplot2). More info at https://cran.r-project.org/web/packages/tidyverse/tidyverse.pdf
library(data.table) #Used to design data tables
library(TruncatedDistributions)
#Packages for visualization of results
library(scales) #Add commas in the plots
library(ggplot2) #Used to plot results from the model. More info at https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf
library(gridExtra) #Allows to include multiple plots in a figure.
library(ggpubr) #Used to stylize plots for publication. More info at https://cran.r-project.org/web/packages/ggpubr/ggpubr.pdf
library(hrbrthemes) #Extra functions for ggplot2
library(RColorBrewer) #Additional color palettes for plotting
library(wacolors) #More color palettes for plotting
#library(caret) #Used to design classification/regression trees. More info at https://topepo.github.io/caret/
#Packages for development of shiny app
library(plotly) #Package used to build interactive plot maps (intended for developing of shiny app). More info at https://plotly.com/r/getting-started/
library(beepr)
library(paletteer)
library(see) #half violinplots
set.seed(100) #Seeding
initial.time=Sys.time() #Stores the time the simulation begins (this is then used to determine how long it takes for the model tu run)

#MODEL OVERVIEW#
#This model describes the transmission dynamics of Salmonella Dublin in a Northeastern US heifer-raising operation. The model considers
#that animales enter the operation at already weaned (56 days old) and return to their farms when they are close to calving (640 days old).
#The model consider different age categories (weaned calves, growing heifers, and pregnant heifers). Each age category differs in their probability
#of becoming clinically ill, S. Dublin-related death rate,  and fecal production
#(relevant for S. Dublin shedding). There are 5 disease states considered in the model: Susceptible (S), Asymptomatic (A), Clinically ill (I), 
#Carrier (C), and Recovered (R). Individuals become infected (A or I) through contact with S. Dublin in the environment (density-dependent indirect transmission).
#For more detail about disease dynamics, check the article draft. Additionally, compartments were separated in raising stages to better estimate the moment 
#at which individuals died due to infection for economic purposes (i.e., knowing the approximate moment in time when they died allow us to estimate the costs of
#raising an individual until it died). There are 8 raising stages: 1 for weaned calves, 4 for growing heifers, and 3 for pregnant heifers. Individuals move from one raising
#stage to the other every 73 days (this is an arbitrary number, it doesn't have any significant in dairy production)

#NOTATION#
#Compartments were coded as follows: State_Age.Category_RaisingStage. 
#For example, A_growing_3 (Asymptomatic growing heifer in raising stage 3)

####PARAMETERS####
n <- 1000 #number of iterations for the Monte-Carlo simulation

duration<-365*2+1 #length of  simulation in days (2 years)

#Heifers
No<- 1000 #Number of heifers in the operation at the start of the simulation.
#Baseline 1000 heifers
#For model CALIBRATION, 148 (154*0.96) heifers + 19 (148*(1/8)) calves = 167. Data for calibration obtained from Nielsen et al. 2013 (maximum herd size is 154, I used this value because we are not considering dairy cows in the calibration process (indicated as Cws in the paper). Herd size was estimated based on Dr. Thomsen who indicated that the ratios of heifers to dairy cows is 0.96 to 1:00.
#For model VALIDATION, 219 heifers (average calculated for dairy farm in Kent et al. 2021). 
#Data for validation was obtained from Kent et al. 2021 and personal communication with one of the co-authors (Dr. Andrea Lear)

asymptomatic_introduced<-1  #Number of infected asymptomatic weaned calves (56-days old) introduced into the herd at the beginning of the simulation
carriers_introduced<-0 #Number of carrier weaned calves (56-days old) introduced into the herd at the beginning of the simulation

#Set the starting values (these values indicate the initial population values; i.e., the number of individuals that would be
#in each compartment at the beginning of the simulation (time=0))

SALL_var <-c(
  ###Raising stage 
  S_weaned_1_init=(No/12) - asymptomatic_introduced - carriers_introduced, #Susceptible weaned heifer
  A_weaned_1_init=asymptomatic_introduced,                                #Asymptomatic weaned heifer
  I_weaned_1_init=0,                                                      #Clinically ill weaned heifer
  C_weaned_1_init=carriers_introduced,                                    #Carrier weaned heifer
  R_weaned_1_init=0,                                                      #Recovered weaned heifer
  
  ###Raising stage 
  S_growing_2_init=(No/12),                                         #Susceptible growing heifer 
  A_growing_2_init=0,                                              #Asymptomatic growing heifer
  I_growing_2_init=0,                                              #Clinically ill growing heifer
  C_growing_2_init=0,                                              #Carrier growing heifer
  R_growing_2_init=0,                                              #Carrier growing heifer
  
  ###Raising stage 
  S_growing_3_init=(No/12),                                         #Susceptible growing heifer
  A_growing_3_init=0,                                              #Asymptomatic growing heifer
  I_growing_3_init=0,                                              #Clinically ill growing heifer
  C_growing_3_init=0,                                              #Carrier growing heifer
  R_growing_3_init=0,                                              #Recovered growing heifer
  
  ###Raising stage 
  S_growing_4_init=(No/12),                                         #...
  A_growing_4_init=0,
  I_growing_4_init=0,
  C_growing_4_init=0,
  R_growing_4_init=0,
  
  ###Raising stage 
  S_growing_5_init=(No/12),                                         #...
  A_growing_5_init=0,
  I_growing_5_init=0,
  C_growing_5_init=0,
  R_growing_5_init=0,
  
  ###Raising stage 
  S_growing_6_init=(No/12),                                         #...
  A_growing_6_init=0,
  I_growing_6_init=0,
  C_growing_6_init=0,
  R_growing_6_init=0,
  
  ###Raising stage
  S_growing_7_init=(No/12),                                         #...
  A_growing_7_init=0,
  I_growing_7_init=0,
  C_growing_7_init=0,
  R_growing_7_init=0,
  
  ###Raising stage 
  S_growing_8_init=(No/12),                                         #...
  A_growing_8_init=0,
  I_growing_8_init=0,
  C_growing_8_init=0,
  R_growing_8_init=0,
  

  ###Raising stage
  S_pregnant_9_init=(No/12), #Susceptible pregnant heifer
  A_pregnant_9_init=0,                                             #Asymptomatic pregnant heifer
  I_pregnant_9_init=0,                #Clinically ill pregnant heifer
  C_pregnant_9_init=0,                     #Carrier pregnant heifer
  R_pregnant_9_init=0,                                             #Recovered pregnant heifer
  
  ###Raising stage 
  S_pregnant_10_init=(No/12), #Susceptible pregnant heifer
  A_pregnant_10_init=0,                                             #Asymptomatic pregnant heifer
  I_pregnant_10_init=0,                #Clinically ill pregnant heifer
  C_pregnant_10_init=0,                     #Carrier pregnant heifer
  R_pregnant_10_init=0,                                             #Recovered pregnant heifer
  
  ###Raising stage 
  S_pregnant_11_init=(No/12), #Susceptible pregnant heifer
  A_pregnant_11_init=0,                                             #Asymptomatic pregnant heifer
  I_pregnant_11_init=0,                #Clinically ill pregnant heifer
  C_pregnant_11_init=0,                     #Carrier pregnant heifer
  R_pregnant_11_init=0,                                             #Recovered pregnant heifer
  
  ###Raising stage 
  S_pregnant_12_init=(No/12),#(No/12), #Susceptible pregnant heifer ##Add 65 here for validation (65 individuals)
  A_pregnant_12_init=0,                                             #Asymptomatic pregnant heifer
  I_pregnant_12_init=0,                #Clinically ill pregnant heifer
  C_pregnant_12_init=0,                     #Carrier pregnant heifer
  R_pregnant_12_init=0,                                             #Recovered pregnant heifer
  
  ###Environment
  E_weaned_1_init=0, #Environment compartment
  
  E_growing_2_init=0, #Environment compartment
  E_growing_3_init=0, #Environment compartment
  E_growing_4_init=0, #Environment compartment
  E_growing_5_init=0, #Environment compartment
  E_growing_6_init=0, #Environment compartment
  E_growing_7_init=0, #Environment compartment
  E_growing_8_init=0, #Environment compartment
  
  E_pregnant_9_init=0, #Environment compartment
  E_pregnant_10_init=0, #Environment compartment
  E_pregnant_11_init=0, #Environment compartment
  E_pregnant_12_init=0, #Environment compartment
  
  Abortions_init=0,
  
  ###Outputs (these are not really compartments but were included to output results for analysis)
  
  D_weaned_1_init=0,                                            #Weaned heifers that died in raising stage 1
  D_growing_2_init=0,                                  #Growing heifers that died in raising stage 2
  D_growing_3_init=0,                                  #Growing heifers that died in raising stage 3
  D_growing_4_init=0,                                  #Growing heifers that died in raising stage 4
  D_growing_5_init=0,                                  #Growing heifers that died in raising stage 4
  D_growing_6_init=0,                                  #Growing heifers that died in raising stage 4
  D_growing_7_init=0,                                  #Growing heifers that died in raising stage 4
  D_growing_8_init=0,                                  #Growing heifers that died in raising stage 4
  
  Carriers_leaving_init=0,                                    #Outputs the number of carriers leaving the operation
  Asymptomatic_leaving_init=0,
  
  Pregnant_sold_init=0,                                       #Number of pregnant heifer sold
  Sacrificed_init=0,                                          #Number of heifers that were sacrificed if sick when leaving the operation
  
  Treatment_weaned_init=0,                                    #Number of calves that became clinically ill (for treatment cost estimation)
  Treatment_growing_and_pregnant_init=0,                      #Number of growing and pregnant heifers that became clinically ill (for treatment cost estimation)
  
  Completed_raising.stage_1_init=0,                           #Number of calves that completed raising stage 1                 
  Completed_raising.stage_2_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_3_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_4_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_5_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_6_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_7_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_8_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_9_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_10_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_11_init=0,                           #Number of growing heifers that completed raising stage 1 
  Completed_raising.stage_12_init=0,                           #Number of growing heifers that completed raising stage 1 
  
  Event_trigger_init=0                                        #This compartment was used as a Event_trigger to activate the root function (explained in more detail below). 
  #In a few words, "Event_trigger" tracks the days in the model. When a multiple of 73 is reached, the root function activates and 
  #animals move from one raising stage to the next (this was done to force individuals to move to the next raising stage
  #every 73 days). 
)

# # ###To extract parameters
# List1 = list()
# List2 = list()
# List3 = list()
# List4 = list()
# List5 = list()
# List6 = list()
# List7 = list()
# List8= list()
# List9= list()
# List10= list()
# List11= list()
# List12= list()
# List13= list()
# List14= list()

###Environmental compartment
#Amount of feces produced by an individual (g day-1). Value changes depending on raising stage.
f1<-7400 #Amount of feces produced by a calf - raising stage 1
f2<-9200 #Amount of feces produced by a growing heifer - raising stage 2
f3<-11000 #Amount of feces produced by a growing heifer - raising stage 3
f4<-12700 #Amount of feces produced by a growing heifer - raising stage 4
f5<-14500 #Amount of feces produced by a pregnant heifer - raising stage 5
f6<-16300 #Amount of feces produced by a pregnant heifer - raising stage 6
f7<-18100 #Amount of feces produced by a pregnant heifer - raising stage 7
f8<-19900 #Amount of feces produced by a pregnant heifer - raising stage 8
f9<-21700 #Amount of feces produced by a pregnant heifer - raising stage 9
f10<-23500 #Amount of feces produced by a pregnant heifer - raising stage 10
f11<-25200 #Amount of feces produced by a pregnant heifer - raising stage 11
f12<-27000 #Amount of feces produced by a pregnant heifer - raising stage 12

#Indirect transmission rate (Salmonella CFU-1 individual-1 day-1)
beta<- 10^-10.0 #baseline: 10^-10

#How many times infection is reduced 
kappa<- 0.005

#In the model, two types of scenarios were considered for cleaning. 1) The operation uses a free-stall system in which the barn alleys are scraped at certain frequency per day.
#2) The operation uses a deep bedding system in which bedding material is applied multiple times per day and cleaning is done only once per day (because bedding material is
#added during the day, then cattle is less exposed to manure with pathogens because manure is constantly covered with new bedding)

#Proportion of feces removed daily (% day-1)
Fk<-0.474 #Baseline: 0.474 (Clean 1x per day)
#% of feces removed every time the skid steer cleans the barn.
#If alley scraper is< included, then cleaning 2x per day=0.737, cleaning 4x per day=0.868, 6x per day=0.912, 8x day=0.934, 12x day= 0.960
#For calibration, we assumed that dairy farms clean 4x per day with an alley scraper=0.868
#For validation, we assumed that dairy farms clean 4x per day=0.868 (Flushing 4x per day).

freq<- 1/7 #Cleaning frequencies considered were 4x, 3x, 2x, 1x per week (baseline = 1x week)
#For deep cleaning, freq=1/30 (clean once per month)

#Rate of S. Dublin removal by cleaning (day-1)
mu<- (-log(1-Fk))*freq #This formula is used to transform the proportion of feces cleaned per day into a rate. 

#Natural mortality rate (days-1)
rho<-0#0.016/365

#Period in the R compartment (days)
D_R<-rtexp(n, rate = 1/140, a = 140, b = Inf)

#Immunity loss rate (day-1)
omega<- 1/D_R

#Period in the A and I compartments
D_I<-rpert(n,min=3, mode=17, max=65)

#Recovery rate from infection (day-1)
gamma<- 1/D_I

#Rate of individuals in A becoming carriers (day-1)
m<- 0.015/D_I  

#Probability of infected weaned calves becoming severely infected (dimensionless)
u1<- rpert(n, min=0.05,mode=0.25,max=0.5)

#Rate of individuals in I becoming carriers (day-1)
w<- 0.18/D_I

#Rate of weaned calves dying if severely infected (day-1) #We assume individuals are treated if sick
d1<- rpert(n, min=0.10,mode=0.2,max=0.5)/D_I #

#Vaccine effect in reducing the probability of death (SD vaccine)
v<- 1-rpert(n,min=1 ,mode = 1, max=1) #If administered, change values for min=0.05 ,mode = 0.2, max=0.8. If not administered, all values should be set to 1.

#Period in the C compartment (days)
D_C<- rpert(n, min=240, mode=365, max = 1095)

#Recovery rate from carrier state (day-1)
eta<- 1/D_C

#Probability of infected growing heifer becoming severely infected (dimensionless)
u2<- rpert(n, min=0.05,mode=0.15,max=0.3)

#Rate of growing heifers dying if severely infected (day-1) #We assume individuals are treated if sick
d2<- rpert(n, min=0.02,mode=0.1,max=0.3)/D_I 

#Probability of infected pregnant heifer becoming severely infected (dimensionless)
u3<- rpert(n, min=0.05,mode=0.1,max=0.3)

#Rate at which pregnant heifers in A and I compartments aborts in late gestation (day-1)
a1<-runif(n, min=0.02,max=0.15)/D_I 

#Rate at which pregnant heifers in C compartment aborts in late gestation (day-1)
a2<-runif(n, min=0.02,max=0.15)/D_C

#Amount of S. Dublin shed in feces by individuals in I (CFU/g)
z<- rpert(n, min=10^2.08, mode=10^2.81, max = 10^4.6) #rpert(n, min=10^3, mode=10^4, max = 10^5.7) 

#Reduction in shedding by individuals in A and C compared to those in I (i.e., carriers and asymptomatic shed 100 times less S. Dublin compared to clinically ill individuals)
s<- 100

#Proportion of feces covered by newly added bedding material (day). Only for deep bedding
covered<- 0 #If deep bedding considered, change to 0.474 (47% feces are covered on average by newly added bedding using the Appendix S1 sheet and a 95% of feces covered per bedding addition)

#Rate of feces being covered by newly added bedding material (day-1). Only for deep bedding.
b<- -log(1-covered)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 1
lambda1<- (f1*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 1
epsilon1 <- (f1*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 2
lambda2<- (f2*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 2
epsilon2 <- (f2*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 3
lambda3<- (f3*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 3
epsilon3 <- (f3*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 4
lambda4<- (f4*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 4
epsilon4 <- (f4*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 5
lambda5<- (f5*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 5
epsilon5 <- (f5*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 5
lambda6<- (f6*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 5
epsilon6 <- (f6*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 5
lambda7<- (f7*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 5
epsilon7 <- (f7*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 5
lambda8<- (f8*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 8
epsilon8 <- (f8*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 9
lambda9<- (f9*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 9
epsilon9 <- (f9*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 10
lambda10<- (f10*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 10
epsilon10 <- (f10*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 11
lambda11<- (f11*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 11
epsilon11 <- (f11*z)

#Pathogen shedding rate by individuals in A and C (CFU-1 day-1 individual) - Raising stage 12
lambda12<- (f12*z)/s

#Pathogen shedding rate by individuals in I (CFU-1 day-1 individual) - Raising stage 12
epsilon12 <- (f12*z)

#The parameter "phase" was used to indicate the season in which the simulation starts (the function for temperature is found below (tempC)
phase<- 6.2 #Summer
#phase<- 33 #Fall
#phase<- 22 #Winter
#phase<- 17.4 #Spring

####ROOT FUNCTION AND EVENTS####
#A root function was used to allow individuals to move from one raising stage to the next. y[71] represents the "compartment" "Event_trigger" (the 71st "compartment" in the model).
#When the "Event_trigger" "compartment" reaches a multiple of 73, the root function is activated and animals move from one raising stage to the next. For instance, when the day 219 is 
#reached in the model, then individuals move from one raising stage to the next.
# Root function
rootfunc <- function(t,y,parms) c(
  y[100]-45, y[100]-(45*2), y[100]-(45*3), y[100]-(45*4), y[100]-(45*5), y[100]-(45*6), y[100]-(45*7), y[100]-(45*8),
  y[100]-(45*9), y[100]-(45*10), y[100]-(45*11), y[100]-(45*12), y[100]-(45*13), y[100]-(45*14), y[100]-(45*15), y[100]-(45*16),
  y[100]-(45*17), y[100]-(45*18), y[100]-(45*19), y[100]-(45*20), y[100]-(45*21), y[100]-(45*22), y[100]-(45*23), y[100]-(45*24),
  y[100]-(45*25)

)

#Events
eventfunc <- function(t,y,parms) { #Events Event_triggered at multiplnes of 73
  
  if(floor(y[100])-45*1==0 | floor(y[100])-45*2==0 | floor(y[100])-45*3==0 | floor(y[100])-45*4==0 | floor(y[100])-45*5==0 |
     floor(y[100])-45*6==0 | floor(y[100])-45*7==0 | floor(y[100])-45*8==0 | floor(y[100])-45*9==0 | floor(y[100])-45*10==0 |
      floor(y[100])-45*11==0 | floor(y[100])-45*12==0 | floor(y[100])-45*13==0 | floor(y[100])-45*14==0 | floor(y[100])-45*15==0 |
        floor(y[100])-45*16==0 | floor(y[100])-45*17==0 | floor(y[100])-45*18==0 | floor(y[100])-45*19==0 | floor(y[100])-45*20==0 |
       floor(y[100])-45*21==0 | floor(y[100])-45*22==0 | floor(y[100])-45*23==0 | floor(y[100])-45*24==0 | floor(y[100])-45*25==0
  )

      {
    
    y[84]<-y[84] + y[56] + y[57] + y[59] + y[60]          #Individuals in the Raising stage 5 return to their farms as pregnant heifers
    y[85]<-y[85] + y[58]                                  #Severely diseased sacrificed
    
    y[88]<- y[88] + y[1] + y[2] + y[3] + y[4] + y[5]      #Records those who left the Raising stage 1
    y[89]<- y[89] + y[6] + y[7] + y[8] + y[9] + y[10]     #Records those who left the Raising stage 2
    y[90]<- y[90] + y[11] + y[12] + y[13] + y[14] + y[15] #Records those who left the Raising stage 3
    y[91]<- y[91] + y[16] + y[17] + y[18] + y[19] + y[20] #Records those who left the Raising stage 4
    y[92]<- y[92] + y[21] + y[22] + y[23] + y[24] + y[25] #Records those who left the Raising stage 5
    y[93]<- y[93] + y[26] + y[27] + y[28] + y[29] + y[30] #Records those who left the Raising stage 6
    y[94]<- y[94] + y[31] + y[32] + y[33] + y[34] + y[35] #Records those who left the Raising stage 7
    y[95]<- y[95] + y[36] + y[37] + y[38] + y[39] + y[40] #Records those who left the Raising stage 8
    y[96]<- y[96] + y[41] + y[42] + y[43] + y[44] + y[45] #Records those who left the Raising stage 9
    y[97]<- y[97] + y[46] + y[47] + y[48] + y[49] + y[50] #Records those who left the Raising stage 10
    y[98]<- y[98] + y[51] + y[52] + y[53] + y[54] + y[55] #Records those who left the Raising stage 11
    y[99]<- y[99] + y[56] + y[57] + y[58] + y[59] + y[60] #Records those who left the Raising stage 12
    
    y[82]<- y[82] + y[59]                                 #Records pregnant cattle that return to their farms as carriers
    y[83]<- y[83] + y[57]                                 #Records pregnant cattle that return to their farms as asymptomatic
    
    y[56] <- y[51]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[57] <- y[52]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[58] <- y[53]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[59] <- y[54]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[60] <- y[55]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[51] <- y[46]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[52] <- y[47]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[53] <- y[48]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[54] <- y[49]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[55] <- y[50]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[46] <- y[41]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[47] <- y[42]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[48] <- y[43]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[49] <- y[44]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[50] <- y[45]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[41] <- y[36]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[42] <- y[37]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[43] <- y[38]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[44] <- y[39]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[45] <- y[40]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[36] <- y[31]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[37] <- y[32]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[38] <- y[33]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[39] <- y[34]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[40] <- y[35]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[31] <- y[26]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[32] <- y[27]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[33] <- y[28]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[34] <- y[29]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[35] <- y[30]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  

    y[26] <- y[21]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[27] <- y[22]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[28] <- y[23]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[29] <- y[24]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[30] <- y[25]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[21] <- y[16]                                        #Susceptible (S) heifers move from Raising stage 4 to Raising stage 5 
    y[22] <- y[17]                                        #Asymptomatic (A) heifers move from Raising stage 4 to Raising stage 5 
    y[23] <- y[18]                                        #Clinically ill (I) heifers move from Raising stage 4 to Raising stage 5  
    y[24] <- y[19]                                        #Carriers (C) heifers move from Raising stage 4 to Raising stage 5  
    y[25] <- y[20]                                        #Recovered (R) heifers move from Raising stage 4 to Raising stage 5  
    
    y[16] <- y[11]                                        #Susceptible (S) heifers move from Raising stage 3 to Raising stage 4 
    y[17] <- y[12]                                        #Asymptomatic (A) heifers move from Raising stage 3 to Raising stage 4 
    y[18] <- y[13]                                        #Clinically ill (I) heifers move from Raising stage 3 to Raising stage 4 
    y[19] <- y[14]                                        #Carriers (C) heifers move from Raising stage 3 to Raising stage 4 
    y[20] <- y[15]                                        #Recovered (R) heifers move from Raising stage 3 to Raising stage 4 
    
    y[11] <- y[6]                                         #Susceptible (S) heifers move from Raising stage 2 to Raising stage 3 
    y[12] <- y[7]                                         #Asymptomatic (A) heifers move from Raising stage 2 to Raising stage 3 
    y[13] <- y[8]                                         #Clinically ill (I) heifers move from Raising stage 2 to Raising stage 3 
    y[14] <- y[9]                                         #Carriers (C) heifers move from Raising stage 2 to Raising stage 3 
    y[15] <- y[10]                                        #Recovered (R) heifers move from Raising stage 2 to Raising stage 3 
    
    y[6] <- y[1]                                          #Susceptible (S) heifers move from Raising stage 1 to Raising stage 2  
    y[7] <- y[2]                                          #Asymptomatic (A) heifers move from Raising stage 1 to Raising stage 2 
    y[8] <- y[3]                                          #Clinically ill (I) heifers move from Raising stage 1 to Raising stage 2 
    y[9] <- y[4]                                          #Carriers (C) heifers move from Raising stage 1 to Raising stage 2 
    y[10] <- y[5]                                         #Recovered (R) heifers move from Raising stage 1 to Raising stage 2
    
    y[1] <- (No/12)                                         #Individuals in Raising stage 1 left. A new batch of Susceptible weaned calves arrive at the operation
    y[2] <- 0                                             #Individuals in Raising stage 1 left and now the compartment is empty.
    y[3] <- 0                                             #Individuals in Raising stage 1 left and now the compartment is empty.
    y[4] <- 0                                             #Individuals in Raising stage 1 left and now the compartment is empty.
    y[5] <- 0                                             #Individuals in Raising stage 1 left and now the compartment is empty. 
  }
  
  
  return(y)
}

####MODEL EQUATIONS#####
SALL_dyn <-function(t,var,par) { 
  
  ## State values ##
  #Weaned heifers
  S_weaned_1=var[1]; #Susceptible 
  A_weaned_1=var[2]; #Asymptomatic
  I_weaned_1=var[3]; #Clinically ill
  C_weaned_1=var[4]; #Carrier
  R_weaned_1=var[5]; #Recovered
  #Growing heifers
  S_growing_2=var[6]; 
  A_growing_2=var[7];
  I_growing_2=var[8];
  C_growing_2=var[9];
  R_growing_2=var[10];
  
  S_growing_3=var[11];
  A_growing_3=var[12];
  I_growing_3=var[13];
  C_growing_3=var[14];
  R_growing_3=var[15];
  
  S_growing_4=var[16];
  A_growing_4=var[17];
  I_growing_4=var[18];
  C_growing_4=var[19];
  R_growing_4=var[20];
  
  S_growing_5=var[21];
  A_growing_5=var[22];
  I_growing_5=var[23];
  C_growing_5=var[24];
  R_growing_5=var[25];
  
  S_growing_6=var[26];
  A_growing_6=var[27];
  I_growing_6=var[28];
  C_growing_6=var[29];
  R_growing_6=var[30];
  
  S_growing_7=var[31];
  A_growing_7=var[32];
  I_growing_7=var[33];
  C_growing_7=var[34];
  R_growing_7=var[35];
  
  S_growing_8=var[36];
  A_growing_8=var[37];
  I_growing_8=var[38];
  C_growing_8=var[39];
  R_growing_8=var[40];
  
  S_pregnant_9=var[41];
  A_pregnant_9=var[42];
  I_pregnant_9=var[43];
  C_pregnant_9=var[44];
  R_pregnant_9=var[45];
  
  S_pregnant_10=var[46];
  A_pregnant_10=var[47];
  I_pregnant_10=var[48];
  C_pregnant_10=var[49];
  R_pregnant_10=var[50];
  
  S_pregnant_11=var[51];
  A_pregnant_11=var[52];
  I_pregnant_11=var[53];
  C_pregnant_11=var[54];
  R_pregnant_11=var[55];
  
  S_pregnant_12=var[56];
  A_pregnant_12=var[57];
  I_pregnant_12=var[58];
  C_pregnant_12=var[59];
  R_pregnant_12=var[60];
  
  E_weaned_1=var[61];
  
  E_growing_2=var[62];
  E_growing_3=var[63];
  E_growing_4=var[64];
  E_growing_5=var[65];
  E_growing_6=var[66];
  E_growing_7=var[67];
  E_growing_8=var[68];
  
  E_pregnant_9=var[69];
  E_pregnant_10=var[70];
  E_pregnant_11=var[71];
  E_pregnant_12=var[72];
  
  Abortions=var[73];

  D_weaned_1=var[74];
  D_growing_2=var[75];
  D_growing_3=var[76];
  D_growing_4=var[77];
  D_growing_5=var[78];
  D_growing_6=var[79];
  D_growing_7=var[80];
  D_growing_8=var[81];
  
  Carriers_leaving[82];
  Asymptomatic_leaving[83];
  
  Pregnant_sold[84];
  Sacrificed[85];
  Treatment_weaned[86];
  Treatment_growing_and_pregnant[87];
  
  Completed_raising.stage_1=var[88];
  Completed_raising.stage_2=var[89];
  Completed_raising.stage_3=var[90];
  Completed_raising.stage_4=var[91];
  Completed_raising.stage_5=var[92];
  Completed_raising.stage_6=var[93];
  Completed_raising.stage_7=var[94];
  Completed_raising.stage_8=var[95];
  Completed_raising.stage_9=var[96];
  Completed_raising.stage_10=var[97];
  Completed_raising.stage_11=var[98];
  Completed_raising.stage_12=var[99];
  
  Event_trigger=var[100]
  
  #Parameters list. Parameters included here are only those that are defined in a distribution. They need to respect the order in which 
  #are presented above and need to be included in this same order in the equations below.
  omega=par[1]; 
  gamma=par[2]; 
  m=par[3];
  u1=par[4];
  w=par[5];
  d1=par[6];
  v=par[7];
  eta=par[8]; 
  u2=par[9];
  d2=par[10];
  u3=par[11];
  a1=par[12];
  a2=par[13];
  lambda1=par[14];  
  epsilon1=par[15]; 
  lambda2=par[16];
  epsilon2=par[17];
  lambda3=par[18]; 
  epsilon3=par[19];
  lambda4=par[20]; 
  epsilon4=par[21];
  lambda5=par[22];
  epsilon5=par[23];
  lambda6=par[24];
  epsilon6=par[25]; 
  lambda7=par[26];
  epsilon7=par[27]; 
  lambda8=par[28];
  epsilon8=par[29]; 
  lambda9=par[30];
  epsilon9=par[31]; 
  lambda10=par[32];
  epsilon10=par[33]; 
  lambda11=par[34];
  epsilon11=par[35]; 
  lambda12=par[36];
  epsilon12=par[37]; 
  
  #Function for temperature oscillation through the year
  tempC= 15*cos((2*pi)*(t)/365+phase)+10 #Temperature goes from about -5 to 25 degrees C. The season in which the simulation
  # will start depends on the value of the "phase" variable.
  # tempC= 10*cos((2*pi)*(t)/365+17.4)+7.7 #Denmark. The model starts in spring (Calibration)
  # tempC= 5.4*cos((2*pi)*(t)/365+9.3)+22.6 #Florida. The model starts in January. Temperatures ranges from 18 to 28 degrees C (Validation)
  
  #Function representing temperature-dependent death rate of Salmonella Dublin in the environment
  phi= ifelse(tempC>0, -1*(-0.009*tempC - 0.455), 1.25) 
  
  ### Derivatives
  #Calves
  dS_weaned_1= - beta*S_weaned_1*E_weaned_1 - beta*kappa*S_weaned_1*(E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_weaned_1 + omega*R_weaned_1 
  dA_weaned_1= - (gamma+ m)*A_weaned_1 - rho*A_weaned_1 + beta*(1-u1)*S_weaned_1*E_weaned_1 + beta*(1-u1)*kappa*S_weaned_1*(E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_weaned_1= - (gamma+ w+d1*(1-v))*I_weaned_1 - rho*I_weaned_1 + beta*u1*S_weaned_1*E_weaned_1 + beta*u1*kappa*S_weaned_1*(E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_weaned_1= - eta*C_weaned_1 - rho*C_weaned_1 +  m*A_weaned_1 +  w*I_weaned_1
  dR_weaned_1= - omega*R_weaned_1 - rho*R_weaned_1 + (gamma)*A_weaned_1 + (gamma)*I_weaned_1 + eta*C_weaned_1 
  
  #Growing heifers
  dS_growing_2= - beta*S_growing_2*E_growing_2 - beta*kappa*S_growing_2*(E_weaned_1+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_2 + omega*R_growing_2 
  dA_growing_2= - (gamma+ m)*A_growing_2 - rho*A_growing_2 + beta*(1-u2)*S_growing_2*E_growing_2 + beta*(1-u2)*kappa*S_growing_2*(E_weaned_1+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_2= - (gamma+ w+d2*(1-v))*I_growing_2 - rho*I_growing_2 + beta*u2*S_growing_2*E_growing_2  + beta*u2*kappa*S_growing_2*(E_weaned_1+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_2= - eta*C_growing_2 - rho*C_growing_2 +  m*A_growing_2 +  w*I_growing_2
  dR_growing_2= - omega*R_growing_2 - rho*R_growing_2 + (gamma)*A_growing_2 + (gamma)*I_growing_2 + eta*C_growing_2 
  
  dS_growing_3= - beta*S_growing_3*E_growing_3 - beta*kappa*S_growing_3*(E_weaned_1+E_growing_2+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_3 + omega*R_growing_3 
  dA_growing_3= - (gamma+ m)*A_growing_3 - rho*A_growing_3 + beta*(1-u2)*S_growing_3*E_growing_3 + beta*(1-u2)*kappa*S_growing_3*(E_weaned_1+E_growing_2+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_3= - (gamma+ w+d2*(1-v))*I_growing_3 - rho*I_growing_3 + beta*u2*S_growing_3*E_growing_3 + beta*u2*kappa*S_growing_3*(E_weaned_1+E_growing_2+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_3= - eta*C_growing_3 - rho*C_growing_3 +  m*A_growing_3 +  w*I_growing_3
  dR_growing_3= - omega*R_growing_3 - rho*R_growing_3 + (gamma)*A_growing_3 + (gamma)*I_growing_3 + eta*C_growing_3 
  
  dS_growing_4= - beta*S_growing_4*E_growing_4 - beta*kappa*S_growing_4*(E_weaned_1+E_growing_2+E_growing_3+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_4 + omega*R_growing_4 
  dA_growing_4= - (gamma+ m)*A_growing_4 - rho*A_growing_4 + beta*(1-u2)*S_growing_4*E_growing_4 + beta*(1-u2)*kappa*S_growing_4*(E_weaned_1+E_growing_2+E_growing_3+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_4= - (gamma+ w+d2*(1-v))*I_growing_4 - rho*I_growing_4 + beta*u2*S_growing_4*E_growing_4 + beta*u2*kappa*S_growing_4*(E_weaned_1+E_growing_2+E_growing_3+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_4= - eta*C_growing_4 - rho*C_growing_4 +  m*A_growing_4 +  w*I_growing_4
  dR_growing_4= - omega*R_growing_4 - rho*R_growing_4 + (gamma)*A_growing_4 + (gamma)*I_growing_4 + eta*C_growing_4 
  
  dS_growing_5= - beta*S_growing_5*E_growing_5 - beta*kappa*S_growing_5*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_5 + omega*R_growing_5 
  dA_growing_5= - (gamma+ m)*A_growing_5 - rho*A_growing_5 + beta*(1-u2)*S_growing_5*E_growing_5 + beta*(1-u2)*kappa*S_growing_5*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_5= - (gamma+ w+d2*(1-v))*I_growing_5 - rho*I_growing_5 + beta*u2*S_growing_5*E_growing_5 + beta*u2*kappa*S_growing_5*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_5= - eta*C_growing_5 - rho*C_growing_5 +  m*A_growing_5 +  w*I_growing_5
  dR_growing_5= - omega*R_growing_5 - rho*R_growing_5 + (gamma)*A_growing_5 + (gamma)*I_growing_5 + eta*C_growing_5 
  
  dS_growing_6= - beta*S_growing_6*E_growing_6 - beta*kappa*S_growing_6*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_6 + omega*R_growing_6 
  dA_growing_6= - (gamma+ m)*A_growing_6 - rho*A_growing_6 + beta*(1-u2)*S_growing_6*E_growing_6 + beta*(1-u2)*kappa*S_growing_6*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_6= - (gamma+ w+d2*(1-v))*I_growing_6 - rho*I_growing_6 + beta*u2*S_growing_6*E_growing_6 + beta*u2*kappa*S_growing_6*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_6= - eta*C_growing_6 - rho*C_growing_6 +  m*A_growing_6 +  w*I_growing_6
  dR_growing_6= - omega*R_growing_6 - rho*R_growing_6 + (gamma)*A_growing_6 + (gamma)*I_growing_6 + eta*C_growing_6 
  
  dS_growing_7= - beta*S_growing_7*E_growing_7 - beta*kappa*S_growing_7*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_7 + omega*R_growing_7 
  dA_growing_7= - (gamma+ m)*A_growing_7 - rho*A_growing_7 + beta*(1-u2)*S_growing_7*E_growing_7 + beta*(1-u2)*kappa*S_growing_7*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_7= - (gamma+ w+d2*(1-v))*I_growing_7 - rho*I_growing_7 + beta*u2*S_growing_7*E_growing_7 + beta*u2*kappa*S_growing_7*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_7= - eta*C_growing_7 - rho*C_growing_7 +  m*A_growing_7 +  w*I_growing_7
  dR_growing_7= - omega*R_growing_7 - rho*R_growing_7 + (gamma)*A_growing_7 + (gamma)*I_growing_7 + eta*C_growing_7 
  
  dS_growing_8= - beta*S_growing_8*E_growing_8 - beta*kappa*S_growing_8*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_growing_8 + omega*R_growing_8 
  dA_growing_8= - (gamma+ m)*A_growing_8 - rho*A_growing_8 + beta*(1-u2)*S_growing_8*E_growing_8 + beta*(1-u2)*kappa*S_growing_8*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_growing_8= - (gamma+ w+d2*(1-v))*I_growing_8 - rho*I_growing_8 + beta*u2*S_growing_8*E_growing_8 + beta*u2*kappa*S_growing_8*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_pregnant_9+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_growing_8= - eta*C_growing_8 - rho*C_growing_8 +  m*A_growing_8 +  w*I_growing_8
  dR_growing_8= - omega*R_growing_8 - rho*R_growing_8 + (gamma)*A_growing_8 + (gamma)*I_growing_8 + eta*C_growing_8 
  
  #Pregnant heifers
  dS_pregnant_9= - beta*S_pregnant_9*E_pregnant_9 - beta*kappa*S_pregnant_9*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_10+E_pregnant_11+E_pregnant_12) - rho*S_pregnant_9 + omega*R_pregnant_9 
  dA_pregnant_9= - (gamma+ m)*A_pregnant_9 - rho*A_pregnant_9 + beta*(1-u3)*S_pregnant_9*E_pregnant_9 + beta*(1-u3)*kappa*S_pregnant_9*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dI_pregnant_9= - (gamma+ w)*I_pregnant_9 - rho*I_pregnant_9 + beta*(u3)*S_pregnant_9*E_pregnant_9 + beta*u3*kappa*S_pregnant_9*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_10+E_pregnant_11+E_pregnant_12)
  dC_pregnant_9= - eta*C_pregnant_9 - rho*C_pregnant_9 +  m*A_pregnant_9 +  w*I_pregnant_9
  dR_pregnant_9= - omega*R_pregnant_9 - rho*R_pregnant_9 + (gamma)*A_pregnant_9 + (gamma)*I_pregnant_9 + eta*C_pregnant_9 
  
  dS_pregnant_10= - beta*S_pregnant_10*E_pregnant_10 - beta*kappa*S_pregnant_10*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_11+E_pregnant_12) - rho*S_pregnant_10 + omega*R_pregnant_10 
  dA_pregnant_10= - (gamma+ m)*A_pregnant_10 - rho*A_pregnant_10 + beta*(1-u3)*S_pregnant_10*E_pregnant_10 + beta*(1-u3)*kappa*S_pregnant_10*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_11+E_pregnant_12)
  dI_pregnant_10= - (gamma+ w)*I_pregnant_10 - rho*I_pregnant_10 + beta*(u3)*S_pregnant_10*E_pregnant_10 + beta*u3*kappa*S_pregnant_10*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_11+E_pregnant_12)
  dC_pregnant_10= - eta*C_pregnant_10 - rho*C_pregnant_10 +  m*A_pregnant_10 +  w*I_pregnant_10
  dR_pregnant_10= - omega*R_pregnant_10 - rho*R_pregnant_10 + (gamma)*A_pregnant_10 + (gamma)*I_pregnant_10 + eta*C_pregnant_10 
  
  dS_pregnant_11= - beta*S_pregnant_11*E_pregnant_11 - beta*kappa*S_pregnant_11*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_12) - rho*S_pregnant_11 + omega*R_pregnant_11 
  dA_pregnant_11= - (gamma+ m)*A_pregnant_11 - rho*A_pregnant_11 + beta*(1-u3)*S_pregnant_11*E_pregnant_11 + beta*(1-u3)*kappa*S_pregnant_11*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_12)
  dI_pregnant_11= - (gamma+ w)*I_pregnant_11 - rho*I_pregnant_11 + beta*(u3)*S_pregnant_11*E_pregnant_11 + beta*u3*kappa*S_pregnant_11*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_12)
  dC_pregnant_11= - eta*C_pregnant_11 - rho*C_pregnant_11 +  m*A_pregnant_11 +  w*I_pregnant_11
  dR_pregnant_11= - omega*R_pregnant_11 - rho*R_pregnant_11 + (gamma)*A_pregnant_11 + (gamma)*I_pregnant_11 + eta*C_pregnant_11 
  
  dS_pregnant_12= - beta*S_pregnant_12*E_pregnant_12 - beta*kappa*S_pregnant_12*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11) - rho*S_pregnant_12 + omega*R_pregnant_12 
  dA_pregnant_12= - (gamma+ m)*A_pregnant_12 - rho*A_pregnant_12 + beta*(1-u3)*S_pregnant_12*E_pregnant_12 + beta*(1-u3)*kappa*S_pregnant_12*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11)
  dI_pregnant_12= - (gamma+ w)*I_pregnant_12 - rho*I_pregnant_12 + beta*(u3)*S_pregnant_12*E_pregnant_12 + beta*u3*kappa*S_pregnant_12*(E_weaned_1+E_growing_2+E_growing_3+E_growing_4+E_growing_5+E_growing_6+E_growing_7+E_growing_8+E_pregnant_9+E_pregnant_10+E_pregnant_11)
  dC_pregnant_12= - eta*C_pregnant_12 - rho*C_pregnant_12 +  m*A_pregnant_12 +  w*I_pregnant_12
  dR_pregnant_12= - omega*R_pregnant_12 - rho*R_pregnant_12 + (gamma)*A_pregnant_12 + (gamma)*I_pregnant_12 + eta*C_pregnant_12 
  
  dE_weaned_1= - phi*E_weaned_1 - mu*E_weaned_1 - b*E_weaned_1 +  A_weaned_1*lambda1 +  I_weaned_1*epsilon1  +  C_weaned_1*lambda1 
  
  dE_growing_2 = - phi*E_growing_2 - mu*E_growing_2 - b*E_growing_2 + A_growing_2*lambda2 +  I_growing_2*epsilon2 +  C_growing_2*lambda2
  dE_growing_3 =  - phi*E_growing_3 - mu*E_growing_3 - b*E_growing_3 + A_growing_3*lambda3 +  I_growing_3*epsilon3 +  C_growing_3*lambda3 
  dE_growing_4 =  - phi*E_growing_4 - mu*E_growing_4 - b*E_growing_4 + A_growing_4*lambda4 +  I_growing_4*epsilon4 +  C_growing_4*lambda4 
  dE_growing_5 =  - phi*E_growing_5 - mu*E_growing_5 - b*E_growing_5 + A_growing_5*lambda5 +  I_growing_5*epsilon5 +  C_growing_5*lambda5 
  dE_growing_6 =  - phi*E_growing_6 - mu*E_growing_6 - b*E_growing_6 + A_growing_6*lambda6 +  I_growing_6*epsilon6 +  C_growing_6*lambda6 
  dE_growing_7 =  - phi*E_growing_7 - mu*E_growing_7 - b*E_growing_7 + A_growing_7*lambda7 +  I_growing_7*epsilon7 +  C_growing_7*lambda7 
  dE_growing_8 =  - phi*E_growing_8 - mu*E_growing_8 - b*E_growing_8 + A_growing_8*lambda8 +  I_growing_8*epsilon8 +  C_growing_8*lambda8 
    
  dE_pregnant_9 = - phi*E_pregnant_9 - mu*E_pregnant_9 - b*E_pregnant_9 + A_pregnant_9*lambda9 +  I_pregnant_9*epsilon9 +  C_pregnant_9*lambda9 
  dE_pregnant_10 = - phi*E_pregnant_10 - mu*E_pregnant_10 - b*E_pregnant_10 + A_pregnant_10*lambda10 +  I_pregnant_10*epsilon10 +  C_pregnant_10*lambda10 
  dE_pregnant_11 = - phi*E_pregnant_11 - mu*E_pregnant_11 - b*E_pregnant_11 + A_pregnant_11*lambda11 +  I_pregnant_11*epsilon11 +  C_pregnant_11*lambda11 
  dE_pregnant_12 = - phi*E_pregnant_12 - mu*E_pregnant_12 - b*E_pregnant_12 + A_pregnant_12*lambda12 +  I_pregnant_12*epsilon12 +  C_pregnant_12*lambda12 
  
  dAbortions= A_pregnant_12*a1 + I_pregnant_12*a1 + C_pregnant_12*a2 
  
  dD_weaned_1= d1*(1-v)*I_weaned_1 
  dD_growing_2= d2*(1-v)*I_growing_2 
  dD_growing_3= d2*(1-v)*I_growing_3 
  dD_growing_4= d2*(1-v)*I_growing_4 
  dD_growing_5= d2*(1-v)*I_growing_5 
  dD_growing_6= d2*(1-v)*I_growing_6 
  dD_growing_7= d2*(1-v)*I_growing_7 
  dD_growing_8= d2*(1-v)*I_growing_8 
  
  dAsymptomatic_leaving=0
  dCarriers_leaving=0
  
  dPregnant_sold= 0
  dSacrificed= 0
  
  dTreatment_weaned = beta*u1*S_weaned_1*E_weaned_1 
  dTreatment_growing_and_pregnant = beta*u2*S_growing_2*E_growing_2 + beta*u2*S_growing_3*E_growing_3 + beta*u2*S_growing_4*E_growing_4 + beta*u2*S_growing_5*E_growing_5 + beta*u2*S_growing_6*E_growing_6 + beta*u2*S_growing_7*E_growing_7 + beta*u2*S_growing_8*E_growing_8 +
    beta*u3*S_pregnant_9*E_pregnant_9+beta*u3*S_pregnant_10*E_pregnant_10+beta*u3*S_pregnant_11*E_pregnant_11+beta*u3*S_pregnant_12*E_pregnant_12
  
  dCompleted_raising.stage_1=0
  dCompleted_raising.stage_2=0
  dCompleted_raising.stage_3=0
  dCompleted_raising.stage_4=0
  dCompleted_raising.stage_5=0
  dCompleted_raising.stage_6=0
  dCompleted_raising.stage_7=0
  dCompleted_raising.stage_8=0
  dCompleted_raising.stage_9=0
  dCompleted_raising.stage_10=0
  dCompleted_raising.stage_11=0
  dCompleted_raising.stage_12=0
  
  
  dEvent_trigger= 1
  
  v = c(dS_weaned_1, dA_weaned_1, dI_weaned_1, dC_weaned_1, dR_weaned_1, 
        dS_growing_2, dA_growing_2, dI_growing_2, dC_growing_2, dR_growing_2, 
        dS_growing_3, dA_growing_3, dI_growing_3, dC_growing_3, dR_growing_3,
        dS_growing_4, dA_growing_4, dI_growing_4, dC_growing_4, dR_growing_4,
        dS_growing_5, dA_growing_5, dI_growing_5, dC_growing_5, dR_growing_5,
        dS_growing_6, dA_growing_6, dI_growing_6, dC_growing_6, dR_growing_6,
        dS_growing_7, dA_growing_7, dI_growing_7, dC_growing_7, dR_growing_7,
        dS_growing_8, dA_growing_8, dI_growing_8, dC_growing_8, dR_growing_8,
        dS_pregnant_9, dA_pregnant_9, dI_pregnant_9, dC_pregnant_9, dR_pregnant_9, 
        dS_pregnant_10, dA_pregnant_10, dI_pregnant_10, dC_pregnant_10, dR_pregnant_10, 
        dS_pregnant_11, dA_pregnant_11, dI_pregnant_11, dC_pregnant_11, dR_pregnant_11, 
        dS_pregnant_12, dA_pregnant_12, dI_pregnant_12, dC_pregnant_12, dR_pregnant_12, 
        
        dE_weaned_1, 
        
        dE_growing_2,
        dE_growing_3,
        dE_growing_4,
        dE_growing_5,
        dE_growing_6,
        dE_growing_7,
        dE_growing_8,
        
        dE_pregnant_9,
        dE_pregnant_10,
        dE_pregnant_11,
        dE_pregnant_12,
        
        dAbortions,
        
        dD_weaned_1, dD_growing_2, dD_growing_3,dD_growing_4,dD_growing_5,dD_growing_6,dD_growing_7,dD_growing_8,
        dCarriers_leaving, dAsymptomatic_leaving,
        dPregnant_sold, dSacrificed, 
        dTreatment_weaned, dTreatment_growing_and_pregnant,
        
        dCompleted_raising.stage_1,
        dCompleted_raising.stage_2,
        dCompleted_raising.stage_3,
        dCompleted_raising.stage_4,
        dCompleted_raising.stage_5,
        dCompleted_raising.stage_6,
        dCompleted_raising.stage_7,
        dCompleted_raising.stage_8,
        dCompleted_raising.stage_9,
        dCompleted_raising.stage_10,
        dCompleted_raising.stage_11,
        dCompleted_raising.stage_12,
        
        dEvent_trigger
  )
  names(v) = c('dS_weaned_1', 'dA_weaned_1', 'dI_weaned_1', 'dC_weaned_1', 'dR_weaned_1', 
               'dS_growing_2', 'dA_growing_2', 'dI_growing_2', 'dC_growing_2', 'dR_growing_2', 
               'dS_growing_3', 'dA_growing_3', 'dI_growing_3', 'dC_growing_3', 'dR_growing_3',
               'dS_growing_4', 'dA_growing_4', 'dI_growing_4', 'dC_growing_4', 'dR_growing_4',
               'dS_growing_5', 'dA_growing_5', 'dI_growing_5', 'dC_growing_5', 'dR_growing_5',
               'dS_growing_6', 'dA_growing_6', 'dI_growing_6', 'dC_growing_6', 'dR_growing_6',
               'dS_growing_7', 'dA_growing_7', 'dI_growing_7', 'dC_growing_7', 'dR_growing_7',
               'dS_growing_8', 'dA_growing_8', 'dI_growing_8', 'dC_growing_8', 'dR_growing_8',
               'dS_pregnant_9', 'dA_pregnant_9', 'dI_pregnant_9', 'dC_pregnant_9', 'dR_pregnant_9', 
               'dS_pregnant_10', 'dA_pregnant_10', 'dI_pregnant_10', 'dC_pregnant_10', 'dR_pregnant_10', 
               'dS_pregnant_11', 'dA_pregnant_11', 'dI_pregnant_11', 'dC_pregnant_11', 'dR_pregnant_11', 
               'dS_pregnant_12', 'dA_pregnant_12', 'dI_pregnant_12', 'dC_pregnant_12', 'dR_pregnant_12', 
               
               'dE_weaned_1', 
               
               'dE_growing_2',
               'dE_growing_3',
               'dE_growing_4',
               'dE_growing_5',
               'dE_growing_6',
               'dE_growing_7',
               'dE_growing_8',
               
               'dE_pregnant_9',
               'dE_pregnant_10',
               'dE_pregnant_11',
               'dE_pregnant_12',
               
               'dAbortions',
               
               'dD_weaned_1', 'dD_growing_2', 'dD_growing_3','dD_growing_4','dD_growing_5','dD_growing_6','dD_growing_7','dD_growing_8',
               'dCarriers_leaving', 'dAsymptomatic_leaving',
               'dPregnant_sold', 'dSacrificed', 
               'dTreatment_weaned', 'dTreatment_growing_and_pregnant',
               
               'dCompleted_raising.stage_1',
               'dCompleted_raising.stage_2',
               'dCompleted_raising.stage_3',
               'dCompleted_raising.stage_4',
               'dCompleted_raising.stage_5',
               'dCompleted_raising.stage_6',
               'dCompleted_raising.stage_7',
               'dCompleted_raising.stage_8',
               'dCompleted_raising.stage_9',
               'dCompleted_raising.stage_10',
               'dCompleted_raising.stage_11',
               'dCompleted_raising.stage_12',
               
               'dEvent_trigger'
  )
  return(list(v))
} 
Sys.time()

#### Run the model
SALL_sol=NULL

S_weaned_1=NULL
A_weaned_1=NULL
I_weaned_1=NULL
C_weaned_1=NULL
R_weaned_1=NULL

S_growing_2=NULL
A_growing_2=NULL
I_growing_2=NULL
C_growing_2=NULL
R_growing_2=NULL

S_growing_3=NULL
A_growing_3=NULL
I_growing_3=NULL
C_growing_3=NULL
R_growing_3=NULL

S_growing_4=NULL
A_growing_4=NULL
I_growing_4=NULL
C_growing_4=NULL
R_growing_4=NULL

S_growing_5=NULL
A_growing_5=NULL
I_growing_5=NULL
C_growing_5=NULL
R_growing_5=NULL

S_growing_6=NULL
A_growing_6=NULL
I_growing_6=NULL
C_growing_6=NULL
R_growing_6=NULL

S_growing_7=NULL
A_growing_7=NULL
I_growing_7=NULL
C_growing_7=NULL
R_growing_7=NULL

S_growing_8=NULL
A_growing_8=NULL
I_growing_8=NULL
C_growing_8=NULL
R_growing_8=NULL

S_pregnant_9=NULL
A_pregnant_9=NULL
I_pregnant_9=NULL
C_pregnant_9=NULL
R_pregnant_9=NULL

S_pregnant_10=NULL
A_pregnant_10=NULL
I_pregnant_10=NULL
C_pregnant_10=NULL
R_pregnant_10=NULL

S_pregnant_11=NULL
A_pregnant_11=NULL
I_pregnant_11=NULL
C_pregnant_11=NULL
R_pregnant_11=NULL

S_pregnant_12=NULL
A_pregnant_12=NULL
I_pregnant_12=NULL
C_pregnant_12=NULL
R_pregnant_12=NULL

E_weaned_1=NULL

E_growing_2=NULL
E_growing_3=NULL
E_growing_4=NULL
E_growing_5=NULL
E_growing_6=NULL
E_growing_7=NULL
E_growing_8=NULL

E_pregnant_9=NULL
E_pregnant_10=NULL
E_pregnant_11=NULL
E_pregnant_12=NULL

Abortions=NULL

D_weaned_1=NULL
D_growing_2=NULL
D_growing_3=NULL
D_growing_4=NULL
D_growing_5=NULL
D_growing_6=NULL
D_growing_7=NULL
D_growing_8=NULL

Carriers_leaving=NULL
Asymptomatic_leaving=NULL

Pregnant_sold=NULL
Sacrificed=NULL

Treatment_weaned=NULL
Treatment_growing_and_pregnant=NULL

Completed_raising.stage_1=NULL
Completed_raising.stage_2=NULL
Completed_raising.stage_3=NULL
Completed_raising.stage_4=NULL
Completed_raising.stage_5=NULL
Completed_raising.stage_6=NULL
Completed_raising.stage_7=NULL
Completed_raising.stage_8=NULL
Completed_raising.stage_9=NULL
Completed_raising.stage_10=NULL
Completed_raising.stage_11=NULL
Completed_raising.stage_12=NULL

Event_trigger=NULL

#####Monte Carlo Simulation
for(i in  1:n){
  SALL_par<- c(
    omega[i], #Rate of immunity loss
    gamma[i], #Recovery rate for A and I
    m[i], #Rate of individuals in A becoming carriers
    u1[i], #Probability of infected weaned calves becoming severely infected
    w[i], #Rate of individuals in I becoming carriers
    d1[i], #Death rate for weaned calves in I
    v[i], #Vaccination effect on mortality in calves  
    eta[i], #Recovery rate C (Heifers)
    u2[i], #Probability of infected growing heifers becoming severely infected
    d2[i], #Death rate for growing heifers in I
    u3[i], #Probability of infected pregnant heifers becoming severely infected
    a1[i],
    a2[i],
    lambda1[i],  #Shedding rate for individuals in A and C (Raising stage 1)
    epsilon1[i], #Shedding rate for individuals in I (Raising stage 1)
    lambda2[i],  #Shedding rate for individuals in A and C (Raising stage 2)
    epsilon2[i], #Shedding rate for individuals in I (Raising stage 2)
    lambda3[i],  #Shedding rate for individuals in A and C (Raising stage 3)
    epsilon3[i], #Shedding rate for individuals in I (Raising stage 3)
    lambda4[i],  #Shedding rate for individuals in A and C (Raising stage 4)
    epsilon4[i], #Shedding rate for individuals in I (Raising stage 4)
    lambda5[i],  #Shedding rate for individuals in A and C (Raising stage 5)
    epsilon5[i], #Shedding rate for individuals in I (Raising stage 5)
    lambda6[i],  #Shedding rate for individuals in A and C (Raising stage 6)
    epsilon6[i], #Shedding rate for individuals in I (Raising stage 6)
    lambda7[i],  #Shedding rate for individuals in A and C (Raising stage 7)
    epsilon7[i], #Shedding rate for individuals in I (Raising stage 7)
    lambda8[i],  #Shedding rate for individuals in A and C (Raising stage 8)
    epsilon8[i], #Shedding rate for individuals in I (Raising stage 8)
    lambda9[i],  #Shedding rate for individuals in A and C (Raising stage 9)
    epsilon9[i], #Shedding rate for individuals in I (Raising stage 9)
    lambda10[i],  #Shedding rate for individuals in A and C (Raising stage 10)
    epsilon10[i], #Shedding rate for individuals in I (Raising stage 10)
    lambda11[i],  #Shedding rate for individuals in A and C (Raising stage 11)
    epsilon11[i], #Shedding rate for individuals in I (Raising stage 11)
    lambda12[i],  #Shedding rate for individuals in A and C (Raising stage 12)
    epsilon12[i] #Shedding rate for individuals in I (Raising stage 12)
  ) 
  
  SALL_dyn(0, SALL_var, SALL_par)
  
  SALL_time <- seq(0,duration,1) 
  
  
  SALL_sol [[i]] <-lsoda(y = SALL_var, times = SALL_time, func = SALL_dyn,
                         parms = SALL_par
                         , rootfun = rootfunc,
                         events = list(func=eventfunc, root=T)
  )
  # List1[[length(List1)+1]] = SALL_par[1]
  # List2[[length(List2)+1]] = SALL_par[2]
  # List3[[length(List3)+1]] = SALL_par[3]
  # List4[[length(List4)+1]] = SALL_par[4]
  # List5[[length(List5)+1]] = SALL_par[5]
  # List6[[length(List6)+1]] = SALL_par[6]
  # List7[[length(List7)+1]] = SALL_par[7]
  # List8[[length(List8)+1]] = SALL_par[8]
  # List9[[length(List9)+1]] = SALL_par[9]
  # List10[[length(List10)+1]] = SALL_par[10]
  # List11[[length(List11)+1]] = SALL_par[11]
  # List12[[length(List12)+1]] = SALL_par[12]
  # List13[[length(List13)+1]] = SALL_par[13]
  # List14[[length(List14)+1]] = SALL_par[14]
  
}

# for (i in 1:length(List1)) {
#   List1[[i]] <- rep(List1[[i]], duration+1)
# }
# List1<- lapply(List1, as.data.frame)
# List1<- lapply(List1, setNames, "omega")
# 
# for (i in 1:length(List2)) {
#   List2[[i]] <- rep(List2[[i]], duration+1)
# }
# List2<- lapply(List2, as.data.frame)
# List2<- lapply(List2, setNames, "gamma")
# 
# 
# for (i in 1:length(List3)) {
#   List3[[i]] <- rep(List3[[i]], duration+1)
# }
# List3<- lapply(List3, as.data.frame)
# List3<- lapply(List3, setNames, "m")
# 
# for (i in 1:length(List4)) {
#   List4[[i]] <- rep(List4[[i]], duration+1)
# }
# List4<- lapply(List4, as.data.frame)
# List4<- lapply(List4, setNames, "u1")
# 
# for (i in 1:length(List5)) {
#   List5[[i]] <- rep(List5[[i]], duration+1)
# }
# List5<- lapply(List5, as.data.frame)
# List5<- lapply(List5, setNames, "w")
# 
# for (i in 1:length(List6)) {
#   List6[[i]] <- rep(List6[[i]], duration+1)
# }
# List6<- lapply(List6, as.data.frame)
# List6<- lapply(List6, setNames, "d1")
# 
# for (i in 1:length(List7)) {
#   List7[[i]] <- rep(List7[[i]], duration+1)
# }
# List7<- lapply(List7, as.data.frame)
# List7<- lapply(List7, setNames, "v")
# 
# for (i in 1:length(List8)) {
#   List8[[i]] <- rep(List8[[i]], duration+1)
# }
# List8<- lapply(List8, as.data.frame)
# List8<- lapply(List8, setNames, "eta")
# 
# for (i in 1:length(List9)) {
#   List9[[i]] <- rep(List9[[i]], duration+1)
# }
# List9<- lapply(List9, as.data.frame)
# List9<- lapply(List9, setNames, "u2")
# 
# for (i in 1:length(List10)) {
#   List10[[i]] <- rep(List10[[i]], duration+1)
# }
# List10<- lapply(List10, as.data.frame)
# List10<- lapply(List10, setNames, "d2")
# 
# for (i in 1:length(List11)) {
#   List11[[i]] <- rep(List11[[i]], duration+1)
# }
# List11<- lapply(List11, as.data.frame)
# List11<- lapply(List11, setNames, "u3")
# 
# for (i in 1:length(List12)) {
#   List12[[i]] <- rep(List12[[i]], duration+1)
# }
# List12<- lapply(List12, as.data.frame)
# List12<- lapply(List12, setNames, "a1")
# 
# for (i in 1:length(List13)) {
#   List13[[i]] <- rep(List13[[i]], duration+1)
# }
# List13<- lapply(List13, as.data.frame)
# List13<- lapply(List13, setNames, "a2")
# 
# for (i in 1:length(List14)) {
#   List14[[i]] <- rep(List14[[i]], duration+1)
# }
# List14<- lapply(List14, as.data.frame)
# List14<- lapply(List14, setNames, "z")
# 
# 
# SALL_sol<- mapply(data.table, SALL_sol, List1,  SIMPLIFY=FALSE)
# rm(List1)
# SALL_sol<- mapply(data.table, SALL_sol, List2,  SIMPLIFY=FALSE)
# rm(List2)
# SALL_sol<- mapply(data.table, SALL_sol, List3,  SIMPLIFY=FALSE)
# rm(List3)
# SALL_sol<- mapply(data.table, SALL_sol, List4,  SIMPLIFY=FALSE)
# rm(List4)
# SALL_sol<- mapply(data.table, SALL_sol, List5,  SIMPLIFY=FALSE)
# rm(List5)
# SALL_sol<- mapply(data.table, SALL_sol, List6,  SIMPLIFY=FALSE)
# rm(List6)
# SALL_sol<- mapply(data.table, SALL_sol, List7,  SIMPLIFY=FALSE)
# rm(List7)
# SALL_sol<- mapply(data.table, SALL_sol, List8,  SIMPLIFY=FALSE)
# rm(List8)
# SALL_sol<- mapply(data.table, SALL_sol, List9,  SIMPLIFY=FALSE)
# rm(List9)
# SALL_sol<- mapply(data.table, SALL_sol, List10,  SIMPLIFY=FALSE)
# rm(List10)
# SALL_sol<- mapply(data.table, SALL_sol, List11,  SIMPLIFY=FALSE)
# rm(List11)
# SALL_sol<- mapply(data.table, SALL_sol, List12,  SIMPLIFY=FALSE)
# rm(List12)
# SALL_sol<- mapply(data.table, SALL_sol, List13,  SIMPLIFY=FALSE)
# rm(List13)
# SALL_sol<- mapply(data.table, SALL_sol, List14,  SIMPLIFY=FALSE)
# rm(List14)


#Extract model parameters
# Parameters<- data.frame(omega, gamma, eta, m, w, u1, u2, u3, d1, d2, a1, a2, z
# ) 

#Check one iteration to see if the model ran correctly
# one.iteration<- as.data.frame(SALL_sol) #Here I am just adding all compartments from one iteration to see how many individuals were in the operation at day 730 (if I get weird numbers it means that something is wrong)
# one.iteration$S_weaned_1_init+one.iteration$A_weaned_1_init+one.iteration$I_weaned_1_init+one.iteration$C_weaned_1_init+one.iteration$R_weaned_1_init+
# one.iteration$S_growing_2_init+one.iteration$A_growing_2_init+one.iteration$I_growing_2_init+one.iteration$C_growing_2_init+one.iteration$R_growing_2_init+
# one.iteration$S_growing_3_init+one.iteration$A_growing_3_init+one.iteration$I_growing_3_init+one.iteration$C_growing_3_init+one.iteration$R_growing_3_init+
# one.iteration$S_growing_4_init+one.iteration$A_growing_4_init+one.iteration$I_growing_4_init+one.iteration$C_growing_4_init+one.iteration$R_growing_4_init+
# one.iteration$S_growing_5_init+one.iteration$A_growing_5_init+one.iteration$I_growing_5_init+one.iteration$C_growing_5_init+one.iteration$R_growing_5_init+
# one.iteration$S_pregnant_5_init+one.iteration$A_pregnant_5_init+one.iteration$I_pregnant_5_init+one.iteration$C_pregnant_5_init+one.iteration$R_pregnant_5_init+
# one.iteration$S_pregnant_7_init+one.iteration$A_pregnant_7_init+one.iteration$I_pregnant_7_init+one.iteration$C_pregnant_7_init+one.iteration$R_pregnant_7_init+
# one.iteration$S_growing_8_init+one.iteration$A_growing_8_init+one.iteration$I_growing_8_init+one.iteration$C_growing_8_init+one.iteration$R_growing_8_init

#####ECONOMIC PARAMETERS#####
###Costs
#Per day raising costs (Data about heifer-raising costs by age was extracted from Karszes and Hill 2020. 
#A formula obtained from heifer-raising costs in FINBIN database (University of Minnesota 2022) 
#was applied as a percentage modifier (range:0 to 100%) to account for the influence of herd size 
#on heifer-raising costs (Hrn))

#Treatment costs for weaned heifers (Trw) and growing heifers (Trg)
Trw<- 79.7 #cost range $79.65-99.75
Trg<- 99.8
#Rendering costs
Ren=125 

###Income
##Daily income (USD/day)
Daily.income= 2.5
#Return from sale of animal parts (USD/carcass) 
Car= 125




#This tell us how long it took for the model to run

