###############################################################
## Accompanying code for paper: 'Is the association between  ##
## childhood leukaemia and population mixing an artefact of  ##
## focusing on 'clusters' of cases? (Berrie et al., 2018)    ##
###############################################################

###################
## Load packages ##
###################

rm(list=ls()); 
library(MASS); library(ggplot2); library(reshape); 
library(extrafont); library(gridExtra); library(grid); library(VGAM); library(Matrix)
library(boot); library(pscl);library(SuppDists);library(distr);library(distrEx);library(stringr)

#############################################
## Graphics settings - colourblind palette ##
#############################################

cbPal   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Dark    <- c("#3333FF", "#CC0000"); Lite <- c("#99CCFF", "#FF9999")

############################################
## Read in the Yorkshire & Humber dataset ##
############################################

## Omitted as dataset not publically available

#############################################################
## GenData Function (adapted from Ruscio & Kaczetow, 2008) ##
#############################################################

GenData <- function(Supp.Data=NULL,n.Fact=0,Max.Trials=5,Init.Mult=1,seed=0,Emp=TRUE,Target.Corr=NULL,N=NULL,k=NULL) {
  ################################################################################################################
  # Initialize variables and (if applicable) set random number seed (step 1) -------------------------------------
  if (Emp) {
    N   <- dim(Supp.Data)[1]                  # Number of cases                                                   
    k   <- dim(Supp.Data)[2] }                # Number of variables                                               
  Data    <- matrix(0,nrow=N,ncol=k)          # Matrix to store the simulated data                                
  Distributions <- matrix(0,nrow=N,ncol=k)    # Matrix to store each variable?s score distribution                
  Iteration <- 0                              # Iteration counter                                                 
  Best.RMSR <- 1                              # Lowest RMSR correlation                                           
  Trials.Without.Improvement <- 0             # Trial counter                                                     
  if (seed != 0) set.seed(seed)               # If user specified a nonzero seed, set it                          
  ################################################################################################################
  # Generate distribution for each variable (step 2) -------------------------------------------------------------
  if (Emp) for (i in 1:k) Distributions[,i] <- sort(sample(Supp.Data[,i],N,replace=TRUE)) else {
    Distributions[,1] <- sort(rnegbin(N,1300,1.6))              # 0 - 14 population
    Distributions[,2] <- sort(rnegbin(N,26,0.7))                # Area              
    Distributions[,3] <- sort(rnegbin(N,500,1.8))              # 'Post' in-migration
    Distributions[,4] <- sort(rnegbin(N,6500,2.0))  }           # Total population  
  # This implementation of GenData bootstraps each variable?s score distribution from a supplied data set.        
  # Users should modify this block of the program, as needed, to generate the desired distribution(s).            
  # For example, to sample from chi-square distributions with 2 df, replace the 2nd line in this block with:      
  #   Distributions[,i] <- sort(rchisq(N,df=2))                                                                   
  # Or, one can drop the loop and use a series of commands that samples variables from specified populations:     
  #   Distributions[,1] <- sort(rnorm(N,0,1))         # Standard normal distribution                              
  #   Distributions[,2] <- sort(runif(N,0,1))         # Uniform distribution ranging from 0 - 1                   
  #   Distributions[,3] <- sort(rlnorm(N,0,1))        # Log-normal distribution, log scale M = 0, SD = 1          
  #   Distributions[,4] <- sort(rexp(N,rate=1))       # Exponential distribution with rate = 1                    
  #   Distributions[,5] <- sort(rpois(N,lambda=4))    # Poisson distribution with lambda = 4                      
  #   Distributions[,6] <- sort(rbinom(N,10,0.25)     # Binominal distribution, size = 10 and p = 0.25            
  #   Distributions[,7] <- sort(rbinom(N,2,0.25)      # Binary distribution with p = 0.25                         
  ################################################################################################################
  # All of the commands shown above draw random samples from specified population distributions. Alternatively,   
  # one can reproduce distributions without sampling error. For example, working with a supplied data set, one can
  # replace the 2nd line in this block with:                                                                      
  #   Distributions[,i] <- Supp.Data[,i]                                                                          
  ################################################################################################################
  # Alternatively, idealized distributions can be reproduced. For example, uniform quantiles can be created and   
  # used to generate data from common distributions:                                                              
  #   Uniform.Quantiles <- seq(from = 0, to = 1, length = (N + 2))[2:(N + 1)] # quantiles 0, 1 dropped            
  #   Distributions[,1] <- qnorm(Uniform.Quantiles,0,1)   # Standard normal distribution                          
  #   Distributions[,2] <- qunif(Uniform.Quantiles,0,1)   # Uniform distribution ranging from 0 to 1              
  #   Distributions[,3] <- qchisq(Uniform.Quantiles,df=2) # Chi-square distribution with 2 df                     
  ################################################################################################################
  # Note that when score distributions are generated from specified populations rather than bootstrapped from a   
  # supplied dataset, the user must provide the target correlation matrix (see the next block). This is true      
  # regardless of whether the distributions incorporate sampling error.                                           
  ################################################################################################################
  # Calculate and store a copy of the target correlation matrix (step 3) -----------------------------------------
  if (Emp) Target.Corr <- cor(Supp.Data) 
  Intermediate.Corr    <- Target.Corr
  # This implementation of GenData calculates the target correlation matrix from a supplied dataset.              
  # Alternatively, the user can modify the program to generate data with user-defined sample size, number of      
  # variables and target correlation matrix by redefining the function as follows:                                
  #   GenData <- function(N,k,Target.Corr,n.Fact=0,Max.Trials=5,Init.Mult=1,seed=0)                               
  # In this case, one would also remove the program lines that calculate N, k, and Target.Corr.                   
  # To generate data in which variables are uncorrelated, one would remove the SsortT function from step 2        
  # and terminate the program before step 3 begins by returning the Distributions object as the dataset.          
  ################################################################################################################
  # If number of latent factors was not specified, determine it through parallel analysis (step 4) ---------------
  if (n.Fact == 0) {
    Eigenvalues.Observed <- eigen(Intermediate.Corr)$values
    Eigenvalues.Random <- matrix(0, nrow = 100, ncol = k)
    Random.Data <- matrix(0, nrow = N, ncol = k)
    for (i in 1:100) {
      for (j in 1:k) Random.Data[,j] <- sample(Distributions[,j], size = N, replace = TRUE)
      Eigenvalues.Random[i,] <- eigen(cor(Random.Data))$values }
    Eigenvalues.Random <- apply(Eigenvalues.Random, 2, mean)        # calculate mean eigenvalue for each factor     
    n.Fact <- max(1, sum(Eigenvalues.Observed > Eigenvalues.Random)) }
  ################################################################################################################
  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------
  Shared.Comp <- matrix(rnorm(N*n.Fact,0,1),nrow=N,ncol=n.Fact)
  Unique.Comp <- matrix(rnorm(N*k,0,1),nrow=N,ncol=k)
  Shared.Load <- matrix(0,nrow=k,ncol=n.Fact)
  Unique.Load <- matrix(0,nrow=k,ncol=1)
  ################################################################################################################
  # Begin loop that ends when specified number of iterations pass without improvement in RMSR correlation --------
  while (Trials.Without.Improvement < Max.Trials) {
    Iteration <- Iteration + 1
    ############################################################################################################
    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) -----------------------
    Fact.Anal <- Factor.Analysis(Intermediate.Corr, Corr.Matrix = TRUE, n.Fact = n.Fact)
    if (n.Fact == 1) Shared.Load[,1] <- Fact.Anal$loadings else 
      Shared.Load <- Fact.Anal$loadings
    Shared.Load[Shared.Load > 1] <- 1
    Shared.Load[Shared.Load < -1] <- -1
    if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
    Shared.Load.sq <- Shared.Load * Shared.Load
    for (i in 1:k)
      if (sum(Shared.Load.sq[i,]) < 1) Unique.Load[i,1] <- (1 - sum(Shared.Load.sq[i,])) else 
        Unique.Load[i,1] <- 0
    Unique.Load <- sqrt(Unique.Load)
    for (i in 1:k)
      Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
    # the %*% operator = matrix multiplication, and the t() function = transpose (both used again in step 13)   
    ############################################################################################################
    # Replace normal with nonnormal distributions (step 9) -----------------------------------------------------
    for (i in 1:k) {
      Data <- Data[sort.list(Data[,i]),]
      Data[,i] <- Distributions[,i] }
    ############################################################################################################
    # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12) ----------
    Reproduced.Corr <- cor(Data)
    Resid.Corr <- Target.Corr - Reproduced.Corr
    RMSR <- sqrt(sum(Resid.Corr[lower.tri(Resid.Corr)]*Resid.Corr[lower.tri(Resid.Corr)])/(0.5*(k*k-k)))
    if (RMSR < Best.RMSR) {
      Best.RMSR <- RMSR
      Best.Corr <- Intermediate.Corr
      Best.Res <- Resid.Corr
      Intermediate.Corr <- Intermediate.Corr + Init.Mult * Resid.Corr
      Trials.Without.Improvement <- 0 } else {
        Trials.Without.Improvement <- Trials.Without.Improvement + 1
        Current.Multiplier <- Init.Mult * .5 ^ Trials.Without.Improvement
        Intermediate.Corr <- Best.Corr + Current.Multiplier * Best.Res }
  } # end of the while loop
  ################################################################################################################
  # Construct the data set with the lowest RMSR correlation (step 13) --------------------------------------------
  Fact.Anal <- Factor.Analysis(Best.Corr, Corr.Matrix = TRUE, n.Fact = n.Fact)
  if (n.Fact == 1) Shared.Load[,1] <- Fact.Anal$loadings else 
    Shared.Load <- Fact.Anal$loadings
  Shared.Load[Shared.Load > 1] <- 1
  Shared.Load[Shared.Load < -1] <- -1
  if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
  Shared.Load.sq <- Shared.Load * Shared.Load
  for (i in 1:k)
    if (sum(Shared.Load.sq[i,]) < 1) Unique.Load[i,1] <- (1 - sum(Shared.Load.sq[i,])) else 
      Unique.Load[i,1] <- 0
  Unique.Load <- sqrt(Unique.Load)
  for (i in 1:k)
    Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
  Data <- apply(Data, 2, scale)       # standardizes each variable in the matrix                                  
  for (i in 1:k) {
    Data <- Data[sort.list(Data[,i]),]
    Data[,i] <- Distributions[,i] }
  ################################################################################################################
  # Report the results and return the simulated data set (step 14) -----------------------------------------------
  Iteration <- Iteration - Max.Trials
  #cat("\nN =",N,", k =",k,",",Iteration,"Iterations,",n.Fact,"Factors, RMSR r =",round(Best.RMSR,3),"\n")
  return(Data) }

Factor.Analysis <- function(Data,Corr.Matrix=FALSE,Max.Iter=50,n.Fact=0) {
  Data <- as.matrix(Data)
  k <- dim(Data)[2]
  if (n.Fact == 0) n.Fact <- k
  if (!Corr.Matrix) Cor.Matrix <- cor(Data) else 
    Cor.Matrix <- Data
  Criterion <- .001
  Old.H2 <- rep(99, k)
  H2 <- rep(0, k)
  Change <- 1
  Iter <- 0
  Factor.Loadings <- matrix(nrow = k, ncol = n.Fact)
  while ((Change >= Criterion) & (Iter < Max.Iter)) {
    Iter <- Iter + 1
    Eig <- eigen(Cor.Matrix)
    L <- sqrt(Eig$values[1:n.Fact])
    for (i in 1:n.Fact)
      Factor.Loadings[,i] <- Eig$vectors[,i] * L[i]
    for (i in 1:k)
      H2[i] <- sum(Factor.Loadings[i,] * Factor.Loadings[i,])
    Change <- max(abs(Old.H2 - H2))
    Old.H2 <- H2
    diag(Cor.Matrix) <- H2 }
  if (n.Fact == k) n.Fact <- sum(Eig$values > 1)
  return(list(loadings = Factor.Loadings[,1:n.Fact], factors = n.Fact)) }

#####################################################
## Function to generate Poisson distributed cases  ##
#####################################################

NullSim   <- function(N, id, x){
  #y      <- rpois(N, x*setPR)
  y       <- rbinom(N, x, setPR)
  Dat     <- data.frame(Id = id, Pop = x, Obs = y)
  PR      <- sum(Dat$Obs)/sum(Dat$Pop)
  Dat     <- cbind(Dat, Exp = Dat$Pop*PR)
  return(Dat)
} 

######################################
## Set-up target correlation matrix ##
######################################

ObsCor  <- matrix(c(1, -0.2961142, 0.8946106, 0.9707670, 
             -0.2961142, 1, -0.3227724, -0.2930421,
             0.8946106, -0.3227724, 1, 0.9178227,
             0.9707670, -0.2930421, 0.9178227, 1), nrow = 4, ncol = 4)


###########################################################################
## Simulate 10000 datasets using GenData function to create correlation  ##
## structure and approximate distributions of observed dataset           ##
## Use NullSim function to generate cases of childhood leukaemia based   ##
## solely on the population of 0-14 year olds. Perform sub-region and    ##
## region-wide methods on simulated datasets and store results.          ##
###########################################################################

##########################################################################################
## Set-up simulation with:                                                              ##
## (1). Seed (2.) Record start time (3). Set 5-year incidence rate of leukaemia         ##
## (4). Set N = number of electoral wards (5). Set k = number of variables to generate: ##
## 0-14 population, area, 'post' in-migration, total population                         ##
## (6). Set Nsim = number of simulations (10000)                                        ##
##########################################################################################

set.seed(1123)
Beg     <- Sys.time()
setPR   <- 0.0002         # Set incidence rate for 5 year period - test against this
N       <- 532            # Number of EWs to generate
k       <- 4              # Number of variables: 0-14 population, area, 'post' in-migration, total population
Nsim    <- 10000

###########################################################################################################
## Create empty vectors to store estimates and p-values from each method performed on the 10000 datasets ##
###########################################################################################################

## Sub-region approach empty vectors
p.binom1 <- p.binom2 <- p.binom3 <- p.binom4 <- p.binom5 <- p.binom6 <- p.binom7 <- NULL
p.binom8 <- p.binom9 <- p.binom10 <- p.binom11 <- p.binom12 <- p.binom13 <- p.binom14 <- p.binom15 <- NULL
p.binom16 <- p.binom17 <- NULL

c.binom1 <- c.binom2 <- c.binom3 <- c.binom4 <- c.binom5 <- c.binom6 <- c.binom7 <- c.binom8 <- NULL
c.binom9 <- c.binom10 <- c.binom11 <- c.binom12 <- c.binom13 <- c.binom14 <- c.binom15 <- NULL
c.binom16 <- c.binom17  <- NULL

## Region-wide approach empty vectors
p.PmD    <- NULL     
p.PmM    <- NULL
p.PmB    <- matrix(NA, nrow = Nsim, ncol = 2)
c.PmD    <- NULL     
c.PmM    <- NULL
c.PmB    <- matrix(NA, nrow = Nsim, ncol = 2)

######################################################################################################################
## Initiate simulation using for loop. Each loop: Generates a dataset using GenData and the set up defined above,   ##
## certain restraints are placed on the data as explained in the paper, each method is performed on the dataset and ##
## results are stored, this is repeated until 10000 datasets have been generated and analysed.                      ##
######################################################################################################################

for (itn in 1:Nsim){
  # Generate 4 variables using GenData, set seed N values apart, target correlation matrix = correlation matric from observed data
  # N = number of electoral wards to generate data for, k = number of variables to generate
  YHsim     <- data.frame(GenData(seed = (itn*N), Emp = FALSE, Target.Corr = ObsCor, N = N, k = k))         
  # Assign column names and IDs to generated variables
  YHsim     <- cbind(Id = 1:N, YHsim); names(YHsim) <- c("Id", "Pop", "Area", "InMig", "Tot_Pop")    
  
  # Replace lower than observed generated values with samples from values between minimum and median observed values
  if (min(YHsim$Tot_Pop) < 450) YHsim$Tot_Pop[which(YHsim$Tot_Pop < 450)] <- sample(450:6000, 1, replace = TRUE)
  if (min(YHsim$Pop) < 70) YHsim$Pop[which(YHsim$Pop < 70)] <- sample(70:1300, 1, replace = TRUE)
  if (min(YHsim$Area) < 0.17)  YHsim$Area[which(YHsim$Area < 0.17)] <- sample(0.17:16, 1, replace = TRUE)

  PrePop   <- YHsim[, 5] - YHsim[, 4]
  if (min(PrePop) < 450)  PrePop[which(PrePop < 450)] <- sample(450:3600, 1, replace = TRUE)
  
  # Calculate 'pre' in-migration proportions
  PreInMig <-  (YHsim[, 4]/PrePop)  
  YHsim    <- cbind(YHsim, PreInMig = PreInMig, PreDen = PrePop/YHsim[, 3])
  
  # Simulate Poisson distributed cases
  OutPois  <- NullSim(length(YHsim$Id), YHsim$Id, YHsim$Pop)  
  SimPois  <- cbind(OutPois, Den = YHsim$PreDen, Mig = YHsim$PreInMig)
  
  #########################
  ## Sub-region strategy ##
  #########################
  
  # Selection 1 - Low population density
  SimTwnDO   <- SimPois[order(SimPois$Den), ]
  SimTwnD2   <- SimTwnDO[0:ceiling((1/2)*length(SimTwnDO[, 1])), ]
  SimTwnD    <- SimPois[sample(SimTwnD2[, 1], 16, replace = FALSE), ]
  
  df1        <- data.frame(Cases = SimTwnD[, 3], NonCases = SimTwnD[, 2] - SimTwnD[, 3])
  p.binom1[itn] <- binom.test(sum(df1[, 1]), sum(df1[, 1]) + sum(df1[, 2]), p = setPR)$p.value
  c.binom1[itn] <- binom.test(sum(df1[, 1]), sum(df1[, 1]) + sum(df1[, 2]), p = setPR)$estimate
  
  # Selection 2 - High inward-migration
  SimTwnMO   <- SimPois[order(-SimPois$Mig), ]
  SimTwnM2   <- SimTwnMO[0:ceiling((1/2)*length(SimTwnMO[, 1])), ]
  SimTwnM    <- SimPois[sample(SimTwnM2[, 1], 16, replace = FALSE), ]
  
  df2        <- data.frame(Cases = SimTwnM[, 3], NonCases = SimTwnM[, 2] - SimTwnM[, 3])
  p.binom2[itn] <- binom.test(sum(df2[, 1]), sum(df2[, 1]) + sum(df2[, 2]), p = setPR)$p.value
  c.binom2[itn] <- binom.test(sum(df2[, 1]), sum(df2[, 1]) + sum(df2[, 2]), p= setPR)$estimate
  
  # Selection 3 - High incidence
  SimObs     <- cbind(SimPois, Inc = SimPois$Obs/SimPois$Exp)
  SimTwnI1   <- SimObs[order(-SimObs$Inc), ]
  SimTwnI2   <- SimTwnI1[0:ceiling((1/2)*length(SimTwnI1[, 1])), ]
  SimTwnI    <- SimObs[sample(SimTwnI2[, 1], 16, replace = FALSE), ]
  
  df3        <- data.frame(Cases = SimTwnI[, 3], NonCases = SimTwnI[, 2] - SimTwnI[, 3])
  p.binom3[itn] <- binom.test(sum(df3[, 1]), sum(df3[, 1]) + sum(df3[, 2]), p = setPR)$p.value
  c.binom3[itn] <- binom.test(sum(df3[, 1]), sum(df3[, 1]) + sum(df3[, 2]), p = setPR)$estimate
  
  # Selection 4 - Low population density and high inward-migration
  OrdD4      <- SimObs[order(SimObs$Den), ]
  SmpD4      <- OrdD4[1:ceiling(0.5*length(OrdD4[, 1])), ]
  OrdM4      <- SmpD4[order(-SmpD4$Mig), ]
  SmpM4      <- OrdM4[1:ceiling(0.5*length(OrdM4[, 1])), ]
  Smp4       <- SimObs[sample(SmpM4[, 1], 16, replace = FALSE), ]
  
  df4      <- data.frame(Cases = Smp4[, 3], NonCases = Smp4[, 2] - Smp4[, 3])
  p.binom4[itn] <- binom.test(sum(df4[, 1]), sum(df4[, 1]) + sum(df4[, 2]), p = setPR)$p.value
  c.binom4[itn] <- binom.test(sum(df4[, 1]), sum(df4[, 1]) + sum(df4[, 2]), p = setPR)$estimate
  
  # Selection 5 - High inward-migration and low population density
  OrdM5      <- SimObs[order(-SimObs$Mig), ]
  SmpM5      <- OrdM5[1:ceiling(0.5*length(OrdM5[, 1])), ]
  OrdD5      <- SmpM5[order(SmpM5$Den), ]
  SmpD5      <- OrdD5[1:ceiling(0.5*length(OrdD5[, 1])), ]
  Smp5       <- SimObs[sample(SmpD5[, 1], 16, replace = FALSE), ]
  
  df5      <- data.frame(Cases = Smp5[, 3], NonCases = Smp5[, 2] - Smp5[, 3])
  p.binom5[itn] <- binom.test(sum(df5[, 1]), sum(df5[, 1]) + sum(df5[, 2]), p = setPR)$p.value
  c.binom5[itn] <- binom.test(sum(df5[, 1]), sum(df5[, 1]) + sum(df5[, 2]), p = setPR)$estimate
  
  # Selection 6 - Low population density and high incidence
  OrdD6      <- SimObs[order(SimObs$Den), ]
  SmpD6      <- OrdD6[1:ceiling(0.5*length(OrdD6[, 1])), ]
  OrdI6      <- SmpD6[order(-SmpD6$Inc), ]
  SmpI6      <- OrdI6[1:ceiling(0.5*length(OrdI6[, 1])), ]
  Smp6       <- SimObs[sample(SmpI6[, 1], 16, replace = FALSE), ]
  
  df6      <- data.frame(Cases = Smp6[, 3], NonCases = Smp6[, 2] - Smp6[, 3])
  p.binom6[itn] <- binom.test(sum(df6[, 1]), sum(df6[, 1]) + sum(df6[, 2]), p = setPR)$p.value
  c.binom6[itn] <- binom.test(sum(df6[, 1]), sum(df6[, 1]) + sum(df6[, 2]), p = setPR)$estimate
  
  # Selection 7 - High incidence and low population density
  OrdI7      <- SimObs[order(-SimObs$Inc), ]
  SmpI7      <- OrdI7[1:ceiling(0.5*length(OrdI7[, 1])), ]
  OrdD7      <- SmpI7[order(SmpI7$Den), ]
  SmpD7      <- OrdD7[1:ceiling(0.5*length(OrdD7[, 1])), ]
  Smp7       <- SimObs[sample(SmpD7[, 1], 16, replace = FALSE), ]
  
  df7      <- data.frame(Cases = Smp7[, 3], NonCases = Smp7[, 2] - Smp7[, 3])
  p.binom7[itn] <- binom.test(sum(df7[, 1]), sum(df7[, 1]) + sum(df7[, 2]), p = setPR)$p.value
  c.binom7[itn] <- binom.test(sum(df7[, 1]), sum(df7[, 1]) + sum(df7[, 2]), p = setPR)$estimate
  
  # Selection 8 - High inward-migration and high incidence
  OrdM8      <- SimObs[order(-SimObs$Mig), ]
  SmpM8      <- OrdM8[1:ceiling(0.5*length(OrdM8[, 1])), ]
  OrdI8      <- SmpM8[order(-SmpM8$Inc), ]
  SmpI8      <- OrdI8[1:ceiling(0.5*length(OrdI8[, 1])), ]
  Smp8       <- SimObs[sample(SmpI8[, 1], 16, replace = FALSE), ]
  
  df8      <- data.frame(Cases = Smp8[, 3], NonCases = Smp8[, 2] - Smp8[, 3])
  p.binom8[itn] <- binom.test(sum(df8[, 1]), sum(df8[, 1]) + sum(df8[, 2]), p = setPR)$p.value
  c.binom8[itn] <- binom.test(sum(df8[, 1]), sum(df8[, 1]) + sum(df8[, 2]), p = setPR)$estimate
  
  # Selection 9 - High incidence and high inward-migration
  OrdI9      <- SimObs[order(-SimObs$Inc), ]
  SmpI9      <- OrdI9[1:ceiling(0.5*length(OrdI9[, 1])), ]
  OrdM9      <- SmpI9[order(-SmpI9$Mig), ]
  SmpM9      <- OrdM9[1:ceiling(0.5*length(OrdM9[, 1])), ]
  Smp9       <- SimObs[sample(SmpM9[, 1], 16, replace = FALSE), ]
  
  df9      <- data.frame(Cases = Smp9[, 3], NonCases = Smp9[, 2] - Smp9[, 3])
  p.binom9[itn] <- binom.test(sum(df9[, 1]), sum(df9[, 1]) + sum(df9[, 2]), p = setPR)$p.value
  c.binom9[itn] <- binom.test(sum(df9[, 1]), sum(df9[, 1]) + sum(df9[, 2]), p= setPR)$estimate
  
  # Selection 10 - Low population density, high inward-migration and high incidence
  OrdD10    <- SimObs[order(SimObs$Den), ]
  SmpD10    <- OrdD10[1:ceiling(0.5*length(OrdD10[, 1])), ]
  OrdM10    <- SmpD10[order(-SmpD10$Mig), ]
  SmpM10    <- OrdM10[1:ceiling(0.5*length(OrdM10[, 1])), ]
  OrdI10    <- SmpM10[order(-SmpM10$Inc), ]
  Smp10     <- OrdI10[1:ceiling(0.5*length(OrdI10[, 1])), ]
  Smp10     <- SimObs[sample(Smp10[, 1], 16, replace = FALSE), ]
  
  df10      <- data.frame(Cases = Smp10[, 3], NonCases = Smp10[, 2] - Smp10[, 3])
  p.binom10[itn] <- binom.test(sum(df10[, 1]), sum(df10[, 1]) + sum(df10[, 2]), p = setPR)$p.value
  c.binom10[itn] <- binom.test(sum(df10[, 1]), sum(df10[, 1]) + sum(df10[, 2]), p = setPR)$estimate
  
  # Selection 11 - Low population density, high incidence and high inward-migration 
  OrdD11    <- SimObs[order(SimObs$Den), ]
  SmpD11    <- OrdD11[1:ceiling(0.5*length(OrdD11[, 1])), ]
  OrdI11    <- SmpD11[order(-SmpD11$Inc), ]
  SmpI11    <- OrdI11[1:ceiling(0.5*length(OrdI11[, 1])), ]
  OrdM11    <- SmpI11[order(-SmpI11$Mig), ]
  Smp11     <- OrdM11[1:ceiling(0.5*length(OrdM11[, 1])), ]
  Smp11    <- SimObs[sample(Smp11[, 1], 16, replace = FALSE), ]
  
  df11      <- data.frame(Cases = Smp11[, 3], NonCases = Smp11[, 2] - Smp11[, 3])
  p.binom11[itn] <- binom.test(sum(df11[, 1]), sum(df11[, 1]) + sum(df11[, 2]), p = setPR)$p.value
  c.binom11[itn] <- binom.test(sum(df11[, 1]), sum(df11[, 1]) + sum(df11[, 2]), p = setPR)$estimate
  
  # Selection 12 - High inward-migration, low population density and high incidence
  OrdM12    <- SimObs[order(-SimObs$Mig), ]
  SmpM12    <- OrdM12[1:ceiling(0.5*length(OrdM12[, 1])), ]
  OrdD12    <- SmpM12[order(SmpM12$Den), ]
  SmpD12    <- OrdD12[1:ceiling(0.5*length(OrdD12[, 1])), ]
  OrdI12    <- SmpD12[order(-SmpD12$Inc), ]
  Smp12     <- OrdI12[1:ceiling(0.5*length(OrdI12[, 1])), ]
  Smp12    <- SimObs[sample(Smp12[, 1], 16, replace = FALSE), ]
  
  df12      <- data.frame(Cases = Smp12[, 3], NonCases = Smp12[, 2] - Smp12[, 3])
  p.binom12[itn] <- binom.test(sum(df12[, 1]), sum(df12[, 1]) + sum(df12[, 2]), p = setPR)$p.value
  c.binom12[itn] <- binom.test(sum(df12[, 1]), sum(df12[, 1]) + sum(df12[, 2]), p = setPR)$estimate
  
  # Selection 13 - High inward-migration, high incidence and low population density
  OrdM13    <- SimObs[order(-SimObs$Mig), ]
  SmpM13    <- OrdM13[1:ceiling(0.5*length(OrdM13[, 1])), ]
  OrdI13    <- SmpM13[order(-SmpM13$Inc), ]
  SmpI13    <- OrdI13[1:ceiling(0.5*length(OrdI13[, 1])), ]
  OrdD13    <- SmpI13[order(SmpI13$Den), ]
  Smp13     <- OrdI13[1:ceiling(0.5*length(OrdD13[, 1])), ]
  Smp13    <- SimObs[sample(Smp13[, 1], 16, replace = FALSE), ]
  
  df13      <- data.frame(Cases = Smp13[, 3], NonCases = Smp13[, 2] - Smp13[, 3])
  p.binom13[itn] <- binom.test(sum(df13[, 1]), sum(df13[, 1]) + sum(df13[, 2]), p = setPR)$p.value
  c.binom13[itn] <- binom.test(sum(df13[, 1]), sum(df13[, 1]) + sum(df13[, 2]), p = setPR)$estimate
  
  # Selection 14 - High incidence, low population density and high inward-migration
  OrdI14    <- SimObs[order(-SimObs$Inc), ]
  SmpI14    <- OrdI14[1:ceiling(0.5*length(OrdI14[, 1])), ]
  OrdD14    <- SmpI14[order(SmpI14$Den), ]
  SmpD14    <- OrdD14[1:ceiling(0.5*length(OrdD14[, 1])), ]
  OrdM14    <- SmpD14[order(-SmpD14$Mig), ]
  Smp14     <- OrdM14[1:ceiling(0.5*length(OrdM14[, 1])), ]
  Smp14    <- SimObs[sample(Smp14[, 1], 16, replace = FALSE), ]
  
  df14      <- data.frame(Cases = Smp14[, 3], NonCases = Smp14[, 2] - Smp14[, 3])
  p.binom14[itn] <- binom.test(sum(df14[, 1]), sum(df14[, 1]) + sum(df14[, 2]), p = setPR)$p.value
  c.binom14[itn] <- binom.test(sum(df14[, 1]), sum(df14[, 1]) + sum(df14[, 2]), p = setPR)$estimate
  
  # Selection 15 - High incidence, high inward-migration and low population density
  OrdI15    <- SimObs[order(-SimObs$Inc), ]
  SmpI15    <- OrdI15[1:ceiling(0.5*length(OrdI15[, 1])), ]
  OrdM15    <- SmpI15[order(-SmpI15$Mig), ]
  SmpM15    <- OrdM15[1:ceiling(0.5*length(OrdM15[, 1])), ]
  OrdD15    <- SmpM15[order(SmpM15$Den), ]
  Smp15    <- OrdD15[1:ceiling(0.5*length(OrdD15[, 1])), ]
  Smp15    <- SimObs[sample(Smp15[, 1], 16, replace = FALSE), ]
  
  df15     <- data.frame(Cases = Smp15[, 3], NonCases = Smp15[, 2] - Smp15[, 3])
  p.binom15[itn] <- binom.test(sum(df15[, 1]), sum(df15[, 1]) + sum(df15[, 2]), p = setPR)$p.value
  c.binom15[itn] <- binom.test(sum(df15[, 1]), sum(df15[, 1]) + sum(df15[, 2]), p = setPR)$estimate
  
  # Selection 16 - Random selection of 16 wards
  Smp16     <- SimObs[sample(SimObs[, 1], 16, replace = FALSE), ]
  
  df16     <- data.frame(Cases = Smp16[, 3], NonCases = Smp16[, 2] - Smp16[, 3])
  p.binom16[itn] <- binom.test(sum(df16[, 1]), sum(df16[, 1]) + sum(df16[, 2]), p = setPR)$p.value
  c.binom16[itn] <- binom.test(sum(df16[, 1]), sum(df16[, 1]) + sum(df16[, 2]), p = setPR)$estimate
  
  # Selection 17 - Incidence less than average
  SimObs     <- cbind(SimPois, Inc = SimPois$Obs/SimPois$Exp)
  SimTwnI17   <- SimObs[order(-SimObs$Inc), ]
  SimTwnI17   <- SimTwnI17[ceiling((1/2)*length(SimTwnI17[, 1])):length(SimTwnI17[, 1]), ]
  SimTwnLI    <- SimObs[sample(SimTwnI17[, 1], 16, replace = FALSE), ]
  
  df17        <- data.frame(Cases = SimTwnLI[, 3], NonCases = SimTwnLI[, 2] - SimTwnLI[, 3])
  p.binom17[itn] <- binom.test(sum(df17[, 1]), sum(df17[, 1]) + sum(df17[, 2]), p = setPR)$p.value
  c.binom17[itn] <- binom.test(sum(df17[, 1]), sum(df17[, 1]) + sum(df17[, 2]), p = setPR)$estimate
  
  ##########################
  ## Region-wide strategy ##
  ##########################
  
  # Create Poisson models
  SimHalf    <- SimPois[sample(SimPois[, 1], 266, replace = FALSE), ]
  PmDtmp     <- glm(Obs ~ offset(log(Pop)) + Den, data = SimHalf, family = poisson(link = log))
  PmMtmp     <- glm(Obs ~ offset(log(Pop)) + Mig, data = SimHalf, family = poisson(link = log))
  PmBtmp     <- glm(Obs ~ offset(log(Pop)) + Den + Mig, data = SimHalf, family = poisson(link = log))
  
  # Store point estimates
  c.PmD[itn]   <- PmDtmp$coefficients[2]
  c.PmM[itn]   <- PmMtmp$coefficients[2]
  c.PmB[itn, 1] <- PmBtmp$coefficients[2]
  c.PmB[itn, 2] <- PmBtmp$coefficients[3]
  
  # Store p-values
  p.PmD[itn]   <- summary(PmDtmp)$coefficients[8]
  p.PmM[itn]   <- summary(PmMtmp)$coefficients[8]
  p.PmB[itn, 1] <- summary(PmBtmp)$coefficients[11]
  p.PmB[itn, 2] <- summary(PmBtmp)$coefficients[12]
}
End   <- Sys.time()

## Find length of time code took to run
End-Beg

######################################################
## Calculate Type I error rate of sub-region method ##
######################################################

Bpdists <- data.frame(BE1 = p.binom1,  BE2 = p.binom2,  BE3 = p.binom3, BE4 = p.binom4, BE5 = p.binom5, BE6 = p.binom6, BE7 = p.binom7, BE8 = p.binom8, 
                      BE9 = p.binom9, BE10 = p.binom10, BE11 = p.binom11, BE12 = p.binom12, BE13 = p.binom13, BE14 = p.binom14, BE15 = p.binom15, BE16 = p.binom16, BE17 = p.binom17)

## Calculate percentages of critical p-values
B.1pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE1 <= x)}))/itn)  
B.2pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE2 <= x)}))/itn)
B.3pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE3 <= x)}))/itn)
B.4pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE4 <= x)}))/itn)
B.5pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE5 <= x)}))/itn)
B.6pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE6 <= x)}))/itn)
B.7pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE7 <= x)}))/itn)
B.8pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE8 <= x)}))/itn)
B.9pPct    <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE9 <= x)}))/itn)
B.10pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE10 <= x)}))/itn)
B.11pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE11 <= x)}))/itn)  
B.12pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE12 <= x)}))/itn)
B.13pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE13 <= x)}))/itn)
B.14pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE14 <= x)}))/itn)
B.15pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE15 <= x)}))/itn)
B.16pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE16 <= x)}))/itn)
B.17pPct   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(Bpdists$BE17 <= x)}))/itn)

B.1pPct; B.2pPct; B.3pPct; B.4pPct; B.5pPct; B.6pPct; B.7pPct; B.8pPct; B.9pPct; B.10pPct; B.11pPct; B.12pPct; B.13pPct; B.14pPct; B.15pPct; B.16pPct; B.17pPct

## Store p-values of sub-region method
write.table(Bpdists, "p-distributions binomial test - all combinations.csv", sep = ",", col.names = TRUE)

#Bpdists <- read.csv("p-distributions binomial test - all combinations .csv", sep = ",", header = TRUE)

#########################################################
## Calculate Type I error rates for region-wide method ##
#########################################################

Pct.p.PmM     <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(p.PmM <= x)}))/itn)
Pct.p.PmD     <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(p.PmD <= x)}))/itn)  
Pct.p.PmB.M   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(p.PmB[, 2] <= x)}))/itn)
Pct.p.PmB.D   <- c(100*do.call(c, lapply(c(0.05, 0.01, 0.001), function(x){sum(p.PmB[, 1] <= x)}))/itn)

Pct.p.PmD; Pct.p.PmM; Pct.p.PmB.D; Pct.p.PmB.M; 

# Store p-values of region-wide method
pdists <- data.frame(PmM = p.PmM, PmD = p.PmD, PmBM = p.PmB[, 2], PmBD = p.PmB[, 1])
write.table(pdists, "p-distributions regression - all combinations.csv", sep = ",", col.names = TRUE)

## Plot the p-values for the region-wide method
pdists <- data.frame("Inward-migration" = pdists[, 1], "Population density" = pdists[, 2], "Inward-migration adj. Population density" = pdists[, 3],
                     "Population density - adj. Inward-migration" = pdists[, 4])
pdists <- stack(pdists)
Xlab <- "p-values"; Ylab <- "Kernel Density"
Mlab    <- "Distribution of p-values for the 'whole region' selection strategy"; windows()
ggplot(pdists, aes(x = values)) + geom_density(aes(group = ind, colour = ind), size = 1.2) + 
  labs(title = Mlab, x = Xlab, y = Ylab, colour = NULL) + 
  scale_colour_manual(values = cbPal, name = "Covariate", label = c("Inward-Migration", "Population Density", "Inward-Migration (adj. Population Density)",
                                                             "Population Density (adj. Inward-Migration)")) +
  theme_bw() + theme(axis.title = element_text(family = "Times New Roman", size = 16), plot.title = element_text(size = 16, family = "Times New Roman", face = "bold"), legend.position = "top",
                   legend.text = element_text(family = "Times New Roman", size = 12), legend.title = element_text(family = "Times New Roman", size = 12, face = "bold"))


## Plot the p-values for the sub-region method
Bpdists <- data.frame("High Inward-Migration" = Bpdists[, 1], "Low Population Density" = Bpdists[, 2], "High In-ward Migration and Low Population Density" = Bpdists[, 3],
                     "High Incidence" = Bpdists[, 4])
Bpdists <- stack(Bpdists)
Xlab <- "p-values"; Ylab <- "Kernel Density"
Mlab    <- "Distribution of p-values for the 'sub-sample' selection strategy"; windows()
ggplot(Bpdists, aes(x = values)) + geom_density(aes(group = ind, colour = ind), size = 1.2) + 
  labs(title = Mlab, x = Xlab, y = Ylab, colour = NULL) + 
  scale_colour_manual(values = cbPal, name = "Selection", label = c("High Inward-Migration", "Low Population Density", "High Inward-Migration and Low Population Density",
                                                            "High Incidence")) +
  theme_bw() + theme(axis.title = element_text(family = "Times New Roman", size = 16), plot.title = element_text(size = 16, family = "Times New Roman", face = "bold"), legend.position = "top",
                   legend.text = element_text(family = "Times New Roman", size = 12), legend.title = element_text(family = "Times New Roman", size = 12,face = "bold"))

################################################################
## Calculation of Monte Carlo Error for the regression models ##
################################################################

var(c.PmM)
var(c.PmD)
var(c.PmB)

##########################################################
## Calculation of Monte Carlo Error for binomial tests ##
#########################################################

var(c.binom1)
var(c.binom2)
var(c.binom3)
var(c.binom4)
var(c.binom5)
var(c.binom6)
var(c.binom7)
var(c.binom8)
var(c.binom9)
var(c.binom10)
var(c.binom11)
var(c.binom12)
var(c.binom13)
var(c.binom14)
var(c.binom15)

#################################################
## Calculate 95% Empircal Confidence Intervals ##
#################################################

## Region-wide method confidence intervals
MCI  <- quantile(sort(c.PmM), c(0.025, 0.5, 0.975))
DCI  <- quantile(sort(c.PmD), c(0.025, 0.5, 0.975))
BDCI <- quantile(sort(c.PmB[, 1]), c(0.025, 0.5, 0.975))
BMCI <- quantile(sort(c.PmB[, 2]), c(0.025, 0.5, 0.975))

DCI; MCI; BDCI; BMCI

## Sub-region method confidence intervals
S1CI <- quantile(sort(c.binom1), c(0.025, 0.5, 0.975))
S2CI <- quantile(sort(c.binom2), c(0.025, 0.5, 0.975))
S3CI <- quantile(sort(c.binom3), c(0.025, 0.5, 0.975))
S4CI <- quantile(sort(c.binom4), c(0.025, 0.5, 0.975))
S5CI <- quantile(sort(c.binom5), c(0.025, 0.5, 0.975))
S6CI <- quantile(sort(c.binom6), c(0.025, 0.5, 0.975))
S7CI <- quantile(sort(c.binom7), c(0.025, 0.5, 0.975))
S8CI <- quantile(sort(c.binom8), c(0.025, 0.5, 0.975))
S9CI <- quantile(sort(c.binom9), c(0.025, 0.5, 0.975))
S10CI <- quantile(sort(c.binom10), c(0.025, 0.5, 0.975))
S11CI <- quantile(sort(c.binom11), c(0.025, 0.5, 0.975))
S12CI <- quantile(sort(c.binom12), c(0.025, 0.5, 0.975))
S13CI <- quantile(sort(c.binom13), c(0.025, 0.5, 0.975))
S14CI <- quantile(sort(c.binom14), c(0.025, 0.5, 0.975))
S15CI <- quantile(sort(c.binom15), c(0.025, 0.5, 0.975))

S1CI; S2CI; S3CI; S4CI; S5CI; S6CI; S7CI; S8CI; S9CI; S10CI; S11CI; S12CI; S13CI; S14CI; S15CI

####################
## Visualisation  ##
####################

##############################################################################
## Region-wide and sub-sample selection analyses performed on observed data ##
##############################################################################

## Omitted as observed dataset not publically available

Estimatesdf  <- data.frame(Estimates =  rbind(S1CI, S2CI, S3CI, S4CI, S5CI, S6CI, S7CI, S8CI, 
                                              S9CI, S10CI, S11CI, S12CI, S13CI, S14CI, S15CI),
  Label = c("S1: D", "S2: M", "S3: I", "S4: DM", "S5: MD", "S6: DI", "S7: ID", "S8: MI", 
                   "S9: IM", "S10: DMI", "S11: DIM", "S12: MDI", "S13: MID", "S14: IDM", "S15: IMD"), 
  Mod = c(rep("Sub-region - Simulated", 15)))

EstimatesdfSim <- Estimatesdf[1:15, ]
EstimatesdfSim$Label <- factor(EstimatesdfSim$Label, as.character(EstimatesdfSim$Label))

x11()
ggplot(EstimatesdfSim, aes(x = Label, y = Estimates.50., group = Mod)) +
  geom_point() +
  geom_errorbar(data = EstimatesdfSim, aes(ymin = Estimates.2.5., ymax = Estimates.97.5.), width = .2)

Estimatesdf$Label  <- factor(Estimatesdf$Label, unique(as.character(Estimatesdf$Label)))

#####################################################
## Figure 3 of paper - showing simulated data only ##
#####################################################

x11()
ggplot(Estimatesdf, aes(x = Label, y = Estimates.50., group = Mod, colour = Mod)) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(data = Estimatesdf, aes(ymin = Estimates.2.5., ymax = Estimates.97.5., colour = Mod), width = .5, position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0.0002, size = 0.2, colour = "black", linetype = "dashed") + theme_bw() +
  scale_y_continuous(name = "Estimate of Childhood Leukemia Incidence") + scale_x_discrete(name = "") +
  annotate("text", label = "Simulated Null", x = 1.0, y = 0.000225, size = 5, colour = "black")  +
  theme(axis.text = element_text(colour = "black", size = 12, family = "Times New Roman"),
        axis.title = element_text(colour = "black", size = 20, face = "bold", family = "Times New Roman"),
        plot.title = element_text(size = 16, face = "bold", family = "Times New Roman"), legend.position = c(0.1, 0.95),
        legend.text = element_text(size = 16, family = "Times New Roman"), axis.text.x = element_text(vjust = 0.5), legend.key = element_blank(), legend.title = element_blank())

#############################################################################
## Create dataframe of coefficients according to specific covariate values ##
#############################################################################

# Log scale
Coef10pc <- exp(0.1*MCI); Coef25pc <- exp(0.25*MCI); Coef50pc <- exp(0.5*MCI)

Coef10pcadj <- exp(0.1*BMCI); Coef25pcadj <- exp(0.25*BMCI); Coef50pcadj <- exp(0.5*BMCI)

Coef100d <- exp(100*DCI); Coef500d <- exp(500*DCI); Coef1000d <- exp(1000*DCI)

Coef100dadj <- exp(100*BDCI); Coef500dadj <- exp(500*BDCI); Coef1000dadj <- exp(1000*BDCI)

# Normal scale
# Coef10pc <- 0.1*MCI; Coef25pc <- 0.25*MCI; Coef50pc <- 0.5*MCI
# 
# Coef10pcadj <- 0.1*BMCI; Coef25pcadj <- 0.25*BMCI; Coef50pcadj <- 0.5*BMCI
# 
# Coef100d <- 100*DCI; Coef500d <- 500*DCI; Coef1000d <- 1000*DCI
# 
# Coef100dadj <- 100*BDCI; Coef500dadj <- 500*BDCI; Coef1000dadj <- 1000*BDCI

##################################################################################################################
## Plot confidence intervals for 500 persons/km^2 increase in population density and 25% increase in population ##
##################################################################################################################

Coeffsdf  <- data.frame(Coeffs = rbind(Coef25pc, Coef25pcadj, Coef500d, Coef500dadj),
                        Label= c("Inward-Migration \n25%",
                                 "Inward-Migration \n(adj. Population Density) - 25%",
                                 "Population Density \n500 persons per km^2",
                                 "Population Density \n(adj. Inward-Migration) \n500 persons per km^2"),
                        Mod=c(rep("Region-wide - Simulated", 4)))

Coeffsdf$Label  <- factor(Coeffsdf$Label, levels = unique(Coeffsdf$Label))

#############################################
## Figure 4 of paper - simulated data only ##
#############################################

x11()
ggplot(Coeffsdf,aes(x = Label, y = Coeffs.50., group = Mod, colour = Mod)) +
  scale_y_continuous(name = "Risk Ratio", breaks = c(0, 1, 2, 4, 6, 8, 10)) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(data = Coeffsdf, aes(ymin = Coeffs.2.5., ymax = Coeffs.97.5., colour = Mod), width = .5, position = position_dodge(width = .5)) +
  geom_hline(yintercept = 1.00, size = 0.2, colour = "black", linetype = "dashed") + theme_bw() +
  scale_x_discrete(name = "") +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(colour = "black", size=12, family = "Times New Roman"), axis.title = element_text(size = 16, family = "Times New Roman", face = "bold"), plot.title = element_text(size = 16, face = "bold", family = "Times New Roman"), legend.position = "top",
        legend.text = element_text(size = 12, family = "Times New Roman")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", label = "Simulated Null", x = 0.45, y = 1.25, size = 5, colour = "black", family = "Times New Roman") +
  coord_trans(y = "log10")

#########################################
## Summary of significant coefficients ##
#########################################

## Summary of proportion of inward-migration coefficient
summary(Sp.PmMSubHi)
summary(Sp.PmMSubLo)
summary(Sp.PmMtmpSub)

## Summary of population density coefficient
summary(Sp.PmDSubHi)
summary(Sp.PmDSubLo)
summary(Sp.PmDtmpSub)

## Summary of proportion of inward-migration coefficient (adjusted for population density)
summary(Sp.PmBMSubHi)
summary(Sp.PmBMSubLo)
summary(Sp.PmBMtmpSub)

## Summary of proportion of inward-migration coefficient (adjusted for population density)
summary(Sp.PmDSubHi)
summary(Sp.PmDSubLo)
summary(Sp.PmBDtmpSub)

###
S.binomS1tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS2tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS3tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS4tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS5tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS6tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS7tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS8tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS9tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS10tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS11tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS12tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS13tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS14tmp    <- matrix(NA, nrow = 10000, ncol = 2)
S.binomS15tmp    <- matrix(NA, nrow = 10000, ncol = 2)

S.binomS1tmp[, 1] <- c.binom1
S.binomS1tmp[, 2] <- p.binom1
S.binomS2tmp[, 1] <- c.binom2
S.binomS2tmp[, 2] <- p.binom2
S.binomS3tmp[, 1] <- c.binom3
S.binomS3tmp[, 2] <- p.binom3
S.binomS4tmp[, 1] <- c.binom4
S.binomS4tmp[, 2] <- p.binom4
S.binomS5tmp[, 1] <- c.binom5
S.binomS5tmp[, 2] <- p.binom5
S.binomS6tmp[, 1] <- c.binom6
S.binomS6tmp[, 2] <- p.binom6
S.binomS7tmp[, 1] <- c.binom7
S.binomS7tmp[, 2] <- p.binom7
S.binomS8tmp[, 1] <- c.binom8
S.binomS8tmp[, 2] <- p.binom8
S.binomS9tmp[, 1] <- c.binom9
S.binomS9tmp[, 2] <- p.binom9
S.binomS10tmp[, 1] <- c.binom10
S.binomS10tmp[, 2] <- p.binom10
S.binomS11tmp[, 1] <- c.binom11
S.binomS11tmp[, 2] <- p.binom11
S.binomS12tmp[, 1] <- c.binom12
S.binomS12tmp[, 2] <- p.binom12
S.binomS13tmp[, 1] <- c.binom13
S.binomS13tmp[, 2] <- p.binom13
S.binomS14tmp[, 1] <- c.binom14
S.binomS14tmp[, 2] <- p.binom14
S.binomS15tmp[, 1] <- c.binom15
S.binomS15tmp[, 2] <- p.binom15

S.binomS1tmpSub <- S.binomS1tmp[which(S.binomS1tmp[, 2] < 0.05), ]
S.binomS1SubHi  <- S.binomS1tmpSub[which(S.binomS1tmpSub[, 1] > 0.0002), ]
S.binomS1SubLo  <- S.binomS1tmpSub[which(S.binomS1tmpSub[, 1] < 0.0002), ]
S.binomS1Sub <- data.frame(cbind(-length(S.binomS1SubLo[, 1])/100, length(S.binomS1SubHi[, 1])/100))

S.binomS2tmpSub <- S.binomS2tmp[which(S.binomS2tmp[, 2] < 0.05), ]
S.binomS2SubHi  <- S.binomS2tmpSub[which(S.binomS2tmpSub[, 1] > 0.0002), ]
S.binomS2SubLo  <- S.binomS2tmpSub[which(S.binomS2tmpSub[, 1] < 0.0002), ]
S.binomS2Sub <- data.frame(cbind(-length(S.binomS2SubLo[, 1])/100, length(S.binomS2SubHi[, 1])/100))

S.binomS3tmpSub <- S.binomS3tmp[which(S.binomS3tmp[, 2] < 0.05), ]
S.binomS3SubHi  <- S.binomS3tmpSub[which(S.binomS3tmpSub[, 1] > 0.0002), ]
S.binomS3SubLo  <- S.binomS3tmpSub[which(S.binomS3tmpSub[, 1] < 0.0002), ]
S.binomS3Sub <- data.frame(cbind(-length(S.binomS3SubLo[, 1])/100, length(S.binomS3SubHi[, 1])/100))

S.binomS4tmpSub <- S.binomS4tmp[which(S.binomS4tmp[, 2] < 0.05), ]
S.binomS4SubHi  <- S.binomS4tmpSub[which(S.binomS4tmpSub[, 1] > 0.0002), ]
S.binomS4SubLo  <- S.binomS4tmpSub[which(S.binomS4tmpSub[, 1] < 0.0002), ]
S.binomS4Sub <- data.frame(cbind(-length(S.binomS4SubLo[1])/100, length(S.binomS4SubHi[, 1])/100))

S.binomS5tmpSub <- S.binomS5tmp[which(S.binomS5tmp[, 2] < 0.05), ]
S.binomS5SubHi  <- S.binomS5tmpSub[which(S.binomS5tmpSub[, 1] > 0.0002), ]
S.binomS5SubLo  <- S.binomS5tmpSub[which(S.binomS5tmpSub[, 1] < 0.0002), ]
S.binomS5Sub <- data.frame(cbind(-length(S.binomS5SubLo[1])/100, length(S.binomS5SubHi[, 1])/100))

S.binomS6tmpSub <- S.binomS6tmp[which(S.binomS6tmp[, 2] < 0.05), ]
S.binomS6SubHi  <- S.binomS6tmpSub[which(S.binomS6tmpSub[, 1] > 0.0002), ]
S.binomS6SubLo  <- S.binomS6tmpSub[which(S.binomS6tmpSub[, 1] < 0.0002), ]
S.binomS6Sub <- data.frame(cbind(-length(S.binomS6SubLo[, 1])/100, length(S.binomS6SubHi[, 1])/100))

S.binomS7tmpSub <- S.binomS7tmp[which(S.binomS7tmp[, 2] < 0.05), ]
S.binomS7SubHi  <- S.binomS7tmpSub[which(S.binomS7tmpSub[, 1] > 0.0002), ]
S.binomS7SubLo  <- S.binomS7tmpSub[which(S.binomS7tmpSub[, 1] < 0.0002), ]
S.binomS7Sub <- data.frame(cbind(-length(S.binomS7SubLo[, 1])/100, length(S.binomS7SubHi[, 1])/100))

S.binomS8tmpSub <- S.binomS8tmp[which(S.binomS8tmp[, 2] < 0.05), ]
S.binomS8SubHi  <- S.binomS8tmpSub[which(S.binomS8tmpSub[, 1] > 0.0002), ]
S.binomS8SubLo  <- S.binomS8tmpSub[which(S.binomS8tmpSub[, 1] < 0.0002), ]
S.binomS8Sub <- data.frame(cbind(-length(S.binomS8SubLo[1])/100, length(S.binomS8SubHi[, 1])/100))

S.binomS9tmpSub <- S.binomS9tmp[which(S.binomS9tmp[, 2] < 0.05), ]
S.binomS9SubHi  <- S.binomS9tmpSub[which(S.binomS9tmpSub[, 1] > 0.0002), ]
S.binomS9SubLo  <- S.binomS9tmpSub[which(S.binomS9tmpSub[, 1] < 0.0002), ]
S.binomS9Sub <- data.frame(cbind(-length(S.binomS9SubLo[, 1])/100, length(S.binomS9SubHi[, 1])/100))

S.binomS10tmpSub <- S.binomS10tmp[which(S.binomS10tmp[, 2] < 0.05), ]
S.binomS10SubHi  <- S.binomS10tmpSub[which(S.binomS10tmpSub[, 1] > 0.0002), ]
S.binomS10SubLo  <- S.binomS10tmpSub[which(S.binomS10tmpSub[, 1] < 0.0002), ]
S.binomS10Sub <- data.frame(cbind(-length(S.binomS10SubLo[, 1])/100, length(S.binomS10SubHi[, 1])/100))

S.binomS11tmpSub <- S.binomS11tmp[which(S.binomS11tmp[, 2] < 0.05), ]
S.binomS11SubHi  <- S.binomS11tmpSub[which(S.binomS11tmpSub[, 1] > 0.0002), ]
S.binomS11SubLo  <- S.binomS11tmpSub[which(S.binomS11tmpSub[, 1] < 0.0002), ]
S.binomS11Sub <- data.frame(cbind(-length(S.binomS11SubLo[, 1])/100, length(S.binomS11SubHi[, 1])/100))

S.binomS12tmpSub <- S.binomS12tmp[which(S.binomS12tmp[, 2] < 0.05), ]
S.binomS12SubHi  <- S.binomS12tmpSub[which(S.binomS12tmpSub[, 1] > 0.0002), ]
S.binomS12SubLo  <- S.binomS12tmpSub[which(S.binomS12tmpSub[, 1] < 0.0002), ]
S.binomS12Sub <- data.frame(cbind(-length(S.binomS12SubLo[, 1])/100, length(S.binomS12SubHi[, 1])/100))

S.binomS13tmpSub <- S.binomS13tmp[which(S.binomS13tmp[, 2] < 0.05), ]
S.binomS13SubHi  <- S.binomS13tmpSub[which(S.binomS13tmpSub[, 1] > 0.0002), ]
S.binomS13SubLo  <- S.binomS13tmpSub[which(S.binomS13tmpSub[, 1] < 0.0002), ]
S.binomS13Sub <- data.frame(cbind(-length(S.binomS13SubLo[, 1])/100, length(S.binomS13SubHi[, 1])/100))

S.binomS14tmpSub <- S.binomS14tmp[which(S.binomS14tmp[, 2] < 0.05), ]
S.binomS14SubHi  <- S.binomS14tmpSub[which(S.binomS14tmpSub[, 1] > 0.0002), ]
S.binomS14SubLo  <- S.binomS14tmpSub[which(S.binomS14tmpSub[, 1] < 0.0002), ]
S.binomS14Sub <- data.frame(cbind(-length(S.binomS14SubLo[, 1])/100, length(S.binomS14SubHi[, 1])/100))

S.binomS15tmpSub <- S.binomS15tmp[which(S.binomS15tmp[, 2] < 0.05), ]
S.binomS15SubHi  <- S.binomS15tmpSub[which(S.binomS15tmpSub[, 1] > 0.0002), ]
S.binomS15SubLo  <- S.binomS15tmpSub[which(S.binomS15tmpSub[, 1] < 0.0002), ]
S.binomS15Sub <- data.frame(cbind(-length(S.binomS15SubLo[, 1])/100, length(S.binomS15SubHi[, 1])/100))

Sdfbinom    <- data.frame(Significant = rbind(S.binomS1Sub, S.binomS2Sub, S.binomS3Sub, S.binomS4Sub, S.binomS5Sub, S.binomS6Sub, S.binomS7Sub, S.binomS8Sub, 
                                              S.binomS9Sub, S.binomS10Sub, S.binomS11Sub, S.binomS12Sub, S.binomS13Sub, S.binomS14Sub, S.binomS15Sub), 
                            Label = c("S1: D", "S2: M", "S3: I", "S4: DM", "S5: MD", "S6: DI", "S7: ID", "S8: MI",
                                    "S9: IM","S10: DMI","S11: DIM","S12: MDI","S13: MID","S14: IDM","S15: IMD"),
                            Mod = rep("Sub-region - Simulated", 15))

Simdf       <- data.frame(rbind(Sdfbinom, Sdf))
Simdf$Label <- factor(Simdf$Label, as.character(Simdf$Label))

Sim <- Simdf

#write.table(Simdf, "Simulated data - significant results.csv", sep = ",", col.names = TRUE)
#Sim <- read.table("Simulated data - significant results.csv", sep = ",", header = TRUE)
Sim$Label <- factor(Sim$Label, as.character(Sim$Label))

Int <- rbind(Sim)
DenU <- -Int[16, 1]
DenL <- -Int[16, 2]

Int[16, 1] <- DenL
Int[16, 2] <- DenU

DenMU <- -Int[18, 1]
DenML <- -Int[18, 2]

Int[18, 1] <- DenML
Int[18, 2] <- DenMU

Int$Label  <- c("S1: D", "S2: M", "S3: I", "S4: DM", "S5: MD", "S6: DI", "S7: ID", "S8: MI",
                "S9: IM", "S10: DMI", "S11: DIM", "S12: MDI", "S13: MID", "S14: IDM", "S15: IMD",
                "Low Population Density", "High Inward-Migration", "Low Population Density (adj. Inward-Migration)",
                "High Inward-Migration (adj. Population Density)")

Int$Label <- factor(Int$Label, unique(as.character(Int$Label)))

############################################
## Figure 2 of Paper - observed data only ##
############################################

x11()
ggplot(Int)+
  geom_bar(aes(Label, Significant.X1, fill = Mod, alpha = Mod, order = Mod), position = "dodge", stat = "identity") +
  geom_bar(aes(Label, Significant.X2, fill = Mod, alpha = Mod, order = Mod), position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0, size = 0.5) + scale_fill_manual(values = c("#08306b", "#a50f15", "#fc9272", "#6baed6")) + theme_bw() + scale_alpha_manual(values = c(0.8, 0.8, 0.8, 0.8)) +
  theme(axis.text = element_text(colour = "black", size = 12, family = "Times New Roman"), axis.title = element_text(colour = "black", size = 20, face = "bold", family = "Times New Roman"),
        plot.title = element_text(size = 16, face = "bold", family = "Times New Roman"), legend.position = c(0.2, 0.95), legend.text = element_text(size = 16, family = "Times New Roman"), axis.text.x = element_text(vjust = 0.5),
        legend.key = element_blank()) + 
  coord_cartesian(ylim = c(-15, 100)) + scale_x_discrete(name = "", labels = function(Method) str_wrap(c("S1: D", "S2: M", "S3: I", "S4: DM", "S5: MD", "S6: DI", "S7: ID", "S8: MI",
                                                                                                 "S9: IM", "S10: DMI", "S11: DIM", "S12: MDI", "S13: MID", "S14: IDM", "S15: IMD", 
                                                                                                "Low Population Density", "High Inward-Migration", "Low Population Density (adj. Inward-Migration)",
                                                                                                "High Inward-Migration (adj. Population Density)"), width=10))+
  scale_y_continuous(name = "Percentage of statistically significant results", breaks = c(-15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)) +
  theme(legend.title = element_blank(),legend.background = element_rect(fill = NA, colour = NA)) + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  annotate("text", label = "Positive coefficients\nabove zero", x = 1.5, y = 10, size = 5, colour = "black") +
  annotate("text", label = "Negative coefficients\nbelow zero", x = 1.5, y = -10, size = 5, colour = "black")
