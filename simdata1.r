# simulated Dail Madsen count data
# assume Gains are result of local population size, apparent survival function of habitat
# assume single stream segment, no confluences

data.fn <- function(nsites, surveys, alpha0, alpha1, cmin, cmax){
  y <- array(dim=c(nsites,surveys))
  X <- sort(runif(n=nsites,min=cmin,max=cmax))
  lam <- exp(alpha0 + alpha1*X)
  N <- rpois(n=nsites, lambda=lam)

  # generate observations for first year
  tmp <- N
  for(i in 1:surveys){
  y[,i] <- rbinom(n=nsites, size=tmp, prob=.4)
  tmp <- tmp-y[,i]
  }
  
return(list(y=y, N=N))
}
  
popdy <- function(N, surveys){
  # generate observations for second year with S & G
  len <- length(N)
  for(i in 1:len){
    S[i] <- sum(rbinom(n=N[i], size=1, prob=.8)) # constant apparent survival across sites
  }
  G <- rbinom(n=N, size=10, prob=.5) # make gains function of site level abundance
  N2 <- S + G
  
  y <- array(dim=c(len, surveys))
  tmp <- N
  for(i in 1:surveys){
    y[,i] <- rbinom(n=nsites, size=tmp, prob=.4)
    tmp <- tmp-y[,i]
  }
  
  
return(list(N2=N2, G=G, S=S, N=N, y=y))
}



set.seed(10)
sim.data <- data.fn(nsites=20,surveys=3,alpha0=1, alpha1=1.1, cmin=-3, cmax=3)
sim.data$y

# simulate gains and survival




cat("
    model{
      # define priors for parameters
      # gamma = recruitment
      gamma ~ dgamma(0.001, 0.001)

      # omega[i] = survival probability of stage i
      # p[i] = detection probability of i;
      # lambda[i] = initial population size of i
      for(i in 1:nstages){
        omega[i] ~ dbeta(1,1)
        p[i] ~ dbeta(1, 1)
        lambda[i] ~ dgamma(0.001, 0.001)
      }
      # phi = transition probability from juveniles to adults
      phi ~ dbeta(1, 1)

      # degine the stage transition matrix
      # TransMat(i_new, i_old) - probability of transitioning to stage i_new from i_old, conditional on survival
      # stage

    }
    
    ")







cat("
model {
  ### priors
    # initial abundance: negative bionomial parameters
    r ~ dunif(0,50)                              
    p ~ dunif(0,1)
    
    # survival
    for(i in 1:nSites){                          
    omega[i] <- 1/(1+exp(-omegaX[i]))    
    omegaX[i] <- alpha + beta1*poolRatioSt[i] + beta2*meanDepthSt[i] +  beta3*tempSt[i] 
    }
    # Jeffery's prior for survival coefficients
    alpha ~ dnorm(0, 0.37); beta1 ~ dnorm(0, 0.37);beta2 ~ dnorm(0, 0.37); beta3 ~ dnorm(0, 0.37) 
    # recruitment
    gamma ~ dunif(0,10)           
    # fixed detection probability based on three-pass depletion  
    Q ~ dunif(0.63, 0.65)                        
    
    ### Dail-Madsen model
    # loop across sites
    for(i in 1:nSites) {
    # Year 1 - initial abundance
    N[i,1] ~ dnegbin(p,r)
    # Detection model
    for(r in 1:nReps){
    y[i,1,r] ~ dbin(Q, N[i,1])
    }
    
    # Year 2
    for(t in 2:nYears) {
    # Estimate survival
    S[i,t-1] ~ dbin(omega[i], N[i,t-1]) 
    }
    }
    
    # Estimate gains: including two sites upstream & downstream
    # Due to locations of tributaries and study area boundaries, this section cannot be completely looped
    # resulting in a lengthy code
    for(t in 2:nYears) {
    # Jefferson Hill Brook
    G[1,t-1] ~ dpois(gamma*(N[1,t-1] + N[2,t-1] + N[3,t-1] + N[73,t-1] + N[74,t-1] + N[75,t-1] + N[76,t-1]))
    G[2,t-1] ~ dpois(gamma*(N[1,t-1] + N[2,t-1] + N[3,t-1] + N[4,t-1] + N[74,t-1] + N[75,t-1]))
    for(i in 3:8){
    G[i,t-1] ~ dpois(gamma*(N[i-2,t-1] + N[i-1,t-1] + N[i,t-1] + N[i+1,t-1] + N[i+2,t-1]))
    }
    G[9,t-1] ~ dpois(gamma*(N[7,t-1] + N[8,t-1] + N[9,t-1] + N[10,t-1] + N[11,t-1] + N[51,t-1]))
    G[10,t-1] ~ dpois(gamma*(N[8,t-1] + N[9,t-1] + N[10,t-1] + N[11,t-1] + N[12,t-1] + N[51,t-1] + N[52,t-1]))  
    G[11,t-1] ~ dpois(gamma*(N[9,t-1] + N[10,t-1] + N[11,t-1] + N[12,t-1] + N[13,t-1] + N[51,t-1] + N[52,t-1])) 
    G[12,t-1] ~ dpois(gamma*(N[10,t-1] + N[11,t-1] + N[12,t-1] + N[13,t-1] + N[14,t-1] + N[51,t-1]))
    for(i in 13:33){
    G[i,t-1] ~ dpois(gamma*(N[i-2,t-1] + N[i-1,t-1] + N[i,t-1] + N[i+1,t-1] + N[i+2,t-1]))
    }
    G[34,t-1] ~ dpois(gamma*(N[32,t-1] + N[33,t-1] + N[34,t-1] + N[35,t-1] + N[36,t-1] + N[59,t-1]))
    G[35,t-1] ~ dpois(gamma*(N[33,t-1] + N[34,t-1] + N[35,t-1] + N[36,t-1] + N[37,t-1] + N[59,t-1] + N[60,t-1]))  
    G[36,t-1] ~ dpois(gamma*(N[34,t-1] + N[35,t-1] + N[36,t-1] + N[37,t-1] + N[38,t-1] + N[59,t-1] + N[60,t-1])) 
    G[37,t-1] ~ dpois(gamma*(N[35,t-1] + N[36,t-1] + N[37,t-1] + N[38,t-1] + N[39,t-1] + N[59,t-1]))
    for(i in 38:46){
    G[i,t-1] ~ dpois(gamma*(N[i-2,t-1] + N[i-1,t-1] + N[i,t-1] + N[i+1,t-1] + N[i+2,t-1]))
    }
    G[47,t-1] ~ dpois(gamma*(N[45,t-1] + N[46,t-1] + N[47,t-1] + N[48,t-1] + N[49,t-1] + N[63,t-1]))
    G[48,t-1] ~ dpois(gamma*(N[46,t-1] + N[47,t-1] + N[48,t-1] + N[49,t-1] + N[50,t-1] + N[63,t-1] + N[64,t-1]))  
    G[49,t-1] ~ dpois(gamma*(N[47,t-1] + N[48,t-1] + N[49,t-1] + N[50,t-1] + N[63,t-1] + N[64,t-1])) 
    G[50,t-1] ~ dpois(gamma*(N[48,t-1] + N[49,t-1] + N[50,t-1] + N[38,t-1] + N[63,t-1]))    
    G[51,t-1] ~ dpois(gamma*(N[9,t-1] + N[10,t-1] + N[11,t-1] + N[12,t-1] + N[51,t-1] + N[52,t-1] + N[53,t-1]))
    G[52,t-1] ~ dpois(gamma*(N[10,t-1] + N[11,t-1] + N[51,t-1] + N[52,t-1] + N[53,t-1] + N[54,t-1]))
    
    for(i in 53:56){
    G[i,t-1] ~ dpois(gamma*(N[i-2,t-1] + N[i-1,t-1] + N[i,t-1] + N[i+1,t-1] + N[i+2,t-1]))
    }
    
    G[57,t-1] ~ dpois(gamma*(N[55,t-1] + N[56,t-1] + N[57,t-1] + N[58,t-1])) 
    G[58,t-1] ~ dpois(gamma*(N[56,t-1] + N[57,t-1] + N[58,t-1]))    
    G[59,t-1] ~ dpois(gamma*(N[34,t-1] + N[35,t-1] + N[36,t-1] + N[37,t-1] + N[59,t-1] + N[60,t-1] + N[61,t-1]))
    G[60,t-1] ~ dpois(gamma*(N[35,t-1] + N[36,t-1] + N[59,t-1] + N[60,t-1] + N[61,t-1] + N[62,t-1]))  
    G[61,t-1] ~ dpois(gamma*(N[59,t-1] + N[60,t-1] + N[61,t-1] + N[62,t-1])) 
    G[62,t-1] ~ dpois(gamma*(N[60,t-1] + N[61,t-1] + N[62,t-1]))    
    G[63,t-1] ~ dpois(gamma*(N[47,t-1] + N[48,t-1] + N[49,t-1] + N[50,t-1] + N[63,t-1] + N[64,t-1] + N[65,t-1]))
    G[64,t-1] ~ dpois(gamma*(N[48,t-1] + N[49,t-1] + N[63,t-1] + N[64,t-1] + N[65,t-1] + N[66,t-1]))  
    G[65,t-1] ~ dpois(gamma*(N[63,t-1] + N[64,t-1] + N[65,t-1] + N[66,t-1])) 
    G[66,t-1] ~ dpois(gamma*(N[64,t-1] + N[65,t-1] + N[66,t-1]))    
    
    # Spruce Brook
    G[67,t-1] ~ dpois(gamma*(N[67,t-1] + N[68,t-1] + N[69,t-1]))
    G[68,t-1] ~ dpois(gamma*(N[67,t-1] + N[68,t-1] + N[69,t-1] + N[70,t-1]))
    for(i in 69:72){
    G[i,t-1] ~ dpois(gamma*(N[i-2,t-1] + N[i-1,t-1] + N[i,t-1] + N[i+1,t-1] + N[i+2,t-1]))
    }
    G[73,t-1] ~ dpois(gamma*(N[1,t-1] + N[71,t-1] + N[72,t-1] + N[73,t-1] + N[74,t-1] + N[75,t-1]))
    G[74,t-1] ~ dpois(gamma*(N[1,t-1] + N[2,t-1] + N[72,t-1] + N[73,t-1] + N[74,t-1] + N[75,t-1] + N[76,t-1]))
    G[75,t-1] ~ dpois(gamma*(N[1,t-1] + N[2,t-1] + N[73,t-1] + N[74,t-1] + N[75,t-1] + N[76,t-1] + N[77,t-1]))
    G[76,t-1] ~ dpois(gamma*(N[1,t-1] + N[74,t-1] + N[75,t-1] + N[76,t-1] + N[77,t-1] + N[78,t-1]))
    
    for(i in 77:144){
    G[i,t-1] ~ dpois(gamma*(N[i-2,t-1] + N[i-1,t-1] + N[i,t-1] + N[i+1,t-1] + N[i+2,t-1]))
    }
    
    G[145,t-1] ~ dpois(gamma*(N[143,t-1] + N[144,t-1] + N[145,t-1] + N[146,t-1] + N[147,t-1] + N[150,t-1]))
    G[146,t-1] ~ dpois(gamma*(N[144,t-1] + N[145,t-1] + N[146,t-1] + N[147,t-1] + N[148,t-1] + N[150,t-1] + N[151,t-1]))
    G[147,t-1] ~ dpois(gamma*(N[145,t-1] + N[146,t-1] + N[147,t-1] + N[148,t-1] + N[149,t-1] + N[150,t-1] + N[151,t-1]))
    G[148,t-1] ~ dpois(gamma*(N[146,t-1] + N[147,t-1] + N[148,t-1] + N[149,t-1] + N[150,t-1]))
    G[149,t-1] ~ dpois(gamma*(N[147,t-1] + N[148,t-1] + N[149,t-1]))
    G[150,t-1] ~ dpois(gamma*(N[145,t-1] + N[146,t-1] + N[147,t-1] + N[148,t-1] + N[150,t-1] + N[151,t-1] + N[152,t-1]))
    G[151,t-1] ~ dpois(gamma*(N[146,t-1] + N[147,t-1] + N[150,t-1] + N[151,t-1] + N[152,t-1]))
    G[152,t-1] ~ dpois(gamma*(N[150,t-1] + N[151,t-1] + N[152,t-1]))  
    } 
    #Sum survival and gain to get total N at each site i in each year t
    for(i in 1:nSites) {
    for(t in 2:nYears){
    N[i,t] <- S[i,t-1] + G[i,t-1]  
    #Detection model
    for(r in 1:nReps){
    y[i,t,r] ~ dbin(Q, N[i,t])
    }    
    }
    } 
    }
    
    ",fill=TRUE,file="mod1.txt")
