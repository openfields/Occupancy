# load library
library(jagsUI)

# specify model file
cat("
model {
  # Define the priors for all parameters
  # gamma = recruitment;
  gamma ~ dgamma(0.001, 0.001)
  # omega[i] = survival probability of stage i;
  # p[i] = detection probability of i;
  # lambda[i]= initial population size of i
  for(i in 1:nStages){
    omega[i] ~ dbeta(1, 1)
    p[i] ~ dbeta(1, 1)
    lambda[i] ~ dgamma(0.001, 0.001)
  }
  # phi = transition probability from juveniles to adults
  phi ~ dbeta(1,1)
  
  # Define the stage transition matrix
  # TransMat(i_new,i_old) - Probability of transitioning to 
  # stage i_new from i_old conditional on survival
  TransMat[1,1] <- omega[1] * (1-phi)
  TransMat[2,1] <- omega[1] * phi
  TransMat[3,1] <- 1-omega[1]
  TransMat[1,2] <- 0
  TransMat[2,2] <- omega[2]
  TransMat[3,2] <- 1-omega[2]
  TransMat[1,3] <- 0
  TransMat[2,3] <- 0
  TransMat[3,3] <- 1
  
  # Define conditional Transition Matrix
  # i.e. probability of transitioning to a cell given that it
  # hasn't previously transitioned
  for(row in 1:(nStages+1)){
    for(col in 1:(nStages+1)){
      TransMatCond[row,col] <- TransMat[row,col] /
        ifelse(TransMat[row,col]==0,1,sum(TransMat[row:(nStages+1),col])
        )
    }}
  
  # Specify the model
  # Loop across all j sites
  for(j in 1:nSites) {
    # Initialize the model in Year 1
    N[1,j,1] ~ dpois(lambda[1])
    N[2,j,1] ~ dpois(lambda[2])
    N[3,j,1] <- 0
    # Estimate stage-specific detection for year one
    for(i in 1:nStages){
      for(k in 1:nReps){
        n[i,j,k,1] ~ dbin(p[i], N[i,j,1])
      }}
    # Specify the model for years 2+
    for(t in 2:nYears) {
      # Recruitment
      G[j,t-1] ~ dpois(gamma*N[2,j,t-1])
      # Specify the transitions between stages
      for(i in 1:(nStages+1)){
        Ncond[i,j,t-1,1] <- N[i,j,t-1]
        T[i,j,t-1,1] ~ dbin(TransMatCond[1,i],Ncond[i,j,t-1,1])
        for(a in 2:(nStages+1)){
          Ncond[i,j,t-1,a] <- N[i,j,t-1]-sum(T[i,j,t-1,1:(a-1)])
          T[i,j,t-1,a] ~ dbin(TransMatCond[a,i],Ncond[i,j,t-1,a]) 
        }
      }
      # Determine annual stage-specific abundances at locations
      for(i in 1:(nStages+1)){
        N[i,j,t] <- sum(T[i,j,t-1,1:nStages]) +
          ifelse(i==1,G[j,t-1],0)
      }
      # Estimate stage-specific detection for subsequent
      for(i in 1:nStages){
        for(k in 1:nReps){
          n[i,j,k,t] ~ dbin(p[i], N[i,j,t])
        }}
    }
  }
  # Sum up the number of individuals in all locations to
  # estimate annual, Ntot
  for (t in 1:nYears){
    for(i in 1:nStages){
      Ntot[i,t] <- sum(N[i,,t])
    }
  }
}
",fill=TRUE, file="model.txt")

# MCMC settings
nc <- 3
ni <- 50000
nb <- 10000
nt <- 1

# compile data
jags.data <- list()

# specify parameters to monitor
parameters <- c("")

# initial values
inits <- function(){}

# run JAGS 
model1 <- jags(jags.data, parameters, "model.txt", n.thin=nt, n.burnin=nb, n.chains=nc, n.iter=ni)

summary(model1)