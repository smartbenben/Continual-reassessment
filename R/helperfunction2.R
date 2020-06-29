#' helper function 2
#'
#' @description stop probability with efficacy
#' @keywords internal
#' export

continue <- function(p, b)
{ #  probability to continue sampling after each patient
  #  Finds probability of each patient outcome inside the upper bounds b
  #  p is prob of an event
  N <- length(b)                # number of patients
  accum <- NULL
  if(p < 0 || p > 1 || N < 2)return(NaN) # check validity of p, N
  new <- old <- c(1, rep(0, N)) # NB: Subscript is one more than # of events
  for(i in 1 : N)               # for each patient accrued...
  {
    for(j in 0 : b[i])          # for all valid events within the bound...
    {
      if(j == 0) new[1] <- old[1] * (1 - p )   # zero events
      else                      # recursive relationship
        new[j + 1] <- old[j + 1] * (1 - p)  +  old[j] * p
    }                           # end j loop: events
    old <- new                  # update
    accum <- c(accum, sum(new)) # accumulate probabilities
  }                             # end i loop: patients
  accum

}


##n1: evaluate after data from n1th patient becomes available
efcont <- function(p, e, f, n1)
  #   probability of completion without crossing either efficacy (e)
  #   or futility (f) boundaries
{ #   p is the probability of a positive treatment response
  N <- min(length(e), length(f)) # maximum sample size
  e <- e[1 : N]                  # trim both bounds to this length
  f <- f[1 : N]

  effstop <- pbinom(q=e[1], size=n1, prob=p, lower.tail=F)                   # prob of stopping for efficacy
  futstop <- pbinom(q=f[1]-1, size=n1, prob=p)

  accum <- NULL
  if(p < 0 || p > 1)return(NaN)
  old <- dbinom(x=0:N, size=n1-1, prob=p)
  new <- c(1, rep(0, N))  # NOTE: subscript # is one more than # of events
  for(i in 1 : N)                # for each patient accrued . . .
  {
    for(j in (f[i] : e[i]))       # valid range of events
    { ##probability of observing j events in the ith patient
      if(j == 0){
        new[1] <- old[1] * (1 - p )
      }else{
        new[j+1] <- old[j+1]*(1-p)+old[j]*p
      }                            # recursive relationship
    }                             # end j loop: number of events
    ##for the ith patient, p of continue to sample is 0 if n events is
    ##below the futility bound
    ##for patient i, any number of events less than the futility bound
    ##cannot get the trial continued
    if(f[i] > 0)new[1 : f[i]] <- 0  # events j = 0,1,...,f(i)-1 are futile
    accum <- c(accum, sum(new))  # probability of continuing to sample
    if(i > 1)                    # probability of stopping . . .
    { ##stop for fulitity if just on the futility bouddary with i-1
      ##patients and observe no event for the next patient
      ##for the previvous n, the number of responses is f[i-1]
      if(f[i] > f[i-1])futstop <- c(futstop,
                                    futstop[i-1] + (1-p) * old[f[i-1] + 1])

      if(f[i] == f[i - 1])futstop <- c(futstop, futstop[i-1])
      ##for the previous n, the number of responses is e[i]
      if(e[i] == e[i - 1])effstop <- c(effstop,
                                       effstop[i-1] + p * old[e[i] + 1])
      if(e[i] > e[i - 1])effstop <- c(effstop, effstop[i-1])
    }
    old <- new                   # update for next patient
  }                              # end i loop on patients
  efcont <- NULL                 # assemble output of function
  efcont$accum <- c(rep(1, n1-1),accum)
  efcont$futstop <- c(rep(0, n1-1), futstop)
  efcont$effstop <- c(rep(0, n1-1), effstop)
  return(efcont)
}


##############################################
##evaluate after data from n1th pt is available
eftcon <- function(p, e, f, b, n1)
{ # probability of continue sampling in continual reassessment
  # evaluate for efficacy (e), futility (f), and toxicity (b)
  # p(2,2) is probabilities of response[2,] and toxicity[,2]
  # find sample size
  N <- min(length(e), length(f), length(b))
  # truncate bounds to match minimum, if needed
  e <- e[1 : N]
  f <- f[1 : N]
  b <- b[1 : N]

  # all subscripts one larger
  old <- new <- matrix(0, (N + 1), (N + 1))
  # zero patients: no responses, no toxicities
  old[1, 1] <- 1

  accum <- NULL

  for (i in 1:n1) # for each patient accrued...
  {
    for(j in 0:i) # for each response
    {
      for (k in 0:b[i]) # for each AE
      { # probability distribution for next patient
        new[j+1, k+1] <- connew(p, old, j, k)
      } # end k loop on toxicities
    }# end j loop on treatment responses
    accum <- c(accum, sum(new))
    old <- new
    new[ , ] <- 0
  }
  ##probability of stopping at n1 th patient
  futstop <- sum(old[1:f[i],0:(b[i]+1)])
  effstop <- sum(old[(e[i]+2):(i+1),0:(b[i]+1)])
  if(f[i]==0)futstop <- 0
  if(e[i]==i)effstop <- 0
  toxstop <- 1-continue(p=sum(p[,2]), b=b)[1:n1]
  ##correct for eff/fut stopping
  tmp <- matrix(0, (N + 1), (N + 1))

  for(j in f[i]:e[i]) # for each response
  {
    for (k in 0:b[i]) # for each AE
    { # probability distribution for next patient
      tmp[j+1, k+1] <- old[j+1,k+1]
    } # end k loop on toxicities
  }# end j loop on treatment responses
  old <- tmp

  #accum[n1] <- sum(old[(f[i]+1):(e[i]+1),0:(b[i]+1)])
  accum[n1] <- sum(old)
  futstop <- c(rep(0,n1-1),futstop)
  effstop <- c(rep(0,n1-1),effstop)
  accum+effstop+futstop+toxstop

  for (i in (n1+1) : N)        # for each patient accrued...
  {
    for(j in f[i] : e[i]) # for each valid number of responses
    {
      for (k in 0 : b[i]) # for valid number of toxicities
      {                   # probability distribution for next patient
        new[j + 1, k + 1] <- connew(p, old, j, k)
      }                   # end k loop on toxicities
    }                     # end j loop on treatment responses

    accum <- c(accum, sum(new))
    ##no response, not toxic
    ## cross futility bound if
    ## the futility bound increase by 1
    ## the number of responses remain the same
    ## sum over all possible toxicity levels
    ## if previuosly below toxicity bound, next
    ## patient have/have no toxicity
    ##if previously on the toxicity bound, next
    ##patient have no toxicity
    if(f[i] > f[i-1])futstop <- c(futstop,
                                  futstop[i-1] +
                                    ifelse(0:b[i-1]<b[i], sum(p[1,]), p[1,1])%*%
                                    old[f[i-1]+1,1:(b[i-1]+1)] )

    if(f[i] == f[i-1])futstop <- c(futstop, futstop[i-1])

    if(e[i] == e[i-1])effstop <- c(effstop,
                                   effstop[i-1] +
                                     ifelse(0:b[i-1]<b[i], sum(p[2,]), p[2,1])%*%
                                     old[e[i-1]+1,1:(b[i-1]+1)] )

    if(e[i] > e[i-1])effstop <- c(effstop, effstop[i-1])

    if(b[i] == b[i-1])toxstop <- c( toxstop,
                                    toxstop[i-1] +
                                      sum(sum(p[,2])*old[(f[i-1]+1):(e[i-1]+1),(b[i-1]+1)]) )
    if(b[i] > b[i-1])toxstop <- c(toxstop, toxstop[i-1])

    old <- new            # update for next patient
    new[ , ] <- 0
  }                       # end i loop on patients
  eftcont <- NULL                 # assemble output of function
  eftcont$accum <- accum
  eftcont$futstop <- futstop
  eftcont$effstop <- effstop
  eftcont$toxstop <- toxstop
  return(eftcont)
}


#N <- 20                #  maximum sample size
#Ptox.rate <- 0.3       #  probability of toxic, adverse event
#Ptox.tail <- 0.9     #  tail area for safety for each patient
#Presp.rate <- 0.15     #  probability of favorable response
#Presp.etail <- 0.955   #  efficacy tail area for each patient
#Presp.ftail <- 0.05   #  futility tail area for each patient
#Bayes <- 25            #  number of Bayesian, virtual patients
#    = alpha + beta in beta prior distribution

#ppair <- matrix(c((1 - Ptox.rate) * (1 - Presp.rate), # paired responses
#                 (1 - Ptox.rate) * Presp.rate,
#                 Ptox.rate * (1 - Presp.rate),
#                Ptox.rate * Presp.rate),
#              2, 2)  # joint prob's of response/toxic

#colnames(ppair) <- c("safe", "toxic")
#row.names(ppair) <- c("no resp", "response")

#b <- bounds(Ptox.tail,  N, Ptox.rate)
#e <- bounds(Presp.etail, N, Presp.rate)
#f <- bounds(Presp.ftail, N, Presp.rate)
#p <- ppair

#n1 <-7
#test <- eftcon(p=ppair,e,f,b=1:N,n1=n1)
#test$futstop+test$effstop+test$toxstop+test$accum

#cbind(eftcon(p=ppair,e,f,b=1:N,n1=n1)$effstop,
#efcont(p=sum(p[2,]),e[n1:N],f[n1:N], n1=n1)$effstop)

#cbind(eftcon(p=ppair,e,f,b=1:N,n1=n1)$futstop,
#     efcont(p=sum(p[2,]),e[n1:N],f[n1:N], n1=n1)$futstop)
