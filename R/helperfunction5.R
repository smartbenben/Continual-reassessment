#'helper function 5
#'@description  efficacy, toxicity and futility stoping bounds
#' @keywords internal
#' export
#'
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
