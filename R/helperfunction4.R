#'helper function 4
#'@description efficacy and futility stoping bounds
#'@importFrom stats pbinom dbinom
#' @keywords internal
#' export
#'
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
