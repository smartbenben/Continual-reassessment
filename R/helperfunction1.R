#'Helper function 1
#'@importFrom stats qbinom
#'@importFrom ExomeDepth qbetabinom
#'
#'@keywords internal
#'export

# # # # # functions# # # # #  # #  #  #  # #
find.bound <- function(tail, N, p, Bayes = NaN)
{  # build general bounds: Bayes = alpha + beta in beta-binomial distribution
  # p is binomial probability or mean of beta-binomial
  if(length(N)>1){
  if(p < 0 || p > 1 || sum(N < 1)>0 ||
     tail < 0 || tail > 1) return(NaN)   # check for validity
  if(is.nan(Bayes))return(qbinom(tail, N, p))  # frequentist bounds
  if(Bayes <= 0) return(NaN)             # more error checking
  qbetabinom(tail, N, p, Bayes)      # beta-binomial tail areas

  }
  else if (length(N)==1){
    if(p < 0 || p > 1 || N < 1 ||
       tail < 0 || tail > 1) return(NaN)   # check for validity
    if(is.nan(Bayes))return(qbinom(tail, 1:N, p))  # frequentist bounds
    if(Bayes <= 0) return(NaN)
    qbetabinom(tail, N, p, Bayes)
  }else{
    return (NaN)
  }
}
##remove tails that generate the same bounds
rm.tails <- function(N, p, tails)
{
  mat <- sapply(tails, function(x)
  {
    find.bound(tail=x, N=N, p=p)
  })
  tails[!duplicated(t(mat))]
}

ess <- function(cont)
{  # Expected sample size for given probabiity of continued sampling
  N <- length(cont)
  ess <-  N * (1 - cont[N])    # full sample size needed
  ess <-  ess + sum((2 : N) * (cont[-1] - cont[-N]))
  ess + cont[1]                # only one observation needed
}

#####################################################################
##to accrual one more patients, the eff/tox outcomes can be ++,+-,-+,--
##need to calculate all 4 probs and add them up
##j,k index starts from 0
connew <- function(p, old, j, k)  #  support function, called by eftcon
{ #  continuation probability for next patient w/j responses and k toxicities
  ## if the previous cohort have j responses and k toxicities
  ## the next patient must have 0 tox, 0 response
  new <- p[1,1] * old[j + 1,  k + 1]            # no resp and not toxic
  if(j > 0)new <- new +  p[2,1] * old[j, k + 1] # response and not toxic
  if(k > 0)new <- new +  p[1,2] * old[j + 1, k] # no resp + toxic
  if(j > 0 && k > 0) new <- new + p[2,2] * old[j, k]  # resp and toxic
  new
}

#####################################################################
jointp <- function(resp, tox, lambda)
{ # joint probability of response and safety with odds ratio lambda
  a <- 1 - lambda                          # quadratic equation
  b <- 1 + (lambda - 1) * (tox + resp)
  c <- -lambda * resp * tox
  det <- b ^ 2 - 4 * a * c
  if(det < 0)return("No valid joint probability")
  if(a == 0)return(resp * tox)
  r1 <- (-b + sqrt(det)) / (2 * a)
  r2 <- (-b - sqrt(det)) / (2 * a)   # usually not the correct root
  return(r1)
}
#####################################################################
prob4 <- function(resp, tox, lambda)
{   # probabilities of four possible patient outcomes
  pout <- rep(0,4)
  # toxicity and positive response
  pout[1] <- jointp(resp, tox, lambda)
  # no toxicity, positive response
  pout[2] <- resp - pout[1]
  # toxicity, negative response
  pout[3] <- tox - pout[1]
  # no toxicity, no response
  pout[4] <- 1 - sum(pout[1:3])
  if(min(pout) <0)return("Invalid probabilities")
  out=cbind(c(pout[4],pout[2]),c(pout[3],pout[1]))
  colnames(out) <- c("safe", "toxic")
  row.names(out) <- c("no resp", "response")
  out
}
