#'Helper function 6
#'@importFrom stats qbinom
#'@importFrom ExomeDepth qbetabinom
#'
#'@keywords internal
#'export

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
