#'helper function 3
#' @keywords internal
#' export
#'
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
