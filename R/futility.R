#' Futility Stoping bound
#'
#' This function generats an futility stoping bound giving the parameters p and b.
#'
#' @param p0 null response rate
#' @param p1 target response rate
#' @param N maximum sample size
#' @param n1 the sample size to start interim monitoring
#' @param alpha alpha
#' @param beta beta
#' @importFrom foreach foreach %do%
#' @export
#' @examples futility(N=44, n1=1, p0=0.6, p1=0.8,alpha=0.05, beta=0.2)

futility <- function(N, n1, p0, p1,alpha, beta)
{
  f_tail <- rm.tails(N, p=p1, seq(0.001, 0.5, 0.001))

  ##design operating characteristics
  design_opers <- sapply(f_tail, function(x){
    e <- rep(0,N)
    f <- find.bound(tail=x, N, p=p1)
    last  <- N #  last chance to cross e
    f <- pmax(f, e[last] - last + (1 : N))

    fun0 <- efcont(p=p0, e, f, n1=n1)
    err1 <- fun0$effstop[N]
    es0 <- ess(1-fun0$accum)
    if(err1< alpha){
    pwr <- efcont(p=p1, e, f, n1=n1)$effstop[N]
    }else{
     pwr <- es0 <- NA
    }
    c(err1, pwr, es0, x)
  })

  alpha_vec <- design_opers[1,]
  pwr_vec <- design_opers[2,]
  es_vec <- design_opers[3,]

  out <- data.frame("alpha"=alpha_vec,
                    "power"=pwr_vec,
                    "EN"=es_vec,
                    "f_tail"=design_opers[4,])

  out <- out[which(pwr_vec>1-beta),]

  if(is.null(dim(out))) warning('no proper design')

  out
}
