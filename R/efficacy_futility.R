#' Efficacy and futility Stoping bound
#'
#' This function generats an efficacy stoping bound giving the parameters p and b.
#'
#' @param p0 null response rate
#' @param p1 target response rate
#' @param N maximum sample size
#' @param n1 the sample size to start interim monitoring
#' @param alpha alpha
#' @param beta beta
#' @importFrom foreach foreach %do%
#' @export
#' @examples two.tail(N=44, n1=1, p0=0.6, p1=0.8,alpha=0.05, beta=0.2)

two.tail <- function(N, n1, p0, p1,alpha,beta)
{
  e_tail <- rm.tails(N, p0, seq(0.5, 0.999, 0.001))

  f_tail <- rm.tails(N, p1, seq(0.001, 0.5, 0.001))

  tails <- expand.grid(e_tail=e_tail, f_tail=f_tail)

  ##design operating characteristics

  design_opers <- apply(as.matrix(tails), 1, function(x){
    e <- find.bound(x[1], N, p0)
    f <- find.bound(x[2], N, p1)
    last  <- min(N, max(which(diff(e) == 0)) + 1) #  last chance to cross e

    if (length(N)>1){
      f <- pmax(f, e[last] - last + (N))
      fun0 <- efcont(p=p0, e, f, n1=n1)
      err1 <- fun0$effstop[N]
      es0 <- ess(1-fun0$accum)
      if(err1< alpha){
        pwr <- exact.stop(nvec=N, e, f, ptrue=p0)$effstop[length(N)]
      }else{
        pwr <- es0 <- NA
      }
      c(err1, pwr, es0, x[1], x[2])
    }else{

    f <- pmax(f, e[last] - last + (1 : N))

    fun0 <- efcont(p=p0, e, f, n1=n1)
    err1 <- fun0$effstop[N]
    es0 <- ess(1-fun0$accum)
    if(err1< alpha){
      pwr <- efcont(p=p1, e, f, n1=n1)$effstop[N]
    }else{
      pwr <- es0 <- NA
    }
    c(err1, pwr, es0, x[1], x[2])
  }
})
  alpha_vec <- design_opers[1,]
  pwr_vec <- design_opers[2,]
  es_vec <- design_opers[3,]

  out <- data.frame("alpha"=alpha_vec,
                    "power"=pwr_vec,
                    "EN"=es_vec,
                    "e_tail"=design_opers[4,],
                    "f_tail"=design_opers[5,])

  out <- out[which(pwr_vec>1-beta),]

  if(is.null(dim(out))) warning('no proper design')

  out
}


#two.tail(N=44, n1=1, p0=0.6, p1=0.8, alpha=0.05,
#         beta=0.2)

#sel.2tail(N=29, n1=1, p0=0.1, p1=0.3, alpha=0.1,
#          beta=0.2, e_tail=seq(0.9,0.999,0.001),
#          f_tail=seq(0,0.1,0.001))
