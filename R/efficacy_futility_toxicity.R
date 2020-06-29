#' Efficacy,toxicity and futility Stoping bound
#'
#' This function generats an efficacy stoping bound giving the parameters p and b.
#'
#' @param p0 null response rate
#' @param p1 target response rate
#' @param N maximum sample size
#' @param n1 the sample size to start interim monitoring
#' @param ptox_l a toxicity rate considered safe and we want to avoid stopping for safety under ptox_l
#' @param ptox_u a toxicity rate above than ptox_u is overly toxic, requring immiediate stop
#' @param lambda eff-tox correlation
#' @param gam probability of mistakenly stopping for toxicity when ptox<=ptox_l
#' @param w weight
#' @param alpha alpha
#' @param beta beta
#' @importFrom foreach foreach %do%
#' @export
#' @examples three.tail(N=44, n1=1, p0=0.6, p1=0.8, ptox_l=0.05, ptox_u=0.35,w=0.5,alpha=0.05, beta=0.20,gam=0.1)

three.tail <- function(N, n1, p0, p1,
                       ptox_l, ptox_u, lambda=1,
                        w,alpha, beta, gam)
{
  ##remove tails that generate the same bound
  t_tail <- rm.tails(N, ptox_u, seq(0.5, 0.999, 0.001))
  ## generate candidate tails
  out <- two.tail(N, n1, p0, p1, alpha, beta)

  i <- NULL
  tails <- foreach(i=1:length(t_tail),.combine="rbind")%do%{
    cbind(out[c("e_tail","f_tail")], rep(t_tail[i], nrow(out)))
  }
  ##safe,nonefficacious
  ppair00 <- prob4(tox=ptox_l, resp=p0, lambda)
  ##safe,efficacious
  ppair01 <- prob4(tox=ptox_l, resp=p1, lambda)
  ##toxic,nonefficacious
  ppair10 <- prob4(tox=ptox_u, resp=p0, lambda)
  ##toxic,efficacious
  ppair11 <- prob4(tox=ptox_u, resp=p1, lambda)

  ##design operating characteristics
  design_opers <- apply(as.matrix(tails), 1, function(x){
    e <- find.bound(tail=x[1], N=N, p=p0)

    f <- find.bound(tail=x[2], N=N, p=p1)
    last  <- min(N, max(which(diff(e) == 0)) + 1) #  last chance to cross e
    f <- pmax(f, e[last] - last + (1 : N))

    b <- find.bound(tail=x[3], N=N, p=ptox_u)

    fun00 <- eftcon(p=ppair00, e, f, b, n1=n1)
    ##expected sample size under h00
    es00 <- ess(1-fun00$accum)
    ##type 1 error rate under h00
    err1 <- fun00$effstop[N]

    if(err1<alpha){
      fun01 <- eftcon(p=ppair01, e, f, b, n1=n1)
      es01 <- ess(1-fun01$accum)
      err2 <- 1-fun01$effstop[N]
      err3 <- max(fun00$toxstop[N],fun01$toxstop[N])
    }else{
      err2 <- err3 <- es01 <- es10 <- es11 <- NA
    }

    if(err1<alpha & err2<beta & err3<gam){

      fun10 <- eftcon(p=ppair10, e, f, b, n1=n1)
      es10 <- ess(1-fun10$accum)

      fun11 <- eftcon(ppair11, e, f, b, n1=n1)
      es11 <- ess(1-fun11$accum)
    }else{
      es10 <- es11 <- NA
    }

    en <- t(w)%*%c(es00, es10, es11)

    c(err1, err2, err3, en, x, es00, es01, es10, es11)
  })

  alpha_vec <- design_opers[1,]
  pwr_vec <- 1-design_opers[2,]
  gam_vec <- design_opers[3,]
  es_vec <- design_opers[4,]

  out <- data.frame("alpha"=alpha_vec,
                    "power"=pwr_vec,
                    "gam"=gam_vec,
                    "EN"=es_vec,
                    "e_tail"=design_opers[5,],
                    "f_tail"=design_opers[6,],
                    "t_tail"=design_opers[7,],
                    "es00"=design_opers[8,],
                    "es01"=design_opers[9,],
                    "es10"=design_opers[10,],
                    "es11"=design_opers[11,])

  out <- out[which(pwr_vec>1-beta & gam_vec<gam), ]

  if(is.null(dim(out))) warning('no proper design')

  out
}


#out <- three.tail(N=41, n1=1, alpha=0.05, beta=0.20,
#gam=0.10,p0=0.2, p1=0.4, ptox_l=0.05, ptox_u=0.35)
#out[which.min(out$EN),]


