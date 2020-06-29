#' Group continul assessment
#'
#'@keywords internal
#'export


#cbind(e,f)

exact.stop <- function(nvec, e, f, ptrue){

  ##count outcomes from 0, that is, 0 responses or more
  new <- estop <- fstop <-
    matrix(0, nrow=length(nvec), ncol=max(nvec)+1)

  ##for each stage
  for(k in seq_along(nvec)){
    ##for each number of responses defined
    for(s in f[k]:e[k]){
      if(k==1){
        new[k, s+1] <- dbinom(s, nvec[1], ptrue)
      }else{
        new[k, s+1] <- old[(1+f[k-1]):(1+e[k-1])]%*%
          dbinom(s-f[k-1]:e[k-1],nvec[k]-nvec[k-1], ptrue)
      }
    }
    old <- new[k,]
  }

  for(k in seq_along(nvec)){
    for(s in (e[k]+1):nvec[k]){
      if(k==1){
        estop[k,s+1] <- dbinom(s, nvec[1], ptrue)
      }else{
        estop[k,s+1] <- new[k-1,(1+f[k-1]):(1+e[k-1])]%*%
          dbinom(s-f[k-1]:e[k-1],nvec[k]-nvec[k-1], ptrue)
      }
    }
    if(e[k]==nvec[k]){estop[k,] <- 0}
  }

  for(k in seq_along(nvec)){
    for(s in 0:(f[k]-1)){
      if(k==1){
        fstop[k,s+1] <- dbinom(s, nvec[1], ptrue)
      }else{
        fstop[k,s+1] <- new[k-1,(1+f[k-1]):(1+e[k-1])]%*%
          dbinom(s-f[k-1]:e[k-1],nvec[k]-nvec[k-1], ptrue)
      }
    }
    if(f[k]==0){fstop[k,] <- 0}
  }

  ef <- NULL
  ef$cont <- rowSums(new)
  ef$effstop <- cumsum(rowSums(estop))
  ef$futstop <- cumsum(rowSums(fstop))

  return(ef)
}

#ef <- exact.stop(nvec=seq(5,30,5),e, f, ptrue=0.1)

#ef$cont+ef$effstop+ef$futstop

#cbind(efcont(p=0.1, e, f, 1)$effstop, ef$effstop)
#cbind(efcont(p=0.1, e, f, 1)$futstop, ef$futstop)
#cbind(efcont(p=0.1, e, f, 1)$accum, ef$cont)
