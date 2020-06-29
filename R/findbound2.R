# # # # # functions# # # # #  # #  #  #  # #
find.bound <- function(tail, N, p, Bayes = NaN)
{  # build general bounds: Bayes = alpha + beta in beta-binomial distribution
  # p is binomial probability or mean of beta-binomial
  if(length(N)>1){
    if(p < 0 || p > 1 || sum(N < 1)>0 ||
       tail < 0 || tail > 1) return(NaN)   # check for validity
    if(is.nan(Bayes))return(qbinom(tail, N, p))  # frequentist bounds
    if(Bayes <= 0) return(NaN)             # more error checking
    c(N[1:(length(N)-1)],qbetabinom(tail, N[length(N)], p, Bayes))      # beta-binomial tail areas

  }
  else if (length(N)==1){
    if(p < 0 || p > 1 || N < 1 ||
       tail < 0 || tail > 1) return(NaN)   # check for validity
    if(is.nan(Bayes))return(qbinom(tail, 1:N, p))  # frequentist bounds
    if(Bayes <= 0) return(NaN)
    c(1:(N-1), qbetabinom(tail, N, p, Bayes))
  }else{
    return (NaN)
  }
}
