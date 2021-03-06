% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efficacy_futility_toxicity.R
\name{three.tail}
\alias{three.tail}
\title{Efficacy,toxicity and futility Stoping bound}
\usage{
three.tail(N, n1, p0, p1, ptox_l, ptox_u, lambda = 1, w, alpha, beta, gam)
}
\arguments{
\item{N}{maximum sample size}

\item{n1}{the sample size to start interim monitoring}

\item{p0}{null response rate}

\item{p1}{target response rate}

\item{ptox_l}{a toxicity rate considered safe and we want to avoid stopping for safety under ptox_l}

\item{ptox_u}{a toxicity rate above than ptox_u is overly toxic, requring immiediate stop}

\item{lambda}{eff-tox correlation}

\item{w}{weight}

\item{alpha}{alpha}

\item{beta}{beta}

\item{gam}{probability of mistakenly stopping for toxicity when ptox<=ptox_l}
}
\description{
This function generats an efficacy stoping bound giving the parameters p and b.
}
\examples{
three.tail(N=44, n1=1, p0=0.6, p1=0.8, ptox_l=0.05, ptox_u=0.35,w=0.5,alpha=0.05, beta=0.20,gam=0.1)
}
