% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_gibbsR.R
\name{cp_gibbsR}
\alias{cp_gibbsR}
\title{A Gibbs sampler using Rcpp}
\usage{
cp_gibbsR(
  length_of_chains = 10000,
  from_point = 1001,
  a = 1,
  b = 1,
  x_range = 10,
  mu_x = NA,
  mu_y = NA
)
}
\arguments{
\item{length_of_chains}{the number of samples}

\item{from_point}{how many number ahead are to abandon}

\item{a, b}{the shape of distribution Y = beta(x+a,n-x+b)}

\item{x_range}{= n the X values uses in distribution X = Binomial(n, y)}

\item{mu_x, mu_y}{the expectation of X and Y}
}
\value{
a random sample of size \code{n}
}
\description{
A Gibbs sampler using Rcpp
}
\examples{
\dontrun{
ZC <- cp_gibbsC()
}
}
