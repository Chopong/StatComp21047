% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cp_gibbsC}
\alias{cp_gibbsC}
\title{A Gibbs sampler using Rcpp}
\usage{
cp_gibbsC(
  length_of_chains = 10000L,
  from_point = 1000L,
  a = 1L,
  b = 1L,
  x_range = 10L,
  mu_x = -1,
  mu_y = -1
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
