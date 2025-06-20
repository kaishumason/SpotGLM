% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_extension_adam.R
\name{glm_fisher}
\alias{glm_fisher}
\title{Fisher Information Matrix for Spot-GLM}
\usage{
glm_fisher(
  X,
  y,
  lambda,
  family,
  beta,
  good_beta,
  dispersion,
  weights = NULL,
  offset = rep(0, length(y))
)
}
\description{
Computes the Fisher Information matrix for a given GLM family and model.
}
\keyword{internal}
