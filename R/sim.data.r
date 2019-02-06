##' @title Simulated binomial data-set over the unit square
##' @description This binomial data-set was simulated by generating a zero-mean stationary Gaussian process over a 35 by 35 grid covering the unit square with Matern correlation sturcture. The parameters used in the simulation are \eqn{\sigma^2 = 0.7}, \eqn{\phi = 0.15}, \eqn{\kappa = 1.5} and \eqn{\tau^2 = 0}. The nugget effect was not included, hence \code{tau2 = 0}.
##' The variables are as follows:
##'
##' \itemize{
##'   \item \code{data} simulated values of the Gaussian process.
##'   \item \code{y} binomial observations.
##'   \item \code{units.m} binomial denominators.
##'   \item \code{geometry} X and Y coordinates.
##' }
##'
##' @docType data
##' @keywords datasets
##' @name sim.data
##' @usage data("sim.data")
##' @format A data frame with 1225 rows and 5 variables
NULL

