\name{tergm-terms}
\alias{tergm-terms}
\alias{btergm-terms}
\alias{tergmterms}
\alias{btergmterms}
\alias{ergm-terms}
\alias{ergmterms}
\title{Temporal dependencies for TERGMs}
\description{Network statistics that span multiple time points.}
\details{
In addition to the ERGM user terms that can be estimated within a single network (see \link[ergm]{ergm-terms}), the \pkg{btergm} package provides additional model terms that can be used within a formula. These additional statistics span multiple time periods and are therefore called "temporal dependencies." Examples include memory terms (i.e., positive autoregression, dyadic stability, edge innovation, or edge loss), delayed reciprocity or mutuality, and time covariates (i.e., functions of time or interactions with time):

\describe{
\item{\code{memory(type = "stability", lag = 1)}}{
Memory terms control for the impact of a previous network on the current network. Four different types of memory terms are available: positive autoregression (\code{type = "autoregression"}) checks whether previous ties are carried over to the current network; dyadic stability (\code{type = "stability"}) checks whether both edges and non-edges are stable between the previous and the current network; edge loss (\code{type = "loss"}) checks whether ties in the previous network have been dissolved and no longer exist in the current network; and edge innovation (\code{type = "innovation"}) checks whether previously unconnected nodes have the tendency to become tied in the current network. The \code{lag} argument accepts integer values and controls whether the comparison is made with the previous network (\code{lag = 1}), the pre-previous network (\code{lag = 2}) etc. Note that as \code{lag} increases, the number of time steps on the dependent variable decreases.
}
\item{\code{delrecip(mutuality = FALSE, lag = 1)}}{
The \code{delrecip} term checks for delayed reciprocity. For example, if node j is tied to node i at t = 1, does this lead to a reciprocation of that tie back from i to j at t = 2? If \code{mutuality = TRUE} is set, this extends not only to ties, but also non-ties. That is, if i is not tied to j at t = 1, will this lead to j not being tied to i at t = 2, in addition to positively reciprocal patterns over time? The \code{lag} argument controls the size of the temporal lag: with \code{lag = 1}, reciprocity over one consecutive time period is checked. Note that as \code{lag} increases, the number of time steps on the dependent variable decreases.
}
\item{\code{timecov(x = NULL, minimum = 1, maximum = NULL, transform = function(t) 1 + (0 * t) + (0 * t^2))}}{
The \code{timecov} model term checks for linear or non-linear time trends with regard to edge formation. Optionally, this can be combined with a covariate to create an interaction effect between a dyadic covariate and time in order to test whether the importance of a covariate increases or decreases over time. In the default case, time is modeled as being constant over time. By tweaking the \code{transform} function, arbitrary functional forms of time can be tested. For example, \code{transform = sqrt} (for a geometrically decreasing time effect), \code{transform = function(x) x^2} (for a geometrically increasing time effect), \code{transform = function(t) t} (for a linear time trend) or polynomial functional forms can be used. For time steps below the \code{minimum} value and above the \code{maximum} value, the time covariate is set to 0. These arguments can be used to create stepping effects, for example to use a value of 0 up to an external event and 1 from that event onwards in order to control for external influences. The model term works in a similar way as the separate \link{timecov} function.
}
}
}
\seealso{
\link{btergm-package} \link{btergm} \link{preprocess} \link{timecov} \link[ergm]{ergm-terms}
}
