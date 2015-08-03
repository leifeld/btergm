\name{preprocess}
\alias{preprocess}
\title{Preprocess lists of network matrices for use with btergm}
\description{Remove NAs and structural zeroes, create lagged covariates, and create memory terms.}
\usage{
preprocess(object, ..., lag = FALSE, covariate = FALSE, 
    memory = c("no", "autoregression", "stability", 
    "innovation", "loss"), na = NA, na.method = "fillmode", 
    structzero = -9, structzero.method = "remove", 
    verbose = FALSE)
}
\details{
For use with the \code{btergm} function, lists of network matrices or covariates can be preprocessed using the \code{preprocess} function. The first object that is provided as an argument is adjusted to all remaining objects that are handed over via the \code{...} argument. The preprocess function can deal with constant or time-varying \code{matrix}, \code{vector} and \code{data.frame} objects. The user can specify whether the object of interest is a dependent network or a covariate and whether NA and structural zero removal should take into account missing data at previous or following time steps (\code{lag = TRUE}).

The preprocessing procedure consists of four steps: (1) remove or impute missing data and/or structural zeroes in each object at each time step; (2) cross-sectional adjustment of matrix dimensions (e.g., if the object of interest has more observations than another object at the current time step, the observations are dropped from the object of interest); (3) backward-looking adjustment of matrix dimensions (if the object of interest is a dependent variable/network and the \code{lag = TRUE} argument is used); and (4) forward-looking adjustment of matrix dimensions (if the object of interest is a covariate and the \code{lag = TRUE} argument is used). The last two steps can be helpful for building lagged or delayed covariates which have to take into account missing data at other time steps as well.

For example, if 26 nodes are present during the first time step, 25 during the second, 23 during the third, and 25 during the fourth time step, the arguments \code{covariate = TRUE} and \code{lag = TRUE} return a list of three objects with 25, 23 and 25 objects, respectively, because the covariate is assumed to lag one step behind the other objects provided after the first object.

If the \code{memory} argument is specified (see below for details on the allowed values), the function will not return a lagged network, but will attempt to create a memory term. A memory term is a dyadic covariate which captures the temporal process.
}
\arguments{
\item{object}{ The object of interest that should be manipulated. The object can be a matrix, vector, data.frame, or a list of several objects of this kind. Usually, this is a list of network matrices (the dependent variable for the \code{btergm} function) or a list of covariates that needs to be adjusted given the missing observations at previous or future time steps. The name of the object is not allowed to match any of the argument names of this function (e.g., handing over an object with the name \code{memory} will cause an error message). }
\item{...}{ One or multiple other objects, separated by commas. These are used for comparison. If data are missing in these objects but not in the object of interest (given as the \code{object} argument), the missing data are also removed from the \code{object} of interest. The names of these objects are not allowed to match any of the argument names of this function (e.g., handing over an object with the name \code{memory} will cause an error message). }
\item{lag}{ Take into account missing or unobserved data at the previous or next time step? Use this argument in conjunction with the \code{covariate} argument. If \code{lag = TRUE} and \code{covariate = TRUE}, the first item in the list is adjusted to the dimensions of the second item, the second item is adjusted to the dimensions of the third item, and so forth, and the last item is omitted because it cannot be used as a covariate anyway. If \code{lag = TRUE} and \code{covariate = FALSE}, the first item is omitted, the second item is adjusted to the dimensions of the first item, the third item is adjusted to the dimensions of the second item, and so forth. If \code{lag = FALSE} (the \code{covariate} argument does not matter in this case), the first item is adjusted to the dimensions of the first item, the second item is adjusted to the dimensions of the second item etc. In all cases, missing data are first removed cross-sectionally across objects. }
\item{covariate}{ Is the object of interest a covariate or a dependent variable/network? This governs whether forward- or backward-looking missing data adjustment is used (when \code{lag = TRUE} is active). See the entry for \code{lag} for more details. }
\item{memory}{ Create a memory term which can be used as an edge covariate. Memory terms describe the temporal stability of the network. Several types of memory terms can be created. The default value \code{memory = "no"} will not create a memory term. \code{memory = "stability"} creates a term that captures dyadic stability between one time step and the next. \code{memory = "autoregression"} creates a lagged dependent network and is equivalent to \code{memory = "no"}. The result is a positive autoregression term, which is the lagged dependent network. This captures the stability of edges (but not of non-edges). \code{memory = "innovation"} measures the extent to which new edges are introduced over time ("edge innovation"). \code{memory = "loss"} captues the extent to which previously existing edges are no longer present at each current time step. The \code{memory} argument only works when \code{lag = TRUE} and \code{covariate = TRUE} are also set. }
\item{na}{ A vector of values that are identified as missing data. These values can be removed or replaced (see the \code{\link{handleMissings}} function). }
\item{na.method}{ What should be done with missing data? Valid values are \code{remove} and \code{fillmode}. See the \code{\link{handleMissings}} function for details. }
\item{structzero}{ A vector of values that are identified as structural zeroes. These values can be removed or replaced (see the \code{\link{handleMissings}} function). }
\item{structzero.method}{ What should be done with structural zeroes? Valid values are \code{remove} and \code{fillmode}. See the \code{\link{handleMissings}} function for details. }
\item{verbose}{ Report the amount of missing data that were removed or replaced in each object? }
}
\examples{
\dontrun{
# This example illustrates the usefulness of the preprocess function.

# first network: nodes a to j present
mat1 <- rbinom(100, 1, 0.1)
mat1 <- matrix(mat1, nrow = 10)  # has 10 nodes
rownames(mat1) <- letters[1:10]
colnames(mat1) <- letters[1:10]

# second network: nodes c to n present
mat2 <- rbinom(144, 1, 0.1)
mat2 <- matrix(mat2, nrow = 12)  # has 12 nodes
rownames(mat2) <- letters[3:14]
colnames(mat2) <- letters[3:14]

# third network: nodes a and d to k present
mat3 <- rbinom(81, 1, 0.1)
mat3 <- matrix(mat3, nrow = 9)  # has 9 nodes
rownames(mat3) <- letters[c(1, 4:11)]
colnames(mat3) <- letters[c(1, 4:11)]

# fourth network: same as second matrix
mat4 <- mat2

networks <- list(mat1, mat2, mat3, mat4)

# btergm without cross-temporal dependencies:
model.1 <- btergm(networks ~ edges + mutual)
summary(model.1)

# When cross-temporal dependencies are specified, the dimensions 
# of the matrices do not match. This would cause a problem for btergm:

\dontrun{
btergm(networks[2:4] ~ edges + mutual + edgecov(networks[1:3])) # ERROR!
}

# This is because the first network in the dependent network and the 
# first network in the lagged covariate are expected to have the same 
# dimensions (and also at the second and third time step, of course).

# Therefore, missing nodes in the covariate (here: {k, l, m, n} at t=1,  
# {a} at t=2, and {c, l, m, n} at t=3) must be removed from the 
# dependent network at t=2, t=3 and t=4 as well:

dep <- preprocess(networks, lag = TRUE, covariate = FALSE)

# This reduces the size of dep from 12 to 8 at t=2, from 9 to 8 at 
# t=3, and from 12 to 8 at t=4, and it removes the first network from 
# the list. Moreover, some nodes are present in the lagged covariate 
# but not in the dependent network (that is, at the next time step). 
# Therefore, node sets {a, b}, {c, l, m, n}, and {a} must be removed 
# from the lagged covariate at t=1, t=2, and t=3, respectively, to make 
# the dimensions compatible:

lag <- preprocess(networks, lag = TRUE, covariate = TRUE)

# To compare the dimensions of the original versus preprocessed 
# dependent networks and covariates, try the following code:

cbind(
    "original_dep" = lapply(networks[2:4], nrow), 
    "original_lag" = lapply(networks[1:3], nrow), 
    "new_dep" = lapply(dep, nrow), 
    "new_lag" = lapply(lag, nrow)
)

# The dependent networks were reduced from 12, 9 and 12 to 8, 8 and 
# 8 nodes, and the lagged networks were reduced from 10, 12 and 9 to 
# 8, 8 and 8 nodes, respectively. The lagged node sets are now 
# compatible. To see this:

cbind(rownames(dep[[1]]), rownames(lag[[1]]))
cbind(rownames(dep[[2]]), rownames(lag[[2]]))
cbind(rownames(dep[[3]]), rownames(lag[[3]]))

# Note, however, that the composition still changes within each list 
# across some of the time steps:

cbind(rownames(dep[[1]]), rownames(dep[[2]]), rownames(dep[[3]]))
cbind(rownames(lag[[1]]), rownames(lag[[2]]), rownames(lag[[3]]))

# We can now use the btergm function on the preprocessed lists:

model.2 <- btergm(dep ~ edges + mutual + edgecov(lag))
summary(model.2)

# The model can now be estimated because the current and lagged networks 
# have the same node sets at each time step. The disadvantage of this 
# approach is that some observations are lost. The advantage, however, 
# is that cross-temporal theories can be tested.

# However, since the node sets still differ across time steps, ROC and 
# PR curves cannot be estimated. This is true because a simulation from 
# nodes {c ... j} cannot be compared to a target network with nodes 
# {d ... k}. Therefore, the following command would compare the wrong 
# sets of nodes to estimate prediction performance:

\dontrun{
gof.2 <- gof(model.2, classicgof = FALSE, rocprgof = TRUE) # PROBLEM!
}

# To solve this problem, the most obvious approach is to estimate the 
# model at earlier time steps and compute the out-of-sample predictive 
# performance only for the last network:

model.3 <- btergm(dep[1:2] ~ edges + mutual + edgecov(lag[1:2]))
gof.3 <- gof(model.3, target = dep[[3]], formula = dep[[3]] ~ edges + 
    mutual + edgecov(lag[[3]]), classicgof = FALSE, rocprgof = TRUE)

# This models time steps 2 and 3 as a function of the lagged network 
# at time steps 1 and 2, uses the resulting coefficients to predict 
# the network at time step 4, and compares network 4 to simulations 
# based on the coefficients from the previous time steps and the 
# lagged network at the third time step. As the matrices within the 
# third list item have identical node sets, predictive performance 
# could be computed. The resulting ROC and PR curves can be plotted 
# as follows:

plot(gof.3, boxplot = FALSE, pr = FALSE, roc.random = TRUE, 
    ylab = "TPR/PPV", xlab = "FPR/TPR", roc.main = "ROC and PR")
plot(gof.3, boxplot = FALSE, roc = FALSE, pr.random = TRUE, 
    rocpr.add = TRUE)
legend("right", legend = c("ROC", "ROC random graph", "PR", 
    "PR random graph"), col = c("#bd0017", "#bd001744", "#5886be", 
    "#5886be44"), lty = 1, lwd = 3)

# For another example with real-world data, see vignette("knecht")
}
}
\seealso{
\link{btergm-package} \link{handleMissings} \link{adjust} \link{timecov}
}
