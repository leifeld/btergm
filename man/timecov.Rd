\name{timecov}
\alias{timecov}
\title{Create interaction terms between covariates and (transformations of) time for TERGMs}
\description{Create interaction terms between covariates and (transformations of) time for TERGMs.}
\usage{
timecov(covariate, minimum = 1, maximum = length(covariate), 
    transform = function(t) 1 + (0 * t) + (0 * t^2), 
    onlytime = FALSE)
}
\details{
The \code{timecov} function takes a list of matrices or networks (a varying dyadic covariate) and creates an interaction term with time. For example, if the rules of network formation are expected to change after some time steps (say, after the third out of six time steps), one can model the effect of the covariate for the first three time points and for the remaining three time points using separate model terms. To achieve this, the covariate matrix at each time step is multiplied by zeros or ones, depending on whether the time step should be incorporated in that model term. In this situation, one would create two dyadic covariates, one where the first three time periods are present and the remaining three are set to zero, and the other one where the first three time periods are set to zero and the remaining ones are present.

Another usage scenario is that time per se may have a polynomial effect (of any shape) on the probability of forming a tie. In this case, the \code{timecov} function can be used to create a covariate list of matrices where the entries of the matrix correspond to polynomial functions of time. For example, if tie formation becomes increasingly likely (as expressed by a linear relationship), a list of matrices with linearly increasing entries over time can be created. Such time effects can also be interacted with other covariates. For example, the further time progresses, the more (or less) the property captured by the covariate becomes important for tie formation.
}
\arguments{
\item{covariate}{ A list of matrices or network objects (a time-varying dyadic covariate). Note that nodal covariates can be manually converted into dyadic covariates to use them with this function. }
\item{minimum}{ The first matrix or network in the covariate list which should be non-zero. All matrices or networks before this time step are filled with zeros. }
\item{maximum}{ The last matrix or network in the covariate list which should be non-zero. All matrices or networks after this time step are filled with zeros. }
\item{transform}{ A function that transforms the time axis. By default, all time steps are weighted equally. For example, if there are five consecutive matrices, all of them are equally important, i.e., time is represented by weights of \code{1-1-1-1-1}. A linearly increasing importance of time can be achieved by using a multiplication with time, for example \code{function(t) = 2 * t} leads to a multiplication of all values in the covariates by the values \code{2-4-6-8-10}, respectively by time step. Quadratic or higher-order polynomials are possible using more complex transformations of time. }
\item{onlytime}{ If \code{onlytime = TRUE} is set, the function returns a list of matrices where all entries within a matrix are identical and reflect the transformation of time. This is useful when a function of time per se should be included in the model as a covariate, e.g., when the tie formation probability increases or decreases over time. If \code{onlytime = FALSE} is set, the aforementioned list of matrices is interacted with the covariate list, i.e., each value in the matrix at time step \code{t} is multiplied by the transformation of time at that time step. }
}
\seealso{
\link{btergm-package} \link{btergm} \link{preprocess}
}
