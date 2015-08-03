\name{plot.btergmgof}
\alias{plot.btergmgof}
\alias{print.btergmgof}
\alias{summary.btergmgof}
\title{Plot or print btergmgof objects}
\description{Plot or print formatted goodness-of-fit statistics and degeneracy checks from btergmgof objects.}
\usage{
\method{plot}{btergmgof}(x, boxplot = TRUE, boxplot.mfrow = TRUE, 
    boxplot.dsp = TRUE, boxplot.esp = TRUE, boxplot.geodist = TRUE, 
    boxplot.degree = TRUE, boxplot.idegree = TRUE, 
    boxplot.odegree = TRUE, boxplot.kstar = TRUE, boxplot.istar = TRUE, 
    boxplot.ostar = TRUE, boxplot.dsp.max = NULL, 
    boxplot.esp.max = NULL, boxplot.geodist.max = NULL, 
    boxplot.degree.max = NULL, boxplot.idegree.max = NULL, 
    boxplot.odegree.max = NULL, boxplot.kstar.max = NULL, 
    boxplot.istar.max = NULL, boxplot.ostar.max = NULL, 
    boxplot.transform = function(x) x, boxplot.border = "darkgray", 
    boxplot.mean.col = "black", boxplot.median.col = "black", 
    boxplot.lwd = 0.8, boxplot.outline = FALSE, boxplot.ylab = "Frequency", 
    boxplot.main = NULL, boxplot.ylim = NULL, roc = TRUE, pr = TRUE, 
    rocpr.add = FALSE, rocpr.avg = c("none", "horizontal", "vertical", 
    "threshold"), rocpr.spread = c("boxplot", "stderror", "stddev"), 
    rocpr.lwd = 3, roc.main = NULL, roc.random = FALSE, roc.col = "#bd0017", 
    roc.random.col = "#bd001744", pr.main = NULL, pr.random = FALSE, 
    pr.col = "#5886be", pr.random.col = "#5886be44", pr.poly = 0, ...)

\method{print}{btergmgof}(x, classicgof = TRUE, 
    rocprgof = TRUE, degeneracy = TRUE, ...)

\method{summary}{btergmgof}(object, classicgof = TRUE, 
    rocprgof = TRUE, degeneracy = TRUE, ...)
}
\details{
These methods plot or print nicely formatted goodness-of-fit statistics from \code{btergmgof} objects to the R console. The typical workflow is to estimate a TERGM using the \link{btergm} function and save it as an object, then use the \link{gof} function to produce a \code{btergmgof} object, and finally use the functions and methods described here to show the output of this resulting object to assess the goodness of fit.

The \code{plot} method plots (1) boxplots of network statistics and draws the observed statistics as lines (classic statnet-like GOF boxplots) and (2) receiver operating characteristics (ROC) and precision recall (PR) curves.

The \code{print} and \code{summary} methods show (1) observed versus simulated network statistics (classic statnet-like GOF tables), (2) the area under the ROC curve (AUC) for each observed network, and (3) degeneracy checks by comparing global statistics of simulated versus observed networks.
}
\arguments{
\item{x}{ A \code{btergmgof} object created by calling the \link{gof.btergm} method. }
\item{boxplot}{ Create a plot where simulated statistics are summarized as boxplots and observed networks as lines (= classic statnet-style GOF plots)? }
\item{boxplot.mfrow}{ Should the classic GOF plots come out separately (\code{boxplot.mfrow = FALSE}), or should all statistics be aligned in a single diagram (\code{boxplot.mfrow = TRUE})? Returning the plots separately can be helpful if the output is redirected to a multipage PDF or TIFF file. }
\item{boxplot.dsp}{ Plot a GOF diagram for dyad-wise shared partners? }
\item{boxplot.esp}{ Plot a GOF diagram for edge-wise shared partners? }
\item{boxplot.geodist}{ Plot a GOF diagram for geodesic distances? }
\item{boxplot.degree}{ Plot a GOF diagram for the degree distribution? }
\item{boxplot.idegree}{ Plot a GOF diagram for the indegree distribution? }
\item{boxplot.odegree}{ Plot a GOF diagram for the outdegree distribution? }
\item{boxplot.kstar}{ Plot a GOF diagram for the kstar distribution? }
\item{boxplot.istar}{ Plot a GOF diagram for the instar distribution? }
\item{boxplot.ostar}{ Plot a GOF diagram for the outstar distribution? }
\item{boxplot.dsp.max}{ Upper bound of dyad-wise shared partners to report. }
\item{boxplot.esp.max}{ Upper bound of edge-wise shared partners to report. }
\item{boxplot.geodist.max}{ Upper bound of geodesic distances to report (excluding \code{Inf} values). }
\item{boxplot.degree.max}{ Upper bound of the degree statistic to report. }
\item{boxplot.idegree.max}{ Upper bound of the indegree statistic to report. }
\item{boxplot.odegree.max}{ Upper bound of the outdegree statistic to report. }
\item{boxplot.kstar.max}{ Upper bound of the kstar statistic to report. }
\item{boxplot.istar.max}{ Upper bound of the istar statistic to report. }
\item{boxplot.ostar.max}{ Upper bound of the ostar statistic to report. }
\item{boxplot.transform}{ A function which transforms the y values used for the boxplots. For example, if some of the values become very large and make the output illegible, \code{boxplot.transform = function(x) x^0.1} or a similar transformation of the values can be used. Note that logarithmic transformations often produce infinite values because \code{log(0) = -Inf}, so one should rather use something like \code{boxplot.transform = function(x) log1p} to avoid infinite values. }
\item{boxplot.border}{ Color of the borders of the boxplots. }
\item{boxplot.mean.col}{ Color of the mean of the observed network statistic. }
\item{boxplot.median.col}{ Color of the median of the observed network statistic. }
\item{boxplot.lwd}{ Scaling factor for the line width of the boxplots and of the line for the observed network statistics. }
\item{boxplot.outline}{ Print outliers in the boxplots? }
\item{boxplot.ylab}{ Label of the y-axis of a GOF plot. }
\item{boxplot.main}{ Main title of a GOF plot. Note: The same main title is used for all resulting diagrams. To obtain diagrams with varying main titles, plot each statistic separately using a new call of the \code{plot.btergm} method. }
\item{boxplot.ylim}{ Vertical limit of the boxplots. Only the maximum value must be provided because the minimum is always 0. Providing a vertical limit may be sensible if multiple models have to be compared side by side. }

\item{roc}{ Should the ROC curve(s) be plotted? }
\item{pr}{ Should the precision-recall curve(s) be plotted? }
\item{rocpr.add}{ Add the ROC and/or PR curve to an existing plot? }
\item{rocpr.avg}{ Averaging pattern for the ROC and PR curve(s) if multiple target time steps were used. Allowed values are \code{"none"} (plot all curves separately), \code{"horizontal"} (horizontal averaging), \code{"vertical"} (vertical averaging), and \code{"threshold"} (threshold (= cutoff) averaging). Note that while threshold averaging is always feasible, vertical and horizontal averaging are not well-defined if the graph cannot be represented as a function x->y and y->x, respectively. More information can be obtained from the help pages of the \pkg{ROCR} package, the functions of which are employed here. }
\item{rocpr.spread}{When multiple target time steps are used and curve averaging is enabled, the variation around the average curve can be visualized as standard error bars (\code{"stderror"}), standard deviation bars (\code{"stddev"}), or by using box plots (\code{"boxplot"}). Note that the function plotCI, which is used internally by the \pkg{ROCR} package to draw error bars, might raise a warning if the spread of the curves at certain positions is 0. More details can be found in the documentation of the \pkg{ROCR} package, the functions of which are employed here. }
\item{rocpr.lwd}{ Line width of the ROC and PR curve(s). }
\item{roc.main}{ Main title of the ROC plot. }
\item{roc.random}{ Draw an additional ROC curve for random graphs with the same tie probability (as a reference value)? }
\item{roc.col}{ Color of the ROC curve(s). By default, a dark red color is used. }
\item{roc.random.col}{ Color of the random graph ROC curve. }
\item{pr.main}{ Main title of the ROC plot. }
\item{pr.random}{ Draw an additional PR curve for random graphs with the same tie probability (as a reference value)? }
\item{pr.col}{ Color of the precision-recall curve(s). By default, a blue color is used. }
\item{pr.random.col}{ Color of the random graph PR curve. }
\item{pr.poly}{ If a value of \code{0} is set, nothing special happens. If a value of \code{1} is set, a straight line is fitted through the PR curve and displayed. Values between \code{2} and \code{9} fit higher-order polynomial curves through the PR curve and display the resulting curve. This argument allows to check whether the imputation of the first precision value in the PR curve yielded a reasonable result (in case the value had to be imputed). }

\item{classicgof}{ Print tables with t-tests comparing observed and simulated networks, as in statnet? }
\item{rocprgof}{ Print tables with the area under the curve (AUC) for the PR and ROC curves at each time step of the observed network(s)? }
\item{degeneracy}{ Print tables with t-tests comparing the global statistics of observed and simulated networks for each observed time step? This is essentially a degeneracy check. }
\item{object}{ A \code{btergmgof} object created by calling the \link{gof.btergm} method. }
\item{...}{ Further arguments. These arguments are passed on to the plot command. This facilitates customization of the plots. }
}
\seealso{
\link{btergm-package} \link{btergm} \link{simulate.btergm} \link[btergm]{gof}
}
