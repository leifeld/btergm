#' Longitudinal classroom friendship network and behavior (Andrea Knecht)
#'
#' Longitudinal classroom friendship network and behavior (Andrea Knecht).
#'
#' The Knecht dataset contains the friendship network of 26 pupils in a Dutch
#' school class measured at four time points along with several demographic and
#' behavioral covariates like age, sex, ethnicity, religion, delinquency,
#' alcohol consumption, primary school co-attendance, and school advice. Some
#' of these covariates are constant while others vary over time.
#'
#' The full dataset (see Knecht 2006 and 2008) contains a large number of
#' classrooms while the dataset presented here is an excerpt based on one
#' single classroom. This excerpt was first used in a tutorial for the software
#' \pkg{Siena} and the corresponding R package \pkg{RSiena} (Snijders, Steglich
#' and van de Bunt 2010). The following description was largely copied from the
#' original data description provided on the homepage of the \pkg{Siena}
#' project (see below for the URL).
#'
#' The data were collected between September 2003 and June 2004 by Andrea
#' Knecht, supervised by Chris Baerveldt, at the Department of Sociology of the
#' University of Utrecht (NL). The entire study is reported in Knecht (2008).
#' The project was funded by the Netherlands Organisation for Scientific
#' Research NWO, grant 401-01-554. The 26 students were followed over their
#' first year at secondary school during which friendship networks as well as
#' other data were assessed at four time points at intervals of three months.
#' There were 17 girls and 9 boys in the class, aged 11--13 at the beginning of
#' the school year. Network data were assessed by asking students to indicate
#' up to twelve classmates which they considered good friends. Delinquency is
#' defined as a rounded average over four types of minor delinquency (stealing,
#' vandalism, graffiti, and fighting), measured in each of the four waves of
#' data collection. The five-point scale ranged from `never' to `more than 10
#' times', and the distribution is highly skewed. In a range of 1--5, the mode
#' was 1 at all four waves, the average rose over time from 1.4 to 2.0, and the
#' value 5 was never observed.
#'
#' @name knecht
#'
#' @aliases knecht friendship demographics primary delinquency alcohol advice
#'
#' @docType data
#'
#' @format
#' Note: the data have to be transformed before they can be used with
#' \pkg{btergm} and related packages (see examples below).
#' \describe{
#' \item{\code{friendship}}{is a list of adjacency matrices at four time
#' points, containing friendship nominations of the column node by the row
#' node. The following values are used: \code{0} = no, \code{1} = yes,
#' \code{NA} = missing, \code{10} = not a member of the classroom (structural
#' zero).}
#'
#' \item{\code{demographics}}{is a data frame with 26 rows (the pupils) and
#' four demographic variables about the pupils:} \itemize{ \item\code{sex}
#' (\code{1} = girl, \code{2} = boy) \item\code{age} (in years)
#' \item\code{ethnicity} (\code{1} = Dutch, \code{2} = other, \code{0} =
#' missing) \item\code{religion} (\code{1} = Christian, \code{2} =
#' non-religious, \code{3} = non-Christian religion, \code{0} = missing) }
#'
#' \item{\code{primary}}{is a 26 x 26 matrix indicating whether two pupils
#' attended the same primary school. \code{0} = no, \code{1} = yes.}
#'
#' \item{\code{delinquency}}{is a data frame with 26 rows (the pupils) and
#' four columns (the four time steps). It contains the rounded average of four
#' items (stealing, vandalizing, fighting, graffiti).  Categories: frequency
#' over last three months, \code{1} = never, \code{2} = once, \code{3} = 2--4
#' times, \code{4} = 5--10 times, \code{5} = more than 10 times; \code{0} =
#' missing.}
#'
#' \item{\code{alcohol}}{is a data frame with 26 rows (the pupils) and 3
#' columns (waves 2, 3, and 4). It contains data on alcohol use (\dQuote{How
#' often did you drink alcohol with friends in the last three months?}).
#' Categories: \code{1} = never, \code{2} = once, \code{3} = 2--4 times,
#' \code{4} = 5--10 times, \code{5} = more than 10 times; \code{0} = missing.}
#'
#' \item{\code{advice}}{is a data frame with one variable, \dQuote{school
#' advice}, the assessment given at the end of primary school about the school
#' capabilities of the pupil (\code{4} = low, \code{8} = high, \code{0} =
#' missing)}
#' }
#'
#' @references
#' Knecht, Andrea (2006): \emph{Networks and Actor Attributes in Early
#' Adolescence} [2003/04]. Utrecht, The Netherlands Research School ICS,
#' Department of Sociology, Utrecht University. (ICS-Codebook no. 61).
#'
#' Knecht, Andrea (2008): \emph{Friendship Selection and Friends' Influence.
#' Dynamics of Networks and Actor Attributes in Early Adolescence}. PhD
#' Dissertation, University of Utrecht. \url{
#' https://dspace.library.uu.nl/handle/1874/25950}.
#'
#' Knecht, Andrea, Tom A. B. Snijders, Chris Baerveldt, Christian E.  G.
#' Steglich, and Werner Raub (2010): Friendship and Delinquency: Selection and
#' Influence Processes in Early Adolescence.  \emph{Social Development} 19(3):
#' 494--514. \doi{10.1111/j.1467-9507.2009.00564.x}.
#'
#' Leifeld, Philip and Skyler J. Cranmer (2019): A Theoretical and Empirical
#' Comparison of the Temporal Exponential Random Graph Model and the Stochastic
#' Actor-Oriented Model. Network Science 7(1): 20--51.
#' \doi{10.1017/nws.2018.26}.
#'
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2018): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1--36.
#' \doi{10.18637/jss.v083.i06}.
#'
#' Snijders, Tom A. B., Christian E. G. Steglich, and Gerhard G. van de Bunt
#' (2010): Introduction to Actor-Based Models for Network Dynamics.
#' \emph{Social Networks} 32: 44--60. \doi{10.1016/j.socnet.2009.02.004}.
#'
#' Steglich, Christian E. G. and Andrea Knecht (2009): Die statistische Analyse
#' dynamischer Netzwerkdaten. In: Stegbauer, Christian and Roger Haeussling
#' (editors), \emph{Handbuch der Netzwerkforschung}, Wiesbaden: Verlag fuer
#' Sozialwissenschaften.
#'
#' @source The data were gathered by Andrea Knecht, as part of her PhD
#' research, building on methods developed by Chris Baerveldt, initiator and
#' supervisor of the project. The project is funded by the Netherlands
#' Organisation for Scientific Research NWO, grant 401-01-554, and is part of
#' the research program "Dynamics of Networks and Behavior" with principle
#' investigator Tom A. B. Snijders.
#'
#' \itemize{
#' \item Complete original data:
#' \url{https://easy.dans.knaw.nl/ui/datasets/id/easy-dataset:48665}
#' \item This excerpt in Siena format:
#' \url{http://www.stats.ox.ac.uk/~snijders/siena/klas12b.zip}
#' \item Siena dataset description:
#' \url{http://www.stats.ox.ac.uk/~snijders/siena/tutorial2010_data.htm}
#' }
#'
#' Permission to redistribute this dataset along with this package was granted
#' by Andrea Knecht on April 17, 2014. Questions about the data or the original
#' study should be directed to her.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' # ====================================================================
#' # The following example was taken from the JSS article about btergm
#' # that is referenced above (Leifeld, Cranmer and Desmarais 2018).
#' # ====================================================================
#'
#' require("texreg")
#' require("sna")
#' require("btergm")
#' require("RSiena")
#' data("knecht")
#'
#' # step 1: make sure the network matrices have node labels
#' for (i in 1:length(friendship)) {
#'   rownames(friendship[[i]]) <- 1:nrow(friendship[[i]])
#'   colnames(friendship[[i]]) <- 1:ncol(friendship[[i]])
#' }
#' rownames(primary) <- rownames(friendship[[1]])
#' colnames(primary) <- colnames(friendship[[1]])
#' sex <- demographics$sex
#' names(sex) <- 1:length(sex)
#'
#' # step 2: imputation of NAs and removal of absent nodes:
#' friendship <- handleMissings(friendship, na = 10, method = "remove")
#' friendship <- handleMissings(friendship, na = NA, method = "fillmode")
#'
#' # step 3: add nodal covariates to the networks
#' for (i in 1:length(friendship)) {
#'   s <- adjust(sex, friendship[[i]])
#'   friendship[[i]] <- network(friendship[[i]])
#'   friendship[[i]] <- set.vertex.attribute(friendship[[i]], "sex", s)
#'   idegsqrt <- sqrt(degree(friendship[[i]], cmode = "indegree"))
#'   friendship[[i]] <- set.vertex.attribute(friendship[[i]],
#'       "idegsqrt", idegsqrt)
#'   odegsqrt <- sqrt(degree(friendship[[i]], cmode = "outdegree"))
#'   friendship[[i]] <- set.vertex.attribute(friendship[[i]],
#'       "odegsqrt", odegsqrt)
#' }
#' sapply(friendship, network.size)
#'
#' # step 4: plot the networks
#' pdf("knecht.pdf")
#' par(mfrow = c(2, 2), mar = c(0, 0, 1, 0))
#' for (i in 1:length(friendship)) {
#'   plot(network(friendship[[i]]), main = paste("t =", i),
#'   usearrows = TRUE, edge.col = "grey50")
#' }
#' dev.off()
#'
#' # step 5: estimate TERGMS without and with temporal dependencies
#' model.2a <- btergm(friendship ~ edges + mutual + ttriple +
#'     transitiveties + ctriple + nodeicov("idegsqrt") +
#'     nodeicov("odegsqrt") + nodeocov("odegsqrt") +
#'     nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
#'     edgecov(primary), R = 100)
#'
#' model.2b <- btergm(friendship ~ edges + mutual + ttriple +
#'     transitiveties + ctriple + nodeicov("idegsqrt") +
#'     nodeicov("odegsqrt") + nodeocov("odegsqrt") +
#'     nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
#'     edgecov(primary) + delrecip + memory(type = "stability"),
#'     R = 100)
#'
#' # step 6: alternatively, estimate via MCMC-MLE:
#' model.2d <- mtergm(friendship ~ edges + mutual + ttriple +
#'     transitiveties + ctriple + nodeicov("idegsqrt") +
#'     nodeicov("odegsqrt") + nodeocov("odegsqrt") +
#'     nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
#'     edgecov(primary) + delrecip + memory(type = "stability"),
#'     control = control.ergm(MCMC.samplesize = 5000, MCMC.interval = 2000))
#'
#' # step 7: GOF assessment with out-of-sample prediction
#' # (note the commentaries and corrections at
#' #  https://doi.org/10.1017/nws.2022.7 and
#' #  https://doi.org/10.1017/nws.2022.6)
#' model.2e <- btergm(friendship[1:3] ~ edges + mutual + ttriple +
#'     transitiveties + ctriple + nodeicov("idegsqrt") +
#'     nodeicov("odegsqrt") + nodeocov("odegsqrt") +
#'     nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
#'     edgecov(primary) + delrecip + memory(type = "stability"),
#'     R = 100)
#'
#' gof.2e <- gof(model.2e, nsim = 100, target = friendship[[4]],
#'     formula = friendship[3:4] ~ edges + mutual + ttriple +
#'     transitiveties + ctriple + nodeicov("idegsqrt") +
#'     nodeicov("odegsqrt") + nodeocov("odegsqrt") +
#'     nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
#'     edgecov(primary) + delrecip + memory(type = "stability"),
#'     coef = coef(model.2b), statistics = c(esp, dsp, geodesic,
#'     deg, triad.undirected, rocpr))
#' pdf("gof-2e.pdf", width = 8, height = 6)
#' plot(gof.2e)
#' dev.off()
#' }
NULL