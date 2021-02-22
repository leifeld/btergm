#' German Toxic Chemicals Policy Network in the 1980s (Volker Schneider)
#' 
#' German Toxic Chemicals Policy Network in the 1980s (Volker Schneider).
#' 
#' The chemnet dataset contains network and attribute data and for the 30 most
#' influential political actors with regard to toxic chemicals regulation in
#' Germany in 1983/1984. While the original dataset contains up to 47 actors,
#' this dataset contains the "complete influence core" of mutually relevant
#' actors. The data are cross-sectional. There are no missing data; the
#' response rate was 100 percent. Volker Schneider (University of Konstanz)
#' collected this dataset for his dissertation (Schneider 1988). The dataset
#' was later re-used for a journal publication on information exchange in
#' policy networks (Leifeld and Schneider 2012).
#' 
#' The chemnet dataset contains network relations on political/strategic and
#' technical/scientific information exchange, influence attribution, and
#' membership in policy committees/forums, as well as nodal attributes on the
#' actor type and opinions about the six most salient issues related to the
#' political process that was leading to a new chemicals law at the time being.
#' 
#' @name chemnet
#' 
#' @aliases chemnet pol scito scifrom infrep committee types intpos
#' 
#' @docType data
#' 
#' @format
#' \describe{
#' \item{\code{pol}}{is a directed 30 x 30 adjancency matrix indicating which
#' row actor sends political/strategic information to which column actor.
#' \code{1} indicates an information exchange tie, and \code{0} indicates the
#' absence of a network tie.}
#' 
#' \item{\code{scito}}{is a directed 30 x 30 adjacency matrix indicating which
#' row actor sends technical/scientific information to which column actor.
#' \code{1} indicates an information exchange tie, and \code{0} indicates the
#' absence of a network tie. In contrast to political/strategic information
#' exchange, two separate survey questions were asked about
#' technical/scientific information exchange: sending information, and
#' receiving information. The two matrices contain the same relation but one
#' time from the sender's perspective and one time from the receiver's
#' perspective. By combining the two matrices, one can create a "confirmed"
#' technical/scientific information exchange relation. The \code{scito} matrix
#' contains ties from the sender's perspective.}
#' 
#' \item{\code{scifrom}}{is a directed 30 x 30 adjacency matrix indicating
#' which row actor receives technical/scientific information from which column
#' actor. \code{1} indicates an information exchange tie, and \code{0}
#' indicates the absence of a network tie. In contrast to political/strategic
#' information exchange, two separate survey questions were asked about
#' technical/scientific information exchange: sending information, and
#' receiving information. The two matrices contain the same relation but one
#' time from the sender's perspective and one time from the receiver's
#' perspective. By combining the two matrices, one can create a "confirmed"
#' technical/scientific information exchange relation. The \code{scifrom}
#' matrix contains ties from the receiver's perspective.}
#' 
#' \item{\code{infrep}}{is a directed 30 x 30 adjancency matrix indicating
#' which row actor deems which column actor "particularly influential".
#' \code{1} indicates such a tie, and \code{0} indicates the absence of an
#' influence attribution tie.}
#' 
#' \item{\code{committee}}{is a 30 x 20 two-mode (bipartite) network matrix
#' indicating which row actor is a member of which policy committee/forum (as
#' indicated by the column labels). \code{1} indicates a membership tie, and
#' \code{0} indicates non-membership.}
#' 
#' \item{\code{types}}{is a one-column data.frame where the \code{type}
#' variable contains the actor type of each node. The following values are
#' possible:}
#' \itemize{
#' \item \code{gov} (government actor, e.g., a federal ministry)
#' \item \code{ig} (interest group)
#' \item \code{io} (international organization)
#' \item \code{par} (political party)
#' \item \code{sci} (scientific organization)
#' }
#' 
#' \item{\code{intpos}}{is a 30 x 6 matrix containing the interest positions of
#' the 30 political actors on the six most salient political issues related to a
#' pending new chemicals law. \code{-1} indicates a negative stance, i.e., the
#' actor rejects the proposal; \code{1} indicates a positive stance, i.e., the
#' actor supports the proposal; and \code{0} indicates a neutral or absent
#' opinion.}
#' }
#' 
#' @references
#' Leifeld, Philip and Volker Schneider (2012): Information Exchange in Policy
#' Networks. \emph{American Journal of Political Science} 53(3): 731--744.
#' \doi{10.1111/j.1540-5907.2011.00580.x}.
#' 
#' Schneider, Volker (1988): \emph{Politiknetzwerke der Chemikalienkontrolle.
#' Eine Analyse einer transnationalen Politikentwicklung}. Walter de Gruyter:
#' Berlin/New York.
#' 
#' Schneider, Volker and Philip Leifeld (2009): Ueberzeugungssysteme,
#' Diskursnetzwerke und politische Kommunikation: Ein zweiter Blick auf die
#' deutsche Chemikalienkontrolle der 1980er Jahre. In: Volker Schneider, Frank
#' Janning, Philip Leifeld and Thomas Malang (editors): \emph{Politiknetzwerke.
#' Modelle, Anwendungen und Visualisierungen}. Pages 139--158. Wiesbaden: VS
#' Verlag fuer Sozialwissenschaften.
#' 
#' @source The data were collected using paper-based questionnaires. The
#' questionnaires were administered in personal interviews (PAPI).  Further
#' information, including the actual survey, data on additional actors, the
#' full names of the policy committees/forums, and the full list of
#' unabbreviated actor names can be found online at
#' \url{http://hdl.handle.net/1902.1/17004} in the replication archive of
#' Leifeld and Schneider (2012).
#' 
#' \itemize{
#' \item Replication archive: \url{http://hdl.handle.net/1902.1/17004}
#' \item AJPS publication: \doi{10.1111/j.1540-5907.2011.00580.x}
#' }
#' 
#' The dataset is publicly available. Questions about the data or the original
#' study should be directed to Volker Schneider
#' <volker.schneider@uni-konstanz.de>, the author of the original study and
#' person who collected the data.
#' 
#' @keywords datasets
#' 
#' @examples
#' \dontrun{
#' # Replication code for Leifeld and Schneider (2012), AJPS.
#' # Note that the estimates can only be reproduced approximately 
#' # due to internal changes in the statnet package.
#' 
#' # preparatory steps
#' library("statnet")
#' library("xergm")
#' library("texreg")
#' seed <- 12345
#' set.seed(seed)
#' data("chemnet")
#' 
#' # create confirmed network relation
#' sci <- scito * t(scifrom)  # equation 1 in the AJPS paper
#' prefsim <- dist(intpos, method = "euclidean")  # equation 2
#' prefsim <- max(prefsim) - prefsim  # equation 3
#' prefsim <- as.matrix(prefsim)
#' committee <- committee %*% t(committee)  # equation 4
#' diag(committee) <- 0 # the diagonal has no meaning
#' types <- types[, 1]  # convert to vector
#' 
#' # create network objects and store attributes
#' nw.pol <- network(pol) # political/stratgic information exchange
#' set.vertex.attribute(nw.pol, "orgtype", types)
#' set.vertex.attribute(nw.pol, "betweenness", 
#'     betweenness(nw.pol)) # centrality
#' 
#' nw.sci <- network(sci) # technical/scientific information exchange
#' set.vertex.attribute(nw.sci, "orgtype", types)
#' set.vertex.attribute(nw.sci, "betweenness", 
#'     betweenness(nw.sci)) # centrality
#' 
#' # ERGM: model 1 in the AJPS paper; only preference similarity
#' model1 <- ergm(nw.pol ~ edges + edgecov(prefsim), 
#'     control = control.ergm(seed = seed))
#' summary(model1)
#' 
#' # ERGM: model 2 in the AJPS paper; complete model
#' model2 <- ergm(nw.pol ~ 
#'     edges + 
#'     edgecov(prefsim) + 
#'     mutual + 
#'     nodemix("orgtype", base = -7) + 
#'     nodeifactor("orgtype", base = -1) + 
#'     nodeofactor("orgtype", base = -5) + 
#'     edgecov(committee) + 
#'     edgecov(nw.sci) + 
#'     edgecov(infrep) + 
#'     gwesp(0.1, fixed = TRUE) + 
#'     gwdsp(0.1, fixed = TRUE), 
#'     control = control.ergm(seed = seed)
#' )
#' summary(model2)
#' 
#' # ERGM: model 3 in the AJPS paper; only preference similarity
#' model3 <- ergm(nw.sci ~ edges + edgecov(prefsim), 
#'     control = control.ergm(seed = seed))
#' summary(model3)
#' 
#' # ERGM: model 4 in the AJPS paper; complete model
#' model4 <- ergm(nw.sci ~ 
#'     edges + 
#'     edgecov(prefsim) + 
#'     mutual + 
#'     nodemix("orgtype", base = -7) + 
#'     nodeifactor("orgtype", base = -1) + 
#'     nodeofactor("orgtype", base = -5) + 
#'     edgecov(committee) + 
#'     edgecov(nw.pol) + 
#'     edgecov(infrep) + 
#'     gwesp(0.1, fixed = TRUE) + 
#'     gwdsp(0.1, fixed = TRUE), 
#'     control = control.ergm(seed = seed)
#' )
#' summary(model4)
#' 
#' # regression table using the texreg package
#' screenreg(list(model1, model2, model3, model4))
#' 
#' # goodness of fit using the btergm package
#' gof2 <- gof(model2, roc = FALSE, pr = FALSE)
#' gof2  # print gof output
#' plot(gof2)  # visual inspection of GOF
#' 
#' gof4 <- gof(model4, roc = FALSE, pr = FALSE)
#' gof4
#' plot(gof4)
#' 
#' # MCMC diagnostics
#' pdf("diagnostics2.pdf")
#' mcmc.diagnostics(model2)
#' dev.off()
#' 
#' pdf("diagnostics4.pdf")
#' mcmc.diagnostics(model4)
#' dev.off()
#' }
NULL