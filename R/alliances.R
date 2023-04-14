#' Longitudinal international defense alliance network, 1981--2000
#'
#' Longitudinal international defense alliance network, 1981--2000.
#'
#' The alliances dataset contains the international defense alliance network
#' among 164 countries, covering the years 1981--2000. In addition to the
#' yearly defense alliance network, it contains data on military capabilities,
#' governing regime type, geographic contiguity and international conflict.
#' This is an excerpt from a dataset that has been used in two published
#' analyses. The full dataset (Cranmer, Desmarais and Menninga 2012; Cranmer,
#' Desmarais and Kirkland 2012) contains a large number of countries and a much
#' longer time series.
#'
#' @name alliances
#' @aliases alliances allyNet contigMat lNet LSP warNet
#' @docType data
#' @format
#' \describe{
#' \item{\code{allyNet}}{is a list of network objects at 20 time points,
#' 1981--2000, containing undirected defense alliance networks. In addition to
#' the alliance ties, each network object contains three vertex attributes.
#' \code{cinc} is the "CINC" or Composite Index of National Capability score
#' (see
#' \url{https://correlatesofwar.org/data-sets/national-material-capabilities/}).
#' \code{polity} is the "polity score" of each country in the respective year.
#' Quoting the online description, "the Polity Score captures this regime
#' authority spectrum on a 21-point scale ranging from -10 (hereditary
#' monarchy) to +10 (consolidated democracy)," (see
#' \url{https://www.systemicpeace.org/polityproject.html}). \code{year} is
#' simply the year recorded as a vertex attribute.}
#'
#' \item{\code{contigMat}}{is a 164 x 164 binary matrix in which a 1 indicates
#' that two countries share a border.}
#'
#' \item{\code{lNet}}{is a list of 20 matrices. Each element is the adjacency
#' matrix from the previous year. This is used to model memory in the ties.}
#'
#' \item{\code{LSP}}{is a list of 20 matrices. Each element is a matrix
#' recording the number of shared partners between countries in the alliance
#' network from the previous year.}
#'
#' \item{\code{warNet}}{is a list of 20 matrices. Each element is a binary
#' matrix that indicates whether two states were in a militarized interstate
#' dispute in the respective year.}
#' }
#'
#' @references Cranmer, Skyler J., Bruce A. Desmarais, and Justin H. Kirkland
#' (2012): Toward a Network Theory of Alliance Formation. \emph{International
#' Interactions} 38(3): 295--324. \doi{10.1080/03050629.2012.677741}.
#'
#' Cranmer, Skyler J., Bruce A. Desmarais, and Elizabeth Menninga (2012):
#' Complex Dependencies in the Alliance Network. \emph{International
#' Interactions} 29(3): 279--313. \doi{10.1177/0738894212443446}.
#'
#' @source The data were gathered by Skyler Cranmer and Bruce Desmarais in the
#' process of writing Cranmer, Desmarais and Menninga (2012) and Cranmer,
#' Desmarais and Kirkland (2012).
#'
#' Permission to redistribute this dataset along with this package was granted
#' by Skyler Cranmer and Bruce Desmarais on December 15, 2015. Questions about
#' the data should be directed to them.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' data("alliances")
#'
#' # btergm formulas look very similar to ERGM formulas.
#' # Note the R argument; usually want R > 1000.
#' # Here it is set to 50 to limit computation time.
#' # First, set the seed for replicability.
#' set.seed(123)
#' model <- btergm(allyNet ~ edges + gwesp(0, fixed = TRUE)
#'     + edgecov(lNet) + edgecov(LSP) + edgecov(warNet)
#'     + nodecov("polity") + nodecov("cinc") + absdiff("polity")
#'     + absdiff("cinc") + edgecov(contigMat) + nodecov("year"),
#'     R = 50)
#'
#' # View estimates and confidence intervals.
#' summary(model)
#'
#' # Evaluate model fit. Simulate 100 networks for each time point.
#' # Calculate edgewise shared partners, degree and geodesic distance
#' # distance distributions.
#' alliance_gof <- gof(model, statistics = c(deg, esp, geodesic))
#'
#' # Plot goodness of fit.
#' plot(alliance_gof)
#' }
NULL