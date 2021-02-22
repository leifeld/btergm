#' Check if a matrix is a one-mode matrix
#' 
#' Check if a matrix is a one-mode matrix.
#' 
#' @param mat A matrix object containing zeros and ones.
#' @return \code{TRUE} if the input matrix \code{mat} represents a one-mode
#'   network and \code{FALSE} otherwise.
#'
#' @noRd
is.mat.onemode <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  } else if (!is.null(rownames(mat)) && !is.null(colnames(mat))
      && any(rownames(mat) != colnames(mat))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Check if a matrix represents a directed network
#' 
#' Check if a matrix represents a directed network.
#' 
#' @param mat A matrix object containing zeros and ones.
#' @return \code{TRUE} if the input matrix \code{mat} represents a directed
#'   network and \code{FALSE} otherwise.
#'
#' @noRd
is.mat.directed <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  } else if (!is.null(rownames(mat)) && !is.null(colnames(mat))
      && any(rownames(mat) != colnames(mat), na.rm = TRUE)) {
    return(FALSE)
  } else {
    if (any(as.matrix(mat) != t(as.matrix(mat)), na.rm = TRUE)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

#' How many NAs are there per row or column?
#' 
#' How many NAs are there per row or column?
#' 
#' @param mat A matrix object containing zeros and ones.
#' @param type Count number of missings in rows (\code{type = "row"}), columns
#'   (\code{type = "col"}), or the sum of both (\code{type = "both"})?
#' @param na Value that represents missing data.
#' 
#' @return A numeric count of missing observations in the matrix.
#' 
#' @noRd
numMissing <- function(mat, type = "both", na = NA) {
  numrow <- apply(as.matrix(mat), 1, function(x) sum(x %in% na))
  numcol <- apply(as.matrix(mat), 2, function(x) sum(x %in% na))
  if (type == "both") {
    return(numrow + numcol)
  } else if (type == "row") {
    return(numrow)
  } else if (type == "col") {
    return(numcol)
  } else {
    stop("Unknown 'type' argument in the 'numMissing' function.")
  }
}

#' Handle missing data in matrices
#' 
#' Process \code{NA} values (= remove nodes with \code{NA}s iteratively).
#' 
#' This function deals with missing data in matrices or network objects used for
#' inferential network analysis. It can either remove missing rows and/or
#' columns iteratively (rows and columns with more \code{NA} values first, then
#' successively rows and columns with fewer \code{NA} entries) or replace
#' missing values by the modal value of the matrix or by \code{0}. The function
#' can return either the manipulated matrix or a matrix with logical values
#' indicating which of the cells should be removed.
#' 
#' @param mat A matrix object.
#' @param na The value that missing data are coded as. Usually \code{NA},
#'   sometimes \code{9} or \code{10}.
#' @param method What should be done with the missing data? If
#'   \code{method = "remove"} is set, the function determines how many missing
#'   entries are in each row and column and iteratively removes rows or columns
#'   with the largest amount of missing data until no missing data are left in
#'   the matrix. If \code{method = "fillmode"} is set, the modal value of the
#'   matrix is identified (usually \code{0} in network matrices) and missing
#'   cells are imputed by filling in this modal value. \code{method = "zero"}
#'   replaces NAs by 0s.
#' @param logical Return a matrix with logical values indicating which cells
#'   should be removed? By default the manipulated matrix is returned.
#' 
#' @return Either a matrix in which missing data were taken care of or a matrix
#'   indicating where missing data are located.
#' 
#' @seealso \code{\link{adjust}}
#' 
#' @importFrom network is.network is.bipartite is.directed network.size
#'   list.vertex.attributes get.vertex.attribute network set.vertex.attribute
#' @export
handleMissings <- function(mat, na = NA, method = "remove", logical = FALSE) {

  # check and convert arguments
  if (is.null(mat)) {
    stop("The 'mat' argument is not valid.")
  } else if ("list" %in% class(mat)) {
    # OK; do nothing; check later in next step
    initialtype <- "list"
  } else if (is.matrix(mat) || is.network(mat) || is.data.frame(mat)) {
    # wrap in list
    initialtype <- class(mat)
    mat <- list(mat)
  } else if (length(mat) > 1) {
    # vector --> wrap in list
    initialtype <- class(mat)
    mat <- list(mat)
  } else if (is.function(mat)) {
    stop(paste("The input object is a function. Did you choose a name",
        "for the input object which already exists as a function in the",
        "workspace?"))
  } else {
    stop("The 'mat' argument is not valid.")
  }

  onemode <- list()  # will indicate whether it is a one- or two-mode network
  directed <- list()  # will indicate whether the network is directed
  attribnames <- list()  # will contain the names of nodal attributes
  attributes <- list()  # will contain the nodal attributes at each time step
  type <- list()  # will indicate the type of data structure at time step i
  for (i in 1:length(mat)) {
    if (is.matrix(mat[[i]])) {
      # check manually if onemode and directed
      onemode[[i]] <- is.mat.onemode(mat[[i]])  # helper function
      directed[[i]] <- is.mat.directed(mat[[i]])  # helper function
      type[[i]] <- "matrix"
    } else if (is.network(mat[[i]])) {
      # save onemode and directed information; save attributes for later use
      if (is.bipartite(mat[[i]])) {
        onemode[[i]] <- FALSE
      } else {
        onemode[[i]] <- TRUE
      }
      if (is.directed(mat[[i]])) {
        directed[[i]] <- TRUE
      } else {
        directed[[i]] <- FALSE
      }
      attribnames[[i]] <- list.vertex.attributes(mat[[i]])
      attrib <- list()  # list of attributes at time i
      if (network.size(mat[[i]]) == 0) {
        attribnames[[i]] <- character()
        attributes[[i]] <- list(character())
      } else {
        for (j in 1:length(attribnames[[i]])) {
          attrib[[j]] <- get.vertex.attribute(mat[[i]], attribnames[[i]][j])
        }
      }
      attributes[[i]] <- attrib
      mat[[i]] <- as.matrix(mat[[i]])
      type[[i]] <- "network"
    } else if (is.data.frame(mat[[i]])) {
      type[[i]] <- "data.frame"
    } else {
      type[[i]] <- class(mat[[i]])
    }
  }

  if (is.null(logical) || !is.logical(logical) || length(logical) > 1) {
    stop("The 'logical' argument should be either TRUE or FALSE.")
  }
  if (is.null(method) || !is.character(method)) {
    stop("The 'method' argument should be a character object.")
  }
  if (length(method) > 1) {
    method <- method[1]
  }

  na.mat <- list()  # will contain matrices indicating which values are NA
  for (i in 1:length(mat)) {
    na.mat[[i]] <- apply(mat[[i]], 1:2, function(x) x %in% NA)
    if (length(mat) == 1) {  # used for reporting later
      time <- ""
    } else {
      time <- paste0("t = ", i, ": ")
    }
    if (is.matrix(mat[[i]])) {
      # matrix objects
      # replace by real NAs, then count NAs
      obs <- length(mat[[i]])
      missing.abs <- length(which(mat[[i]] %in% na))
      missing.perc <- round(100 * missing.abs / obs, digits = 2)

      # do the actual work
      if (method == "fillmode") {
        # fill with modal value (often 0 but not always)
        nwunique <- unique(as.numeric(mat[[i]][!mat[[i]] %in% na]))
        nwmode <- nwunique[which.max(tabulate(match(mat[[i]][!mat[[i]] %in%
            na], nwunique)))]
        mat[[i]][mat[[i]] %in% na] <- nwmode
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ",
            missing.abs, " ties) were replaced by the mode (", nwmode,
            ") because they were NA."))
      } else if (method == "zero") {
        # impute 0 when NA
        mat[[i]][mat[[i]] %in% na] <- 0
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ",
            missing.abs, " ties) were replaced by 0 because they were NA."))
      } else if (method == "remove") {
        # remove rows and columns with NA values iteratively
        rowLabels <- rownames(mat[[i]])
        colLabels <- colnames(mat[[i]])
        if (onemode[[i]] == TRUE) {
          while(sum(numMissing(mat[[i]], na = na)) > 0) {
            indices <- which(numMissing(mat[[i]], na = na) ==
                max(numMissing(mat[[i]], na = na)))
            mat[[i]] <- mat[[i]][-indices, -indices]
            rowLabels <- rowLabels[-indices]
            colLabels <- colLabels[-indices]
            na.mat[[i]][indices, ] <- TRUE
            na.mat[[i]][, indices] <- TRUE
            if ("network" %in% type[[i]]) {
              if (length(attribnames[[i]]) > 0) {
                for (j in 1:length(attribnames[[i]])) {
                  attributes[[i]][[j]] <- attributes[[i]][[j]][-indices]
                }
              }
            }
          }
        } else {
          while(sum(numMissing(mat[[i]], type = "row", na = na)) +
              sum(numMissing(mat[[i]], type = "col", na = na)) > 0) {
            rowNAs <- numMissing(mat[[i]], type = "row", na = na)
            colNAs <- numMissing(mat[[i]], type = "col", na = na)
            maxNA <- max(c(rowNAs, colNAs))
            if (length(which(rowNAs == maxNA)) > 0) {
              indices <- which(rowNAs == maxNA)
              mat[[i]] <- mat[[i]][-indices, ]
              rowLabels <- rowLabels[-indices]
              na.mat[[i]][indices, ] <- TRUE
              if ("network" %in% type[[i]]) {
                if (length(attribnames[[i]]) > 0) {
                  for (j in 1:length(attribnames[[i]])) {
                    attributes[[i]][[j]] <- attributes[[i]][[j]][-indices]
                  }
                }
              }
            } else if (length(which(colNAs == maxNA)) > 0) {
              indices <- which(colNAs == maxNA)
              mat[[i]] <- mat[[i]][, -indices]
              colLabels <- colLabels[-indices]
              na.mat[[i]][, indices] <- TRUE
              # in bipartite networks, attributes for rows and columns are
              # saved in a single vector consecutively
              indices.bip <- nrow(mat[[i]]) + indices
              if ("network" %in% type[[i]]) {
                if (length(attribnames[[i]]) > 0) {
                  for (j in 1:length(attribnames[[i]])) {
                    attributes[[i]][[j]] <- attributes[[i]][[j]][-indices.bip]
                  }
                }
              }
            }
          }
        }
        rownames(mat[[i]]) <- rowLabels
        colnames(mat[[i]]) <- colLabels
        removed.abs <- obs - length(mat[[i]])
        removed.perc <- round(100 * removed.abs / obs, digits = 2)
        if (is.nan(removed.perc)) {
          removed.perc <- 0
        }
        if (is.nan(missing.perc)) {
          missing.perc <- 0
        }
        message(paste0("t = ", i, ": ", removed.perc, "% of the data (= ",
            removed.abs, " ties) were dropped due to ", missing.perc, "% (= ",
            missing.abs, ") missing ties."))
      } else {
        stop("Method not supported.")
      }

      # convert back into network if initial item was a network
      if ("network" %in% type[[i]]) {
        bip <- (onemode[[i]] == FALSE)
        mat[[i]] <- network(mat[[i]], directed = directed[[i]],
            bipartite = bip)
        if (length(attribnames[[i]]) > 0) {
          for (j in 1:length(attribnames[[i]])) {
            mat[[i]] <- set.vertex.attribute(mat[[i]], attribnames[[i]][j],
                attributes[[i]][[j]])
          }
        }
      }
    } else if (is.data.frame(mat[[i]])) {
      # data.frame objects
      # replace by real NAs, then count NAs
      for (j in 1:nrow(mat[[i]])) {
        for (k in 1:ncol(mat[[i]])) {
          if (mat[[i]][j, k] %in% na) {
            mat[[i]][j, k] <- NA
          }
        }
      }
      obs <- nrow(mat[[i]]) * ncol(mat[[i]])
      missing.abs <- length(which(is.na(mat[[i]])))
      missing.perc <- round(100 * missing.abs / obs, digits = 2)

      # do the actual work
      if (method == "fillmode") {
        # fill with modal value (often 0 but not always)
        for (j in 1:ncol(mat[[i]])) {
          if (is.numeric(mat[[i]][, j])) {
            nwunique <- unique(as.numeric(mat[[i]][, j]))
            nwmode <- nwunique[which.max(tabulate(match(mat[[i]][, j],
                nwunique)))]
            mat[[i]][, j][is.na(mat[[i]][, j])] <- nwmode
          }
        }
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ",
            missing.abs, " ties) were replaced by the mode in the respective ",
            "column because they were NA."))
      } else if (method == "zero") {
        # impute 0 when NA
        for (j in 1:ncol(mat[[i]])) {
          mat[[i]][, j][is.na(mat[[i]][, j])] <- 0
        }
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ",
            missing.abs, " elements) were replaced by 0 because they were NA."))
      } else if (method == "remove") {
        # remove rows with NA values
        before <- nrow(mat[[i]])
        mat[[i]] <- mat[[i]][complete.cases(mat[[i]]), ]
        after <- nrow(mat[[i]])
        removed <- before - after
        rem.perc <- 100 * (1 - after / before)
        message(paste0("t = ", i, ": ", removed, " rows (", rem.perc,
            "% of all rows) were removed due to missing elements."))
      } else {
        stop("Method not supported.")
      }
    } else if (length(mat[[i]]) > 1) {
      # vectors of arbitrary content
      mat[[i]][mat[[i]] %in% na] <- NA
      obs <- length(mat[[i]])
      missing.abs <- length(which(is.na(mat[[i]])))
      missing.perc <- round(100 * missing.abs / obs, digits = 2)

      # do the actual work
      if (method == "fillmode") {
        # fill with modal value (often 0 but not always)
        if (!is.numeric(mat[[i]])) {
          stop("'fillmode' is only compatible with numeric objects.")
        }
        nwunique <- unique(as.numeric(mat[[i]]))
        nwmode <- nwunique[which.max(tabulate(match(mat[[i]], nwunique)))]
        mat[[i]][is.na(mat[[i]])] <- nwmode
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ",
            missing.abs, " ties) were replaced by the mode in the respective ",
            "column because they were NA."))
      } else if (method == "zero") {
        # impute 0 when NA
        mat[[i]][is.na(mat[[i]])] <- 0
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ",
            missing.abs, " ties) were replaced by 0 because they were NA."))
      } else if (method == "remove") {
        # remove NA values
        mat[[i]] <- mat[[i]][!is.na(mat[[i]])]
        message(paste0(time, missing.perc, "% of the data (= ",
            missing.abs, " elements) were removed because they were NA."))
      } else {
        stop("Method not supported.")
      }
    }
  }

  if (logical == TRUE) {
    if (length(na.mat) == 1 && !"list" %in% initialtype) {
      return(na.mat[[1]])
    } else {
      return(na.mat)
    }
  } else {
    if (length(mat) == 1 && !"list" %in% initialtype) {
      return(mat[[1]])
    } else {
      return(mat)
    }
  }
}

#' Adjust the dimensions of a source object to the dimensions of a target object
#' 
#' Adjust the dimensions of a source object to the dimensions of a target
#' object.
#' 
#' An adjacency matrix (the \code{source} matrix) is compared to another
#' adjacency matrix (the \code{target} matrix) by matching the row or column
#' labels. If the target matrix contains rows/columns which are not present in
#' the source matrix, new rows and columns with the corresponding labels and
#' \code{NA} values in the cells are inserted into the source matrix. If the
#' source matrix contains rows/columns which are not present in the target
#' matrix, these rows and columns are removed from the source matrix. In
#' addition to adjacency matrices, two-mode matrices, network objects (also with
#' vertex attributes), and vectors are supported.
#' 
#' Note that it is not necessary to use this function to preprocess any data
#' before estimating a TERGM. The estimation functions in the \pkg{btergm}
#' package call this function repeatedly to mutually adjust all data as needed.
#' 
#' @param source A matrix, network, list or data.frame object or a vector which
#'   should be adjusted.
#' @param target A matrix, network, list or data.frame object or a vector to
#'   which the source object is compared with regard to its labels.
#' @param remove Should rows and columns that are not present in the target
#'   object be removed?
#' @param add Should rows and columns that are present in the target object but
#'   not in the source object be added to the source object?
#' @param value The value to be inserted if a new row or column is added. By
#'   default, new cells are filled with \code{NA} values, but other sensible
#'   values may include \code{-Inf} or \code{0}.
#' @param returnlabels Return a list of added and removed row and column labels
#'   rather than the actual matrix, vector, or network object?
#' 
#' @examples
#' # create sociomatrix a with 13 vertices a to m
#' vertices <- letters[1:13]
#' a <- matrix(rbinom(length(vertices)^2, 1, 0.1), nrow = length(vertices))
#' rownames(a) <- colnames(a) <- vertices
#' 
#' # create matrix b with the same vertices except f and k, but additional n
#' vertices <- c(vertices[-c(6, 11)], "n")
#' b <- matrix(rbinom(length(vertices)^2, 1, 0.1), nrow = length(vertices))
#' rownames(b) <- colnames(b) <- vertices
#' 
#' # check dimensions
#' dim(a)  # 13 x 13
#' dim(b)  # 12 x 12
#' 
#' # adjust a to b: add n and fill up with NAs; remove f and k
#' adjust(a, b, add = TRUE, remove = TRUE)
#' 
#' \dontrun{
#' # more complex example with additional attributes stored in the network
#' # object; convert a to network object with additional vertex and network
#' # attributes
#' nw <- network(a)
#' vertices <- letters[1:13]
#' nwattrib1 <- matrix(rbinom(length(vertices)^2, 1, 0.1),
#'                     nrow = length(vertices))
#' nwattrib2 <- nwattrib1
#' rownames(nwattrib1) <- colnames(nwattrib1) <- vertices
#' set.network.attribute(nw, "nwattrib1", nwattrib1)
#' set.network.attribute(nw, "nwattrib2", nwattrib2)
#' set.vertex.attribute(nw, "vattrib", 1:length(vertices))
#' 
#' # check presence of the two attributes
#' list.network.attributes(nw)  # nwattrib1 and nwattrib2 are listed
#' get.network.attribute(nw, "nwattrib1")  # returns sociomatrix with labels
#' get.network.attribute(nw, "nwattrib2")  # returns sociomatrix without labels
#' list.vertex.attributes(nw)  # vattrib is listed
#' get.vertex.attribute(nw, "vattrib")  # returns numeric vector 1:13
#' 
#' # adjust the network including the two attributes
#' nw.adjusted <- adjust(nw, b, add = TRUE, remove = TRUE)
#' as.matrix(nw.adjusted)  # note that the order of nodes may have changed
#' get.network.attribute(nw.adjusted, "nwattrib1")  # returns adjusted matrix
#' get.network.attribute(nw.adjusted, "nwattrib2")  # returns adjusted matrix
#' get.vertex.attribute(nw.adjusted, "vattrib")  # returns adjusted vector
#' }
#' 
#' @seealso \code{\link{handleMissings}}
#' 
#' @importFrom network is.network list.vertex.attributes get.vertex.attribute
#'   is.bipartite is.directed list.network.attributes get.network.attribute
#'   set.vertex.attribute network set.network.attribute
#' @export
adjust <- function(source, target, remove = TRUE, add = TRUE, value = NA,
    returnlabels = FALSE) {

  # make sure the source is a list
  if (is.null(source)) {
    stop("The 'source' argument was not recognized.")
  } else if (is.matrix(source)) {
    # wrap in list
    sources <- list()
    sources[[1]] <- source
    sources.initialtype <- "matrix"
  } else if (is.network(source)) {
    # wrap in list
    sources <- list()
    sources[[1]] <- source
    sources.initialtype <- "network"
  } else if ("list" %in% class(source)) {
    # rename
    sources <- source
    sources.initialtype <- "list"
  } else if (is.vector(source)) {
    # vector of some type; wrap in list
    sources <- list()
    sources[[1]] <- source
    sources.initialtype <- "vector"
  } else {
    stop(paste("Source data type not supported. Supported types are 'matrix',",
        "'network', and 'list' objects and vectors."))
  }

  # make sure the target is a list
  if (is.null(target)) {
    stop("The 'target' argument was not recognized.")
  } else if (is.matrix(target)) {
    # wrap in list
    targets <- list()
    targets[[1]] <- target
    targets.initialtype <- "matrix"
  } else if (is.network(target)) {
    # wrap in list
    targets <- list()
    targets[[1]] <- target
    targets.initialtype <- "network"
  } else if ("list" %in% class(target)) {
    # rename
    targets <- target
    targets.initialtype <- "list"
  } else if (is.vector(target)) {
    # vector of some type; wrap in list
    targets <- list()
    targets[[1]] <- target
    targets.initialtype <- "vector"
  } else {
    stop(paste("Target data type not supported. Supported types are 'matrix',",
        "'network', and 'list' objects and vectors."))
  }

  # make sure that both lists (sources and targets) have the same length
  if (length(sources) == length(targets)) {
    # OK; do nothing
  } else if (length(sources) == 1) {
    for (i in 2:length(targets)) {
      sources[[i]] <- sources[[1]]
    }
  } else if (length(targets) == 1) {
    for (i in 2:length(sources)) {
      targets[[i]] <- targets[[1]]
    }
  } else {
    stop("Different numbers of sources and targets were provided.")
  }

  # convert each item if necessary and save nodal attributes
  sources.attribnames <- list()  # names of additional vertex attributes
  sources.attributes <- list()  # additional vertex attributes
  sources.types <- list()  # matrix, network etc.
  sources.onemode <- list()  # is the source network a one-mode network?
  sources.directed <- list()  # is the source network directed?
  sources.matrixnames <- list()  # names of additional matrices
  sources.matrices <- list()  # additional matrices stored in the source network
  targets.attribnames <- list()  # names of additional vertex attributes
  targets.attributes <- list()  # additional vertex attributes
  targets.types <- list()  # matrix, network etc.
  targets.onemode <- list()  # is the target network a one-mode network?
  targets.directed <- list()  # is the source network directed?
  for (i in 1:length(sources)) {
    sources.types[[i]] <- class(sources[[i]])
    if (is.network(sources[[i]])) {
      # save source attributes and other meta information in list
      sources.attribnames[[i]] <- list.vertex.attributes(sources[[i]])
      attributes <- list()
      if (!is.null(sources.attribnames[[i]]) &&
            length(sources.attribnames[[i]]) > 0) {
        for (j in 1:length(sources.attribnames[[i]])) {
          attributes[[j]] <- get.vertex.attribute(sources[[i]],
              sources.attribnames[[i]][j])
        }
      }
      sources.attributes[[i]] <- attributes
      sources.onemode[[i]] <- !is.bipartite(sources[[i]])
      sources.directed[[i]] <- is.directed(sources[[i]])

      # network attributes (= other matrices)
      temp <- list.network.attributes(sources[[i]])
      temp <- temp[!temp %in% c("bipartite", "directed", "hyper", "loops",
          "mnext", "multiple", "n")]
      if (length(temp) > 0) {
        for (j in length(temp):1) {
          cl <- class(get.network.attribute(sources[[i]], temp[j]))
          if (!"network" %in% cl && !"matrix" %in% cl && !"Matrix" %in% cl) {
            temp <- temp[-j]
          }
        }
      }
      sources.matrixnames[[i]] <- temp
      matrices <- list()
      if (!is.null(sources.matrixnames[[i]]) &&
            length(sources.matrixnames[[i]]) > 0) {
        for (j in 1:length(sources.matrixnames[[i]])) {
          matrices[[j]] <- get.network.attribute(sources[[i]],
              sources.matrixnames[[i]][j])
        }
      }
      sources.matrices[[i]] <- matrices
      rm(temp)

      sources[[i]] <- as.matrix(sources[[i]])  # convert to matrix
    } else if (is.matrix(sources[[i]])) {
      sources.onemode[[i]] <- is.mat.onemode(sources[[i]])
      sources.directed[[i]] <- is.mat.directed(sources[[i]])
    } else {
      sources[[i]] <- as.matrix(sources[[i]], ncol = 1)
    }

    targets.types[[i]] <- class(targets[[i]])
    if (is.network(targets[[i]])) {
      # save target attributes and other meta information in list
      targets.attribnames[[i]] <- list.vertex.attributes(targets[[i]])
      attributes <- list()
      if (!is.null(targets.attribnames[[i]]) &&
            length(targets.attribnames[[i]]) > 0) {
        for (j in 1:length(targets.attribnames[[i]])) {
          attributes[[j]] <- get.vertex.attribute(targets[[i]],
              targets.attribnames[[i]][j])
        }
      }
      targets.attributes[[i]] <- attributes
      targets.onemode[[i]] <- !is.bipartite(targets[[i]])
      targets.directed[[i]] <- is.directed(targets[[i]])
      targets[[i]] <- as.matrix(targets[[i]])  # convert to matrix
    } else if (is.matrix(targets[[i]])) {
      targets.onemode[[i]] <- is.mat.onemode(targets[[i]])
      targets.directed[[i]] <- is.mat.directed(targets[[i]])
    } else {
      targets[[i]] <- as.matrix(targets[[i]], ncol = 1)
    }
  }

  # impute row or column labels if only one of them is present
  for (i in 1:length(sources)) {
    if (is.null(rownames(sources[[i]])) && !is.null(colnames(sources[[i]])) &&
        nrow(sources[[i]]) == ncol(sources[[i]])) {
      rownames(sources[[i]]) <- colnames(sources[[i]])
    }
    if (is.null(colnames(sources[[i]])) && !is.null(rownames(sources[[i]])) &&
        nrow(sources[[i]]) == ncol(sources[[i]])) {
      colnames(sources[[i]]) <- rownames(sources[[i]])
    }
    if (is.null(rownames(targets[[i]])) && !is.null(colnames(targets[[i]])) &&
        nrow(targets[[i]]) == ncol(targets[[i]])) {
      rownames(targets[[i]]) <- colnames(targets[[i]])
    }
    if (is.null(colnames(targets[[i]])) && !is.null(rownames(targets[[i]])) &&
        nrow(targets[[i]]) == ncol(targets[[i]])) {
      colnames(targets[[i]]) <- rownames(targets[[i]])
    }
  }

  # throw error if there are duplicate names (first sources, then targets)
  for (i in 1:length(sources)) {
    if ("matrix" %in% class(sources[[i]]) || "data.frame" %in% class(sources[[i]])) {
      # row names
      if (!is.null(rownames(sources[[i]]))) {
        test.actual <- nrow(sources[[i]])
        test.unique <- length(unique(rownames(sources[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif,
              " duplicate source row names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif,
              " duplicate source row name."))
        }
      }
      # column names
      if (!is.null(colnames(sources[[i]]))) {
        test.actual <- ncol(sources[[i]])
        test.unique <- length(unique(colnames(sources[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif,
              " duplicate source column names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif,
              " duplicate source column name."))
        }
      }
    } else {
      # vector names
      if (!is.null(names(sources[[i]]))) {
        test.actual <- length(sources[[i]])
        test.unique <- length(unique(names(sources[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif,
              " duplicate source names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif,
              " duplicate source name."))
        }
      }
    }
  }
  for (i in 1:length(targets)) {
    if ("matrix" %in% class(targets[[i]]) || "data.frame" %in% class(targets[[i]])) {
      # row names
      if (!is.null(rownames(targets[[i]]))) {
        test.actual <- nrow(targets[[i]])
        test.unique <- length(unique(rownames(targets[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif,
              " duplicate target row names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif,
              " duplicate target row name."))
        }
      }
      # column names
      if (!is.null(colnames(targets[[i]]))) {
        test.actual <- ncol(targets[[i]])
        test.unique <- length(unique(colnames(targets[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif,
              " duplicate target column names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif,
              " duplicate target column name."))
        }
      }
    } else {
      # vector names
      if (!is.null(names(targets[[i]]))) {
        test.actual <- length(targets[[i]])
        test.unique <- length(unique(names(targets[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif,
              " duplicate target names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif,
              " duplicate target name."))
        }
      }
    }
  }

  # add original labels to saved network attributes (= matrices) if necessary
  for (i in 1:length(sources)) {
    if ("network" %in% sources.types[[i]] && !is.null(sources.matrices[[i]])
        && length(sources.matrices[[i]]) > 0) {
      for (j in 1:length(sources.matrices[[i]])) {
        if (nrow(as.matrix(sources.matrices[[i]][[j]])) !=
            nrow(as.matrix(sources[[i]])) ||
            ncol(as.matrix(sources.matrices[[i]][[j]])) !=
            ncol(as.matrix(sources[[i]]))) {
          warning(paste("Network attribute", sources.matrixnames[[i]][j],
              "does not have the same dimensions as the source network at",
              "time step", i, "."))
        }
        if (is.network(sources.matrices[[i]][[j]])) {
          if (sources.onemode[[i]] == TRUE) {
            sources.matrices[[i]][[j]] <- set.vertex.attribute(
                sources.matrices[[i]][[j]], "vertex.names",
                rownames(as.matrix(sources[[i]])))
          } else {
            sources.matrices[[i]][[j]] <- set.vertex.attribute(
                sources.matrices[[i]][[j]], "vertex.names",
                c(rownames(as.matrix(sources[[i]])),
                colnames(as.matrix(sources[[i]]))))
          }
        } else {
          rownames(sources.matrices[[i]][[j]]) <-
              rownames(as.matrix(sources[[i]]))
          colnames(sources.matrices[[i]][[j]]) <-
              colnames(as.matrix(sources[[i]]))
        }
      }
    }
  }

  # go through sources and targets and do the actual adjustment
  for (i in 1:length(sources)) {
    if (!is.vector(sources[[i]]) && !is.matrix(sources[[i]]) && !is.network(sources[[i]])) {
      stop(paste("Source item", i, "is not a matrix, network, or vector."))
    }
    if (!is.vector(targets[[i]]) && !is.matrix(targets[[i]]) && !is.network(targets[[i]])) {
      stop(paste("Target item", i, "is not a matrix, network, or vector."))
    }

    # add
    add.row.labels <- character()
    add.col.labels <- character()
    if (add == TRUE) {
      # compile source and target row and column labels
      nr <- nrow(sources[[i]])  # save for later use
      source.row.labels <- rownames(sources[[i]])
      if (!"matrix" %in% sources.types[[i]] && !"network" %in% sources.types[[i]]) {
        source.col.labels <- rownames(sources[[i]])
      } else {
        source.col.labels <- colnames(sources[[i]])
      }
      if ("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]]) {
        if (is.null(source.row.labels)) {
          stop(paste0("The source at t = ", i,
              " does not contain any row labels."))
        }
        if (is.null(source.col.labels)) {
          stop(paste0("The source at t = ", i,
              " does not contain any column labels."))
        }
      }

      target.row.labels <- rownames(targets[[i]])
      if (!"matrix" %in% targets.types[[i]] && !"network" %in% targets.types[[i]]) {
        target.col.labels <- rownames(targets[[i]])
      } else {
        target.col.labels <- colnames(targets[[i]])
      }
      if (is.null(target.row.labels)) {
        stop(paste0("The target at t = ", i,
            " does not contain any row labels."))
      }
      if ("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]]) {
        if (is.null(target.col.labels)) {
          stop(paste0("The target at t = ", i,
              " does not contain any column labels."))
        }
      }

      add.row.indices <- which(!target.row.labels %in% source.row.labels)
      add.row.labels <- target.row.labels[add.row.indices]
      add.col.indices <- which(!target.col.labels %in% source.col.labels)
      add.col.labels <- target.col.labels[add.col.indices]

      # adjust rows
      if (length(add.row.indices) > 0) {
        for (j in 1:length(add.row.indices)) {
          insert <- rep(value, ncol(sources[[i]]))
          part1 <- sources[[i]][0:(add.row.indices[j] - 1), ]
          if (!is.matrix(part1)) {
            if ("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]]) {
              part1 <- matrix(part1, nrow = 1)
            } else {
              part1 <- matrix(part1, ncol = 1)
            }
          }
          rownames(part1) <- rownames(sources[[i]])[0:(add.row.indices[j] - 1)]
          if (add.row.indices[j] <= nrow(sources[[i]])) {
            part2 <- sources[[i]][add.row.indices[j]:nrow(sources[[i]]), ]
          } else {
            part2 <- matrix(ncol = ncol(sources[[i]]), nrow = 0)
          }
          if (!is.matrix(part2)) {
            part2 <- matrix(part2, nrow = 1)
          }
          if (nrow(part2) > 0) {
            rownames(part2) <- rownames(sources[[i]])[add.row.indices[j]:
                nrow(sources[[i]])]
            sources[[i]] <- rbind(part1, insert, part2)
          } else {
            sources[[i]] <- rbind(part1, insert)
          }
          rownames(sources[[i]])[add.row.indices[j]] <- add.row.labels[j]

          # adjust nodal attributes (in the one-mode case)
          if ("network" %in% sources.types[[i]] && sources.onemode[[i]] == TRUE) {
            for (k in 1:length(sources.attributes[[i]])) {
              at1 <- sources.attributes[[i]][[k]][0:(add.row.indices[j] - 1)]
              at2 <- sources.attributes[[i]][[k]][add.row.indices[j]:length(
                  sources.attributes[[i]][[k]])]
              if (sources.attribnames[[i]][k] == "vertex.names") {
                sources.attributes[[i]][[k]] <- c(at1, add.row.labels[j], at2)
              } else if (sources.attribnames[[i]][k] == "na") {
                sources.attributes[[i]][[k]] <- c(at1, TRUE, at2)
              } else {
                sources.attributes[[i]][[k]] <- c(at1, value, at2)
              }
            }
          }
        }
      }

      # adjust columns
      if (length(add.col.indices) > 0 && ("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]])) {
        for (j in 1:length(add.col.indices)) {
          insert <- rep(value, nrow(sources[[i]]))
          part1 <- sources[[i]][, 0:(add.col.indices[j] - 1)]
          if (!is.matrix(part1)) {
            part1 <- matrix(part1, ncol = 1)
          }
          colnames(part1) <- colnames(sources[[i]])[0:(add.col.indices[j] - 1)]
          if (add.col.indices[j] <= ncol(sources[[i]])) {
            part2 <- sources[[i]][, add.col.indices[j]:ncol(sources[[i]])]
          } else {  # if last column, add empty column as second part
            part2 <- matrix(nrow = nrow(sources[[i]]), ncol = 0)
          }
          if (!is.matrix(part2)) {
            part2 <- matrix(part2, ncol = 1)
          }
          if (ncol(part2) > 0) {
            colnames(part2) <- colnames(sources[[i]])[add.col.indices[j]:
                ncol(sources[[i]])]
            sources[[i]] <- cbind(part1, insert, part2)
          } else {
            sources[[i]] <- cbind(part1, insert)
          }
          colnames(sources[[i]])[add.col.indices[j]] <- add.col.labels[j]
        }
      }

      # adjust nodal attributes for two-mode networks
      if ("network" %in% sources.types[[i]] && sources.onemode[[i]] == FALSE) {
        add.col.indices <- sapply(add.col.indices, function(x) x + nr)
        combined.indices <- c(add.row.indices, add.col.indices)
        for (j in 1:length(sources.attributes[[i]])) {
          if (length(combined.indices) > 0) {
            for (k in 1:length(combined.indices)) {
              at1 <- sources.attributes[[i]][[j]][0:(combined.indices[k] - 1)]
              at2 <- sources.attributes[[i]][[j]][combined.indices[k]:length(
                  sources.attributes[[i]][[j]])]
              if (sources.attribnames[[i]][j] == "vertex.names") {
                sources.attributes[[i]][[j]] <- c(at1, add.col.labels[j], at2)
              } else if (sources.attribnames[[i]][j] == "na") {
                sources.attributes[[i]][[j]] <- c(at1, TRUE, at2)
              } else {
                sources.attributes[[i]][[j]] <- c(at1, value, at2)
              }
            }
          }
        }
      }
    }

    removed.rows <- character()
    removed.columns <- character()
    if (remove == TRUE) {
      # compile source and target row and column labels
      nr <- nrow(sources[[i]])  # save for later use
      source.row.labels <- rownames(sources[[i]])
      if (!("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]])) {
        source.col.labels <- rownames(sources[[i]])
      } else {
        source.col.labels <- colnames(sources[[i]])
      }
      if ("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]]) {
        if (nr == 0) {
          stop(paste0("The source at t = ", i, " has no rows."))
        }
        if (is.null(source.row.labels)) {
          stop(paste0("The source at t = ", i,
              " does not contain any row labels."))
        }
        if (is.null(source.col.labels)) {
          stop(paste0("The source at t = ", i,
              " does not contain any column labels."))
        }
      }

      target.row.labels <- rownames(targets[[i]])
      if (!("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]])) {
        target.col.labels <- rownames(targets[[i]])
      } else {
        target.col.labels <- colnames(targets[[i]])
      }
      if ("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]]) {
        if (is.null(target.row.labels)) {
          stop(paste0("The target at t = ", i,
              " does not contain any row labels."))
        }
        if (is.null(target.col.labels)) {
          stop(paste0("The target at t = ", i,
              " does not contain any column labels."))
        }
      }

      # remove
      source.row.labels <- rownames(sources[[i]])
      source.col.labels <- colnames(sources[[i]])
      target.row.labels <- rownames(targets[[i]])
      target.col.labels <- colnames(targets[[i]])
      keep.row.indices <- which(source.row.labels %in% target.row.labels)
      if (("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]]) &&
          ("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]])) {
        keep.col.indices <- which(source.col.labels %in% target.col.labels)
      } else if (("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]])
          && !("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]])) {
        # target is a vector -> keep all columns of source if not onemode
        if (sources.onemode[[i]] == TRUE) {  # columns same as rows
          keep.col.indices <- keep.row.indices
        } else {
          keep.col.indices <- 1:ncol(sources[[i]])
        }
      } else {
        keep.col.indices <- 1
      }
      removed.rows <- which(!1:nrow(as.matrix(sources[[i]])) %in%
          keep.row.indices)
      removed.columns <- which(!1:ncol(as.matrix(sources[[i]])) %in%
          keep.col.indices)

      sources[[i]] <- as.matrix(sources[[i]][keep.row.indices,
          keep.col.indices])
      if ("network" %in% sources.types[[i]]) {
        if (sources.onemode[[i]] == TRUE) {
          for (j in 1:length(sources.attributes[[i]])) {
            sources.attributes[[i]][[j]] <- sources.attributes[[i]][[j]][
                keep.row.indices]
          }
        } else {
          keep.col.indices <- sapply(keep.col.indices, function(x) x + nr)
          combined.indices <- c(keep.row.indices, keep.col.indices)
          for (j in 1:length(sources.attributes[[i]])) {
            sources.attributes[[i]][[j]] <- sources.attributes[[i]][[j]][
                combined.indices]
          }
        }
      }
    }

    # sort source (and attributes) according to row and column names of target
#    if (length(sources.attributes[[i]]) > 0) {
#      for (j in 1:length(sources.attributes[[i]])) {
#        if (!is.null(sources.attributes[[i]][[j]]) &&
#            length(sources.attributes[[i]][[j]]) > 0) {
#          if (sources.onemode[[i]] == TRUE) {
#            names(sources.attributes[[i]][[j]]) <- rownames(sources[[i]])
#            sources.attributes[[i]][[j]] <-
#                sources.attributes[[i]][[j]][rownames(sources[[i]])]
#          } else {
#            names(sources.attributes[[i]][[j]]) <- c(rownames(sources[[i]]),
#                rownames(sources[[i]]))
#            sources.attributes[[i]][[j]] <-
#                c(sources.attributes[[i]][[j]][rownames(sources[[i]])],
#                sources.attributes[[i]][[j]][colnames(sources[[i]])])
#          }
#        }
#      }
#    }
#
    if (("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]]) &&
        ("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]]) &&
        nrow(sources[[i]]) == nrow(targets[[i]]) &&
        ncol(sources[[i]]) == ncol(targets[[i]])) {
      sources[[i]] <- sources[[i]][rownames(targets[[i]]),
          colnames(targets[[i]])]
    } else if (("matrix" %in% sources.types[[i]] || "network" %in% sources.types[[i]]) &&
        !("matrix" %in% targets.types[[i]] || "network" %in% targets.types[[i]]) &&
        nrow(sources[[i]]) == nrow(targets[[i]])) {
      sources[[i]] <- sources[[i]][rownames(targets[[i]]),
          rownames(targets[[i]])]
    } else if (length(sources[[i]]) == nrow(targets[[i]])) {
      # source is a vector, irrespective of the target
      sources[[i]] <- sources[[i]][rownames(targets[[i]]), ]
    } else if (add == FALSE && (nrow(sources[[i]]) < nrow(targets[[i]]) ||
        any(rownames(sources[[i]]) != rownames(targets[[i]])))) {
    }

    # convert back into network
    if ("network" %in% sources.types[[i]]) {
      sources[[i]] <- network(sources[[i]], directed = sources.directed[[i]],
          bipartite = !sources.onemode[[i]])
      for (j in 1:length(sources.attribnames[[i]])) {
        sources[[i]] <- set.vertex.attribute(sources[[i]],
            sources.attribnames[[i]][j], sources.attributes[[i]][[j]])
      }
    }

    # convert vectors back from one-column matrices to vectors
    if (!"matrix" %in% sources.types[[i]] && !"network" %in% sources.types[[i]] &&
        is.matrix(sources[[i]]) && ncol(sources[[i]]) == 1) {
      sources[[i]] <- sources[[i]][, 1]
    }

    # return added and removed labels instead of actual objects
    if (returnlabels == TRUE) {
      sources[[i]] <- list()
      sources[[i]]$removed.row <- removed.rows
      sources[[i]]$removed.col <- removed.columns
      sources[[i]]$added.row <- add.row.labels
      sources[[i]]$added.col <- add.col.labels
    }
  }

  # adjust network attributes (= matrices) recursively and add back in
  for (i in 1:length(sources)) {
    if ("network" %in% sources.types[[i]] && !is.null(sources.matrixnames[[i]])
        && length(sources.matrixnames[[i]]) > 0) {
      for (j in 1:length(sources.matrixnames[[i]])) {
        mat <- adjust(source = sources.matrices[[i]][[j]],
            target = sources[[i]], add = add, remove = remove, value = value)
        set.network.attribute(sources[[i]], sources.matrixnames[[i]][j], mat)
      }
    }
  }

  if ("list" %in% sources.initialtype) {
    return(sources)
  } else {
    return(sources[[1]])
  }
}