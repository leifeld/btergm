\name{handleMissings}
\alias{handleMissings}
\title{Handle missing data in matrices.}
\description{Handle missing data in matrices.}
\usage{
handleMissings(mat, na = NA, method = "remove", logical = FALSE)
}
\details{
This function deals with missing data in matrices or network objects used for inferential network analysis. It can either remove missing rows and/or columns iteratively (rows and columns with more NA values first, then successively rows and columns with fewer NA entries) or replace missing values by the modal value of the matrix or by 0. The function can return either the manipulated matrix or a matrix with logical values indicating which of the cells should be removed.
}
\arguments{
\item{mat}{ A matrix object. }
\item{na}{ The value that missing data are coded as. Usually \code{NA}, sometimes \code{9} or \code{10}. }
\item{method}{ What should be done with the missing data? If \code{method = "remove"} is set, the function determines how many missing entries are in each row and column and iteratively removes rows or columns with the largest amount of missing data until no missing data are left in the matrix. If \code{method = "fillmode"} is set, the modal value of the matrix is identified (usually \code{0} in network matrices) and missing cells are imputed by filling in this modal value. \code{method = "zero"} replaces NAs by 0s. }
\item{logical}{ Return a matrix with logical values indicating which cells should be removed? By default the manipulated matrix is returned. }
}
\seealso{
\link{btergm-package} \link{adjust} \link{preprocess}
}
