\name{adjust}
\alias{adjust}
\title{Adjust the dimensions of a matrix to the dimensions of another matrix}
\description{Adjust the dimensions of a matrix to the dimensions of another matrix.}
\usage{
adjust(source, target, remove = TRUE, add = TRUE, value = NA, 
    returnlabels = FALSE)
}
\details{
An adjacency matrix (the \code{source} matrix) is compared to another adjacency matrix (the \code{target} matrix) by matching the row or column labels. If the target matrix contains rows/columns which are not present in the source matrix, new rows and columns with the corresponding labels and NA values in the cells are inserted into the source matrix. If the source matrix contains rows/columns which are not present in the target matrix, these rows and columns are removed from the source matrix. In addition to adjacency matrices, two-mode matrices, network objects (also with vertex attributes), and vectors are supported.
}
\arguments{
\item{source}{ A matrix, network, list or data.frame object or a vector which should be adjusted. }
\item{target}{ A matrix, network, list or data.frame object or a vector to which the source object is compared with regard to its labels. }
\item{remove}{ Should rows and columns that are not present in the target object be removed? }
\item{add}{ Should rows and columns that are present in the target object but not in the source object be added to the source object? }
\item{value}{ The value to be inserted if a new row or column is added. By default, new cells are filled with \code{NA} values, but other sensible values may include \code{-Inf} or \code{0}. }
\item{returnlabels}{ Return a list of added and removed row and column labels rather than the actual matrix, vector, or network object? }
}
\seealso{
\link{btergm-package} \link{handleMissings} \link{preprocess}
}
