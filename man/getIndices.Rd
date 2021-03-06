% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/GatingHierarchy-Methods.R
\docType{methods}
\name{getIndices,GatingHierarchy,character-method}
\alias{getIndices}
\alias{getIndices,GatingHierarchy,character-method}
\title{Get the membership indices for each event with respect to a particular gate in a GatingHierarchy}
\usage{
\S4method{getIndices}{GatingHierarchy,character}(obj, y)
}
\arguments{
\item{obj}{A \code{GatingHierarchy} representing a sample.}

\item{y}{A \code{character} giving the name or full(/partial) gating path of the population / node of interest.}
}
\value{
A logical vector of length equal to the number of events in the FCS file that determines whether each event is or is not included in the current gate.
}
\description{
Returns a logical vector that describes whether each event in a sample is included or excluded by this gate.
}
\details{
Returns a logical vector that describes whether each event in the data file is included in the given gate of this \code{GatingHierarchy}. The indices are for all events in the file, and do not reflect the population counts relative to the parent but relative to the root. To get population frequencies relative to the parent one cross-tabulate the  indices of \code{y} with the indices of its parent.
}
\note{
Generally you should not need to use \code{getIndices} but the more convenient methods \code{getProp} and \code{getPopStats} which return population frequencies relative to the parent node.
The indices returned reference all events in the file and are not directly suitable for computing population statistics, unless subsets are taken with respect to the parent populations.
}
\examples{
\dontrun{
    #G is a gating hierarchy
    #Return the indices for population 5 (topological sort)
    getIndices(G,getNodes(G,tsort=TRUE)[5]);
}
}
\seealso{
\code{\link{getPopStats}}
}

