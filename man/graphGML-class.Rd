% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/read.gatingML.cytobank.R
\docType{class}
\name{graphGML-class}
\alias{graphGML-class}
\title{A graph object returned by 'read.gatingML.cytobank' function.}
\description{
Each node corresponds to a population(or GateSet) defined in gatingML file. The actual gate object (both global and tailored gates) is
associated with each node as nodeData. Compensation and transformations are stored in graphData slot.
}
\details{
The class simply extends the graphNEL class and exists for the purpose of method dispatching.
}

