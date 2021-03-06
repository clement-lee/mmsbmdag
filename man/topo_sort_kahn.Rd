% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topo_order.R
\name{topo_sort_kahn}
\alias{topo_sort_kahn}
\title{Return a sorted vector of nodes id}
\usage{
topo_sort_kahn(g, random = FALSE)
}
\arguments{
\item{g}{An igraph graph object}

\item{random}{Boolean, whether the order of selected nodes is randomised in the process}
}
\description{
Return a sorted vector of nodes id
}
\examples{
g0 <- igraph::graph_from_data_frame(data.frame(a = c("1", "2"), b = c("2", "3")))
topo_sort_kahn(g0)
}
