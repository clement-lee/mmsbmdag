% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_by_membership}
\alias{plot_by_membership}
\title{Plot adjacency matrix by membership}
\usage{
plot_by_membership(df, edgelist, size = 0.1, alpha = 0.5)
}
\arguments{
\item{df}{a data frame of group membership, with variables id, id_num and group}

\item{edgelist}{a data frame of edges, with varaibles citing, cited, id_num, and group}

\item{size}{size of dots in ggplot2::geom_point()}

\item{alpha}{alpha of dots in ggplot2::geom_point()}
}
\description{
Plot adjacency matrix by membership
}
