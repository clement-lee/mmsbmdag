#' Extract node list of a graph object.
#'
#' @param g An igraph graph object.
#' @return An |N| by 1 tibble, where |N| is the number of nodes in g.
#' @examples
#' g0 <- igraph::graph_from_data_frame(data.frame(a = c("1", "2"), b = c("2", "3")))
#' graph_to_nodelist(g0)
#' @importFrom igraph degree
#' @importFrom tibble as_tibble
#' @export
graph_to_nodelist <- function(g) {
    v <- igraph::degree(g)
    df0 <- data.frame(id = names(v), value = v, stringAsFactors = FALSE)
    tibble::as_tibble(df0)
}

#' Extract edge list of a graph object.
#'
#' @param g an igraph graph object.
#' @return An |E| by 2 tibble, where |E| is the number of edges in g.
#' @importFrom igraph as_data_frame
#' @importFrom tibble as_tibble
#' @examples
#' g0 <- igraph::graph_from_data_frame(data.frame(a = c("1", "2"), b = c("2", "3")))
#' graph_to_edgelist(g0)
#' @export
graph_to_edgelist <- function(g) {
    df0 <- igraph::as_data_frame(g)
    tibble::as_tibble(df0)
}

#' From NULL to character NA
#'
#' @param x vector of characters
#' @export
null_to_character_na <- function(x) {
    if (is.null(x)) as.character(NA) else x
}

#' From NULL to integer NA
#'
#' @param x vector of integers
#' @export
null_to_integer_na <- function(x) {
    if (is.null(x)) as.integer(NA) else x
}

#' Extract id of a group
#'
#' @param i integer index representing the group
#' @param df.com output of detect_arrange
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
names_in_group <- function(i, df.com) {
    dplyr::filter(df.com, .data$group == i)$id
}
