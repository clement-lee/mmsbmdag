#' Detect community for a graph
#'
#' @param g an igraph object
#' @param f a community detection function
#' @return a data fram with group memberships and associated numerical id
#' @importFrom tibble data_frame
#' @importFrom dplyr left_join count rename arrange mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
detect_arrange <- function(g, f) {
    ## detect community for graph g & return df with group memberships and associated numerical id
    l <- f(g)
    df0 <- tibble::data_frame(id = l$names, group = l$membership %>% as.integer)
    dplyr::left_join(df0, dplyr::count(df0, .data$group), "group") %>%
        dplyr::rename(group.size = .data$n) %>%
        dplyr::arrange(.data$group, .data$id) %>%
        dplyr::mutate(id_num = seq_along(.data$id))
}

#' Compute list of descendants of id & loop indeices, if any
#'
#' @param id a character string which is a vertex in the graph
#' @param edgelist a data frame representing the edgelist of the graph
#' @param max.level maximum level of descendants to compute; default to 50
#' @return a list of descendants and loop indices
#' @importFrom tibble data_frame
#' @importFrom dplyr semi_join select distinct anti_join bind_rows
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
list_descendants <- function(id, edgelist, max.level = 50L) {
    l0 <- list()
    df0 <- data_frame(citing = vector("character")) # all descendants
    df1 <- data_frame(citing = id) # current descendants
    i <- 1L
    loop_i <- vector("integer")
    while (i <= max.level && nrow(df1) > 0L) {
        ## find ith gen. descendants
        df1 <- edgelist %>%
            dplyr::semi_join(df1, c("cited" = "citing")) %>%
            dplyr::select(.data$citing) %>%
            dplyr::distinct(.data$citing) %>%
            dplyr::anti_join(df0, c("citing"))
        ## check with id to detect loops
        if (id %in% df1$citing) {
            loop_i <- c(loop_i, i)
        }
        ## save & incr.
        df0 <- dplyr::bind_rows(df0, df1)
        l0[[i]] <- df1
        i <- i + 1L
    }
    list(descendants = l0, loop_indices = loop_i)
}
