#' Summary of MCMC
#'
#' @param m matrix; each row is one MCMC iteration
#' @importFrom tibble data_frame
#' @importFrom magrittr %>%
mcmc_summary <- function(m) {
    ## ess <- m %>% coda::mcmc() %>% coda::effectiveSize
    mean <- m %>% colMeans
    data_frame(
        id = colnames(m),
        ## ess = ess,
        mean = mean
    )
}

#' Summary of a list of MCMC output
#'
#' @param l list; essentially an MCMC output
#' @importFrom dplyr mutate
#' @importFrom magrittr %>% %<>% extract
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom stringr str_sub
#' @export
mcmc_summary_all <- function(l) {
    ## list of mcmc summaries
    l0 <- l %>%
        magrittr::extract(c("alpha", "C", "D"))
    if (!is.null(l$i)) {
        l0$i <- l$i
    }
    l0 %<>% purrr::map(mcmc_summary)
    l0$C %<>%
        dplyr::mutate(row = .data$id %>% stringr::str_sub(-2L, -2L) %>% as.integer(),
                      col = .data$id %>% stringr::str_sub(-1L, -1L) %>% as.integer())
    l0$D %<>%
        dplyr::mutate(group = .data$id %>% stringr::str_sub(-1L, -1L),
                      id_nogroup = .data$id %>% stringr::str_sub(1L, -3L))
    l0
}

#' Assign clusters according to membership probabilities, for plot_D_heatmap()
#'
#' @param df data frame with 3 variables: group, id_nogroup, mean
#' @param misc character strings representing the miscellaneous groups
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr count filter select group_by top_n ungroup summarise left_join transmute desc
#' @importFrom rlang .data
#' @export
assign_cluster <- function(df, misc) {
    if (is.integer(misc)) {
        misc %<>% as.character()
    }
    K <- df %>% dplyr::count(.data$group) %>% nrow
    df0 <- df %>%
        dplyr::filter(!(.data$group %in% misc)) %>%
        dplyr::select(.data$id_nogroup, .data$group, .data$mean) %>%
        dplyr::group_by(.data$id_nogroup) %>%
        dplyr::top_n(1L, .data$mean) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$mean > 1.0 / K)
    df1 <- df %>%
        dplyr::filter(.data$group %in% misc) %>%
        dplyr::select(.data$id_nogroup, .data$group, .data$mean) %>%
        dplyr::group_by(.data$id_nogroup) %>%
        dplyr::top_n(1L, .data$mean) %>%
        dplyr::ungroup()
    df2 <- df %>%
        dplyr::group_by(.data$id_nogroup) %>%
        dplyr::summarise(entropy = -sum(.data$mean * log(.data$mean)))
    dplyr::left_join(df1, df0, "id_nogroup") %>%
        dplyr::transmute(.data$id_nogroup,
                        group = ifelse(is.na(.data$group.y), .data$group.x, .data$group.y),
                        mean = ifelse(is.na(.data$mean.y), .data$mean.x, .data$mean.y)) %>%
            dplyr::arrange(dplyr::desc(.data$group), .data$mean) %>%
            dplyr::mutate(y = seq_along(.data$id_nogroup)) %>%
            dplyr::left_join(df2, "id_nogroup")
}

#' Attach colours and sizes to vertices
#'
#' @param df data frame with memberships, id_nogroup as id, & entropy
#' @param g igraph object
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom tibble data_frame
#' @importFrom igraph V
#' @export
graph_colour_size <- function(df, g) {
    df0 <- tibble::data_frame(id_nogroup = names(igraph::V(g))) %>%
        dplyr::left_join(df, "id_nogroup")
    igraph::V(g)$color <- df0$group
    igraph::V(g)$size <- sqrt(df0$entropy) * 8.0 # manual scaling
    g
}

#' Data frame of position in order
#'
#' @param m the actual chain
#' @importFrom magrittr %>%
#' @importFrom tibble data_frame
#' @export
df_i <- function(m) {
    tibble::data_frame(i = m %>% as.vector,
                       id = rep(m %>% colnames, each = m %>% nrow))
}
