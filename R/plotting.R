#' Plot adjacency matrix by membership
#'
#' @param df a data frame of group membership, with variables id, id_num and group
#' @param edgelist a data frame of edges, with varaibles citing, cited, id_num, and group
#' @param size size of dots in ggplot2::geom_point()
#' @param alpha alpha of dots in ggplot2::geom_point()
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr left_join mutate
#' @importFrom ggplot2 ggplot geom_point labs scale_y_reverse scale_x_continuous theme_bw theme aes
#' @export
plot_by_membership <- function(df, edgelist, size = 0.1, alpha = 0.5) {
    edgelist %>%
        dplyr::left_join(df, c("citing" = "id")) %>%
        dplyr::left_join(df, c("cited" = "id")) %>%
        dplyr::mutate(group = ifelse(.data$group.x == .data$group.y, as.character(.data$group.x), "0")) %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(.data$id_num.y, .data$id_num.x, colour = .data$group), size = size, alpha = alpha) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::scale_y_reverse(breaks = NULL, expand = c(0.01, 0.0)) +
        ggplot2::scale_x_continuous(breaks = NULL, expand = c(0.01, 0.0)) +
        ggplot2::theme_bw(20) +
        ggplot2::theme(legend.position = "none")
}

#' Plot an igraph object with specific arguments
#'
#' @param g an igraph object
#' @param seed a number for initial seed, interpreted as a positive integer
#' @param ... any other arguments passed to graphics::plot()
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom graphics par plot
#' @importFrom igraph layout_with_fr
#' @export
plot_network <- function(g, seed = 110L, ...) {
    graphics::par(mar = rep(0.0, 4L)) # rm white space
    set.seed(seed)
    g %>%
        graphics::plot(layout = igraph::layout_with_fr, 
                       edge.width = 0.5,
                       edge.arrow.size = 0.1,
                       vertex.label = NA,
                       margin = rep(0.0, 4L), # rm space
                       ...)
    ## force directed layouts in igraph:
    ## 1) layout_with_drl
    ## 2) layout_with_fr
    ## 3) layout_with_gem
    ## 4) layout_with_graphopt
}

#' Raster plot of group-to-group probabilities
#'
#' @param df a data frame with variables mean, row, and col
#' @param size size of text in ggplot2::geom_text(); default to 10
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_x_continuous scale_y_reverse scale_fill_gradient2 theme_void theme labs aes element_blank
#' @export
plot_C_raster <- function(df, size = 10) {
    df %>% ggplot2::ggplot() +
        ggplot2::geom_tile(ggplot2::aes(x = col, y = row, fill = mean)) +
        ggplot2::geom_text(ggplot2::aes(x = col, y = row, label = round(mean, 3)), size = size) +
        ggplot2::scale_x_continuous(breaks = NULL) +
        ggplot2::scale_y_reverse(breaks = NULL) +
        ggplot2::scale_fill_gradient2(midpoint = 0.5, low = "#ece2f0", mid = "#a6bddb", high = "#1c9099") +
        ggplot2::theme_void(20) +
        ggplot2::theme(legend.title = ggplot2::element_blank()) +
        ggplot2::labs(x = NULL, y = NULL)
}

#' Plot the membership heatmaps of all nodes
#'
#' @param df0 an output of mcmc_summary_all()
#' @param df1 an output of assign_cluster()
#' @importFrom dplyr left_join select
#' @importFrom ggplot2 ggplot geom_raster aes scale_y_continuous scale_fill_gradient2 labs coord_cartesian theme_minimal theme element_blank
#' @export
plot_D_heatmap <- function(df0, df1) {
    df0 %>% 
        dplyr::left_join(df1 %>% dplyr::select(.data$id_nogroup, .data$y), "id_nogroup") %>%
        ggplot2::ggplot() +
        ggplot2::geom_raster(ggplot2::aes(.data$group, .data$y, fill = .data$mean)) +
        ggplot2::scale_y_continuous(breaks = NULL, labels = df1$id_nogroup) + # chg. breaks to df1$y if necessary
        ggplot2::scale_fill_gradient2(midpoint = 0.5, low = "#ece2f0", mid = "#a6bddb", high = "#1c9099") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::theme_minimal(20) +
        ggplot2::theme(legend.title = ggplot2::element_blank())
}

#' Plot the positions of all nodes in topo. order
#'
#' @param df an output of mcmc_summary_all()
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarise mutate left_join
#' @importFrom ggplot2 ggplot geom_bin2d aes scale_y_continuous scale_fill_gradient2 theme_minimal theme element_blank labs
#' @export
plot_i_heatmap <- function(df) {
    df0 <- df %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarise(mean = mean(.data$i)) %>%
        dplyr::mutate(y = rank(.data$mean))
    df %>%
        dplyr::left_join(df0, "id") %>%
        ggplot2::ggplot() +
        ggplot2::geom_bin2d(ggplot2::aes(.data$i, .data$y), binwidth = 1) +
        ggplot2::scale_y_continuous(breaks = NULL) + # breaks = df0$y, labels = df0$id
        ggplot2::scale_fill_gradient2(low = "#ece2f0", mid = "#a6bddb", high = "#1c9099") +
        ggplot2::theme_minimal(20) +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       legend.position = "none",
                       panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::labs(x = "Position", y = NULL)
}

#' Plot the (mean) positions of nodes against year
#'
#' @param df0 an output of mcmc_summary_all()
#' @param df1 df0.link in draft.Rnw
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarise left_join
#' @importFrom ggplot2 ggplot geom_point aes theme_minimal labs
#' @export
plot_i_vs_year <- function(df0, df1) {
    df0 %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarise(mean = mean(.data$i)) %>%
        dplyr::left_join(df1, "id") %>%
        ggplot2::ggplot() + 
        ggplot2::geom_point(ggplot2::aes(.data$year, .data$mean)) +
        ggplot2::theme_minimal(20) +
        ggplot2::labs(x = "Year of publication", 
                      y = "Mean position")
}

#' From ggs to a list of gg objects
#'
#' @param df a data frame & ggs object
#' @param density boolean (default to TRUE); should the density be plotted instead of the histogram
#' @importFrom ggmcmc ggs_traceplot ggs_autocorrelation
#' @importFrom ggplot2 theme_bw ggplot geom_histogram aes facet_wrap 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
ggs_list <- function(df, density = TRUE) {
    gg0.tcp <- df %>%
        ggmcmc::ggs_traceplot() +
        ggplot2::theme_bw()
    if (density) {
        gg0.den <- df %>%
            ggmcmc::ggs_density() +
            ggplot2::theme_bw()
    } else {
        gg0.den <- df %>%
            ggplot2:: ggplot() +
            ggplot2::geom_histogram(ggplot2::aes(.data$value), binwidth = 1L) +
            ggplot2::facet_wrap(~Parameter, ncol = 1L) +
            ggplot2::theme_bw()
    }
    gg0.acf <- df %>%
        ggmcmc::ggs_autocorrelation() +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~Parameter, ncol = 1L, scales = "free_y")
    list(gg0.tcp, gg0.den, gg0.acf)
}

#' From group-to-group probabilities to ggs object
#'
#' @param m_C matrix of group-to-group probabilities
#' @param col column to select
#' @importFrom glue glue
#' @importFrom coda mcmc
#' @importFrom ggmcmc ggs
#' @importFrom dplyr select contains
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @export
matrix_to_ggs_C <- function(m_C, col = 1L) {
    m_C %>%
        tibble::as_tibble() %>%
        dplyr::select(dplyr::contains(glue::glue("C_{col}"))) %>%
        coda::mcmc() %>%
        ggmcmc::ggs()
}

#' From memberships to ggs object
#'
#' @param m_D matrix of membership probabilities
#' @param prefix prefix or id to select from
#' @importFrom stringr str_sub
#' @importFrom coda mcmc
#' @importFrom ggmcmc ggs
#' @importFrom dplyr select contains
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @export
matrix_to_ggs_D <- function(m_D, prefix = "abfx08") {
    df0 <- m_D %>%
        tibble::as_tibble() %>%
        dplyr::select(dplyr::contains(prefix))
    colnames(df0) <- paste0("Group ", colnames(df0) %>% stringr::str_sub(-1L))
    df0 %>% coda::mcmc() %>%
        ggmcmc::ggs()
}
