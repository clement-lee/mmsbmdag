#' Wrapper of MCMC algorithm for mixed membership stochastic block model (MMSBM)
#'
#' @param g a graph object, assumed to be a DAG
#' @param alpha_init initial value of alpha
#' @param s_alpha_init initial value of proposal standard deviation for alpha's Metropolis step
#' @param a,b hyperparameters for alpha's prior
#' @param K number of clusters, fixed
#' @param n_swap number of swaps in order per iteration
#' @param N,thin,burnin,print_freq MCMC quantities
#' @param seed random seed to initiate MCMC
#' @param f MCMC algo in C++; currently supports rgs_mmsbm()
#' @return a list of data frames and matrices representing the MCMC output
#' @importFrom igraph as_adjacency_matrix
#' @importFrom tibble data_frame
#' @importFrom dplyr left_join mutate
#' @importFrom magrittr set_colnames extract2 %>%
#' @importFrom tidyr unnest
#' @importFrom Rcpp sourceCpp
#' @useDynLib mmsbmdag
#' @export
gs_mmsbm <- function(g, alpha_init, s_alpha_init, a, b, K, n_swap, N, thin, burnin, print_freq, seed, f) {
    ## wrapper of rgs_mmsbm() & cgs_mmsbm() in C++
    ## g assumed to be a DAG
    ## f is either C++ Gibbs sampler function
    Y <- g %>%
        igraph::as_adjacency_matrix() %>%
        as.matrix() # orig. order
    df.names <- tibble::data_frame(id = colnames(Y))
    df.nodes <- dplyr::left_join(df.names, topo_sort_kahn(g), "id") # orig. order
    o_init <- order(df.nodes$id_num) - 1L
    Y_star <- reorder_once(Y, o_init) # upper tri.
    ##    Y_star %>% as("dgCMatrix") %>% tril %>% sum # check Y_star is upper tri.; expected 0
    A <- B <- matrix(1.0, K, K) # simply assume uniform for now
    set.seed(seed)
    t0 <- system.time({
        fit0 <- f(Y, o_init, alpha_init, s_alpha_init, a, b, A, B, n_swap, N, thin, burnin, print_freq)
    })
    ## post-process
    fit0$time <- t0
    colnames(fit0$alpha) <- "alpha"
    fit0$C <- fit0$C %>%
        magrittr::set_colnames(
            expand.grid(j = seq(K), i = seq(K)) %>%
            dplyr::mutate(name = glue("C_{i}{j}")) %>%
            magrittr::extract2("name")
        )
    fit0$D <- fit0$D %>%
        magrittr::set_colnames(
            df.nodes %>%
            dplyr::mutate(group = list(seq(K))) %>%
            tidyr::unnest() %>%
            dplyr::transmute(name = glue("{id}_{group}")) %>%
            magrittr::extract2("name")
            )
    colnames(fit0$o) <- df.nodes$id
    colnames(fit0$i) <- df.nodes$id
    fit0
}
