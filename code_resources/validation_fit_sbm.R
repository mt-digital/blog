################################### Install and/or load STRAND
library(devtools)
library(STRAND)
library(checkmate)
library(igraph)

library(purrrlyr)

# install_github('ctross/STRAND')

# Helper functions for STRAND data formatting.
groups_df = function(group) { return (data.frame(Group = group)) }
groupsf_df = function(group) { return (data.frame(Group = factor(group))) }

sbm_make_min_maj = function(metapop_size = 100, min_frac = 0.05, 
                            mean_degree = 20,
                            homophily_min = 0.0,
                            homophily_maj = 0.0) {

  # If N vertices, then N(N-1) possible edges in network.
  n_possible_edges = metapop_size * (metapop_size - 1)

  # The probability of an edge is the desired mean degree by possible edges count,
  # but this next line is apparently not right, see note and correction below...
  # prob_edge = (mean_degree*metapop_size) / n_possible_edges
  # XXX: for some reason adding a 2 here gets closer to desired mean degree.
  prob_edge = (2 * mean_degree * metapop_size) / n_possible_edges

  # Distribute total probability of an edge according to asymmetric homophilies.
  prob_min_teach_min = prob_edge * ((1 + homophily_min) / 2.0)
  prob_maj_teach_min = prob_edge * ((1 - homophily_min) / 2.0)

  prob_maj_teach_maj = prob_edge * ((1 + homophily_maj) / 2.0)
  prob_min_teach_maj = prob_edge * ((1 - homophily_maj) / 2.0)

  # Initialize block matrix with four entries for each dyad possibility.
  B = matrix(NA, nrow = 2, ncol = 2) 

  # Helper function for calculating log-odds from edge probabilities.
  log_odds = function(prob) { return (log(prob / (1 - prob))) }

  # Fill block matrix for stochastic block model.
  B[1, 1] = log_odds(prob_min_teach_min); B[1, 2] = log_odds(prob_min_teach_maj)
  B[2, 1] = log_odds(prob_maj_teach_min); B[2, 2] = log_odds(prob_maj_teach_maj)

  # Set the number of grouping variables to 1. 
  V = 1

  # Set group membership to be either Minority or Majority after calculating sizes.
  min_size = min_frac * metapop_size
  maj_size = metapop_size - min_size
  group = c(rep(1, min_size), rep(2, maj_size))

  # Generate network using stochastic block model.
  sim_dat = simulate_sbm_network(N_id = metapop_size, 
                                 B = list(B = B), 
                                 V = V, 
                                 groups = groups_df(group),
                                 individual_predictors = NULL,
                                 dyadic_predictors = NULL,
                                 individual_effects = NULL,
                                 dyadic_effects = NULL,        
                                 mode="bernoulli"
                                )

  return (sim_dat)
}


fit_sbm_to_minmaj_net = function(minmaj_net) {
  # 
}


# Default parameter settings for simulated min-maj data with block-only SBM fit.
fit_sbm_to_sbm_min_maj = function(sim_minmaj_dat, iter_warmup = 20, 
                                  iter_sampling = 20, chains = 4, 
                                  seed = NULL) {
  
  strand_minmaj_dat = make_strand_data(self_report = list(sim_minmaj_dat$network),
                                       block_covariates = groupsf_df(sim_minmaj_dat$group_ids$Group),
                                       individual_covariates = NULL,
                                       dyadic_covariates = NULL)

  fit <-
    fit_block_model(
      data = sim_strand_dat,
      block_regression = ~ Group,
      focal_regression = ~ 1,
      target_regression = ~ 1,
      dyad_regression = ~ 1,
      mode="mcmc",
      stan_mcmc_parameters = 
        list(seed = seed, chains = chains, parallel_chains = chains, 
             refresh = 1, iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
             max_treedepth = NULL, adapt_delta = NULL)
      )
  
  return (fit)

}


plot_network = function(net, min_frac) {
   
  # net = graph_from_adjacency_matrix(sim_dat$network, mode = c("directed"))

  # group_vec = sim_dat$group_ids$Group
  metapop_size = length(net)
  N_min = round(min_frac * metapop_size)
  N_maj = metapop_size - N_min

  group_vec = c(rep(1, N_min), rep(2, N_maj))

  V(net)$color = c("turquoise4", "goldenrod3")[group_vec]

  coords = matrix(NA, nrow = length(group_vec), ncol = 2)

  N_min = sum(group_vec == 1)
  N_maj = sum(group_vec == 2)

  N = length(group_vec)
  coords = matrix(0.0, nrow = N, ncol = 2)

  coords[group_vec == 1, 1] = -7.0
  coords[group_vec == 2, 1] = 7.0
  
  coords[group_vec == 1] = 
    coords[group_vec == 1] + cbind(rnorm(N_min, sd = .75), rnorm(N_min, sd = 1.5))

  coords[group_vec == 2] = 
    coords[group_vec == 2] + cbind(rnorm(N_maj, sd = 1.25), rnorm(N_maj, sd = 1.25))

  edges_df = as_data_frame(net)
  edge_color_options = c(adjustcolor("turquoise4", alpha.f = 0.8),
                         adjustcolor("pink", alpha.f = 0.215),
                         adjustcolor("darkgreen", alpha.f = 0.5), 
                         adjustcolor("goldenrod3", alpha.f = 0.5)
                         )

  edge_color_lookup = function(row) {
    color_code = NA
    if (group_vec[row$from] == 1 & group_vec[row$to] == 1) {
      color_code = 1
    } else if (group_vec[row$from] == 1 & group_vec[row$to] == 2) {
      color_code = 2
    } else if (group_vec[row$from] == 2 & group_vec[row$to] == 1) {
      color_code = 3
    } else {
      color_code = 4
    }

    return (edge_color_options[color_code])
  }

  edge_color = by_row(as_data_frame(net), edge_color_lookup, .collate = "rows")$.out

  p = plot(net, edge.arrow.size = 0.5, edge.width = 1.5, 
       edge.curved = 0.4, edge.color = edge_color, 
       xlim = c(-9, 11), ylim = c(-1.0, 1.0), vertex.label=NA, 
       vertex.size = 35, layout = coords, rescale = FALSE
      )

  return (p)
}


# Make a minority-majority network using construction from Turner, Reynolds, and
# Jones.
make_min_maj_net = function(metapop_size = 200, min_frac = 0.05, 
                            mean_degree = 10,
                            homophily_min = 0.0,
                            homophily_maj = 0.0) {

  n_edges = mean_degree * metapop_size
  E = n_edges

  E_min = round(n_edges * min_frac)
  E_maj = E - E_min

  E_min_teach_min = round(E_min * ((1 + homophily_min) / 2.0))
  E_maj_teach_min = E_min - E_min_teach_min

  E_maj_teach_maj = round(E_maj * ((1 + homophily_maj) / 2.0))
  E_min_teach_maj = E_maj - E_maj_teach_maj
  
  min_size = round(min_frac * metapop_size)
  maj_size = metapop_size - min_size
  group = c(rep(1, min_size), rep(2, maj_size))
  
  network <- make_graph(c(), n = metapop_size, directed = TRUE)

  minority_ids <- 1:min_size
  majority_ids <- (min_size + 1):metapop_size

  for (xx in 1:E_min_teach_min) {

    # Generate a new edge that may already exist in the social network.
    resample = TRUE
    while (resample) {
      new_edge = sample(minority_ids, 2)
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  for (xx in 1:E_maj_teach_maj) {
    
    resample = TRUE
    while (resample) {
      new_edge = sample(majority_ids, 2)
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  for (xx in 1:E_min_teach_maj) {
    
    resample = TRUE
    while (resample) {
      new_edge = c(sample(minority_ids, 1), sample(majority_ids, 1))
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  for (xx in 1:E_maj_teach_min) {
    
    resample = TRUE
    while (resample) {
      new_edge = c(sample(majority_ids, 1), sample(minority_ids, 1))
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  return (network)
}


visualize_minmaj_collection = function(metapop_size = 200, min_frac = 0.05,
                                       mean_degree = 10,
                                       h_min_maj_pairs = 
                                         rbind(c(0.0, 0.0),
                                               c(-0.4, 0.7),
                                               c(0.7, 0.7)),
                                       savedir = 
                                         file.path("figures", "minmaj_networks")
                                       ) {

  if (!dir.exists(savedir))
    dir.create(savedir, recursive = TRUE)

  for (pair_idx in 1:nrow(h_min_maj_pairs)) {

    h_min = h_min_maj_pairs[pair_idx, 1]
    h_maj = h_min_maj_pairs[pair_idx, 2]

    net = make_min_maj_net(metapop_size, min_frac, mean_degree, h_min, h_maj)

    savefile = file.path(
      savedir, paste0("h_min=", h_min, "_", "h_maj=", h_maj, ".pdf")
    )

    pdf(savefile, 7.5, 5)
      plot_network(net, min_frac)
    dev.off()

  }

  
}
