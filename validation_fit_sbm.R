################################### Install and/or load STRAND
library(devtools)
install_github('ctross/STRAND')
library(STRAND)

######################################################### Create groups. 2 here, but it scales to larger K
library(igraph)
# set.seed(1)
V = 1            # One group variable
G = 2            # Two categories in this variable
N_id = 100       # Number of people

group = sample(1:2, N_id, replace=TRUE, prob=c(0.8, 0.2)) # Get group IDS


######################################################################## First let minority have low teaching influence to group 2
B = matrix(NA, nrow=G, ncol=G) # Tie weights on log-odds scale

B[1,1] = -2   # Group 1 to group 1
B[2,2] = -2.5   # Group 2 to group 2

B[1,2] = -5    # Group 1 teaching group 2
B[2,1] = -3   # Group 2 teaching group 1

A = simulate_sbm_network(N_id = N_id, B=list(B=B), V=V, groups = data.frame(group=factor(group)),
                         individual_predictor=matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1), 
                         individual_effects=matrix(c(0, 0),ncol=1, nrow=2),
                         mode="bernoulli"
)

Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids$group]

plot(Net, edge.arrow.size =0.3, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)

# dat = make_strand_data(outcome)

# fit = fit_block_model( data=model_dat,
#                        block_regression = ~ Ethnicity,
#                        focal_regression = ~ Age * NoFood,
#                        target_regression = ~ Age * NoFood,
#                        dyad_regression = ~ Relatedness + Friends * SameSex,
#                        mode="mcmc",
#                        stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, 
#                        refresh = 1, iter_warmup = 100, iter_sampling = 100,
#                        max_treedepth = NULL, adapt_delta = NULL)
#                        )



