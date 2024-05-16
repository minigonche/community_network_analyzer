# Script for simulations
library("MCMCpack")
library("SPRING")
library("igraph")

setwd("~/Dropbox/Projects/TartuU/community_network_analyzer")

# Parameters
num_otus = 100 # Number of OTUs
num_samples = 200 # Number of Samples
min_depth = 1000
max_depth = 2000

max_value <- round(0.01*max_depth)

stab = 1000 # Stability

alpha <- runif(num_otus, min = 0, max = 1)
alpha = alpha*stab
alpha[alpha < 1] = 1


# Sample a random vector from a Dirichlet distribution
depth <- runif(num_samples, min = min_depth, max = max_depth)
samples = round(rdirichlet(num_samples, alpha)*depth)

# Adds Error
error_random_matrix <- matrix(sample(0:max_value, num_samples * num_otus, replace = TRUE), nrow = num_samples, ncol = num_otus)
samples = samples+error_random_matrix

fit.spring <- SPRING(samples, Rmethod = "approx", quantitative = TRUE, 
                     lambdaseq = "data-specific", nlambda = 50, rep.num = 50, verbose = FALSE, seed = seed)

adj_matrix= as.matrix(fit.spring$fit$est$path[[opt.K]])
# Create a graph object from the adjacency matrix
graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# Find the optimal clustering
communities<-cluster_fast_greedy(graph)##

modularity<-modularity(communities,membership(communities))

# Print the maximum modularity
print(modularity)
nrow(data.frame(sizes(communities)))




