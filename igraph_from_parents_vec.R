library(igraph)
library(stringr)



vals <- str_split(str_replace_all('0  0  0  0  0  0  0  0  0  0  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  6
7  8  9 10 11 12 13 14 15 16 17 18 19 20  1  2  3  4  5  6  7  8  9 10 11 12
13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38
39 40', pattern = '\n', replacement = ' '),' ')[[1]]
parent_vec <- integer()
for(val in vals){
  if(grepl(pattern = '\\S', x = val)){
    parent_vec <- c(parent_vec, as.integer(val))
  }
}


# test_parent_mat <- matrix(data = 0, nrow = )

# as.dendrogram(parent_vec)

adj_mat <- matrix(data = 0, nrow = length(parent_vec), ncol = length(parent_vec))
edge_from <- c()
edge_to <- c()
vertex_sizes = rep(1, length(parent_vec)) # so that even cells that aren't a parent have size 1
for(i in seq_len(length(parent_vec))){
  if(parent_vec[i] > 0){
    adj_mat[parent_vec[i], i] <- 1
    edge_from <- c(edge_from, parent_vec[i])
    edge_to <- c(edge_to, i)
    vertex_sizes[parent_vec[i]] <- vertex_sizes[parent_vec[i]] + 1
  }
}
edge_df <- data.frame(cbind(edge_from, edge_to))
colnames(edge_df) <- c('from', 'to')

lineage <- graph_from_adjacency_matrix(adj_mat)
plot(lineage, edge.arrow.size = 0.1, vertex.size = vertex_sizes)


init_pop_size <- 10
num_generations <- 10
big_vals <- rep(0, init_pop_size)
for(gen in 1:num_generations){
  big_vals <- append(big_vals, seq(1, init_pop_size*2^(gen-1)))
}
big_adj_mat <- matrix(data = 0, nrow = length(big_vals), ncol = length(big_vals))
big_edge_from <- c()
big_edge_to <- c()
big_vertex_sizes = rep(1, length(big_vals))
for(i in seq_len(length(big_vals))){
  if(big_vals[i] > 0){
    big_adj_mat[big_vals[i], i] <- 1
    big_vertex_sizes[big_vals[i]] <- big_vertex_sizes[big_vals[i]] + 1
    big_edge_from <- c(big_edge_from, big_vals[i])
    big_edge_to <- c(big_edge_to, i)
  }
}
big_edge_df <- data.frame(cbind(big_edge_from, big_edge_to))
colnames(big_edge_df) <- c('from', 'to')
big_lineage <- graph_from_adjacency_matrix(big_adj_mat)
plot(big_lineage, edge.arrow.size = 0.1, vertex.size = big_vertex_sizes, 
     vertex.label = NA)

