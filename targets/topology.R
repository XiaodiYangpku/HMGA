library('igraph')
argv <- commandArgs(TRUE)
input_fname <- argv[1]
output_fname <- argv[2]
data <- read.table(input_fname, head = TRUE, sep = "\t")
g <- graph.data.frame(data, directed=FALSE)

dg_value <- degree(g)
bw_value <- betweenness(g, directed = F,normalized = T) # use normalization
cl_value <- closeness(g, normalized = T)  
tr_value <- transitivity(g,type="local") # clustering coefficient
pg_value <- page_rank(g)$vector
ec_value <- eccentricity(g)
ei_value <- eigen_centrality(g)$vector

tr_value[is.na(tr_value)] <- 0 #  substitute NA with 0

forumG <- data.frame(V(g)$name, dg_value, bw_value, cl_value, tr_value, pg_value, ec_value, ei_value)
write.table(forumG, file=output_fname, quote=F, sep=" ", row.names=F, col.names = c('Uniprot', 'Degree', 'Betweenness', 'Closeness', 'Transitivity', 'Pagerank', 'Eccentricity', 'Eigenvector'))






