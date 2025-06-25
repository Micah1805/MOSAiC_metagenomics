################################################################################
####################Create a Gene-Cluster network###############################
################################################################################

####Housekeeping
library(igraph)                                                                 #Network analysis
library(RCy3)                                                                   #Communication with Cytoscape
library(readr)                                                                  #Read in data
library(tidyverse)


####Create a table containing genes, clusters, time, depth and sample

#Load Gene_Cluster table
gene.DETECTION.COVERAGE_no_ice <- read.delim("//wsl.localhost/Ubuntu/home/micah/Metagenomics_MOSAiC/mosaic-gene-clusters/02_GENE_DET_COV_SUMMARY/gene-DETECTION-COVERAGE_no_ice.txt", header=FALSE)

colnames(gene.DETECTION.COVERAGE_no_ice)[1] <- "Gene_id"
colnames(gene.DETECTION.COVERAGE_no_ice)[2] <- "Cluster_representative "
colnames(gene.DETECTION.COVERAGE_no_ice)[3] <- "Detection"
colnames(gene.DETECTION.COVERAGE_no_ice)[4] <- "Coverage"
colnames(gene.DETECTION.COVERAGE_no_ice)[5] <- "Sample"

gene.DETECTION.COVERAGE_no_ice <- gene.DETECTION.COVERAGE_no_ice %>%
  relocate("Sample")

str(gene.DETECTION.COVERAGE_no_ice)

#Load metadata
Genes_Depth_Time <- read.csv("//wsl.localhost/Ubuntu/home/micah/Metagenomics_MOSAiC/mosaic-gene-clusters/02_GENE_DET_COV_SUMMARY/Subsets/Genes_Depth_Time.csv", header=FALSE, sep=";")

colnames(Genes_Depth_Time)[1] <- "Sample"
colnames(Genes_Depth_Time)[2] <- "Depth"
colnames(Genes_Depth_Time)[3] <- "Date"

#Create a table containing Genes, clusters, depth and time 
Gene_Cluster_table  <- inner_join(gene.DETECTION.COVERAGE_no_ice,Genes_Depth_Time,by="Sample")

Gene_Cluster_table <- Gene_Cluster_table %>% 
  relocate("Gene_id","Cluster_representative.","Sample","Depth","Date")

Gene_Cluster_table <- Gene_Cluster_table[1:5]

write.csv(Gene_Cluster_table , 
          "~/AWI/MOSAiC/Networks/Data_ananalysis/Gene_Cluster_table.CSV",
          row.names = FALSE)

#Gene_Cluster_table <- read.csv("~/AWI/MOSAiC/Networks/Data_ananalysis/Gene_Cluster_table.CSV")
Subset_mosaic_water <- Gene_Cluster_table[1:1000,1:5]


####Create a network
#Network showing distribution of genes in clusters
# Only use unique gnes
gene_nodes <- Subset_mosaic_water %>%
  distinct(Gene_id, .keep_all = TRUE) %>%
  mutate(
    id = Gene_id,
    type = "gene"
  ) %>%
  select(id, type, Sample, Date, Depth, Cluster_representative.)

# Cluster-nodes
cluster_nodes <- gene_nodes %>%
  distinct(Cluster_representative.) %>%
  mutate(
    id = paste0("CLUSTER_", Cluster_representative.),
    type = "cluster"
  ) %>%
  select(id, type)

# Merge nodes
nodes <- bind_rows(gene_nodes %>% select(id, type, Sample, Date, Depth, Cluster_representative.),
                   cluster_nodes) %>%
  distinct(id, .keep_all = TRUE)

# Edges: Gene → Cluster
edges <- gene_nodes %>%
  mutate(
    source = id,
    target = paste0("CLUSTER_", Cluster_representative.)
  ) %>%
  select(source, target) %>%
  distinct()

# Calculate network
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Send network to Cytoscape
cytoscapePing()
createNetworkFromIgraph(g, "Gene-Cluster-Netzwerk")

# Transfer attributes (sample_id, region etc.) into Cytoscape
loadTableData(nodes, data.key.column = "id", table = "node")





####Network showing impact of time and water depth on genes 

#Gene_Cluster_table <- read.csv("~/AWI/MOSAiC/Networks/Data_ananalysis/Gene_Cluster_table.CSV")
Subset_mosaic_water <- Gene_Cluster_table[1:1000,1:5]

#Define cluster states (space-time specific)
#Each cluster in a sample at a specific time/space becomes an independent node

Subset_mosaic_water$cluster_state <- paste(Subset_mosaic_water$Cluster_representative., Subset_mosaic_water$Date, 
                                           Subset_mosaic_water$Sample, 
                                           Subset_mosaic_water$Depth, 
                                           sep = "_")

# Optional: Sortierung nach Zeit (für die Kreispositionierung später)
#Subset_mosaic_water$Date <- factor(Subset_mosaic_water$Date, levels = sort(unique(Subset_mosaic_water$Date)))  # explizite Zeitachsen-Reihenfolge


cluster_states <- unique(Subset_mosaic_water[, c("Cluster_representative.", "Date", "cluster_state")])


# Sort after cluster and time
cluster_states <- cluster_states[order(cluster_states$Cluster_representative., cluster_states$Date), ]

# Edges between times following on eachother per cluster
time_edges <- data.frame(from = character(), to = character(), interaction = character())

library(dplyr)

cluster_groups <- split(cluster_states, cluster_states$Cluster_representative.)

for(cl in names(cluster_groups)) {
  grp <- cluster_groups[[cl]]
  if(nrow(grp) > 1) {
    froms <- head(grp$cluster_state, -1)
    tos <- tail(grp$cluster_state, -1)
    time_edges <- rbind(time_edges,
                        data.frame(from = froms, to = tos, interaction = "temporal_followup"))
  }
}



# Gene → cluster-state (Affiliation)
edges_gc <- data.frame(
  from = Subset_mosaic_water$Gene_id,
  to = Subset_mosaic_water$cluster_state,
  interaction = "member_of"
)

# Gene → Sample (Origin)
edges_gp <- data.frame(
  from = Subset_mosaic_water$Gene_id,
  to = Subset_mosaic_water$Sample,
  interaction = "originates_from"
)

# Merge edges
edges_all <- rbind(edges_gc, 
                   edges_gp,
                   time_edges)


#Define nodes and sort to type
nodes <- unique(c(edges_all$from, edges_all$to))
node_df <- data.frame(
  id = nodes,
  type = ifelse(nodes %in% Subset_mosaic_water$Gene_id, "gene",
                ifelse(nodes %in% Subset_mosaic_water$cluster_state, "cluster_state",
                       ifelse(nodes %in% Subset_mosaic_water$Sample, "sample", "unknown")))
)



# Build network
g <- graph_from_data_frame(edges_all, vertices = node_df, directed = TRUE)

# Send to Cytoscape
createNetworkFromIgraph(g, "Cluster-Cyrcle_time")

# Export node attributes
loadTableData(node_df, data.key.column = "id", table = "node")















# Erstelle eindeutige Knoten pro Gen-Probe
nodes <- Subset_mosaic_water %>%
  mutate(
    node_id = paste(Gene_id, Sample, sep = "_")
  ) %>%
  select(node_id, Gene_id, Cluster_representative., Sample, Date, Depth)

# Kanten: Optional – verbinde gleiche Gene in verschiedenen Proben
edges <- nodes %>%
  select(Gene_id, node_id) %>%
  inner_join(., ., by = "Gene_id") %>%
  filter(node_id.x != node_id.y) 

%>%
  distinct() %>%
  rename(source = node_id.x, target = node_id.y)

# Netzwerk erstellen
g <- graph_from_data_frame(edges, vertices = nodes %>% rename(id = node_id), directed = FALSE)

# Mit Cytoscape verbinden
cytoscapePing()
createNetworkFromIgraph(g, "Gene-Probe-Cluster Network")

# Node-Attribute übertragen
loadTableData(nodes, data.key.column = "node_id", table = "node")



