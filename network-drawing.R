################################################################################
##### How to draw networks in R.
##### Ecological Synthesis Lab (SintECO).
##### https://marcomellolab.wordpress.com.
##### Authors: Marco Mello & Renata Muylaert.
##### E-mail: marmello@usp.br. 
##### See README for further info.
##### https://github.com/marmello77/network-drawing/blob/master/README.md
################################################################################


################################################################################
##### SUMMARY
################################################################################


# 1. Get ready
# 2. Monolayer matrix in bipartite
# 3. Monolayer bipartite graph in bipartite
# 4. Monolayer energy-minimization graph in igraph
# 5. Monolayer energy-minimization graph with multiple vertex types in igraph
# 6. Multilayer energy-minimization graph in igraph


################################################################################
##### 1. Get ready
################################################################################


# Warning: There is no single magic way to draw all kinds of networks. 
# There are several network drawing algorithms implemented in different
# R packages and stand-alone software. Study their logic and algorithms, 
# see some papers in which they were used. Think it through, and only 
# then decide which drawing method to use in your study. For guidelines
# on which drawing algorithm to choose, see the suggested readings.

# Set the working directory automatically to the source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Remove all previous objects
rm(list= ls())

# Load the required packages
library(bipartite)
library(igraph)
library(reshape2)
library(tidyverse)

#Import the data in matrix format, ready to be used with bipartite
net1_bi <- as.matrix(read.delim("data/net1.txt", row.names=1))
net2_bi <- as.matrix(read.delim("data/net2.txt", row.names=1))
net3an_bi <- as.matrix(read.delim("data/net3an.txt", row.names=1))
net3mu_bi <- as.matrix(read.delim("data/net3mu.txt", row.names=1))

#Inspect the bipartite networks
class(net1_bi)
net1_bi
dim(net1_bi)
min(net1_bi)
max(net1_bi)

class(net2_bi)
net2_bi
dim(net2_bi)
min(net2_bi)
max(net2_bi)

class(net3an_bi)
net3an_bi
dim(net3an_bi)
min(net3an_bi)
max(net3an_bi)

class(net3mu_bi)
net3mu_bi
dim(net3mu_bi)
min(net3mu_bi)
max(net3mu_bi)

#Create versions of all networks formatted for igraph
net1_ig <- graph_from_incidence_matrix(net1_bi, directed = F, weighted = TRUE)
net2_ig <- graph_from_incidence_matrix(net2_bi, directed = F, weighted = TRUE) 
net3an_ig <- graph_from_incidence_matrix(net3an_bi, directed = F, weighted = TRUE) 
net3mu_ig <- graph_from_incidence_matrix(net3mu_bi, directed = F, weighted = TRUE) 

#Inspect the igraph networks
class(net1_ig)
net1_ig
E(net1_ig)
V(net1_ig)$name

class(net2_ig)
net2_ig
E(net2_ig)
V(net2_ig)$name

class(net3an_ig)
net3an_ig
E(net3an_ig)
V(net3an_ig)$name

class(net3mu_ig)
net3mu_ig
E(net3mu_ig)
V(net3mu_ig)$name


################################################################################
##### 2. Monolayer matrix in bipartite
################################################################################

# Let's start by drawing a matrix. Use net1_bi as an example. Take a look at it:
net1_bi

# Plot the network as a matrix and export the output as a PNG image
png(filename= "figures/net1_bipartite_matrix.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

visweb(net1_bi, #the network you want to plot
       type = "nested", #choose the drawing method
       textsize = 2, #choose vertex label size
       square = "interaction") #choose how the cells will be filled
# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 3. Monolayer bipartite graph in bipartite
################################################################################


# Let's continue with net1_bi Do you remember how it looks?
net1_bi

# Plot the network as a bipartite graph and export the output as a PNG image
png(filename= "figures/net1_bipartite_graph.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 2000, width= 3000) #set the dimensions in the chosen unit

plotweb(net1_bi, #the network you want to plot
        method ="cca", #set the drawing method
        col.low ="darkgreen", #set the color of the row vertices
        bor.col.low = "white", #set the color of the row vertices' borders
        col.high ="gold", #set the color of the column vertices
        bor.col.high = "white", #set the color of the column vertices' borders
        col.interaction ="grey90", #set edge color
        bor.col.interaction = "white", #set link border color
        text.rot ="90", #set the rotation of vertex labels
        labsize = 2) #set vertex label size
# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 4. Monolayer energy-minimization graph in igraph
################################################################################


# We'll continue with net1
# But this time we'll use its igraph version: net1_ig
# Take a look at it. Its heading is very informative, check it's meaning
# in this package's help
net1_ig

# Plot the network as an energy--minimization graph and
# export the output as a PNG image
png(filename= "figures/net1_igraph_graph.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plot(net1_ig, #the network you want to plot
     layout=layout_nicely, #choose a drawing layout for the graph
     vertex.shape = "circle", #choose a vertex shape
     vertex.size = 12, #set vertex size
     vertex.color = c("darkgreen", "gold"), #set vertex colors
     vertex.frame.color = NA, #set vertex border color
     vertex.label.cex = 0.7, #set vertex label size
     vertex.label.color = "white", #set vertex label color
     edge.width = E(net1_ig)$weight/100, #set edge width and adjust it
     edge.color = adjustcolor("grey", alpha = 0.5), #set edge color and adjust
     edge.curved=0.3) #set edge curvature
# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 5. Monolayer energy-minimization graph with multiple vertex types in igraph
################################################################################


# In this example, we are going to work with the network net2_ig
# It's a network composed of three taxonomic groups: bats, hawkmoths, and plants
# We want to differentiate bats and hawkmoths in the graph
# To do this, we have to add info on vertex types
# First, let's remember what the network looks like
net2_ig

# The first 4 rows represent bats, the other rows, hawkmoths
V(net2_ig)$name[1:4]

#Create a new vertex attribute with the taxonomic groups
V(net2_ig)$set = c(rep("Bats", 4),
                   rep("Moths", nrow(net2_bi)-4),
                   rep("Plants", ncol(net2_bi))
                   )

#Set vertex colors by taxonomic group
V(net2_ig)$color = V(net2_ig)$set
V(net2_ig)$color = gsub("Bats","gold",V(net2_ig)$color)
V(net2_ig)$color = gsub("Moths","purple",V(net2_ig)$color)
V(net2_ig)$color = gsub("Plants","darkgreen",V(net2_ig)$color)

# Plot the network as an energy--minimization graph and
# export the output as a PNG image
png(filename= "figures/net1_igraph_graph_taxon.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plot(net2_ig, #the network you want to plot
     layout = layout_nicely, #choose a drawing layout for the graph
     vertex.shape = "circle", #choose a vertex shape
     vertex.size = 12, #set vertex size
     vertex.color = V(net2_ig)$color, #set vertex colors
     vertex.frame.color = NA, #set vertex border color
     vertex.label.cex = 0.7, #set vertex label size
     vertex.label.color = "white", #set vertex label color
     edge.width = E(net1_ig)$weight/50, #set edge width and adjust it
     edge.color = adjustcolor("grey", alpha = 0.5), #set edge color and adjust
     edge.curved = 0.3) #set edge curvature
# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 6. Multilayer energy-minimization graph in igraph
################################################################################


# Now let's go for something more complex: a multilayer network

# The best way to work with multilayer networks is by using vertex and edge
# lists. However, in most cases, ecologists usually have incidence matrices
# to begin the analysis. So let's assume this starting point.

# We are going to work with a multilayer network composed of two layers:
# one with antagonistic interactions and the other with mutualistic 
# interactions. These layers are represented by two matrices with equal
# dimensions and exactly the same label order for the rows and columns.

# Let's take a look at those matrices, which we've already imported. See that
# they are like reflections in a mirror, but with different values in their
# cells.

# Antagonistic matrix:
net3an_bi

# Mutualistic matrix:
net3mu_bi

# Fortunately, we've already converted them into igraph format. Do you remember
# their attributes?
net3an_ig
net3mu_ig

# First, let's create a new attributed with info on interaction types
E(net3an_ig)$layer <- "antagonistic"
E(net3mu_ig)$layer <- "mutualistic"

# Check the layers
E(net3an_ig)$layer
E(net3mu_ig)$layer

# Now let's inform which vertices are animals and plants
V(net3an_ig)$taxon = ifelse(V(net3an_ig)$type == FALSE, "animals", "plants")
V(net3mu_ig)$taxon = ifelse(V(net3mu_ig)$type == FALSE, "animals", "plants")

# Check the taxonomic groups
V(net3an_ig)$taxon
V(net3mu_ig)$taxon

# Now let's export the edge lists (links)
net3an_links <- as.data.frame(as_edgelist(net3an_ig, names = TRUE))
net3mu_links <- as.data.frame(as_edgelist(net3mu_ig, names = TRUE))

# Check the edge lists
net3an_links
net3mu_links

# Add info on the layers
net3an_links$layer <- "antagonistic"
net3mu_links$layer <- "mutualistic"

# Give the columns informative names
colnames(net3an_links) <- c("animals", "plants", "layer")
colnames(net3mu_links) <- c("animals", "plants", "layer")

# Bring back the information on edge weights
net3an_links$weight <- E(net3an_ig)$weight
net3mu_links$weight <- E(net3mu_ig)$weight

# Now create vertex lists
net3an_nodes <- data.frame(node=character(length(V(net3an_ig)$name)), 
                           taxon=character(length(V(net3an_ig)$taxon)), 
                           stringsAsFactors=FALSE) 
net3an_nodes$node <- V(net3an_ig)$name
net3an_nodes$taxon <- V(net3an_ig)$taxon

net3mu_nodes <- data.frame(node=character(length(V(net3mu_ig)$name)), 
                           taxon=character(length(V(net3mu_ig)$taxon)), 
                           stringsAsFactors=FALSE) 
net3mu_nodes$node <- V(net3mu_ig)$name
net3mu_nodes$taxon <- V(net3mu_ig)$taxon

# Check the vertex lists
net3an_nodes
net3mu_nodes

# Now let's merge the two layers
net3_links = rbind(net3an_links, net3mu_links)
net3_nodes = rbind(net3an_nodes, net3mu_nodes)
net3_nodes = unique(net3_nodes)


##### ATTENTION: #####
# If your data are already organized as vertex and edge lists, you can begin
# right here, instead of having to do all previous operations.
######################


# Create a new multilayer network formatted for igraph
net3_multi <- graph_from_data_frame(d=net3_links, vertices=net3_nodes, directed=F) 
class(net3_multi)
net3_multi

# Add information on the bipartite structure
V(net3_multi)$type <- ifelse(V(net3_multi)$taxon == "animals", "animals", "plants")

# Check the network's attributes again to check if it's undirected, named, 
# weighted, and bipartite. Check also its number of vertices and edges
net3_multi

# As net2, this network does also contain 3 taxonomic groups: marsupials, 
# rodents, and plants. So let's recover this information now.
# Create a new vertex attribute with the taxonomic groups
V(net3_multi)$taxon2 = c(c("Rodents", "Rodents", "Marsupials", "Marsupials",
                        "Marsupials", "Rodents", "Rodents", "Marsupials",
                        "Rodents"),
                   rep("Plants", ncol(net3an_bi))
                   )

# Check the taxonomic groups
V(net3_multi)$taxon2

# Set vertex colors by taxonomic group
V(net3_multi)$color = V(net3_multi)$taxon2
V(net3_multi)$color = gsub("Marsupials","gold",V(net3_multi)$color)
V(net3_multi)$color = gsub("Rodents","purple",V(net3_multi)$color)
V(net3_multi)$color = gsub("Plants","darkgreen",V(net3_multi)$color)

# Check vertex colors
V(net3_multi)$color

# Set edge colors by layer
E(net3_multi)$color = E(net3_multi)$layer
E(net3_multi)$color = gsub("antagonistic", "red",E(net3_multi)$color)
E(net3_multi)$color = gsub("mutualistic","blue",E(net3_multi)$color)

# Check edge colors
E(net3_multi)$color

# Adjust a function to set edge curvatures, so the edges do not overlap
curves = curve_multiple(net3_multi)

# Plot the network as an energy--minimization graph and
# export the output as a PNG image
png(filename= "figures/net1_igraph_graph_multi.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plot(net3_multi, #the network you want to plot
     layout = layout_nicely, #choose a drawing layout for the graph
     vertex.shape = "circle", #choose a vertex shape
     vertex.size = 8, #set vertex size
     vertex.color = V(net3_multi)$color, #set vertex colors
     vertex.frame.color = NA, #set vertex border color
     vertex.label.cex = 0.5, #set vertex label size
     vertex.label.color = "white", #set vertex label color
     edge.width = E(net1_ig)$weight/100, #set edge width and adjust it
     edge.color = adjustcolor(E(net3_multi)$color, alpha = 0.3), #set edge color
     edge.curved = curves) #set edge curvature
# You can set many other arguments. Check this function's help

dev.off()


#################################### END #######################################
