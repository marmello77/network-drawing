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
# 2. Monolayer matrix
# 3. Monolayer bipartite graph
# 4. Monolayer energy-minimization graph
# 5. Monolayer energy-minimization graph multiple vertex types
# 6. Multilayer energy-minimization graph multiple vertex types
# 7. Monolayer matrix with modules
# 8. Monolayer energy-minimization graph with modules from igraph
# 9. Monolayer energy-minimization graph with modules from bipartite
# 10. Monolayer energy-minimization graph with modules +multiple vertex types


################################################################################
##### 1. Get ready
################################################################################


# Set the working directory automatically to the source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Remove all previous objects
rm(list= ls())

# Load the required packages
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(tidyverse)


# Warning: 

# There is no single magic way to draw all kinds of networks. 
# There are several network drawing algorithms implemented in different
# R packages and stand-alone software. Study their logic and algorithms, 
# see some papers in which they were used. Think it through, and only 
# then decide which drawing method to use in your study. For guidelines
# on which drawing algorithm to choose, see the suggested readings.

# The elements of a network may be called nodes or vertices. Depending on the
# software or package used, you'll see one term or the other.

# The relationships between the elements of a network may be called links or
# edges. Depending on the software or package used, you'll see one term or
# the other.


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
net1_ig <- graph_from_incidence_matrix(net1_bi, 
                                       directed = F, weighted = TRUE)
net2_ig <- graph_from_incidence_matrix(net2_bi, 
                                       directed = F, weighted = TRUE) 
net3an_ig <- graph_from_incidence_matrix(net3an_bi, 
                                         directed = F, weighted = TRUE) 
net3mu_ig <- graph_from_incidence_matrix(net3mu_bi, 
                                         directed = F, weighted = TRUE) 

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
##### 2. Monolayer matrix
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
##### 3. Monolayer bipartite graph
################################################################################


# Let's continue with net1_bi. Do you remember how it looks?
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
##### 4. Monolayer energy-minimization graph
################################################################################


# We'll continue using net1
# But this time we'll use its igraph version: net1_ig
# Take a look at it. Its heading is very informative. Check the meaning of all
# those codes in this package's help
net1_ig

# Plot the network as an energy-minimization graph and
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
##### 5. Monolayer energy-minimization graph multiple vertex types
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

# Plot the network as an energy-minimization graph and
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
##### 6. Multilayer energy-minimization graph
################################################################################


# Now let's go for something more complex: a multilayer network

# The best way to work with multilayer networks is by using vertex and edge
# lists. However, in most cases, ecologists usually have incidence matrices
# to begin the analysis. 

# If your data are organized as edge and vertex lists, you life will be much
# easier. See another script focused on multilayer networks:
# https://github.com/marmello77/multilayer-networks

# Otherwise, proceed to the next step in this script. You'll need to transform
# your matrices.

# OK, we are going to work with a multilayer network composed of two layers:
# one with antagonistic interactions and the other with mutualistic 
# interactions. These layers are represented by two matrices with equal
# dimensions and exactly the same label order for rows and columns.

# Let's take a look at the matrices, which we've already imported. See that
# they are like reflections in a mirror, but with different values in their
# cells.

# Antagonistic matrix:
net3an_bi

# Mutualistic matrix:
net3mu_bi

#Transform these matrices into a combined edge list
net3list <- bind_rows(
    as.data.frame(as.table(net3an_bi)),
    as.data.frame(as.table(net3mu_bi)),
    .id = "layer") %>%
    
    filter(Freq != 0) %>%
    select(
        from = Var1,
        to = Var2,
        layer,
        Freq)

#Check the data
head(net3list)

#Give the columns informative names
colnames(net3list) <- c("animals", "plants", "layer", "weight")

#Check the data
head(net3list)

#Transform the combined edge list into an igraph object
net3_multi <- graph_from_data_frame(net3list, directed = FALSE)

#Check the multilayer network
net3_multi
V(net3_multi)
E(net3_multi)
attributes(V(net3_multi))
attributes(E(net3_multi))
V(net3_multi)$type
V(net3_multi)$name
E(net3_multi)$layer
E(net3_multi)$weight

# Create a new edge attribute with info on interaction types (layers)
E(net3_multi)$layer <- ifelse(E(net3_multi)$layer == "1",
                              "antagonistic",
                              "mutualistic")

# Check the layers
E(net3_multi)$layer

# Add information on the bipartite structure by assigning vertex classes
V(net3_multi)$type = c(rep(0, nrow(net3an_bi)), 
                       rep(1, ncol(net3an_bi)))

# Check the network's attributes and the vertex classes
net3_multi
V(net3_multi)$type

# Like net2, this network also contains 3 taxonomic groups. In this case,
# marsupials, rodents, and plants. So let's recover this information now.
# Create a new vertex attribute with the taxonomic groups
V(net3_multi)$taxon = c(c("Rodents", "Rodents", "Marsupials", "Marsupials",
                           "Marsupials", "Rodents", "Rodents", "Marsupials",
                           "Rodents"),
                         rep("Plants", ncol(net3an_bi)))

# Check the taxonomic groups
V(net3_multi)$taxon

# Set vertex colors by taxonomic group
V(net3_multi)$color = V(net3_multi)$taxon
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


################################################################################
##### 7. Monolayer matrix with modules
################################################################################


# What if you wanted to represent the modular structure of a network? You can do
# it in different ways. The first is by plotting it as a matrix.

# Let's continue with net1_bi. Do you remember how it looks?
net1_bi

# First you need to analyze the network's modularity. This involves running a
# modularity algorithm of your choice to calculate the degree of modularity Q 
# and identify which vertices belong to which modules.

#Calculate modularity
net1_bi.mod <- computeModules(net1_bi,
                              method = "Beckett") #DIRTLPAwb+ algorithm

#Check the modularity information
net1_bi.mod

# Check module membership
listModuleInformation(net1_bi.mod)

# Plot the network as a matrix and export the output as a PNG image
png(filename= "figures/net1_bipartite_matrix_modules.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plotModuleWeb(net1_bi.mod,
              weighted = T,
              labsize = 1,
              displayAlabels = T,
              displayBlabels = T,
              xlabel = "Bees",
              ylabel = "Flowers",
              square.border	= "white")

# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 8. Monolayer energy-minimization graph with modules from igraph
################################################################################


# We'll use net1 formatted for igraph
net1_ig

# First you need to analyze the network's modularity. This involves running a
# modularity algorithm of your choice to calculate the degree of modularity Q 
# and identify which vertices belong to which modules.

# You can run a modularity analysis in different R packages. There are several
# algorithms implemented in igraph. And there are algorithms implemented in
# bipartite, which are custom-made for weighted, bipartite networks, such as
# net1, our focus here.

# Here we'll use a very popular algorithm, which can be calculate for different
# kinds of networks, both unipartite or bipartite, binary or weighted: the
# Louvain algorithm (Blondel 2008, JStatMech). This algorithm is very fast
# and stable, so it's the primary choice for very large networks.

# Calculate modularity in igraph
net1_ig.mod_ig <- cluster_louvain(net1_ig) #Louvain algorithm

# Check the modularity information
net1_ig.mod_ig

# Check module membership
net1_ig.mod_ig$membership

# We'll plot the network's modules calculated with igraph.
# That's pretty straightforward, as you'll see. You can use the output produced
# before as a plotting input. This way, you can plot the modules with different
# node colors and also draw some clouds around them.

# Plot the network as an energy-minimization graph and
# export the output as a PNG image
png(filename= "figures/net1_igraph_graph_modules_ig.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plot(net1_ig.mod_ig, #use the modularity output as the first input
     net1_ig, #the network you want to plot
     layout = layout_nicely, #choose a drawing layout for the graph
     vertex.shape = "circle", #choose a vertex shape
     vertex.size = 12, #set vertex size
     vertex.frame.color = NA, #set vertex border color
     vertex.label.cex = 0.7, #set vertex label size
     vertex.label.color = "white", #set vertex label color
     edge.width = E(net1_ig)$weight/100, #set edge width and adjust it
     edge.curved = 0.3) #set edge curvature

# In this plot, black edges represent within-module connections, while red edges
# represent between-module connections.

# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 9. Monolayer energy-minimization graph with modules from bipartite
################################################################################


# We'll continue using net1
# We need both its bipartite and igraph versions. Take another look at them
# to remember their structure
net1_bi
net1_ig

# First you need to analyze the network's modularity. This involves running a
# modularity algorithm of your choice to calculate the degree of modularity Q 
# and identify which vertices belong to which modules.

# You can run a modularity analysis in different R packages. There are several
# algorithms implemented in igraph. And there are algorithms implemented in
# bipartite, which are custom-made for weighted, bipartite networks, such as
# net1, our focus here.

# Here we'll use a relatively new algorithm  made for weighted, bipartite
# networks: DIRTLPAwb+ (Beckett 2016, RoyalScoOpSci). This is currently the
# best choice for working with this kind of network, the most common among
# studies on ecological interactions.

#Calculate modularity in bipartite
net1_bi.mod_bi <- computeModules(net1_bi,
                                 method = "Beckett") #DIRTLPAwb+ algorithm

# Check the modularity information
net1_bi.mod_bi

# Check module membership
listModuleInformation(net1_bi.mod_bi)

# How much do the Q-scores differ between Louvain and DIRTLPAwb+? Are there
# differences in the outputs produced by these two algorithms? Let's take a
# look for the sake of curiosity.
net1_bi.mod_ig$modularity
net1_bi.mod_bi@likelihood

# Let's get back to the point.
# Plotting the output of a modularity analysis produced in bipartite together
# with an igraph object is not so straightforward. The plotting function of
# igraph accepts only lists to represent module clouds, not vectors of data
# frames. Therefore, we need to do some data transformation first. 

# Save module membership as a data frame
net1_bi.modules.bi <- module2constraints(net1_bi.mod_bi)
net1_bi.modules.bi
net1_bi.df.bi <- data.frame(c(rownames(net1_bi), colnames(net1_bi)),
                            net1_bi.modules.bi) 
net1_bi.df.bi
colnames(net1_bi.df.bi) <- c("vertices", "modules")
net1_bi.df.bi

# Save module membership as a list to be inputted in igraph
net1_bi.list.bi <- split(net1_bi.df.bi$vertices, net1_bi.df.bi$modules)
net1_bi.list.bi

##Set node and cloud colors by modularity
net1_ig.colors <- rainbow(length(net1_bi.list.bi), alpha = 1.0, s = 1, v = 0.8)
V(net1_ig)$color <- net1_ig.colors[net1_bi.df.bi$modules]
net1_ig.clouds = net1_ig.colors

# Plot the network as an energy-minimization graph and
# export the output as a PNG image
png(filename= "figures/net1_igraph_graph_modules_bi.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plot(net1_ig, #the network you want to plot
     layout=layout_nicely, #choose a drawing layout for the graph
     col = V(net1_ig)$color, #set node color
     mark.groups = net1_bi.list.bi, #set cloud color
     mark.border = net1_ig.clouds, #set cloud border color as null
     mark.col = adjustcolor(net1_ig.clouds, alpha = 0.2), #set cloud color
     vertex.shape = "circle", #choose a vertex shape
     vertex.size = 12, #set vertex size
     vertex.frame.color = NA, #set vertex border color
     vertex.label.cex = 0.7, #set vertex label size
     vertex.label.color = "white", #set vertex label color
     edge.width = E(net1_ig)$weight/100, #set edge width and adjust it
     edge.curved=0.3) #set edge curvature
# You can set many other arguments. Check this function's help

dev.off()


################################################################################
##### 10. Monolayer energy-minimization graph with modules +multiple vertex types
################################################################################


# We'll use net2 formatted for igraph
net2_ig

# Calculate modularity in igraph
net2_ig.mod_ig <- cluster_louvain(net2_ig) #Louvain algorithm

# Check the modularity information
net2_ig.mod_ig

# Check module membership
net2_ig.mod_ig$membership

##Set node and cloud colors by modularity
net2_ig.colors <- rainbow(length(unique(net2_ig.mod_ig$membership)),
                          alpha = 1.0, s = 1, v = 0.8)
V(net2_ig)$color <- net2_ig.colors[net2_ig.mod_ig$membership]
net2_ig.clouds = net2_ig.colors

#Create a new vertex attribute with the taxonomic groups
V(net2_ig)$set = c(rep("Bats", 4),
                   rep("Moths", nrow(net2_bi)-4),
                   rep("Plants", ncol(net2_bi))
                   )

# As colors will represent modules, we'll use shape to represent vertex types.
# The default shapes are not so nice, so let's add a customized shape by using
# a drawing function. This is a modification of the "mytriangle" function
# provided in the igraph manual.

#Import the "diamond" vertex shape
source("MyDiamond.R")

#Set vertex shapes
V(net2_ig)$shape <- V(net2_ig)$set
V(net2_ig)$shape <- gsub("Bats","diamond",V(net2_ig)$shape)
V(net2_ig)$shape <- gsub("Moths","square",V(net2_ig)$shape)
V(net2_ig)$shape <- gsub("Plants","circle",V(net2_ig)$shape)

# We'll plot the network's modules calculated with igraph.
# That's pretty straightforward, as you'll see. You can use the output produced
# before as a plotting input. This way, you can plot the modules with different
# node colors and also draw some clouds around them.

# Plot the network as an energy-minimization graph and
# export the output as a PNG image
png(filename= "figures/net2_igraph_graph_modules_ig.png", #name the file
    units = "px", #set the units to pixels
    res= 300, #set the resolution in dpi (most journals require at least 300)
    height= 3000, width= 3000) #set the dimensions in the chosen unit

plot(net2_ig.mod_ig, #use the modularity output as the first input
     net2_ig, #the network you want to plot
     layout=layout_nicely, #choose a drawing layout for the graph
     col = V(net2_ig)$color, #set node color
     vertex.shape = V(net2_ig)$shape, #choose a vertex shape
     vertex.size = 9, #set vertex size
     vertex.label.cex = 0.5, #set vertex label size
     vertex.label.color = "white", #set vertex label color
     vertex.frame.color = NA, #set vertex border color
     edge.width = E(net2_ig)$weight/2, #set edge width and adjust it
     edge.curved=0.3) #set edge curvature

# In this plot, black edges represent within-module connections, while red edges
# represent between-module connections.

# You can set many other arguments. Check this function's help

dev.off()
