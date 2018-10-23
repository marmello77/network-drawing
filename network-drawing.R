############################################################################
#                                                                          # 
#                SCRIPT FOR DRAWING NETWORKS IN R                          #
#                                                                          # 
############################################################################

##### Ecological Synthesis Lab (SintECO)
##### https://marcomellolab.wordpress.com
##### Author: Marco Mello
##### E-mail: marmello@gmail.com 
##### Script: Script for drawing networks in R
##### How to cite: Mello MAR. 2017. Script for drawing networks in R. Available
##### at https://marcomellolab.wordpress.com.
##### Published on April 25th, 2017 (English version).
##### Updated on November 22nd, 2017 (English version).
##### Run in R 3.3.3 (2017-03-06) -- "Another Canoe"

##### Disclaimer: You may use this script freely for non-comercial
##### purposes at your own risk. We assume no responsibility or
##### liability for the use of this software, convey no license
##### or title under any patent, copyright, or mask work right
##### to the product. We reserve the right to make changes in
##### the software without notification. 
##### We also make no representation or warranty that such
##### application will be suitable for the specified use without
##### further testing or modification. If this script helps you
##### produce any academic work (paper, book, chapter, 
##### dissertation etc.), please acknowledge the authors and
##### cite the source.


#############################################################


####SUMMARY####

# 1. Package bipartite
# 2. Package igraph
# 3. Suggested readings


#############################################################


# PREPARATION

# Load the packages

library(bipartite)
library(igraph)
library(reshape2)
library(ggplot2)
library(rstudioapi)


# Set the working directory automatically to the source file location 

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )


# Warning 1: this script works both with binary and weighted networks.

# Warning 2: this script was designed for two-mode (bipartite) networks,
# but some functions work also for one-mode (unipartite) networks.

# Warning 3: there is no single magic way to draw all kinds of network.
# There are several network drawing methods implemented in R packages
# and stand-alone software. Study their logic and algorithms,
# see some papers in which they were used, think it through,
# and only then decide which drawing method to use in your study.


#############################################################


####1.PACKAGE BIPARTITE####

# Drawing mode: bipartite graph

# Create the object to be analyzed.
# It should be a two-mode matrix formatted for bipartite.
# The format is actually very simple:
# A tab-delimited TXT matrix with row and column labels.
# See the example: a network (net1) analyzed by 
# Bezerra et al. 2009: http://dx.doi.org/10.1111/j.1365-2656.2009.01567.x.

net1 <- read.table("net1.txt", head=TRUE)

# Check the object

net1

#Draw and export the image

png(filename= "net1_bipartitegraph.png", 
    
    #Set image resolution in dpi
    res= 300, 
    
    #Set image dimensions
    height= 2500, width= 3000)

plotweb(net1, 
        
        #Set the drawing method. 
        method="cca", 
        
        #Set the color of the row nodes
        col.low="lightblue",
        
        #Set the color of the column nodes
        col.high="pink", 
        
        #Set the link color
        col.interaction="grey90", 
        
        #Set the rotation of node labels
        text.rot="90", 
        
        #Set the size of node labels
        labsize=1)
dev.off()

# Tip: There are many other parameters that can be set to customize the drawing.
# Explore them!


##################


# Drawing mode: bipartite matrix

# Use the same object as before

net1

# Draw and export the image

png(filename= "net1_bipartite_matrix.png", 
    
    # Set image resolution in dpi
    res= 300, 
    
    # Set image dimensions
    height= 2500, width= 3000)

visweb(net1, 
       
       #Set the drawing mode.
       type="diagonal",
       
       #Set the cell fills.
       square="interaction")
dev.off()

# Tip: There are many other parameters that can be set to customize the drawing.
# Explore them!



#############################################################



####2.PACKAGE IGRAPH####


# Transform the previous bipartite object into an igraph object

net2 = cbind.data.frame(reference=row.names(net1),net1)
str(net2)
netlist = melt(net2, na.rm = T)
colnames(netlist) = c("rows", "columns", "weight")
netlist[,1]=as.character(paste(netlist[,1]))
netlist[,2]=as.character(paste(netlist[,2]))
netlist2 <- subset(netlist, weight > 0)
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)
links = netlist2
net3 <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
class(net3)

#Check the main properties of the igraph object
V(net3)
E(net3)

# Draw and export the image

png(filename= "net1_igraph.png",
    
    # Set image resolution in dpi
    res= 300, 
    
    # Set image dimensions
    height= 2500, width= 3000)

plot.igraph(net3,
            
            # Set the drawing mode.
            # This package contains several drawing methods; try them!
            layout=layout_nicely,
            
            # Set node shapes
            vertex.shape = "circle",
            
            # Set node sizes
            vertex.size = 12,
            
            # Set link width proportional to link weights
            # You can transform the values, if they are too different or too large
            edge.width = E(net3)$weight/100,
            
            # Set node colors
            vertex.color = c("green", "yellow"),
            
            # Set link colors
            edge.color = "lightblue",
            
            # Set link curvature from 0 to 1
            edge.curved=0.3
            )
dev.off()

# Tip: There are many other parameters that can be set to customize the drawing.
# Explore them!



#############################################################



####3. SUGGESTED READINGS####


#Barabasi, A.L. (2016) Network Science, 1st ed. Cambridge 
 #University Press, Cambridge. Available at:
 #http://barabasi.com/networksciencebook/.

#Bascompte, J. & Jordano, P. (2014) Mutualistic Networks, 1st ed.
 #Princeton University Press, Princeton.

#Mello MAR, Muylaert RL, Pinheiro RBP & Félix GMF. 2016. 
 #Guia para análise de redes ecológicas. Edição dos autores, 
 #Belo Horizonte. 112 p. ISBN-13: 978-85-921757-0-2.
 #Available at: www.marcomello.org

#Ognyanova K. 2017. Static and dynamic network visualization with R.
 #Available at: http://kateto.net/network-visualization.

#Pocock, M. J. O., D. M. Evans, C. Fontaine, M. Harvey, R. Julliard,
 #Ó. McLaughlin, J. Silvertown, A. Tamaddoni-Nezhad, P. C. L. White,
 #and D. A. Bohan. 2016. The Visualisation of Ecological Networks, 
 #and Their Use as a Tool for Engagement, Advocacy and Management. Pages 41–85.


