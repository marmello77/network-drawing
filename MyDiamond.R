MyDiamond <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], 
          bg=vertex.color, fg=vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("diamond", clip=shapes("circle")$clip,
          plot=MyDiamond, parameters=list(vertex.frame.color="white",
                                           vertex.frame.width=1))