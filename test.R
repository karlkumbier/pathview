library(pathview)
library(tidyverse)
library(visNetwork)
load('test.Rdata')

p <- pathview_graphviz(
    gene.data=gene.data,
    pathway.id=pathway.id,
    species="hsa",
    out.suffix="image",
    kegg.dir="kegg_xml",
    low=low,
    mid=mid,
    high=high,
    kegg.native=FALSE,
    split.group=TRUE
)


node.data <- p$nd 
nodes <- do.call(cbind, p$nodeAttrs) %>%
  as.data.frame() %>%
  dplyr::rename(id=label) %>%
  mutate(shape=unlist(node.data$shape[id])) %>%
  mutate(label=unlist(node.data$labels[id])) %>%
  mutate(label=str_replace_all(label, ' ', '\n')) %>%
  mutate(genes=sapply(node.data$kegg.names[id], str_c, collapse=', ')) %>%
  mutate(title=str_c(label, ':', genes)) %>%
  mutate(label=str_sub(label, 1, 10)) %>%
  mutate(height=unlist(node.data$height[id])) %>%
  mutate(width=unlist(node.data$width[id])) %>%
  mutate(x=node.data$x[id] * 2) %>%
  mutate(y=node.data$y[id] * 2)

edge.ids <- names(p$eAttrs$label)
edges <- do.call(cbind, p$eAttrs) %>%
  as.data.frame() %>%
  mutate(to=str_remove_all(edge.ids, '^.*~')) %>%
  mutate(from=str_remove_all(edge.ids, '~.*$')) 

visNetwork(nodes, edges, width="1201px", height="1200px") %>%
  visNodes(fixed=TRUE)

