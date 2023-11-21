library(pathview)
load('test.Rdata')

pathview(
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

keggview.graph(
    gene.data=gene.data,
    pathway.id=pathway.id,
    species="hsa",
    out.suffix="image",
    kegg.dir="kegg_xml",
    low=low,
    mid=mid,
    high=high,
)

