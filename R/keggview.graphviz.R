keggview.graphviz <- function(
                           plot.data.gene=NULL,
                           plot.data.cpd=NULL,
                           cols.ts.gene=NULL,
                           cols.ts.cpd=NULL,
                           node.data,
                           path.graph,
                           pathway.name,
                           out.suffix="pathview",
                           pdf.size=c(7,7),
                           multi.state=TRUE,
                           same.layer=TRUE,
                           match.data=TRUE,
                           rankdir=c("LR","TB")[1],
                           is.signal=TRUE,
                           split.group=F,
                           afactor=1,
                           text.width=15, #k
                           cex=0.5,
                           map.cpdname=FALSE, #k
                           cpd.lab.offset=1.0,
                           discrete=list(gene=FALSE, cpd=FALSE),
                           limit=list(gene=1, cpd=1),
                           bins=list(gene=10, cpd=10),
                           both.dirs=list(gene=T, cpd=T),
                           low = list(gene = "green", cpd = "blue"),
                           mid = list(gene = "gray", cpd = "gray"),
                           high = list(gene = "red", cpd = "yellow"),
                           na.col="transparent",
                           new.signature=TRUE,
                           plot.col.key=TRUE,
                           key.align="x",
                           key.pos="topright",
                           sign.pos="bottomright",#g
                           ...){

  gR1=path.graph

  #group nodes mapping and merge
  grp.idx=node.data$size>1
  if(sum(grp.idx)>0 & !split.group){
    sub2grp=cbind(unlist(node.data$component[grp.idx], use.names=F), rep(names(grp.idx)[grp.idx], node.data$size[grp.idx]))
    du.idx=duplicated(sub2grp[,1])
    if(sum(du.idx)>0){
      du.rn=sub2grp[,1] %in% sub2grp[du.idx,1]
      message("Warning: reconcile groups sharing member nodes!")
      print(sub2grp[du.rn,])
      du.grps=sub2grp[du.idx,]
      rn=which(du.idx)
      for(r in rn){
        comps=node.data$component[[sub2grp[r,2]]]
        comps=comps[comps!=sub2grp[r,1]]
        node.data$component[[sub2grp[r,2]]]=comps
        node.data$size[sub2grp[r,2]]=node.data$size[sub2grp[r,2]]-1
      }
      sub2grp=sub2grp[!du.idx,]
    }
    rownames(sub2grp)=sub2grp[,1]
  } else sub2grp=NULL

  if(sum(grp.idx)>0 & !split.group){
    for(gn in names(grp.idx)[grp.idx]){
      gR1=combineKEGGnodes(node.data$component[[gn]], gR1, gn)
    }
  } else if(split.group){
    gR1=subGraph(nodes(gR1)[node.data$size==1], gR1)
  }
  nNames=nodes(gR1)
  nSizes=node.data$size[nNames]

  #unconnected nodes processing
  deg=degree(gR1)
  deg=deg$inDegree+deg$outDegree
  if(is.signal & sum(deg<1)>0){
    gR2=subKEGGgraph(nNames[deg>0], gR1)
    nNames=nNames[deg>0]
    nSizes=nSizes[deg>0]
    if(!is.null(sub2grp)){
      #  sub.idx=!sub2grp[,1] %in% names(deg[deg<1])
      sub.idx=sub2grp[,1] %in% nNames |sub2grp[,2] %in% nNames
    } else sub.idx=0
  } else {
    gR2=gR1
    if(!is.null(sub2grp)){
      sub.idx=rep(T, nrow(sub2grp))
    } else sub.idx=0
  }

  if(length(nNames)<2){
    msg=sprintf("%s not rendered, 0 or 1 connected nodes!\nTry \"kegg.native=T\" instead!", pathway.name)
    message("Note: ", msg)
    return(list())
  }

  #give up the KEGG positions, use graphviz layout

  #general attributes
  attrs=list()
  attrs$graph$rankdir="LR"
  attrs$node <- list(fixedsize=FALSE)

  #node attributes
  ntype=node.data$type[nNames]
  cpd.idx=ntype=="compound"
  map.idx=ntype=="map"
  rect.idx=!(cpd.idx|map.idx)
  nAttr=list()
  nAttr$label=rep('', length(nNames))
  shapes=node.data$shape[nNames]
  if(any(cpd.idx)) shapes[cpd.idx]="ellipse"
  if(any(map.idx)) shapes[map.idx]="plaintext"
  nAttr$shape=shapes
  nAttr$height=.75*17/46*nSizes*afactor
  nAttr$width=rep(.75, length(nNames))*afactor

  if(any(cpd.idx)){
    nAttr$height[cpd.idx]=nAttr$height[cpd.idx]*1.5
    nAttr$width[cpd.idx]=nAttr$width[cpd.idx]*1.5
  }
  
  if(any(map.idx)){
    nAttr$height[map.idx]=nAttr$height[map.idx]*1.5
    nAttr$width[map.idx]=nAttr$width[map.idx]*2
  }
  
  nAttr<- lapply(nAttr, function(x) {names(x) <- nNames})

  na.col=colorpanel2(1, low=na.col, high=na.col)
  fillcol=rep(na.col, length(nNames))
  names(fillcol)=nNames

  #edge attributes
  subdisplay <- subtypeDisplay(gR2)
  if(length(subdisplay)<1) {
    eAttrs=list() 
  } else{
    na.rn <- apply(subdisplay, 2, function(x) sum(is.na(x))==7)
    if(sum(na.rn)>0) subdisplay[,na.rn] <- KEGGEdgeSubtype[KEGGEdgeSubtype[,1]=="others", rownames(subdisplay)]
    eLabel <- subdisplay["label", ]
    eCol <- subdisplay["color", ]
    eTextCol <- subdisplay["fontcolor", ]
    eLty <- subdisplay["style", ]
    eArrowhead <- subdisplay["arrowhead", ]
    
    if (ncol(subdisplay) == 1) {
      tmp <- colnames(subdisplay)[1]
      names(eLabel) <- names(eCol) <- names(eTextCol) <- tmp
      names(eLty) <- names(eArrowhead) <- tmp
    }
    
    eAttrs <- list(lty=eLty, col=eCol, textCol=eTextCol, label=eLabel, arrowhead=eArrowhead)
  }


  layoutType=ifelse(is.signal, "dot", "neato")
  return(list(attrs=attrs, nodeAttrs=nAttr, eAttrs=eAttrs, layout=layoutType, nd=node.data))
}

