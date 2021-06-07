

ustplot <- function(x){
  n <- 10
  df <- as.data.frame(x)
  id <- df$ID[1:n]
  des <- df$Description[1:n]
  glist <- geneInCategory(x)[id]
  names(glist) <- des
  g <- unique(unlist(glist))
  
  dat <- matrix(0, nrow=length(g), ncol=length(id))
  rownames(dat) <- g
  for (i in 1:length(id)) {
    dat[glist[[i]], i] <- 1
  }
  colnames(dat) <- des
  dat <- as.data.frame(dat)
  dat$Name <- row.names(dat)
  upset(dat, nsets = n, mb.ratio = c(0.5, 0.5),order.by = c("freq"), decreasing = c(TRUE))
}



enrichDotplot <- function(object, x = "geneRatio", color = "p.adjust",
                             showCategory=20, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x", decreasing=TRUE) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  
  #nr <- unique(object@result[,x])
  #if(length(nr)<20){
    df <- object@result %>% arrange(p.adjust) %>% head(showCategory)
    df$GeneRatio <- parse_ratio(df$GeneRatio)
  #}else{
  #  df <- fortify(object, showCategory = showCategory, split=split)
  #}
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ylab("Description") + ggtitle(title) + theme_dose(font.size) + scale_size(range=c(3, 8))
  
}

#enrichDotplot(go, title = 'nm')



enrichBarplot <- function(height, x="Count", color='p.adjust', showCategory=20, font.size=12, title="", ...) {
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  object <- height
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  
  #nr <- unique(object@result[,x])
  #if(length(nr)<20){
    df <- object@result %>% arrange(p.adjust) %>% head(showCategory)
    df$Description <- factor(df$Description, levels=rev(unique(df$Description)))
  #}else{
  #  df <- fortify(object, showCategory=showCategory, by=x, ...)
  #}
  
  
  if(colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
  } else {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = "Description")) +
      theme_dose(font.size) +
      theme(legend.position="none")
  }
  p + geom_bar(stat = "identity") + coord_flip() + ggtitle(title) + xlab("Description") + ylab("GeneRatio")
}

#enrichBarplot(go, title = 'nm')





## enrichBarplot <- function(reactome, title){
##   library(ggplot2)
##   reactome@result$GR <- as.numeric(str_split(reactome@result$GeneRatio,'/',simplify = TRUE)[,1])
##   reactdf <- reactome@result
##   reactdf <- reactdf %>% head(min(nrow(reactome@result),20)) 
##   reactdf <-arrange(reactdf, -row(reactdf)[,1])
##   terms <- factor(row(reactdf)[,1], labels = reactdf$Description)
##   
##   p <- ggplot(reactdf,aes(y=reactdf$GR, x=terms, fill=reactdf$p.adjust))
##   p1 <- p + geom_bar(stat="identity") + coord_flip()
##   p1 <- p1 + scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
##   p1 <- p1 + labs(x="Description",
##                   y="GeneRatio",
##                   title=title)
##   p1 <- p1 + scale_y_continuous(expand = c(0,0.025))+theme_linedraw(base_size=16) 
##   return(p1)
## }
## 
## enrichBarplot(go, title = 'nm')
## 
## enrichDotplot <- function(reactome, title){
##   library(ggplot2)
##   grt <- str_split(reactome@result$GeneRatio,'/',simplify = TRUE)
##   reactome@result$GR <- as.numeric(grt[,1])/as.numeric(grt[,2])
##   reactdf <- reactome@result
##   reactdf <- reactdf %>% head(min(nrow(reactome@result),20)) %>% arrange(GR, -p.adjust)
##   terms <- factor(as.integer(rownames(reactdf)), labels = reactdf$Description)
##   
##   p <- ggplot(reactdf,aes(y=reactdf$GR, x=terms, size=Count, color=p.adjust))
##   p1 <- p + geom_point(stat="identity") + coord_flip()
##   p1 <- p1 + scale_color_continuous(low="red", high="blue", name = 'p.adjust', guide=guide_colorbar(reverse=TRUE))
##   p1 <- p1 + theme(legend.key = element_blank())
##   p1 <- p1 + labs(x="Description",
##                   y="GeneRatio",
##                   title=title)
##   p1 <- p1 + theme_linedraw(base_size=16) #+ scale_y_continuous(expand = c(0,0))
##   return(p1)
## }
## 




