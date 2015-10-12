# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#this was from vince buffalo's book
#swl<-function (x) {
#  sapply(strsplit(sub("([\\d+MYX]+):(\\d+)-(\\d+)", "\\1;;\\2;;\\3",x, perl=TRUE),";;"),"[[",2)
#}

swl<- function(x) {as.numeric(sapply(strsplit(as.character(x),":|-"),function(x) x[2]))}


narrowPeakToGRanges<-function(file) {
  x <- read.table(file,stringsAsFactors=F)
  gr <-GRanges(seqnames=x$V1, ranges = IRanges(start=x$V2, end=x$V3),
               strand="*", score=x$V5, e=x$V7,summit=x$V10)
  return(gr)
}

setClass("fileset", representation( filename="character",count="numeric"))

countFileset<-function (x) {
  for (i in 1:length(x@filename)) {
  x@count[i]<-countBam(x@filename[i],param=ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE),))$records
  }
  return(x)
}

reds <- colorRampPalette(brewer.pal(n = 9, "Reds"))

p_param <- PileupParam(distinguish_nucleotides=FALSE,
                       distinguish_strands=FALSE,
                       min_nucleotide_depth=0)


tornado <- function(gr,dataset,pad=250,ord=1,window=5,color="red2") {

  #if gr is not uniform width, find the center of the gr and set width to 1
  if (length(gr) > 1 & var(width(gr))!=0) {
    gr<-GRanges(seqnames = seqnames(gr),
                IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                               (width(gr)+1)/2),width=1))
  }
  #set the color to whatever the user wants.
  reds<-colorRampPalette(c('snow', color))

  #this tells Rsamtools what and where we need data out of the file
  sb_param = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr+pad,
                          flag = scanBamFlag(isUnmappedQuery = FALSE))

  #scan the bam files

  files <- gsub("\\.bam","",sapply(strsplit(dataset@filename,"\\/"), function (x) x[length(x)]))

  resL<-list()

    for (i in files) {
    resL[[i]]<-pileup(dataset@filename[which(files==i)],
                      scanBamParam=sb_param, pileupParam=p_param)
    #this normalizes to read depth
    resL[[i]]$count<-ceiling(resL[[i]]$count*10e6/dataset@count[which(files==i)])
    }
  #print(paste0("Bam file read in:  ",t1))

  #swl<- function(x) {as.numeric(strsplit(as.character(x),":|-")[[1]][2])}

  #figure out order for which label
  which_factor_levels <- resL %>%
    bind_rows(.id="genotype") %>%
    dplyr::filter(genotype==files[ord]) %>%
    group_by(which_label) %>%
    summarize(s = sum(count)) %>%
    mutate(which_label = as.character(which_label)) %>%
    mutate(startpos = swl(which_label)) %>%
    ungroup() %>%
    arrange(s) %>%
    #select(which_label)
    select(which_label,startpos)

    # Make a named vector of startpos
    spHash<-which_factor_levels$startpos
    names(spHash)<-which_factor_levels$which_label
  #which_factor_levels<-data.frame(level=which_factor_levels)[,1],startpos=swl(which_factor_levels)[,1]))
  #which_factor_levels$startpos<-swl(which_factor_levels)


    #Add an option to set order to 0 and use the rev() of gr order for the heatmap


  suppressWarnings(resL %>%
                     bind_rows(.id="genotype") %>%
                     mutate(genotype=factor(genotype)) %>%
                     mutate(which_label = factor(which_label,levels=which_factor_levels$which_label)) %>%
                     mutate(startpos=plyr::round_any(pos-swl(which_label)-pad,window)) %>%
                     #mutate(startpos=plyr::round_any(pos-spHash[as.character(which_label)]-pad,window)) %>%
                     group_by(genotype,which_label,startpos) %>%
                     summarize(count=sum(count)) %>%
                     mutate(count=log2(count)) %>%
                     ungroup() %>%
                     ggplot(aes(x=startpos,y=which_label,fill=count)) +
                     geom_tile() +
                     scale_fill_gradient2(low = reds(256)[1],
                                          mid = reds(256)[64],
                                          high = reds(256)[256],
                                          midpoint = 0, space = "Lab",
                                          name = "Log2(cpm)") +
                     ylab("Genomic Regions") + xlab("Position Relative to Peak Center") +
                     facet_grid(. ~ genotype,scales="free_y") +
                     theme_bw() + theme(axis.text.y=element_blank(),
                                        axis.ticks.y=element_blank(),
                                        panel.grid.major=element_blank())
  )

# print(paste0("Dplyr:  ",t3))
}
