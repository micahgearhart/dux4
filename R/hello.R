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
wls<- function(x) { paste0(seqnames(x),":",start(x),"-",end(x))}

strstrip<- function(s,i) sapply(strsplit(s,"_"),function(x) x[i])

label_wrap <- function(variable, value) {
  lapply(strwrap(as.character(value), width=25, simplify=FALSE),
         paste, collapse="\n")
}

#http://stackoverflow.com/questions/17319487/median-and-quartile-on-violin-plots-in-ggplot2
median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out)
}


narrowPeakToGRanges<-function(file) {
  x <- read.table(file,stringsAsFactors=F)
  gr <-GRanges(seqnames=x$V1, ranges = IRanges(start=x$V2, end=x$V3),
               strand="*", score=x$V5, e=x$V7,summit=x$V10)
  return(gr)
}

center<-function(gr) {
 GRanges(seqnames = seqnames(gr),
              IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                             (width(gr)+1)/2),width=1))
}

setClass("fileset", representation( filename="character",count="numeric",labels="character"))

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
####################
# TORNADO FUNCTION #
####################

tornado <- function(gr,dataset,pad=250,ord=1,window=5,color="red2") {

  #if gr is not uniform width, find the center of the gr and set width to 1
  if (length(gr) > 1 & var(width(gr))!=0) {
    gr<-GRanges(seqnames = seqnames(gr),
                IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                               (width(gr)+1)/2),width=1))
  }

  #add which_label to gr
  gr2<-gr+pad
  gr2$which_label<-paste0(seqnames(gr2),":",start(gr2),"-",end(gr2))


  #set the color to whatever the user wants.
  reds<-colorRampPalette(c('snow', color))

  #this tells Rsamtools what and where we need data out of the file
  sb_param = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
                          flag = scanBamFlag(isUnmappedQuery = FALSE))

  #scan the bam files

  files <- gsub("\\.bam","",sapply(strsplit(dataset@filename,"\\/"), function (x) x[length(x)]))
  names(files)<-dataset@labels

  #plot labeller
  #plot_labeller <- function(variable,value){ return(names(which(files==value))) }

  resL<-list()

    for (i in files) {
    resL[[i]]<-pileup(dataset@filename[which(files==i)],
                      scanBamParam=sb_param, pileupParam=p_param)
    #this normalizes to read depth
    resL[[i]]$count<-ceiling(resL[[i]]$count*10e6/dataset@count[which(files==i)])
    }
  #print(paste0("Bam file read in:  ",t1))

  #swl<- function(x) {as.numeric(strsplit(as.character(x),":|-")[[1]][2])}

  #if not ordering use the first file to do the ordering
  filter_file <- ifelse(ord==0,files[1],files[ord])


    which_factor_levels<-rev(as.character(gr2$which_label))

 if (ord==0) {
    spHash<-swl(which_factor_levels)
    names(spHash)<-which_factor_levels

 } else {
   #print("ORDER IS NOT 0")
   temp2 <- resL %>%
     bind_rows(.id="genotype") %>%
     dplyr::filter(genotype==filter_file) %>%
     group_by(which_label) %>%
     summarize(s = sum(count)) %>%
     arrange(s) %>%
     mutate(which_label = as.character(which_label)) %>%
     dplyr::select(which_label)

   temp1<-which_factor_levels[!which_factor_levels %in% temp2$which_label]
   which_factor_levels<-c(temp1,temp2$which_label)

   spHash<-swl(which_factor_levels)
   names(spHash)<-which_factor_levels


}

  suppressWarnings(  resL %>%
                     bind_rows(.id="genotype") %>%
                     mutate(genotype=factor(genotype,levels=files,labels=names(files))) %>%
                     #mutate(genotype=factor(genotype,levels=c(1,2),labels=c("green","blue"))) %>%
                     mutate(which_label = factor(which_label,levels=names(spHash))) %>%
                     #mutate(which_label = factor(which_label,levels=which_factor_levels$which_label)) %>%
                     mutate(startpos=0.001*plyr::round_any(pos-swl(which_label)-pad - 0.5*width(gr[1]),window)) %>%
                     #mutate(startpos=plyr::round_any(pos-spHash[as.character(which_label)]-pad,window)) %>%
                     group_by(genotype,which_label,startpos) %>%
                     summarize(count=sum(count)) %>%
                     #dplyr::filter(startpos != max(startpos) ) %>%
                     mutate(count=log2(count)) %>%
                     ungroup() %>%
                     ggplot(aes(x=startpos,y=which_label,fill=count)) +
                     geom_tile() +
                     scale_fill_gradient2(
                       low="white",mid=reds(256)[16],high=color,midpoint=0,
                       #low = reds(256)[1],
                         #                 mid = reds(256)[64],
                        #                  high = reds(256)[256],
                       #                   midpoint = 0, space = "Lab",
                                          name = "Log2(cpm)") +

                       ylab("Genomic Regions") + xlab("Position Relative to Peak Center (kb)") +
                     facet_grid(. ~ genotype,scales="free_y",labeller=label_wrap) +
                     theme_bw() + theme(axis.text.y=element_blank(),
                                        axis.ticks.y=element_blank(),
                                        panel.grid.major=element_blank(),
                                        panel.grid.minor=element_blank())
  )


}





######################
# Twister Function #
######################


twister <- function(gr,dataset,pad=250,ord=1,window=5,color="red2",ya=c(0,20)) {

  #if gr is not uniform width, find the center of the gr and set width to 1
  if (length(gr) > 1 & var(width(gr))!=0) {
    gr<-GRanges(seqnames = seqnames(gr),
                IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                               (width(gr)+1)/2),width=1))
  }

  #add which_label to gr
  gr2<-gr+pad
  gr2$which_label<-paste0(seqnames(gr2),":",start(gr2),"-",end(gr2))


  #set the color to whatever the user wants.
  reds<-colorRampPalette(c('snow', color))

  #this tells Rsamtools what and where we need data out of the file
  sb_param = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
                          flag = scanBamFlag(isUnmappedQuery = FALSE))

  #scan the bam files

  files <- gsub("\\.bam","",sapply(strsplit(dataset@filename,"\\/"), function (x) x[length(x)]))
  names(files)<-dataset@labels

  #plot labeller
  #plot_labeller <- function(variable,value){ return(names(which(files==value))) }

  resL<-list()

  for (i in files) {
    resL[[i]]<-pileup(dataset@filename[which(files==i)],
                      scanBamParam=sb_param, pileupParam=p_param)
    #this normalizes to read depth
    resL[[i]]$count<-ceiling(resL[[i]]$count*10e6/dataset@count[which(files==i)])
  }
  #print(paste0("Bam file read in:  ",t1))

  #swl<- function(x) {as.numeric(strsplit(as.character(x),":|-")[[1]][2])}

  #if not ordering use the first file to do the ordering
  filter_file <- ifelse(ord==0,files[1],files[ord])


  which_factor_levels<-rev(as.character(gr2$which_label))
  spHash<-swl(which_factor_levels)
  names(spHash)<-which_factor_levels


  suppressWarnings( resL %>%
                       bind_rows(.id="genotype") %>%
                       mutate(genotype=factor(genotype,levels=files,labels=names(files))) %>%
                       #mutate(genotype=factor(genotype,levels=c(1,2),labels=c("green","blue"))) %>%
                       mutate(which_label = factor(which_label,levels=names(spHash))) %>%
                       #mutate(which_label = factor(which_label,levels=which_factor_levels$which_label)) %>%
                       #mutate(startpos=plyr::round_any(pos-swl(which_label)-pad,window)) %>%
                       mutate(startpos=0.001*plyr::round_any(pos-swl(which_label)-pad - 0.5*width(gr[1]),window)) %>%
                       #mutate(startpos=plyr::round_any(pos-spHash[as.character(which_label)]-pad,window)) %>%
                       group_by(genotype,startpos) %>%
                       summarize(count=sum(count)) %>%
                       dplyr::filter(startpos != max(startpos) ) %>%
                       #  dplyr::filter(count> 4) %>%
                       mutate(count=log2(count)) %>%
                       dplyr::filter(abs(startpos)!=pad) %>% #remove artifacts at the end.
                       ungroup() %>%
                       ggplot(aes(x=startpos,y=count)) +
                       geom_line(aes(color=genotype)) +
                    #   scale_fill_gradient2(
                    #     low="white",mid=reds(256)[16],high=color,midpoint=4,
                    #     #low = reds(256)[1],
                    #     #                 mid = reds(256)[64],
                    #     #                  high = reds(256)[256],
                    #     #                   midpoint = 0, space = "Lab",
                    #    name = "Log2(cpm)") +

                       ylab("Log2(cpm)") + xlab("Position Relative to Peak Center (kb)") +
		       ylim(ya) +
                     #  facet_grid(. ~ genotype,scales="free_y",labeller=label_wrap) +
                       theme_bw() + theme(panel.grid.major=element_blank(),
                                          panel.grid.minor=element_blank())

  )

}


## DUX4 Binding Site Analysis

findTFBS<-function(gr,motif=dux4.pwm,score="85%") {
  seq<-as.character(seqnames(gr))
  gr.seq<-getSequence(gr,genome=hg19)

  #forward Strand
  (hits_f<-matchPWM(motif,gr.seq[[1]],score,with.score=T))
  #hits_fGR<-GRanges(seqnames=seq,IRanges(start=start(gr)+start(hits_f)-1,
  #                                       width=width(hits_f)), strand="+",score=elementMetadata(hits_f)$score)
  #reverse Strand
  (hits_r<-matchPWM(motif,reverseComplement(gr.seq[[1]]),score,with.score=T))
  #hits_rGR<-GRanges(seqnames=seq,IRanges(start=start(gr)+width(gr)-start(hits_r)-width(hits_r),
  #                                       width=width(hits_r)),strand="-",score=elementMetadata(hits_r)$score)
  #x<-c(hits_fGR,hits_rGR)
  #  if (length(x) > 0) {return(x)} else {return(NULL)}
  return (length(hits_f)+length(hits_r))
}



#####################
# NORMALIZE BIGWIG  #
#####################

normalizeBigWig<-function(fs) {
  for (i in 1:length(fs@labels)) {

    system2("/usr/ngs/bedtools2/bin/bedtools",
            paste0("genomecov -scale ",1e6/fs@count[i],
                   " -split -bg -ibam ", fs@filename[i]),
            stdout=gsub(".bam",".bedGraph",basename(fs@filename[i])),
            stderr="testerr.out")

    system2("/usr/ngs/bin/bedGraphToBigWig",
            paste0(gsub(".bam",".bedGraph",basename(fs@filename[i])),
                   " /usr/ngs/hg19.chromsize ",gsub(".bam",".bigWig",basename(fs@filename[i]))),
            stdout="test2.out",stderr="testerr2.out")

  }
}
