#!/usr/bin/env Rscript
suppressMessages({
    library(optparse)
    library(data.table)
    library(ShortRead)
})

.dustySplit <- function(x, batch=100000) {
     unlist(unname(lapply(split(x, ceiling(seq_along(x) / batch)), dustyScore)))
}

.dustScore <- function(vir) {
    ## maximim possible dusty score for a given sequence length
    max.low <- data.table(
      len=4:1000,
      mlow=dustyScore(DNAStringSet(sapply(4:1000, function(x) paste0(rep("A", x), collapse="")))))
    setkey(max.low, len)
    if(nrow(vir)>0) {
        ## 1. split dustyScore into chunks to save some memory?
        tmp.tbl <- vir[,.(
          qname,
          gi,
          len=nchar(seq),
          raw.low=.dustySplit(DNAStringSet(seq), 100000),
          low=0
        )]
        setkey(tmp.tbl, len)
        tmp.tbl[max.low, low:=raw.low / mlow]
    } else {
        tmp.tbl <- data.table(qname=character(0),
                              gi=character(0),
                              len=integer(0),
                              raw.low=numeric(0),
                              low=numeric(0)
                              )
    }
    return(tmp.tbl)
}

.importReads <- function(bam.fn) {
    idxs <- fread(paste(system2("samtools", sprintf("idxstats %s", bam.fn), stdout=TRUE), collapse="\n"),
                  select=c(1,3))[grep("gi|chrEBV",V1)][V3>0]
    sels <- sprintf("%s:1-1000000000", idxs$V1)
    if (length(sels)>0) {
        vir <- rbindlist(lapply(sels, function(fn) {
            fread(paste(system2("samtools", sprintf("view %s \"%s\" | cut -f1,3,10 | head -n 10000", bam.fn, fn), stdout=TRUE),
                        collapse="\n"), header=FALSE, showProgress=FALSE)
        }))
        setnames(vir, c("qname", "gi", "seq"))
    } else {
        vir <- data.table(qname=character(0), gi=character(0), seq=character(0))
    }
    return(vir)
}


.processRun <- function(bam.fn, low.cut){
    vir <- .importReads(bam.fn)
    dst <- .dustScore(vir)
    res <- dst[,.(n.total=.N, n.good=nrow(.SD[low<low.cut])), by=gi]
    setkey(res, gi)
    return(res)
}


option_list = list(
)

parser = optparse::OptionParser("virus-call.R",
  description=c("\n"),
  epilogue=c(
    "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
    "Michigan Center for Translational Pathology (c) 2016\n"),
  option_list=option_list)


opt <- optparse::parse_args(parser, positional_arguments=TRUE)

bam.fn <- opt$args[1]
out.fn <- opt$args[2]

res <- .processRun(bam.fn, 0.05)
write.table(res, out.fn, sep="\t", row.names=FALSE, quote=FALSE)
