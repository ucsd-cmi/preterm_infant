#Stephen Wandro
#7/2/2020

#Decontam Rscript

#Usage
#Rscript run_decontam.R table.biom metadata.tsv ouput.tsv

#Arg1: biom abundance table
#Arg2: metadata. col1:sample name, col2:TRUE/FALSE is blank
#Arg3: output tsv

library(decontam)
library(biomformat)

args <- commandArgs(trailingOnly = TRUE)
#args[[1]] <- "16S_table_with_blanks.biom"
#args[[2]] <- "blank_metadata.tsv"



#Import table
dat <- read_biom(args[[1]])
#Convert to df
dat <- t( as.matrix(biom_data(dat)) )
#Import metadata
md <- read.table(args[[2]], sep='\t', header=TRUE)

  #Order samples
dat <- dat[order(match(rownames(dat),md[,1])),]

#Confirm same number of samples
if (nrow(md) != nrow(dat)){
  stop("Differing number of samples between data and metadata.")
}

#Confirm same order
if (any(md[,1] != rownames(dat))){
  stop("Not all sample names identical between data and metadata")
}

#Run Decontam prevalence
#Threshhold set at default 0.1, but can set custom threshholds with output table.
prev.contam.results <- isContaminant(as.matrix(dat), method="prevalence", neg=md[,2], threshold=0.1)

#Write result
write.table(x = prev.contam.results, file= args[[3]], sep = '\t')
