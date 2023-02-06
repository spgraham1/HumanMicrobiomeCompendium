#read in compendium data
#rows are samples, columns are taxa
library(stringr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("please input your clr-transformed data, desired sample size, desired repetition number,region of interest, and metadata")
} else if (length(args) < 7) {
  #default output file
  region <- args[4]
  abbr <- paste(unlist(str_extract_all(region,"([:upper:])")),collapse="") #get region abbreviation
  sampleSize <- args[2]
  args[6] = paste0("results/",abbr,"_n",sampleSize,"_min.txt")
  medianOutfile <<- paste0("results/",abbr,"_n",sampleSize,"_median.txt")
  allDistsFile <<-paste0("results/",abbr,"_n",sampleSize,"_all.txt") 
  print(paste0("will write to output file ",args[6]))
}

if (length(args) == 7) {
  threads <<- args[7]
  print("threads")
  print(threads)
} else {
  threads <<- 32
}
print(args[1])
d2s <<- read.table(file=args[1],sep="\t",row.names=1,header=TRUE)
sampleSize <- args[2]
reps <<- args[3]
region <<- args[4]
metadata <- args[5]
metadata <<- read.table(file=metadata,sep="\t")
minOutfile <<- args[6]
print(minOutfile)
print(medianOutfile)

if (file.exists(allDistsFile)) {
  file.remove(allDistsFile)
}

library(compositions)
library(parallel)

euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))
}


getDist <- function(sample1,sample2) {
  return(euclidean(d2s[sample1,],d2s[sample2,]))
}

###################################################################
# this function will return a list of n+1 samples, where the last #
# sample in the list is the "new" sample from the chosen region   #
###################################################################

getSample  <- function(n) {
  print("get sample")
  numSamples <- nrow(d2s)
  metadata$chosen <- FALSE #set flag to false (no sample has been selected yet)
  available <- which(metadata$region == region) 
  newSample <- sample(available,1,replace=FALSE) #select a sample from the region of interest
  metadata$chosen[newSample] = TRUE
  available <- which(metadata$chosen == FALSE) 
  randomSample <- sample(available,n,replace=FALSE) #select n samples from anywhere, as long as they haven't already been chosen
  sample <- append(randomSample,newSample) 
  return(sample)
}

###################################################################
# find the smallest distance between the new sample and all other #
# samples chosen. return the min and median of these distances    #
###################################################################
getMinDist <- function(n) {
  numSamples <- nrow(d2s)
  sample <- getSample(n) #
  if (length(sample) != n+1) {
    stop("incorrect sample size selected")
  }
  newSample <- sample[length(sample)]
  randomSample <- sample[1:(length(sample)-1)]
  allDists <- mclapply(randomSample,getDist,sample2=newSample,mc.preschedule = TRUE,mc.cores=threads) #compare each microbiome in original sample to the new sample
  allDists <- unlist(allDists)
  write(allDists,file=allDistsFile,append=TRUE,sep="\t") #save all the calculated distances to another file
  minDist <- min(allDists)
  medianDist <- median(allDists)
  return(c(minDist,medianDist))
}

###################################################################
# run the minimum distance calculation for the specified number   #
# of repetitions                                                  #
###################################################################
getRepetions <- function(n) {
  sampleSize <- c(rep(n,reps)) #apply function will iterate over this vector
  sampleDists <- sapply(sampleSize,getMinDist) 
  if (ncol(sampleDists) != reps) {
    print("incorrect dimensions of sampleDists")
    print(paste0("Reps: ",reps))
    print(paste0("Dimensions: ",dim(sampleDists)))
  }
  return(sampleDists)
}

###################################################################
# calculate the minimum distance between 1 sample from the region #
# of interest and "sampleSize" other samples, repeated "reps"     #
# times                                                           #
###################################################################

oneSampleSize <- function(sampleSize) {
  if (!(identical(metadata$sample,rownames(d2s)))) {
    stop("samples out of order, cannot proceed")
  }
  reps <- as.numeric(reps)
  sampleSize <- as.numeric(sampleSize)
  output <- getRepetions(sampleSize)
  minOutput <- output[1,]
  medianOutput <- output[2,]
  write.table(minOutput,file=minOutfile,sep="\t")
  write.table(medianOutput,file=medianOutfile,sep="\t")
}

oneSampleSize(sampleSize)
warnings()
