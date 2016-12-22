# General Tips
#
# Use whitespace freely in comments and code as it makes things easier to read
#
# Break down a complex task by using well placed function or two
#

#accumulation curve graph

rm(list=ls())

library(ggplot2)
library(dplyr) # useful for data manipulation
library(splitstackshape)
library(tidyr)

#reads in MQ data
#
# Don't use the row.names=6 argument because you put those row.names back into a column later anyway 
#
# Use na.strings="0" as this will properly mark missing values as NA rather than 0. 
# MaxQuant unfortunately uses 0 to denote missing values but this is generally a terrible 
# idea as it can be mistaken for a genuine 0 reading which is different from a missing reading.
# By using NA it will also be easier to extract subsets of interest
#
# MQ_data = read.table("proteinGroups_PNGaseFSL.txt",sep="\t",header=TRUE, row.names = 6)
#
# 
MQ_data = read.table("proteinGroups_PNGaseFSL.txt",sep="\t",header=TRUE,na.strings = "0")


#takes out iBAQ values
#
# This should work however I provide a safer and more convenient way to do this below.
# I recommend using this approach as it eliminates issues due to accidental miscounting of column numbers
# The code below uses the pipe %>% operator which allows functions to be chained without using tons of nested brackets
#
# iBAQdata = MQ_data[c(68,69,70,71,72,73,74,75)]

iBAQdata = MQ_data %>% select(Protein.IDs,contains("iBAQ"))

# Commented this code as it is not needed due to using NA's and not using the row.names argument
#
#subsets the data to get rid of proteins with a value of '0' abundance
#iBAQdata = subset(iBAQdata, !(iBAQ.1_160309183003 == 0 & iBAQ.2 == 0 & iBAQ.3 == 0 & iBAQ.4 == 0 & iBAQ.1P == 0 & iBAQ.2P == 0 & iBAQ.3P == 0 & iBAQ.4P == 0))
#adds a row names columns
#iBAQdata$names = row.names(iBAQdata)

# Before subsetting we need to cSplit the data so that protein groups are split into multiple rows. 
# We need to perform the analysis at this level (ie individual proteins not proteingroups) because 
# it is likely that protein group membership will differ between samples and we want to know which 
# proteins with the group are in common between samples.
#
iBAQ_by_protein <- cSplit(iBAQdata,"Protein.IDs",sep=";",direction = "long")

# This places everything into tall format.  The sample is encoded in a column now.
# Firstly this means we can get rid of all the NA values in on step
# The NA removal is chained using %>%
#
iBAQ_by_protein_tall <- iBAQ_by_protein %>% 
    gather(sample,iBAQ_value,-Protein.IDs) %>% 
    filter(!is.na(iBAQ_value)) %>% # Removes NAs
    filter(!grepl("CON_",Protein.IDs)) # Removes contaminants



# The next benefit is that it is now relatively easy to get the list of proteins found in any particular sample
#
# For example
# iBAQ_by_protein_tall %>% filter(sample=="iBAQ.2") %>% select(Protein.IDs)
# 
# We capture this idea in a function. 
# Note that the function allows multiple samples at once which is very useful. By passing two samples
# eg c("iBAQ.2","iBAQ.1","iBAQ.3") we will get the list of proteins covered by all three samples
#
distinct_proteins_in_samples <- function(samples){
  iBAQ_by_protein_tall %>% 
    filter(sample %in% samples) %>% 
    select(Protein.IDs) %>% 
    distinct()
}


# This function that compares proteins present in two lists of samples
# It returns the new proteins present in list2 that are not in list1
#
new_proteins <- function(samples_1,samples_2){
  prots_samples1 <- distinct_proteins_in_samples(samples_1)
  prots_samples2 <- distinct_proteins_in_samples(samples_2)
  
  setdiff(prots_samples2,prots_samples1)
}

# Now we can generate an accumulation curve given a set of sample in the desired order
#
accumulation <- function(samples){
  n1 <- nrow(distinct_proteins_in_samples(samples[1]))
  n_prots <- sapply(2:length(samples),function(i){
    nrow(new_proteins(samples[1:(i-1)],samples[i]))
  })
  c(n1,n_prots)
}

# Finally we need to be able to generate accumulation curves for various orderings
# This can be done using the combinat package

library(combinat)

# Get names of non-pngase-F samples (for example)
non_pngasef <- colnames(iBAQ_by_protein)[c(3,5,7,9)]
pngasef <- colnames(iBAQ_by_protein)[c(2,4,6,8)]

# All possible permutations are generated using permn like this
#
# permn(non_pngasef)
#

# Using sapply we can easily get the accumulation curves of all permutations
#
acc_allNP <- sapply(permn(non_pngasef),accumulation)
acc_allP <- sapply(permn(pngasef),accumulation)

# This just needs to be transposed and have row and column labels added and it is ready for plotting
# We let ggplot to the averaging

plot_dataNP <- as.data.frame(t(apply(acc_allNP,2,cumsum))) %>% gather(sample_num,total_proteins)
plot_dataP <- as.data.frame(t(apply(acc_allP,2,cumsum))) %>% gather(sample_num,total_proteins)

plot_dataNP2 <- aggregate(total_proteins~sample_num,plot_dataNP,mean)
plot_dataP2 <- aggregate(total_proteins~sample_num,plot_dataP,mean)


ggplot(plot_data,aes(x=sample_num,y=total_proteins)) + geom_point()

ggplot() + 
  geom_point(data = plot_dataNP2, aes(x=sample_num,y=total_proteins, color = "red")) +
  geom_point(data = plot_dataP2, aes(x=sample_num,y=total_proteins, color = "blue")) + 
  scale_color_manual(labels = c("Non-PNGase F", "PNGase F") , values = c("blue", "red")) +
  xlab('sample_num') + 
  ylab ('total_proteins')


