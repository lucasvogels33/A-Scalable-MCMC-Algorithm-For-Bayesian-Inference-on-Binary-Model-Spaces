#############
#Step 1: download the excel file containing the data
#       - Go to http://rstats.immgen.org/DataPage/
#       - Search for GSE15907
#       - Click "GSE15907 Normalized Data" (wait for an hour :), this is a big file)
#       - Save excel file in the same folder as this R script

library(readxl) #for loading excel file to R
library(bioRad) #for checking NaN values in dataframe
library(tidyverse) #for ordering dataframe
library(huge) #for Gaussianization of the data 
library(stringr) #count appearances of a character in a string

##################################################################
#####load data, logtranform and rank by variance #################
##################################################################

#load data
setwd("C:/this_folder")
data_orig <- read_excel("data.xlsx")

#save data_frame
data = data_orig

#test for NA or NaN
sum(is.na(data)) # test for NA (if 0, then there are no NA's)
sum(is.na.data.frame(data)) #test for NaN (if 0, then there are no NaN's)

#save names of genes (variables) and cells (obserations)
gene_ID = data$ProbeSetID #save the IDs of the 24922 genes
gene_names = data$GeneSymbol #save the names of the 24922 genes

#remove non-number columns (we want our matrix to contain only data, for the log transform)
data$ProbeSetID = NULL 
data$GeneSymbol = NULL

#save the names of the cells
cell_names = colnames(data) #save the names of the 653 cells

#check whether dimensions, max, min make sense
dim(data)[2] #n should be 653
dim(data)[1] #p should be 24922
max(data)
min(data)

#observe heteroskedacity
var_vec = apply(data,1,function(x) var(x))
mean_vec = apply(data,1,function(x) mean(x))
plot(NA,xlim=c(0,10000),ylim=c(0,20000000),xlab="mean of the variable",ylab="variance of the variable",main="heteroskedacity plot")
points(x=mean_vec,y=var_vec)

#log 2 transformation to avoid heteroskedacity
data = log(data,2)

#observe that heteroskedacity disappeared after log 2 transformation
var_vec = apply(data,1,function(x) var(x))
mean_vec = apply(data,1,function(x) mean(x))
plot(NA,xlim=c(0,15),ylim=c(0,15),xlab="mean of the variable",ylab="variance of the variable",main="heteroskedacity plot")
points(x=mean_vec,y=var_vec,col = rgb(0, 0, 0, alpha = 0.03),cex = 0.2,asp=1,pch=16)

#Add a column with the row variance
var_vec = apply(data,1,function(x) var(x))
data = cbind(data,variance = var_vec)

#make histogram of variable variances
hist(data$variance,breaks=seq(0,13,0.25),xlab="variance")

#add removed columns back
data = cbind(data,gene_ID = gene_ID) #add the gene IDs 
data = cbind(data,gene_name = gene_names) #add the gene names

#order the rows by descending variance (top row has the highest variance)
#data_ordered = data %>% arrange(desc(variance))  # arrange in descending order
data_ordered = data[order(-data$variance),]  # arrange in descending order

#########################################################################
############ Select 2.5% most variable genes#############################
#########################################################################

#remove non-number columns from the data 
data_ordered$gene_name = NULL 
data_ordered$gene_ID = NULL 
data_ordered$variance = NULL 

#select the top 2.5% variables
p = nrow(data_ordered)
select = round(0.025*p)
data_select = data_ordered[1:select,] 
cell_names = colnames(data_select)

#########################################################################
############ Center each observation within each cell type#############################
#########################################################################

#confirm that every cell_names has exactly one hashtag in it. If the for loop does not print anything, it means every cell name has exactly one hashtag
for (i in 1:length(cell_names)){
  if (str_count(cell_names[i],"#")!=1){
    print(cell_names[i])
  }
}

#confirm that the hashtag is always the penultimate character
for (i in 1:length(cell_names)){
  len = str_length(cell_names[i])
  pos = unlist(gregexpr("#",cell_names[i]))[1]
  if ((len - pos) !=1){
    print(cell_names[i])
  }
}

#find all cell types (the cell type is the part that comes in front of the hashtag)
#make a dataframe
cell_types_all = sub("\\#.*$", "", cell_names) #retrieve the cell_type (the string before the hastag )
cell_types = unique(cell_types_all)
cell_df = data.frame(
  id = 1:653,
  name = cell_names,
  type = cell_types_all
)

#count how often a cell type appears only once
which(table(cell_types_all)==1) #answer: only twice does a cell type appear once

#center each gene within each cell type
p = nrow(data_select)
data_select_old = data_select
for (cell_type in cell_types){
  indices = which(cell_df$type == cell_type)
  for (gene in 1:p){
    average = mean(as.matrix(data_select[gene,indices]))
    data_select[gene,indices] = data_select[gene,indices] - average
  }
}

#check whether the centering went correctly, by means of three sanity checks
#1. Is the mean of each row equal to zero? Answer: yes it is
for (i in 1:p){
  if (mean(as.matrix(data_select[i,]))>0.0001){
    print(i)
  }
}
#2. What are the maximum and minimum values
max(data_select)
min(data_select)
#3. Are the columns that are alone in their cell type a mean of zero
for (cell_type in cell_types_all){
  indices = which(cell_types_all==cell_type)
  if (length(indices)==1){
    print(indices)
  }
}
data_select[1:5,455]
data_select[1:5,456]

#What are the variances? Answers: variances fluctuate between 0.02 and 0.42 with a mean of 0.07
var_vec = c()
for (i in 1:p){
  var_vec = c(var_vec,var(as.matrix(data_select[i,])))
}

#transpose data (both normalisation (next step) and algorithms require n (rows) x p (columns) data )
data_select = t(data_select)

#verify that dimension are n x p = 653 x 623
nrow(data_select) #n should be 653
ncol(data_select) #p should be 623

#observe data is more or less normally distributed
hist(data_select[,1]) 
hist(data_select[,200]) 

#We still gaussianize the data, so that the variance is equal to one everywhere
data_normal = huge.npn(x=data_select,npn.func="shrinkage") 

#observe data is normal
hist(data_normal[,1]) #observe data is normal
hist(data_normal[,2]) #observe data is normal

#save cleaned data in Rdata file
output = list(data = data_normal)
filename = paste0("cleaned_data.Rdata")
save(output,file=filename)