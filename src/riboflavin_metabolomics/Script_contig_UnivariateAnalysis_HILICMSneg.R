###################################### Updates, install and call packages 
install.packages("ggsci")
install.packages("matrixStats")
install.packages("ggrepel")
install.packages("tidyverse")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(matrixStats)
  library(ggrepel)
  library(VennDiagram)
})

Sys.time()
rm(list=ls())

###################################### Input and process files:
Directory<- c("C:/Users/selegato/Desktop/R projects/Projects/Antwerp_Caroline_BacterialCulture_HILICMSneg_targeted/")
setwd(Directory)

### 1) Import data
# We need to import the feature table and metadata
ft_url <- paste0(Directory, "/Quant_Table_20250703.csv")
# metadata
md_url <-  paste0(Directory, "/metadata_20250703.csv")


# read data
ft <- read.csv(ft_url, header = T, check.names = F, row.names = 1)
md <- read.csv(md_url, header = T, check.names = F, row.names = 1)

### 2)  Remove noise values (bellow 30000)
ft[ft < 30000] <- 0.1

# dimension of the data
dim(ft)
dim(md)

head(ft)
head(md)

### 3) Bring feature table and metadata in the correct format:
# how many files in the metadata are also present in the feature table
table(rownames(md) %in% colnames(ft))
# which file names in the metadata are not in the feature table?
setdiff(rownames(md),colnames(ft))
md <- md[rownames(md) %in% colnames(ft),]
dim(md)


### 4) transpose ft
ft <- t(ft)
head(ft)
dim(ft)


### 5) In order to perform a univariate as described below, it is important that the filenames
# in our metadata are identical as well as in the same order as the filenames in 
# our feature table. Let's make sure this is true, using the below code (this 
# should return TRUE).
identical(rownames(ft),rownames(md))
# put the rows in the feature table and metadata in the same order
ft <- ft[match(rownames(md),rownames(ft)),]
identical(rownames(ft),rownames(md))


### 6) Merge feature table and metadata to one dataframe
ft<-data.frame(ft)
DataI <- cbind.data.frame(md,ft)
dim(DataI)

### 7) Filter data
Data<-DataI %>%
  filter(ATTRIBUTE_RELATED=="Related") %>%
  filter(ATTRIBUTE_STRAIN!="Waterblank")%>%
  filter(ATTRIBUTE_STRAIN!="V370")%>%
  filter(ATTRIBUTE_MEDIA!="MRS-R")%>%
  filter(ATTRIBUTE_PRODUCER!="COMMERCIAL_STRAINS_CONSUMER")%>%
  filter(ATTRIBUTE_PRODUCER!="BAD_GROWERS_Streptococcaceae")%>%
  filter(ATTRIBUTE_RELATED!="Commercial Strains")










###################################### Create the Function: InsideLevels
InsideLevels <- function(metatable){
  lev <- c()
  typ<-c()
  for(i in 1:ncol(metatable)){
    x <- levels(droplevels(as.factor(metatable[,i])))
    if(is.double(metatable[,i])==T){x=round(as.double(x),2)}
    x <-toString(x)
    lev <- rbind(lev,x)
    
    y <- class(metatable[,i])
    typ <- rbind(typ,y)
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)),ATTRIBUTES=colnames(metatable),LEVELS=lev,TYPE=typ,row.names=NULL)
  return(out)
}



###################################### T-test

### 1)  Select the metadata you want to use to compare:
nmeta<-ncol(md)+1
InsideLevels(Data[,1:nmeta-1])
if(exists("subset_data")==T){data <-subset_data}else{data <-Data[,1:nmeta-1]}
InsideLevels(data)

ColumnfromMetadata <- 14
Levels_Cdtn <- levels(droplevels(as.factor(data[,ColumnfromMetadata[1]])))
Levels_Cdtn 


### 2)   Separate the data based on these two conditions:
First<-Levels_Cdtn[1]
Second<-Levels_Cdtn[2]

DataFirst<-Data %>%
  filter(Data[ColumnfromMetadata]==First)

DataFirst<-DataFirst[,nmeta:ncol(DataFirst)]

DataSecond<-Data %>%
  filter(Data[ColumnfromMetadata]==Second)
DataSecond<-DataSecond[,nmeta:ncol(DataSecond)]

### 3)   Run t-test 
# Run t-test for the first contig
# Here, you can select if the test will be paired/unpaired and one/two sided
ttest<-NULL
for (i in 1:ncol(DataFirst)){
  aa<-DataFirst[,i]
  bb<-DataSecond[,i]
  temp<-mean(aa)/mean(bb)
  
  if (temp==1){
    x<-1
    ttest <- cbind(ttest, x)
  }
  else {
    x <- wilcox.test(aa, bb, paired = FALSE)
    ttest <- cbind(ttest, x)
  }}

# Create table with statistical values
temp<- data.frame(ttest)
temp1<-(sapply(temp[1,], as.numeric))
output_ttest <- data.frame(colnames(DataFirst[1:ncol(DataFirst)]))
colnames(output_ttest)[1] <- "stats_Metabolites"
output_ttest["p"] <- (sapply(temp[3,], as.numeric))
output_ttest["p_fdr"] <- p.adjust(output_ttest$p,method="fdr")
output_ttest["significant"] <- ifelse(output_ttest$p_fdr<0.0001,"Significant","Nonsignificant")


# Calculate fold change for the first contig
FCm=NULL
mUP=NULL
mDOWN=NULL
for (i in 1:ncol(DataFirst)){
  mUP[i]<-mean(DataFirst[,i])
  mDOWN[i]<-mean(DataSecond[,i])
  xx <- mUP[i]/mDOWN[i]
  FCm <- cbind(FCm, xx)
}
colnames(FCm)<-colnames(Data[nmeta:ncol(Data)])


# Call the condition where the fold change is higher/grouping for the first contig
Grouping=NULL
for(i in 1:ncol(FCm)) { 
  if (FCm[i]>1) {
    Grouping[i]<- First}
  else{
    Grouping[i] <- Second
  }
}


# Extract FC values to t-test output list for the first contig
output_ttest["FC"] <- (sapply(FCm, as.numeric))
output_ttest["log2FC"] <- log(abs(output_ttest$FC),base=2)
output_ttest["Group"] <- (sapply(Grouping, as.character))

output_ttest_final<-output_ttest


# Loop for the analysis of every contig
for (j in 15:ncol(md)){
      nmeta<-ncol(md)+1
      InsideLevels(Data[,1:nmeta-1])
      if(exists("subset_data")==T){data <-subset_data}else{data <-Data[,1:nmeta-1]}
      InsideLevels(data)
      
      ColumnfromMetadata <- j
      Levels_Cdtn <- levels(droplevels(as.factor(data[,ColumnfromMetadata[1]])))
      Levels_Cdtn 
      
      First<-Levels_Cdtn[1]
      Second<-Levels_Cdtn[2]
      
      DataFirst<-Data %>%
        filter(Data[ColumnfromMetadata]==First)
      
      DataFirst<-DataFirst[,nmeta:ncol(DataFirst)]
      
      DataSecond<-Data %>%
        filter(Data[ColumnfromMetadata]==Second)
      DataSecond<-DataSecond[,nmeta:ncol(DataSecond)]
      
      ttest<-NULL
      for (i in 1:ncol(DataFirst)){
        aa<-DataFirst[,i]
        bb<-DataSecond[,i]
        temp<-mean(aa)/mean(bb)
        
      x <- wilcox.test(aa, bb, paired = FALSE)
      ttest <- cbind(ttest, x)
      }
      
      temp<- data.frame(ttest)
      temp1<-(sapply(temp[1,], as.numeric))
      output_ttest <- data.frame(colnames(DataFirst[1:ncol(DataFirst)]))
      colnames(output_ttest)[1] <- "stats_Metabolites"
      output_ttest["p"] <- (sapply(temp[3,], as.numeric))
      output_ttest["p_fdr"] <- p.adjust(output_ttest$p,method="fdr")
      output_ttest["significant"] <- ifelse(output_ttest$p_fdr<0.0001,"Significant","Nonsignificant")
      
      FCm=NULL
      mUP=NULL
      mDOWN=NULL
      for (i in 1:ncol(DataFirst)){
        mUP[i]<-mean(DataFirst[,i])
        mDOWN[i]<-mean(DataSecond[,i])
        xx <- mUP[i]/mDOWN[i]
        FCm <- cbind(FCm, xx)
      }
      colnames(FCm)<-colnames(Data[nmeta:ncol(Data)])
      
      Grouping=NULL
      for(i in 1:ncol(FCm)) { 
        if (FCm[i]>1) {
          Grouping[i]<- First}
        else{
          Grouping[i] <- Second
  }
}


output_ttest["FC"] <- (sapply(FCm, as.numeric))
output_ttest["log2FC"] <- log(abs(output_ttest$FC),base=2)
output_ttest["Group"] <- (sapply(Grouping, as.character))

output_ttest<-output_ttest[,2:ncol(output_ttest)]
output_ttest_final<-cbind(output_ttest_final, output_ttest)

} 


### 4) Save t-test table
write.csv(output_ttest_final, file.path(Directory,"output_ttest_final_MRS.csv"),row.names =TRUE)


### 5) Get contig names
Levels_Cdtn=NULL
j=14

for (j in 14:ncol(md)){
  oo <-colnames(Data[j])
  Levels_Cdtn <-rbind(Levels_Cdtn, oo)
}

write.csv(Levels_Cdtn, file.path(Directory,"ContigNames.csv"),row.names =TRUE)
