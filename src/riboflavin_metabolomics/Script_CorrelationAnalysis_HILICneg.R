###################################### Updates, install and call packages 
# Install packages
install.packages("Hmisc")
install.packages("pheatmap")

# Load necessary libraries
library(pheatmap)
library(Hmisc)
library(dplyr)
library(ggplot2)

Sys.time()
rm(list=ls())

###################################### Input and process files:
Directory<- c("C:/Users/selegato/Desktop/R projects/Projects/Antwerp_Caroline_BacterialCulture_HILICMSneg_targeted/")
setwd(Directory)

# We need to import the feature table and metadata
ft_url <- paste0(Directory, "/Quant_Table_20250814.csv")
# metadata
md_url <-  paste0(Directory, "/metadata_20250703.csv")


# read data
ft <- read.csv(ft_url, header = T, check.names = F, row.names = 1)
md <- read.csv(md_url, header = T, check.names = F, row.names = 1)

# create results folder
dirs <- dir(path=paste(getwd(), sep=""), full.names=TRUE, recursive=TRUE)
folders <- unique(dirname(dirs))
files <- list.files(folders, full.names=TRUE)
files_1 <- basename((files))
files_2 <- dirname((files))
# Creating a Result folder
dir.create(path=paste(files_2[[1]], "/MetaboliteCorrelation_Results", sep=""), showWarnings = TRUE)
fName <-paste(files_2[[1]], "/MetaboliteCorrelation_Results", sep="")


###################################### Data processing
# Remove noise values (bellow 30000)
ft[ft < 30000] <- 0.1

dim(ft)# dimension of the data
dim(md)# dimension of the data
head(ft)# how it looks
head(md)# how it looks

# Bring feature table and metadata in the correct format:
# how many files in the metadata are also present in the feature table
table(rownames(md) %in% colnames(ft))
# which file names in the metadata are not in the feature table?
setdiff(rownames(md),colnames(ft))
md <- md[rownames(md) %in% colnames(ft),]
dim(md)

# transpose ft
ft <- t(ft)
head(ft)
dim(ft)

# It is important that the filenames in our metadata are identical as well 
# as in the same order as the filenames in  our feature table. Let's make 
# sure this is true, using the below code (this should return TRUE).
identical(rownames(ft),rownames(md))
# put the rows in the feature table and metadata in the same order
ft <- ft[match(rownames(md),rownames(ft)),]
identical(rownames(ft),rownames(md))

# Merge feature table and metadata to one dataframe
ft<-data.frame(ft)
DataI <- cbind.data.frame(md,ft)
dim(DataI)

# Filter samples based on metadata
Data<-DataI %>%
  filter(ATTRIBUTE_RELATED=="Related") %>%
  filter(ATTRIBUTE_STRAIN!="Waterblank")%>%
  filter(ATTRIBUTE_STRAIN!="V370")%>%
  filter(ATTRIBUTE_MEDIA!="MRS-R")%>%
  filter(ATTRIBUTE_PRODUCER!="COMMERCIAL_STRAINS_CONSUMER")%>%
  filter(ATTRIBUTE_PRODUCER!="BAD_GROWERS_Streptococcaceae")%>%
  filter(ATTRIBUTE_RELATED!="Commercial Strains")%>%
  filter(ATTRIBUTE_STRAIN=="V339")



###################################### Compute correlations and p-values
# rcorr gives correlation matrix and p-values
nmeta<-ncol(md)+1
data_t<-Data[,nmeta:ncol(Data)]

corr_results <- rcorr(as.matrix(data_t), type = "pearson") # or spearman
cor_matrix <- corr_results$r # Correlation coefficients
cor_matrix[is.na(cor_matrix )] <- 0
write.csv(cor_matrix, file.path(fName,"cor_matrix.csv"),row.names =TRUE)


p_matrix <- corr_results$P     # P-values
p_matrix [is.na(p_matrix )] <- 0
write.csv(p_matrix, file.path(fName,"p_matrix.csv"),row.names =TRUE)


###################################### Create matrix of significance labels
sig_labels <- ifelse(p_matrix < 0.05, "*", "")
diag(sig_labels) <- ""  # Remove diagonal labels



###################################### Filter the metabolites you want to see on the heatmap
# Define which metabolites to keep
selected_metabolites1 <- c( "Riboflavine_375.13045_2.269",
                            "RiboflavineCl_411.1075_2.5",
                            "X4.Aminobenzoic.acid_136.03988_2.752",
                            "X5..2.oxoethylideneamino..6.D.ribitylaminouracil..5.OE.RU._315.0941_2.5",
                            "X5..2.oxopropylideneamino..6.D.ribitylaminouracil..5.OP.RU._329.10975_2.5",
                            "X5.6.Dimethylbenzimidazole_131.07349_1.98",
                            "X5.amino.6..D.ribitylamino.uracil..5.A.RU._275.09918_1.5",
                            "X6.7.Dimethy8.ribityllumazine..R6.7.diMe._325.1148_4.13",
                            "X7.Hydroxy.6.methy8..1.D.ribityl.lumazine..R6.me.7.OH._327.0941_5",
                            "D.pantothenic.acid_218.10287_4.717",
                            "Hypoxanthine_135.03071_2.553",
                            "Indole.3.carboxaldehyde_144.04496_1.422",
                            "Indole.3.lactic.acid_204.06606_3.254",
                            "Methionine_148.04325_5.34",
                            "Quinolinic.acid_166.014_11.22",
                            "Riboflavin.5.Monophosphate_477.07873_10.377",
                            "Xanthine_151.02563_2.537",
                            "Lactate_90.03169_3.78")


selected_metabolites2 <- c( "Riboflavine_375.13045_2.269",
                            "RiboflavineCl_411.1075_2.5",
                            "X4.Aminobenzoic.acid_136.03988_2.752",
                            "X5..2.oxoethylideneamino..6.D.ribitylaminouracil..5.OE.RU._315.0941_2.5",
                            "X5..2.oxopropylideneamino..6.D.ribitylaminouracil..5.OP.RU._329.10975_2.5",
                            "X5.6.Dimethylbenzimidazole_131.07349_1.98",
                            "X5.amino.6..D.ribitylamino.uracil..5.A.RU._275.09918_1.5",
                            "X6.7.Dimethy8.ribityllumazine..R6.7.diMe._325.1148_4.13",
                            "X7.Hydroxy.6.methy8..1.D.ribityl.lumazine..R6.me.7.OH._327.0941_5",
                            "D.pantothenic.acid_218.10287_4.717",
                            "Hypoxanthine_135.03071_2.553",
                            "Indole.3.carboxaldehyde_144.04496_1.422",
                            "Indole.3.lactic.acid_204.06606_3.254",
                            "Methionine_148.04325_5.34",
                            "Quinolinic.acid_166.014_11.22",
                            "Riboflavin.5.Monophosphate_477.07873_10.377",
                            "Xanthine_151.02563_2.537",
                            "Lactate_90.03169_3.78")
                            
                            

only_in_A <- setdiff(selected_metabolites2, colnames(cor_matrix))
only_in_B <- setdiff(selected_metabolites1, colnames(cor_matrix))

# Subset both rows and columns to maintain square matrix
sub_cor_matrix <- cor_matrix[selected_metabolites1,selected_metabolites2]
sub_p_matrix <- p_matrix[selected_metabolites1,selected_metabolites2]
sig_labels2 <- ifelse(sub_p_matrix < 0.05, "*", "")
#diag(sig_labels) <- ""  # Remove diagonal labels
#sub_cor_matrix <- cor_matrix[selected_metabolites1,]

# Subset significance labels too (if using them)
sub_sig_labels <- sig_labels[selected_metabolites1,selected_metabolites2]
#sub_sig_labels <- sig_labels[selected_metabolites1,]


###################################### Plot heatmap with significance
# pheatmap allows annotation with display_numbers
# Plot heatmap for the subset
heatmapplot<-pheatmap(sub_cor_matrix, 
         display_numbers = sig_labels2,
         number_color = "black",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         breaks = seq(-1, 1, length.out = 51),
         main = "Subset Correlation Heatmap",
         fontsize = 8,
         fontsize_number = 10,
         cluster_rows = F,
         cluster_cols = TRUE)
heatmapplot

#Save
namesave=paste0("/Heatmap_Correlation.png")
diir=paste0(fName,namesave)
ggsave(heatmapplot, file=diir, width = 41, height = 20, units = "cm", bg = "white")

namesave=paste0("/Heatmap_Correlation.svg")
diir=paste0(fName,namesave)
ggsave(heatmapplot, file=diir, width = 41, height = 20, units = "cm", bg = "white")
