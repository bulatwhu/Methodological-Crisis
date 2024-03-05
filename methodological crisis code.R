library(limma)
library(dplyr)
library(GEOquery)
library(sva)
library(R.utils)
library(tibble)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

packageVersion("limma")
packageVersion("sva")
packageVersion("GEOquery")
packageVersion("EnhancedVolcano")
packageVersion("pheatmap")


setwd("~/GEO/Meta analysis/GSE98770")
#GSE98770 (Agilent)
GSE98770 <- getGEO("GSE98770", destdir = '.')
GSE98770 <- GSE98770[[1]]

#phenodata
P98770 <- pData(GSE98770)
P98770$AAD <- P98770$source_name_ch1
P98770 <- P98770 %>%
  mutate(AAD = ifelse(AAD == "Intima-media ascending aorta, ataad patient", "1", 
                        ifelse(AAD == "Intima-media ascending aorta, donor", "0", AAD)))
P98770$Dataset <- 'GSE98770'
P98770 <- P98770[, c(2, 39, 40)]

#SKIP getGEOSuppFiles("GSE98770", destdir = '.') 
#SKIP untar("GSE98770/GSE98770_RAW.tar", exdir = ".")

GPL14550 <- getGEO("GPL14550", destdir = ".")
GPL14550 <- Table(GPL14550)

#limma
SDRF <- pData(phenoData(GSE98770))[,c(2,10,34)]
SDRF <- SDRF %>%
  mutate(supplementary_file = gsub("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2611nnn/", "", supplementary_file)) %>%
  mutate (characteristics_ch1 = gsub ("tissue: ", "", characteristics_ch1)) %>%
  mutate(supplementary_file = gsub("GSM2611536_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611537_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611538_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611539_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611540_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611541_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611542_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611542_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611543_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611544_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611545_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("GSM2611546_", "_", supplementary_file)) %>%
  mutate(supplementary_file = gsub("/suppl/", "", supplementary_file))

intens98770 <- read.maimages(SDRF[,"supplementary_file"], source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")
dim(intens98770)


#boxplot before
fe2 <- intens98770$E
colnames(fe2) <- substr(colnames(fe2), 1, 10)
boxplot(log2(fe2), range=0,main = "GSE98770 (raw data)", cex.axis=0.8,ylab="log2 (intensity)")
boxplot(fe2, main = "GSE98770 (raw data)", ylab="intensity")
summary(fe2)

#BC and normalization
n98770 <- backgroundCorrect(intens98770, method="normexp")
n98770 <- normalizeBetweenArrays(n98770, method="quantile")

#box plot after
fe <- n98770$E
colnames(fe) <- substr(colnames(fe), 1, 10)
summary(fe)
boxplot(fe, range=0, main = "GSE98770 after BC, N, and log2",cex.axis=0.8, ylab="intensity")


#gene filtering
Control <- n98770$genes$ControlType==1L
IsExpr <- rowSums(n98770$other$gIsWellAboveBG > 0) >= 6
table(IsExpr)
f98770 <- n98770[!Control & IsExpr, ]
#gene annotation
data_matrix98770 <- cbind(f98770$E, f98770$genes)
data_matrix98770 <- cbind(data_matrix98770, GPL14550[match(data_matrix98770$ProbeName, GPL14550$ID),1:8])
# Remove rows with empty GENE_SYMBOL
empty_rows <- which(data_matrix98770$GENE_SYMBOL == "")
data_matrix98770 <- data_matrix98770[-empty_rows, ]
data_matrix98770 <- data_matrix98770[, -c(12:22, 24)]
#remove duplicate gene symbol rows, unspecific probes and get a mean value for each cell 
data_matrix98770 <- data_matrix98770 %>%
  group_by(GENE_SYMBOL) %>%
  summarise(across(everything(), \(x) median(x, na.rm = TRUE))) %>%
  ungroup()
colnames(data_matrix98770) <- substr(colnames(data_matrix98770), 1, 11)

boxplot(data_matrix98770[,-1], range=0, main = "Normalized GSE98770 after filtering",cex.axis=0.8, ylab="intensity")


#Identification of DEGs 98770
design98770 <- model.matrix(~ P98770$AAD)
fit98770 <- lmFit(data_matrix98770,design98770)
fit98770 <- eBayes(fit98770,trend=TRUE,robust=TRUE)

diff_98770 <- topTable(fit98770, coef = 2, adjust = "fdr", number = Inf)

logFC = 1

expression98770=ifelse(diff_98770$adj.P.Val>0.05,'stable',
                  ifelse( diff_98770$logFC > logFC,'upregulated',
                          ifelse(diff_98770$logFC < -logFC,'downregulated','stable') ))

expression98770 <- mutate(diff_98770, expression98770)
head(expression98770)
table(expression98770$expression)

DEGs_98770 <- expression98770[expression98770$expression != "stable", ]

#VOLCANO_PLOT
EnhancedVolcano(expression98770,
  lab = expression98770$GENE_SYMBOL,  
  x = 'logFC',             #Log2 fold change on X
  y = 'adj.P.Val',           #-log10(p-value) on Y
  title = "GSE98770",
  pCutoff = 0.05,         
  FCcutoff = logFC,     
  pointSize = 3,           
  xlim = c(-6, 6),         # Limit for X
  ylim = c(0, 6),         # Limit for Y
  legendLabels = c("NS", "|log2FC|>1", "adj. P-value<0.05", "|log2FC|>1, adj. P-value<0.05"), 
  legendPosition = "top",  
  legendLabSize = 12       
)



#GSE52093 (Illumina)
setwd("~/GEO/Meta analysis/GSE52093")

GSE52093 <- getGEO("GSE52093", destdir = '.')
GSE52093 <- GSE52093[[1]]
GPL10558 <- getGEO("GPL10558", destdir = ".")
GPL10558 <- Table(GPL10558)

#phenodata
P52093 <- pData(GSE52093)
P52093$AAD <- P52093$source_name_ch1
P52093 <- P52093 %>%
  mutate(AAD = ifelse(AAD == "Aorta dissected", "1", 
                        ifelse(AAD == "Aorta normal", "0", AAD)))
P52093$Dataset <- 'GSE52093'
P52093 <- P52093[, c(2, 36, 37)]


#SKIP getGEOSuppFiles("GSE52093", destdir = '.') 
#SKIP untar("GSE52093/GSE52093_RAW.tar", exdir = ".")

#read files
intens52093 <- read.ilmn("GSE52093/GSE52093_non-normalized.txt.gz", probeid="ID_REF", expr = "SAMPLE")

#boxplot before normalization
fe4 <- as.data.frame(intens52093$E)
colnames(fe4)[1:12] <- c("GSM1259275", "GSM1259276", "GSM1259277","GSM1259278","GSM1259279",
                                "GSM1259280","GSM1259281","GSM1259282","GSM1259283", "GSM1259284",
                                "GSM1259285","GSM1259286")
boxplot(log2(fe4), main="GSE52093 (raw data)", range=0, cex.axis = 0.8, ylab="log2 (intensity)")
summary(fe4)

#normalization
n52093 <- neqc(intens52093)
fe3 <- as.data.frame(n52093$E)
colnames(fe3)[1:12] <- c("GSM1259275", "GSM1259276", "GSM1259277","GSM1259278","GSM1259279",
                                 "GSM1259280","GSM1259281","GSM1259282","GSM1259283", "GSM1259284",
                                 "GSM1259285","GSM1259286")
boxplot(fe3, main="GSE52093 after BC, N, and log2", range=0, cex.axis = 0.8, ylab="intensity")

#filtering
expressed <- rowSums(n52093$other$Detection < 0.05) >= 6
n52093 <- n52093[expressed,]
n52093 <- as.data.frame(n52093)

#gene annotation
symbol = GPL10558[,c(1,13)]
symbol$probe_id <- as.character(symbol$ID)
symbol <- symbol[,-1]
colnames(symbol)[1] <- 'GENE_SYMBOL'
symbol <- symbol[!grepl("///", symbol$GENE_SYMBOL),]
n52093$probe_id <- as.character(rownames(n52093))
data_matrix52093 <- inner_join(symbol, n52093, by = "probe_id")
data_matrix52093 <- data_matrix52093[,-2]

#remove duplicate gene symbol rows, unspecific probes and get a mean value for each cell 
f52093 <- data_matrix52093 %>%
  group_by(GENE_SYMBOL) %>%
  summarise(across(everything(), \(x) median(x, na.rm = TRUE))) %>%
  ungroup()
f52093 <- f52093[-1,]

#boxplot after norm
colnames(f52093)[2:13] <- c("GSM1259275", "GSM1259276", "GSM1259277","GSM1259278","GSM1259279",
                         "GSM1259280","GSM1259281","GSM1259282","GSM1259283", "GSM1259284",
                         "GSM1259285","GSM1259286")
boxplot(f52093[,-1], main="Normalized GSE52093 after filtering", range=0, cex.axis = 0.8, ylab="intensity")


#DEGs
design52093 <- model.matrix(~ P52093$AAD)
fit52093 <- lmFit(f52093,design52093)
fit52093 <- eBayes(fit52093,trend=TRUE,robust=TRUE)
diff_52093 <- topTable(fit52093, coef = 2, adjust = "fdr", number = Inf)

boxplot(f52093[,-1], range=0, main = "Normalized GSE52093 after filtering",cex.axis=0.8, ylab="intensity")

expression52093=ifelse(diff_52093$adj.P.Val>0.05,'stable',
                       ifelse( diff_52093$logFC > logFC,'upregulated',
                               ifelse(diff_52093$logFC < -logFC,'downregulated','stable') ))

expression52093 <- mutate(diff_52093, expression52093)
head(expression52093)
table(expression52093$expression)

DEGs_52093 <- expression52093[expression52093$expression != "stable", ]

#VOLCANO_PLOT
EnhancedVolcano(expression52093,
                lab = expression52093$GENE_SYMBOL,  
                x = 'logFC',             # Log2 fold change on X
                y = 'adj.P.Val',           # -log10(p-value) on Y
                title = "GSE52093",
                pCutoff = 0.05,          
                FCcutoff = 1,    
                pointSize = 3,           
                xlim = c(-5, 5),         # Limit for X
                ylim = c(0, 5),         # Limit for Y
                legendLabels = c("NS", "|log2FC|>1", "adj. P-value<0.05", "|log2FC|>1, adj. P-value<0.05"), 
                legendPosition = "top", 
                legendLabSize = 12)


#setwd("~/GEO/Meta analysis/Merged")
#merge datasets f52093 and data_matrix98770
merged_exprs <- merge(f52093, data_matrix98770, by = "GENE_SYMBOL")
row.names(merged_exprs) <- merged_exprs$GENE_SYMBOL
merged_exprs <- merged_exprs[,-1]
x <- merged_exprs

#combined pData
merged_pData <- rbind(P52093, P98770)

#check for batch effects
pca_merged <- prcomp(t(merged_exprs), scale. = TRUE)
pca1_var <- round(summary(pca_merged)$importance[2, 1], 5)
pca2_var <- round(summary(pca_merged)$importance[2, 2], 5)

pca_merged <- as.data.frame(pca_merged$x[, 1:2])  
pca_merged$Dataset <- merged_pData$Dataset 
pca_merged$AAD <- merged_pData$AAD
ggplot(pca_merged, aes(PC1, PC2, color = Dataset, shape = AAD)) + geom_point(size = 3) +
  geom_point() +
  labs(title = "Merged GSE52093&GSE98770 before adjusting for batch effects", x = "PC1 (68.50%)", y = "PC2 (7.13%)")

boxplot(merged_exprs, range=0, main = "Merged GSE52093&GSE98770 before adjusting for batch effects",cex.axis=0.5, ylab="intensity")

# Remove batch effects (combat)
merged_pData$Dataset <- as.factor(merged_pData$Dataset)
combat_x <- ComBat(dat = x, batch = merged_pData$Dataset)

# check for batch effect after combat
pca_combat <- prcomp(t(combat_x), scale. = TRUE) 
pca1_var_combat <- round(summary(pca_combat)$importance[2, 1], 5)
pca2_var_combat <- round(summary(pca_combat)$importance[2, 2], 5)

pca_combat <- as.data.frame(pca_combat$x[, 1:2])  
pca_combat$Dataset <- merged_pData$Dataset
pca_combat$AAD <- merged_pData$AAD

ggplot(pca_combat, aes(PC1, PC2, color = Dataset, shape = AAD)) + geom_point(size = 3) +
  geom_point() +
  labs(title = "Merged GSE52093&GSE98770 after adjusting for batch effects", x = "PC1 (22.21%)", y = "PC2 (15.04%)")

boxplot(combat_x, range=0, main = "Merged GSE52093&GSE98770 after adjusting for batch effects",cex.axis=0.5, ylab="intensity")


#DEGs
design_merge <- model.matrix(~ merged_pData$AAD)
fit_merge <- lmFit(combat_x,design_merge)
fit_merge <- eBayes(fit_merge,trend=TRUE,robust=TRUE)
diff_merge <- topTable(fit_merge, coef = 2, adjust = "fdr", number = Inf)

expression_merge=ifelse(diff_merge$adj.P.Val>0.05,'stable',
                       ifelse( diff_merge$logFC > logFC,'upregulated',
                               ifelse(diff_merge$logFC < -logFC,'downregulated','stable') ))

expression_merge <- mutate(diff_merge, expression_merge)
head(expression_merge)
table(expression_merge$expression)

DEGs_merge <- expression_merge[expression_merge$expression != "stable", ]

#for volcano and heatmaps
em <- rownames_to_column(expression_merge, var = "GENE_SYMBOL") #DEGs

#Volcano plot
EnhancedVolcano(em,
                lab = em$GENE_SYMBOL,  
                x = 'logFC',             # Log2 fold change on X
                y = 'adj.P.Val',           # -log10(p-value) on Y
                title = "Merged GSE52093&GSE98770",
                pCutoff = 0.05,         
                FCcutoff = logFC,     
                pointSize = 3,           # Size of points
                xlim = c(-5, 5),         # Limit for X
                ylim = c(0, 10),         # Limit for Y
                legendLabels = c("NS", "|log2FC|>1", "adj. P-value<0.05", "|log2FC|>1, adj. P-value<0.05"), 
                legendPosition = "top",  
                legendLabSize = 12)


#HEATMAP
y <- rownames_to_column(as.data.frame(combat_x), var = "GENE_SYMBOL") #gene expression
eh <- rownames_to_column(DEGs_merge, var = "GENE_SYMBOL")

# Merge by the row names
merged_data <- merge(eh, y, by = "GENE_SYMBOL", all.x = TRUE)
row.names(merged_data) <- merged_data$GENE_SYMBOL
merged_data <- merged_data[order(merged_data$adj.P.Val),]

#Select and include the top 50 DEGs
top_50_DEGs <- merged_data$GENE_SYMBOL[1:50]

non_matching_symbols <- setdiff(top_50_DEGs, rownames(merged_data))
if (length(non_matching_symbols) > 0) {
  print(paste("Gene symbols not found in merged_data:", paste(non_matching_symbols, collapse = ", ")))
} else {
  top_50_DEGs_data <- merged_data[top_50_DEGs, ]
}

top_50_DEGs_data_matrix <- as.matrix(top_50_DEGs_data[,c(9:31)])

#create a single row (phenodata)
pheno_row <- t(merged_pData[, "AAD", drop = FALSE])
#match column names
merged_data <- merged_data[, c(9:31)]
colnames(pheno_row) <- colnames(merged_data)
# create groups
groups <- as.character(pheno_row[1, ])
# create a data file with top50 DEGs as rown and all samples from the merged dataset as columns
AAD_samples <- sample(colnames(top_50_DEGs_data_matrix)[groups == "1"])
control_samples <- sample(colnames(top_50_DEGs_data_matrix)[groups == "0"])
selected_samples <- c(AAD_samples, control_samples)
selected_data <- top_50_DEGs_data_matrix[, selected_samples]
#create a vector
genes <- names(tail(sort(apply(selected_data,1,sd)),50))
genesexpr <- selected_data[genes, ]
# annotation
annotation <- data.frame(group = t(pheno_row))

# Create the heatmap (sorted by similarity)
color_gradient <- colorRampPalette(c("red", "white", "blue"))(100)
heatmap_plot <- pheatmap(genesexpr,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         annotation_col = annotation,
                         main = "Merged GSE52093&GSE98770",
                         scale = "row",
                         color = color_gradient,
                         cluster_cols = FALSE)  # Disable clustering of columns (samples)  

