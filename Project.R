####################################################################
# Name: Kush Patel
# Project
# cohort: GDC TCGA Cervical Cancer (CESC) - HTSeq - Counts
# Exploring Grade Gene Expression (G1 vs G3) Analysis in GDC TCGA Cervical Cancer (CESC) Using RNA-Seq Data
# https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Cervical%20Cancer%20(CESC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
####################################################################

# Libraries
library(dplyr)
library(ggplot2)
library(UCSCXenaTools) # needed to retrieve data
library(edgeR) # needed for processing, such as TMM
library(limma) # needed to find DE probes

##############################################################################
# 1) Retrieve Data
##############################################################################

data(XenaData)

#GDC TCGA Cervical Cancer (CESC)

# limit to desired cohort
cesc <- XenaData %>% filter(XenaCohorts == 'GDC TCGA Cervical Cancer (CESC)')

# Get the phenotype / clinical data
cli_query = cesc %>%
  filter(Label == "Phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%     # generate the query
  XenaDownload()      # download the data

# prepare (load) the data into R
cesc_pheno <- XenaPrepare(cli_query)

# Get the RNA-seq data, including the "probe map"
cli_query <- cesc %>% filter(Label == 'HTSeq - Counts') %>%
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload(download_probeMap = TRUE)

# prepare (load) the data into R
cesc_counts <- XenaPrepare(cli_query)

########################################################
# (2) Data pre-processing: we need to do a fair amount of 
# filtering and re-arranging to work with the data, so
# that the expression and phenotype data are aligned
########################################################

# First, let's use more manageable names
#   - X: expression data, with probes as row names
#   - probeMap: the probeMap
#   - Y: the pheno/clinical data

# for X, we need to set the rownames and remove the probe column
# from the data matrix
X <- data.frame(cesc_counts$TCGA.CESC.htseq_counts.tsv.gz)
rownames(X) <- X$Ensembl_ID
X <- X[,-1]  # remove the probe name column

# probeMap = probe names
probeMap <- cesc_counts$gencode.v22.annotation.gene.probeMap

# Y = pheno data
Y <- cesc_pheno
# The expression and clinical data need to match; currently, e.g.,
# The first column of the expression data does not correspond
# to the first row of the pheno data; the data is also not
# in a consistent format (one has '.' and the other has '-')

# compare sample names between X and Y; they do not match, and are 
# not even in the same format!
print(colnames(X)[1])
print(Y$submitter_id.samples[1])

# 'change '.' to '-' so sample ID format is consistent
colnames(X) <- gsub('\\.', '-', colnames(X))

# Note that the sample ID is a barcode that has a special meaning:
#   https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# In particular, the 4th section describes the 'Sample' which is 
#   either tumor (01 - 09) or normal (10-19). For details see:
#   https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

# Keep only the '01A' tumor samples
g <- grep('01A$', colnames(X))
X <- X[,g]

# We still need to match the expression data with the clinical data
# Let's do that by first finding the samples that are common
# between the expression and clinical data. We can use 
# intersect(a,b) to return a vector containing the elements common
# to vectors 'a' and 'b'

common_samples <- intersect(colnames(X), Y$submitter_id.samples)

# we then use match(x, t) to get a vector of indices. The value
# x[i] is the index of 't' containing the i^th value of 'x'

mx <- match(common_samples, colnames(X))
my <- match(common_samples, Y$submitter_id.samples)

X <- X[,mx]
Y <- Y[my,]

# Make sure that the samples match -- if they don't, this will produce an error
stopifnot(all(colnames(X) == Y$submitter_id.samples))

#Total number of samples and probes before processing
#"Total number of samples"
print(ncol(X))
#"Total number of probes"
print(nrow(X))


#############################################
# Setup step 3: Process the expression data
#############################################

# convert from log2(count + 1) to count data
X <- round(2**X - 1)

# remove genes with low counts
dge <- DGEList(counts=X)
keep <- filterByExpr(dge,min.prop = .10 )
dge <- dge[keep,,keep.lib.sizes=FALSE]

# apply TMM normalization, which computes the normalization 
# factors. The actual normalization is done in a later step
dge <- calcNormFactors(dge, method = "TMM")

# Calculate the log CPM values, using the normalization factors;
# 3 counts are added to each observation to prevent log 0 values
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

#Total number of samples and probes after processing
#"Total number of samples"
print(ncol(logCPM))
#Total number of probes"
print(nrow(logCPM))

#print(colnames(X)[1])
#print(Y$submitter_id.samples[1])


#############################################
#5 Box plot for the first 10 samples
#############################################

boxplot(logCPM[, 1:10], 
        main = "Boxplot of Normalized Log-CPM Values for First 10 Samples",
        ylab = "Log-CPM Values")

#############################################
#6 Done
#7 Extract the column for the grade
#############################################

grade <- Y$neoplasm_histologic_grade
grade[is.na(grade)] <- "Unknown"
design <- model.matrix(~-1+grade)
colnames(design) <- c("G1", "G2", "G3", "G4", "GX", "Unknown")
#design
columns_to_keep <- colnames(design)[!(colnames(design) %in% c("G2", "G4", "GX", "Unknown"))]
design <- design[, columns_to_keep]
head(design)

#############################################
#8 Use Limma to find probes across groups using fdr of 10%
#############################################

#fit the linear model to each row of the expression matrix
fit <- lmFit(logCPM, design)
#(G3 - G1)
contrast.matrix <- makeContrasts(G3 - G1, levels = design)
#fit model based on contrasts (e.g., G3 - G1)
fit <- contrasts.fit(fit, contrast.matrix)
#apply the 'eBayes' step to calculate moderated t statistics
fit.de <- eBayes(fit, trend = TRUE)
#FDR < %100
#FDR < %10 Resulted in 0 Values
tt.1 <- topTable(fit.de, sort.by = "p", p.value = 1, number = 30)
nrow(tt.1)
#Top 30 (All of them)
tt.1

#############################################
#9 Top probe boxplot - LOOK OVER
#############################################

probe <- rownames(tt.1)[1]
m <- match(probe, rownames(logCPM))

df <- data.frame(expr = logCPM[m,], grade = grade)

# filter out 'unknown' samples
df <- df[!(df$grade %in% c("G2", "G4", "GX","Unknown")), ]

# convert from logFC to FC #
logFC <- tt.1$logFC[1]
2**logFC

## visualize ##
FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of ", probe, ", ", FC, ", FDR = 100%")

ggplot(df, aes(x = grade, y = expr, fill = grade)) + geom_boxplot() +
  ylab("log2 expression") + ggtitle(main) +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_classic() + theme(legend.position = "none")

#############################################
#10 Df of all genes 
#############################################
# gene names, probe names, logFC, and adjusted p-values
probe_names <- rownames(tt.1)
m <- match(probe_names, probeMap$id)
gene_names <- probeMap$gene[m]
logFc_Values <- tt.1$logFC
p_Values <- tt.1$adj.P.Val

#DF
result_df <- data.frame(
  GeneNames = gene_names,
  ProbeNames = probe_names,
  LogFC = logFc_Values,
  AdjPValues = p_Values
)

head(result_df, 5)

#############################################
#11 Heat Map - Recheck
#############################################

m <- match(result_df$ProbeNames, rownames(logCPM))
expr <- logCPM[m,]
rownames(expr) <- result_df$GeneNames
# create a color range consisting of 200 values between yellow and blue
col.heat <- colorRampPalette(c("yellow", "blue"))(200)

# set colors for grade
col.grade <- as.integer(as.factor(!(grade %in% c("G2", "G4", "GX","Unknown"))))
col.grade <- c("magenta", "lightgreen")[col.grade]

# Generate the heatmap
heatmap(expr, ColSideColors = col.grade, col = col.heat, scale = "none")

#############################################
#12 DAVID
#############################################

#genes <- unique(result_df$GeneNames) # get unique set of genes (removes duplicates)
#genes
#write.table(genes, row.names = FALSE, quote = FALSE, file = "top_genes.txt")

