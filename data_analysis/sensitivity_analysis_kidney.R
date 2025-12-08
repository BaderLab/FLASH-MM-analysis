# Load required libraries
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library('FLASHMM')

# Load data
datafile <- "~/FLASH-MM/Kidneydata/data/Kidney_raw_counts_decontaminated.rds"
dat <- readRDS(file = datafile)
counts <- GetAssayData(dat, layer = 'counts', assay = 'RNA')
coldata <- dat@meta.data
stopifnot(all(rownames(coldata) == colnames(counts)))

# Filtering parameters to evaluate
param_sets <- list(
  stringent = list(minCells = 32, minCellsgrp = 15, minCounts = 2^7, maxCounts = 2^19, mincpc = 0.01),
  default = list(minCells = 16, minCellsgrp = 10, minCounts = 2^6, maxCounts = 2^20, mincpc = 0.005)
  #lenient = list(minCells = 8, minCellsgrp = 5, minCounts = 2^5, maxCounts = 2^21, mincpc = 0.0025),
)

for (tag in names(param_sets)) {
  cat("\nRunning filtering and LMM for:", tag, "\n")
  pars <- param_sets[[tag]]
  
  ## Filter cells
  nFeature <- colSums(counts > 0)
  libsize <- colSums(counts)
  nCellsType <- table(coldata$Cell_Types_Broad)
  hist(libsize)
  abline(v = 2^9, col = "red", lwd = 2, lty = 3)
  abline(v = 2^16, col = "red", lwd = 2, lty = 3)
  dev.off()
  
  j <- (nFeature >= 100) & (libsize >= 2^9) & (libsize <= 2^16)
  j <- j & (coldata$Cell_Types_Broad %in% names(nCellsType)[nCellsType >= 20])
  counts_filt <- counts[, j] ### filtering cells is fixed
  coldata_filt <- coldata[j, ]
  stopifnot(all(colnames(counts_filt) == rownames(coldata_filt)))
  
  ## Filter genes
  nCells <- rowSums(counts_filt >= 1)
  nCounts <- rowSums(counts_filt)
  cpc <- nCounts / ncol(counts_filt) ### counts per cell
  
  nCellsgrp <- do.call(cbind, 
                       tapply(1:ncol(counts_filt), as.factor(coldata_filt$sex), 
                              function(j) rowSums(counts_filt[, j, drop = FALSE] >= 1)))
  
  index <- (nCells >= pars$minCells) & 
    (rowSums(nCellsgrp >= pars$minCellsgrp) >= ncol(nCellsgrp)) &
    (nCounts >= pars$minCounts) & (nCounts <= pars$maxCounts) &
    (cpc > pars$mincpc)
  
  counts_use <- counts_filt[index, ]
  coldata_use <- coldata_filt
  stopifnot(all(colnames(counts_use) == rownames(coldata_use)))
  
  dim(counts_use) # 16201 27550
  ## Run LMM
  Y <- round(counts_use)
  loglib <- log(colSums(Y))
  ngrp <- 6
  sizegrp <- ceiling(nrow(Y)/ngrp)
  for (i in 1:ngrp) {
    j <- (1+(i-1)*sizegrp):(min(nrow(Y), i*sizegrp))
    Y[j, ] <- log2(1 + Y[j, ])
  }
  
  X <- model.matrix(~ loglib + Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata_use)
  x = colnames(X)
  x <- gsub("Cell_Types_Broad", "", x)
  x <- gsub("sex", "", x)
  x <- gsub("\\+", "p", x)
  x <- gsub("-", "_", x)
  x <- gsub(" ", ".", x)
  colnames(X) = x
  #colnames(X) <- gsub("Cell_Types_Broad|sex|\\+|\\-| ", c("", "", "p", "_", "."), colnames(X), perl=TRUE)
  Z <- model.matrix(~ 0 + sampleID, data = coldata_use)
  d <- ncol(Z)
  epsilon <- 1e-5
  max.iter <- 100
  
  t1 <- Sys.time()
  fit <- lmmfit(Y, X, Z, d=d, max.iter=max.iter, epsilon=epsilon)
  t2 <- Sys.time()
  rtlmm <- difftime(t2, t1, units = "secs")
  
  outname <- paste0("~/FLASHMM-analysis/Results_data/sensitivity/kidney-counts-lmm-", tag, ".RData")
  save(fit, epsilon, max.iter, rtlmm, file = outname)
  cat("Saved results to:", outname, "\n")
}


load("~/FLASHMM-analysis/Results_data/sensitivity//kidney-counts-lmm-default.RData")   
fit_default <- fit
load("~/FLASHMM-analysis/Results_data/sensitivity//kidney-counts-lmm-lenient.RData")   
fit_lenient <- fit
load("~/FLASHMM-analysis/Results_data/sensitivity//kidney-counts-lmm-stringent.RData") 
fit_stringent <- fit


genes_default <- rownames(t(fit_default$coef))
genes_lenient <- rownames(t(fit_lenient$coef)) 
genes_stringent <- rownames(t(fit_stringent$coef))

common_genes <- Reduce(intersect, list(genes_default, genes_lenient, genes_stringent))
length(common_genes)


term <- "CNT:Male"

b_default <- t(fit_default$coef)[common_genes, term] 
b_lenient <- t(fit_lenient$coef)[common_genes, term]
b_stringent <- t(fit_stringent$coef)[common_genes, term]

#### correlation between the coefficients for CNT:Male based on variable thresholds
cor(b_default, b_lenient)
cor(b_default, b_stringent)
cor(b_lenient, b_stringent)

# > cor(b_default, b_lenient)
# [1] 0.9999977
# > cor(b_default, b_stringent)
# [1] 0.9999956
# > cor(b_lenient, b_stringent)
# [1] 0.9999976

dev.off()
par(mfrow = c(1, 2))  # 1 row, 2 columns

# Plot 1: Default vs Lenient
plot(b_default, b_lenient, pch = 16, cex = 0.5, col = "#00000044",
     xlab = "Effect size (Default thresholds)",
     ylab = "Effect size (Lenient thresholds)",
     main = paste(term, "\n(Default vs. Lenient filtering)"))
abline(0, 1, col = "red")

# Plot 2: Default vs Stringent
plot(b_default, b_stringent, pch = 16, cex = 0.5, col = "#00000044",
     xlab = "Effect size (Default thresholds)",
     ylab = "Effect size (Stringent thresholds)",
     main = paste(term, "\n(Default vs. Stringent filtering)"))
abline(0, 1, col = "red")

par(mfrow = c(1, 1))  # reset to default layout

