## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options(digits=3)

## ----setup,message=FALSE,warning=FALSE----------------------------------------
library(langevitour)
library(airway)               # airway dataset
library(edgeR)                # RPM calculation
library(limma)                # makeContrasts
library(MASS)                 # ginv generalized matrix inverse
library(GPArotation)          # Bentler rotation
library(EnsDb.Hsapiens.v86)   # Gene names
library(ggplot2)
library(dplyr)
library(tibble)

data(airway)

treatment <- colData(airway)$dex == "trt"
cell <- factor(c(1,1,2,2,3,3,4,4))
design <- model.matrix(~ 0 + cell + treatment)

dge <- airway |>
    assay("counts") |>
    DGEList() |>
    calcNormFactors()

# Convert to log2 Reads Per Million.
# prior.count=5 applies some moderation for counts near zero.
rpms <- cpm(dge, log=TRUE, prior.count=5)

# Only show variable genes (mostly for speed)
keep <- apply(rpms,1,sd) >= 0.5
table(keep)
y <- rpms[keep,,drop=F]

# Use shorter sample names
colnames(y) <- paste0(ifelse(treatment,"T","U"), cell)

# Get friendly gene names
symbols <- 
    AnnotationDbi::select(EnsDb.Hsapiens.v86, keys=rownames(y), keytype="GENEID", columns="SYMBOL") |>
    deframe()
name <- symbols[rownames(y)]
name[is.na(name)] <- rownames(y)[is.na(name)]

# Colors for the samples
colors <- ifelse(treatment,"#f00","#080")

## -----------------------------------------------------------------------------
y[1:5,]

## ----contrasts----------------------------------------------------------------
coefficient_estimator <- MASS::ginv(design)

contrasts <- makeContrasts(
        average=(cell1+cell2+cell3+cell4)/4+treatmentTRUE/2,
        treatment=treatmentTRUE,
        "cell1 vs others" = cell1-(cell2+cell3+cell4)/3,
        "cell2 vs others" = cell2-(cell1+cell3+cell4)/3,
        "cell3 vs others" = cell3-(cell1+cell2+cell4)/3,
        "cell4 vs others" = cell4-(cell1+cell2+cell3)/3,
        levels=design)

contrastAxes <- t(coefficient_estimator) %*% contrasts

contrastAxes

## ----components---------------------------------------------------------------
y_centered <- sweep(y, 1, rowMeans(y), "-")

pca <- prcomp(y_centered, scale=FALSE, rank=4)
plot(pca)
pcaAxes <- pca$rotation

bentlerAxes <- pca$rotation %*% bentlerT(pca$x)$Th
colnames(bentlerAxes) <- paste0("Bentler",seq_len(ncol(bentlerAxes)))

## ----plot---------------------------------------------------------------------
langevitour(
    y, scale=15, axisColor=colors, name=name, 
    extraAxes=cbind(contrastAxes, pcaAxes, bentlerAxes))

