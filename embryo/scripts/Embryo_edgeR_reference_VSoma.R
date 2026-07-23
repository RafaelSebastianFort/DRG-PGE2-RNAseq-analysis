# Differential expression analysis: SARTools/edgeR
# Reference group: VSoma

rm(list=ls())

workDir <- "results/embryo_edgeR/reference_VSoma"
projectName <- "Embryo_edgeR_reference_VSoma"
author <- "RAFAEL_FORT"

targetFile <- "../../../metadata/targets_embryo.txt"
rawDir <- "../../../data/counts_embryo"

featuresToRemove <- c("alignment_not_unique",
                      "ambiguous", "no_feature",
                      "not_aligned", "too_low_aQual")

varInt <- "group"
condRef <- "VSoma"
batch <- NULL

alpha <- 0.05
pAdjustMethod <- "BH"

cpmCutoff <- 1
gene.selection <- "pairwise"
normalizationMethod <- "TMM"

colors <- c("dodgerblue","firebrick1",
            "MediumVioletRed","SpringGreen")

forceCairoGraph <- FALSE

dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
setwd(workDir)

library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

save.image(file=paste0(projectName, ".RData"))

writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=getwd(), projectName=projectName, author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, cpmCutoff=cpmCutoff,
                  colors=colors, gene.selection=gene.selection, normalizationMethod=normalizationMethod)
