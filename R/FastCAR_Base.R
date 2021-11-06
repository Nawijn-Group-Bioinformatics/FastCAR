###############################################################################
# FastCAR package
# Marijn Berg
# m.berg@umcg.nl
###############################################################################
# FastCAR removes ambient RNA from 10X data by looking at the profile of the
# empty droplets and removing per gene the highest value found in the ambient
# RNA if it's found to contaminate more than a certain threshold of droplets
###############################################################################
# TODO
#
###############################################################################
library(Matrix)
library(Seurat)
library(qlcMatrix)
library(pheatmap)
library(ggplot2)
library(gridExtra)

###############################################################################

# had to make this function to efficiently modify sparse matrices on a per gene basis
# A dgCMatrix object has the following elements that matter for this function
# i: a sequence of the row locations of each non-zero element
# x: the non-zero values in the matrix
# p: the where in 'i' and 'x' a new column starts
# Dimnames: The names of the rows and columns
remove.background = function(geneCellMatrix, ambientRNAprofile){
  # do some input checks first
  if(!("dgCMatrix" %in% class(geneCellMatrix))){
    cat("geneCellMatrix should be a sparse matrix of the \"dgCMatrix\" class")
    stop()
  }


  # Here is the actual functionality
  for(gene in names(ambientRNAprofile[ambientRNAprofile > 0])){
    # Determine the locations where the gene is not zero, therefore referenced in i
    iLocs = which(geneCellMatrix@Dimnames[[1]] == gene)
    # Determine the location of the actual values,
    xLocs = which(geneCellMatrix@i == iLocs-1) # -1 because of 0 and 1 based counting systems
    # Remove the contaminating RNA
    geneCellMatrix@x[xLocs] = geneCellMatrix@x[xLocs] - ambientRNAprofile[gene]
  }
  # correct for instances where the expression was corrected to below zero
  geneCellMatrix@x[geneCellMatrix@x < 0] = 0
  # remove the zeroes and return the corrected matrix
  return(drop0(geneCellMatrix, is.Csparse = TRUE))
}
##############

determine.background.to.remove = function(fullCellMatrix, emptyDropletCutoff, contaminationChanceCutoff){

  # determines the highest expression value found for every gene in the droplets that we're sure don't contain cells
  backGroundMax   = as.vector(qlcMatrix::rowMax(fullCellMatrix[,Matrix::colSums(fullCellMatrix) < emptyDropletCutoff]))
  names(backGroundMax) = rownames(fullCellMatrix)

  # droplets that are empty but not unused barcodes, unused barcodes have zero reads assigned to them.
  nEmpty = table((Matrix::colSums(fullCellMatrix) < emptyDropletCutoff) &(Matrix::colSums(fullCellMatrix) > 0))[2]
  # rowSum on a logical statement returns the number of TRUE occurences
  occurences = Matrix::rowSums(fullCellMatrix[,Matrix::colSums(fullCellMatrix) < emptyDropletCutoff] !=0)

  #probability of a background read of a gene ending up in a cell
  probabiltyCellContaminationPerGene = occurences / nEmpty

  # if the probablity of a gene contaminating a cell is too low we set the value to zero so it doesn't get corrected
  backGroundMax[probabiltyCellContaminationPerGene < contaminationChanceCutoff] = 0
  return(backGroundMax)
}

##############

read.cell.matrix = function(cellFolderLocation){
  cellMatrix = Read10X(cellFolderLocation)
  return(cellMatrix)
}

##############

read.full.matrix = function(fullFolderLocation){
  fullMatrix = Read10X(fullFolderLocation)
  return(fullMatrix)
}

###############################################################################
getExpressionThreshold = function(gene, expMat, percentile){
  orderedExpression = expMat[gene, order(expMat[gene,], decreasing = TRUE)]
  CS = cumsum(orderedExpression)
  return(orderedExpression[which(CS/max(CS) > percentile)[1]])
}

###############################################################################
celltypeSpecificityScore = function(gene, expMat){
  CS = cumsum(expMat[gene, order(expMat[gene,], decreasing = TRUE)])
  return((sum(CS/max(CS))/ncol(expMat))  )
}

###############################################################################
describe.correction.effect = function (allExpression, cellExpression, startPos, stopPos, byLength, contaminationChanceCutoff){

  # Make somewhere to store all the data that needs to be returned to the user
  ambientScoreProfileOverview = data.frame(row.names = rownames(cellExpression))

  # do a quick first run to see which genes get corrected at the highest setting
  ambientProfile = determine.background.to.remove(allExpression, stopPos, contaminationChanceCutoff)
  genelist = names(ambientProfile[ambientProfile > 0])

  print(paste0("Calculating cell expression score for ", length(genelist), " genes"))
  ctsScores = vector(mode = "numeric", length = nrow(cellExpression))
  names(ctsScores) = rownames(cellExpression)

  for(gene in genelist){
    ctsScores[gene] = celltypeSpecificityScore(gene, cellExpression)
  }

  # loop over every threshold to test
  # Starts at the highest value so
  for(emptyDropletCutoff in seq(from = startPos, to = stopPos, by = byLength)){
    ambientProfile = determine.background.to.remove(allExpression, emptyDropletCutoff, contaminationChanceCutoff)

    print(paste0("Profiling at cutoff ", emptyDropletCutoff))

    ambientScoreProfile = data.frame(ambientProfile, ctsScores)
    #ambientScoreProfile = ambientScoreProfile[ambientScoreProfile$ctsScores > 0.85, ]
    ambientScoreProfile$stillOverAmbient = 0
    ambientScoreProfile$belowCellexpression = 0

    genesToCheck = names(ambientProfile[ambientProfile > 0])
    if(exists("overAmbientGenes")){
      genesToCheck = overAmbientGenes
    }

    print(paste0("Calculating profiles for ", length(genesToCheck), " genes"))

    for(gene in genesToCheck){
      expThresh = getExpressionThreshold(gene, cellExpression, 0.95)

      if(emptyDropletCutoff == startPos){
        ambientScoreProfile[gene, "belowCellexpression"] = table(cellExpression[gene,] > 0  & cellExpression[gene,] < expThresh)["TRUE"]
      }
      ambientScoreProfile[gene, "stillOverAmbient"] =  table(cellExpression[gene,] > ambientScoreProfile[gene, "ambientProfile"]  & cellExpression[gene,] < expThresh)["TRUE"]
    }

    ambientScoreProfile[is.na(ambientScoreProfile)] = 0
    ambientScoreProfile$contaminationChance = ambientScoreProfile$stillOverAmbient / ambientScoreProfile$belowCellexpression
    ambientScoreProfile[is.na(ambientScoreProfile)] = 0

    # Genes that have already been completely removed don't need to be checked at higher resolution
    overAmbientGenes = rownames(ambientScoreProfile[ambientScoreProfile$stillOverAmbient > 0,])

    ambientScoreProfile[genelist,"AmbientCorrection"]

    ambientScoreProfileOverview[names(ctsScores), "ctsScores"] = ctsScores
    if(emptyDropletCutoff == startPos){
      ambientScoreProfileOverview[rownames(ambientScoreProfile), "belowCellexpression"] = ambientScoreProfile$belowCellexpression
    }
    ambientScoreProfileOverview[rownames(ambientScoreProfile), paste0("stillOverAmbient", as.character(emptyDropletCutoff))] = ambientScoreProfile$stillOverAmbient
    ambientScoreProfileOverview[rownames(ambientScoreProfile), paste0("AmbientCorrection", as.character(emptyDropletCutoff))] = ambientScoreProfile$ambientProfile
  }
  ambientScoreProfileOverview = ambientScoreProfileOverview[!is.na(ambientScoreProfileOverview$ctsScores),]
  ambientScoreProfileOverview[is.na(ambientScoreProfileOverview)] = 0

  ambientScoreProfileOverview[,paste0("Threshold_", seq(from = startPos, to = stopPos, by = byLength))] = ambientScoreProfileOverview[,paste0("stillOverAmbient", as.character(seq(from = startPos, to = stopPos, by = byLength)))] / ambientScoreProfileOverview$belowCellexpression
  ambientScoreProfileOverview[is.na(ambientScoreProfileOverview)] = 0

  ambientScoreProfileOverview[,paste0("contaminationChance", as.character(seq(from = startPos, to = stopPos, by = byLength)))] = ambientScoreProfileOverview[,paste0("stillOverAmbient", as.character(seq(from = startPos, to = stopPos, by = byLength)))] / ambientScoreProfileOverview$belowCellexpression

  return(ambientScoreProfileOverview)
}


##############
# Turns out that cellranger output looks different from WriteMM output and Read10X can't read the latter
# TODO
# Commented out until I find a fix

# write.corrected.matrix = function(correctedMatrix, folderLocation, correctedSignal){
#   dir.create(folderLocation)
#
#   writeMM(obj = correctedMatrix, file=paste(folderLocation, "/matrix.mtx", sep = ""))
#
#   # save genes and cells names
#   write(x = rownames(correctedMatrix), file = paste(folderLocation, "/genes.tsv", sep = ""))
#   write(x = colnames(correctedMatrix), file = paste(folderLocation, "/barcodes.tsv", sep = ""))
#
#   correctedSignal = correctedSignal[correctedSignal > 0]
#   write.table(list(correctedSignal), file = paste(folderLocation, "/genesCorrectedFor.csv", sep = ""), row.names = TRUE, col.names = FALSE)
#
# }

# describe the number of genes identified in the background
# and the number of genes failing the contamination chance threshold
#
describe.ambient.RNA.sequence = function(fullCellMatrix, start, stop, by, contaminationChanceCutoff){
  cutoffValue = seq(start, stop, by)
  genesInBackground  = vector(mode = "numeric", length = length(seq(start, stop, by)))
  genesContaminating = vector(mode = "numeric", length = length(seq(start, stop, by)))
  nEmptyDroplets     = vector(mode = "numeric", length = length(seq(start, stop, by)))

  ambientDescriptions = data.frame(nEmptyDroplets, genesInBackground, genesContaminating, cutoffValue)
  rownames(ambientDescriptions) = seq(start, stop, by)
  for(emptyCutoff in seq(start, stop, by)){
    nEmpty = table((Matrix::colSums(fullCellMatrix) < emptyCutoff) &(Matrix::colSums(fullCellMatrix) > 0))[2]

    occurences = Matrix::rowSums(fullCellMatrix[,Matrix::colSums(fullCellMatrix) < emptyCutoff] !=0)

    #probability of a background read of a gene ending up in a cell
    probabiltyCellContaminationPerGene = occurences / nEmpty
    nFailingThreshold = sum(probabiltyCellContaminationPerGene > contaminationChanceCutoff)

    nGenes = sum(occurences != 0)
    ambientDescriptions[as.character(emptyCutoff), c(1,2,3)] = c(nEmpty ,nGenes, nFailingThreshold)
  }
  return(ambientDescriptions)
}


plot.correction.effect.chance = function(correctionProfile){
  pheatmap(correctionProfile[correctionProfile[,2] > 0, colnames(correctionProfile)[grep("contaminationChance", colnames(correctionProfile))]],
           cluster_cols = FALSE,
           treeheight_row = 0,
           main = "Chance of affecting DE analyses")
}

plot.correction.effect.removal = function(correctionProfile){
  pheatmap(correctionProfile[(correctionProfile[,2] > 0) ,colnames(correctionProfile)[grep("AmbientCorrection", colnames(correctionProfile))]],
           cluster_cols = FALSE, treeheight_row = 0,
           main = "Counts removed from each cell")
}


plot.ambient.profile = function(ambientProfile){

  p1 = ggplot(ambientProfile, aes(x=cutoffValue, y=genesInBackground)) + geom_point()


  p2= ggplot(ambientProfile, aes(x=cutoffValue, y=genesContaminating)) + geom_point()

  p3 = ggplot(ambientProfile, aes(x=cutoffValue, y=nEmptyDroplets)) + geom_point()

  grid.arrange(p1, p2, p3, nrow = 3)

}

# I noticed that the number of genes removed tends to even out over time
# Test whether the point where this first happens is a good empty cutoff point
recommend.empty.cutoff = function(ambientProfile){
  highestNumberOfGenes = max(ambientProfile[,3])
  firstOccurence = match(highestNumberOfGenes, ambientProfile[,3])
  return(as.numeric(rownames(ambientProfile[firstOccurence,])))
}










