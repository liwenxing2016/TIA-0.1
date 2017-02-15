source("I:/TEST/Protocol/AnalysisProtocol.R");

# Required Parameters
codePath <- "I:/TEST/Protocol"
dataPath <- "I:/TEST/Protocol/Results/5. Normalized Dataset/Normalized Dataset.txt";
samplePhenotype <- "I:/TEST/Protocol/Sample Phenotype.txt";
outPath <- "I:/TEST/Protocol";

platformPath <- "I:/Platforms";
datasetsInfo <- data.frame(
	Study=c("GSE36076", "GSE51725"),
	Platform=c("GPL570-13270.txt", "GPL570-13270.txt")
	);
compare <- "Tissue";

integrityTargetPath <- "I:/Thomson Reuters Integrity/Target and Pathways";
diseaseName <- "Gastrointestinal Stromal Tumor (GIST)";

# Note: function pathview is unstable
runAnalysis(analysisType=3, codePath=codePath, dataPath=dataPath, samplePhenotype=samplePhenotype, outPath=outPath, 
	platformPath=platformPath, datasetsInfo=datasetsInfo, compare=compare, runDiffExpr=TRUE, runGO=TRUE, runKEGG=TRUE, runTF=TRUE, 
	runTarget=TRUE, integrityTargetPath=integrityTargetPath, diseaseName=diseaseName);
