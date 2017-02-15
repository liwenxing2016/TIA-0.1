#==============================================================================#
#  An Integrate Analysis Protocol for Multi-platform Gene Expression Datasets  #
#    --from Raw Data to Published Results                                      #
#                                                    Designed by: Wen-Xing Li  #
#                                                    Last Update: 2016-07-12   #
#==============================================================================#

library(oligo);
library(limma);
library(RankProd);
library(pathview);
library(pheatmap);
library(GO.db);
library(KEGG.db);
library(GOstats);
library(org.Hs.eg.db);
library(org.Mm.eg.db);
library(org.Rn.eg.db);

runAnalysis <- function(
	analysisType=1,
	studyID=NULL,
	codePath=NULL,
	dataPath=NULL,
	samplePhenotype=NULL,
	outPath=NULL,
	platformPath=NULL,
	datasetsInfo=NULL,
	isAnnotated=FALSE,
	isMerged=FALSE,
	isIntegrate=FALSE,
	isRenormalized=FALSE,
	mergeMethod="mean",
	repeatSymbol=1,
	filterData=FALSE,
	runTotal=TRUE,
	compare=NA,
	compareType=1,
	runDiffExpr=TRUE,
	DEGs_method="eBayes",
	FC_cutoff=2.0,
	Pval_cutoff=0.05,
	adjPval_cutoff=NULL,
	runGO=TRUE,
	runGOBP=TRUE,
	runGOCC=TRUE,
	runGOMF=TRUE,
	runKEGG=TRUE,
	runTF=TRUE,
	runTarget=TRUE,
	integrityTargetPath=NA,
	diseaseName=NA
	) {
	pData <- read.table(file=samplePhenotype, header=TRUE, stringsAsFactors=FALSE, sep="\t");
	rownames(pData) <- pData$SampleID;
	reportPath <- paste(outPath, "Report", sep="/");
	if (!file.exists(reportPath)) {
		dir.create(reportPath);
	}
	resultsPath <- paste(outPath, "Results", sep="/");
	if (!file.exists(resultsPath)) {
		dir.create(resultsPath);
	}
	if (analysisType == 1) {
		celFile_to_expr_matrix(dataPath, studyID, resultsPath);
		integrate_multi_datasets(analysisType, datasetsInfo, platformPath, samplePhenotype, resultsPath, 
			isAnnotated, isMerged, isIntegrate, isRenormalized, mergeMethod, repeatSymbol);
		if (isIntegrate == FALSE) {
			cat("Drawing integrated datasets value distribution\n");
			filename <- paste(outPath, "Results/4. Integrated Dataset/Integrated Dataset.txt", sep="/");
			eset <- read.table(file=filename, header=TRUE, sep="\t");
			filename <- paste(outPath, "Results/4. Integrated Dataset/Integrated Dataset.png", sep="/");
			valueDistribution(eset, pData, groupName="Study", main="Integrated Datasets Value Distribution", filename=filename);
			rm(eset);
		}
		if (isRenormalized == FALSE) {
			cat("Drawing re-normalized datasets value distribution\n");
			filename <- paste(outPath, "Results/5. Normalized Dataset/Normalized Dataset.txt", sep="/");
			eset <- read.table(file=filename, header=TRUE, sep="\t");
			filename <- paste(outPath, "Results/5. Normalized Dataset/Normalized Dataset.png", sep="/");
			valueDistribution(eset, pData, groupName="Study", main="Re-normalized Datasets Value Distribution", filename=filename);
			rm(eset);
		}
		filename <- paste(outPath, "Results/5. Normalized Dataset/Normalized Dataset.txt", sep="/");
		eset <- read.table(file=filename, header=TRUE, sep="\t");
	}
	if (analysisType == 2) {
		integrate_multi_datasets(analysisType, datasetsInfo, platformPath, samplePhenotype, resultsPath, 
			isAnnotated, isMerged, isIntegrate, isRenormalized, mergeMethod, repeatSymbol);
		if (isIntegrate == FALSE) {
			cat("Drawing integrated datasets value distribution\n");
			filename <- paste(outPath, "Results/4. Integrated Dataset/Integrated Dataset.txt", sep="/");
			eset <- read.table(file=filename, header=TRUE, sep="\t");
			filename <- paste(outPath, "Results/4. Integrated Dataset/Integrated Dataset.png", sep="/");
			valueDistribution(eset, pData, groupName="Study", main="Integrated Datasets Value Distribution", filename=filename);
			rm(eset);
		}
		if (isRenormalized == FALSE) {
			cat("Drawing re-normalized datasets value distribution\n");
			filename <- paste(outPath, "Results/5. Normalized Dataset/Normalized Dataset.txt", sep="/");
			eset <- read.table(file=filename, header=TRUE, sep="\t");
			filename <- paste(outPath, "Results/5. Normalized Dataset/Normalized Dataset.png", sep="/");
			valueDistribution(eset, pData, groupName="Study", main="Re-normalized Datasets Value Distribution", filename=filename);
			rm(eset);
		}
		filename <- paste(outPath, "Results/5. Normalized Dataset/Normalized Dataset.txt", sep="/");
		eset <- read.table(file=filename, header=TRUE, sep="\t");
	}
	if (analysisType == 3) {
		eset <- read.table(file=dataPath, header=TRUE, sep="\t");
	}
	if (filterData == TRUE) {
		pData <- extract_datasets(pData, study, species, sex, age, tissue, tissueType, stage);
		eset <- eset[,pData$SampleID];
	}
	diffExprPath <- paste(outPath, "Results/6. Differential Expression Analysis", sep="/");
	if (!file.exists(diffExprPath)) {
		dir.create(diffExprPath);
	}
	GOPath <- paste(outPath, "Results/7. DEGs GO Enrichment", sep="/");
	if (!file.exists(GOPath)) {
		dir.create(GOPath);
	}
	KEGGPath <- paste(outPath, "Results/8. DEGs KEGG Enrichment", sep="/");
	if (!file.exists(KEGGPath)) {
		dir.create(KEGGPath);
	}
	TFPath <- paste(outPath, "Results/9. Transcription factor Analysis", sep="/");
	if (!file.exists(TFPath)) {
		dir.create(TFPath);
	}
	targetPath <- paste(outPath, "Results/10. Disease Target Analysis", sep="/");
	if (!file.exists(targetPath)) {
		dir.create(targetPath);
	}
	if (runTotal == TRUE) {
		groupName <- "Total";
		analysis_flow(eset, pData, runDiffExpr, DEGs_method, FC_cutoff, Pval_cutoff, adjPval_cutoff, 
			runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, integrityTargetPath,
			diseaseName, groupName, diffExprPath, GOPath, KEGGPath, TFPath, targetPath, codePath);
	}
	if (!is.na(compare)) {
		compare_id <- 0;
		for (i in c(1:length(colnames(pData)))) {
			if (colnames(pData)[i] == compare & compare != "Category") {
				compare_id <- i;
				break;
			}
		}
		if (compare_id == 0) {
			return;
		}
		if (compareType == 1) {
			groupNames <- as.character(unique(pData[,compare_id]));
			groupNames <- sort(groupNames[!is.na(groupNames)]);
			for (i in c(1:length(groupNames))) {
				select_pData <- pData[!is.na(pData[,compare_id]) & pData[,compare_id]==groupNames[i],];
				select_eset <- eset[,select_pData$SampleID];
				analysis_flow(select_eset, select_pData, runDiffExpr, DEGs_method, FC_cutoff, Pval_cutoff, adjPval_cutoff, 
					runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, integrityTargetPath, diseaseName, 
					groupNames[i], diffExprPath, GOPath, KEGGPath, TFPath, targetPath, codePath);
			}
		}
		if (compareType == 2) {
			groupNames <- setdiff(as.character(unique(pData[,compare_id])), "Control");
			groupNames <- sort(groupNames[!is.na(groupNames)]);
			for (i in c(1:length(groupNames))) {
				select_pData <- pData[!is.na(pData[,compare_id]) & (pData[,compare_id]==groupNames[i] | pData[,compare_id]=="Control"),];
				select_eset <- eset[,select_pData$SampleID];
				analysis_flow(select_eset, select_pData, runDiffExpr, DEGs_method, FC_cutoff, Pval_cutoff, adjPval_cutoff, 
					runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, integrityTargetPath, diseaseName, 
					groupNames[i], diffExprPath, GOPath, KEGGPath, TFPath, targetPath, codePath);
			}
		}
		if (compareType == 3) {
			studies <- as.character(unique(pData[,"Study"]));
			for (i in c(1:length(studies))) {
				pData_1 <- pData[pData$Study==studies[i],];
				groupNames <- setdiff(as.character(unique(pData_1[,compare_id])), "Control");
				groupNames <- sort(groupNames[!is.na(groupNames)]);
				for (j in c(1:length(groupNames))) {
					pData_2 <- pData_1[!is.na(pData_1[,compare_id]) & (pData_1[,compare_id]==groupNames[j] | pData_1[,compare_id]=="Control"),];
					select_eset <- eset[,pData_2$SampleID];
					analysis_flow(select_eset, pData_2, runDiffExpr, DEGs_method, FC_cutoff, Pval_cutoff, adjPval_cutoff, 
						runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, integrityTargetPath, diseaseName, 
						paste(studies[i], groupNames[j]), diffExprPath, GOPath, KEGGPath, TFPath, targetPath, codePath);
				}
			}
		}
		group_compare(eset, pData, runDiffExpr, runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, 
			compare, compareType, diffExprPath, GOPath, KEGGPath, TFPath, targetPath);
	}
}

analysis_flow <- function(eset, pData, runDiffExpr, DEGs_method, FC_cutoff, Pval_cutoff, adjPval_cutoff, 
	runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, integrityTargetPath, diseaseName, 
	groupName, diffExprPath, GOPath, KEGGPath, TFPath, targetPath, codePath) {
	cat(paste("Run analysis in \"", groupName, "\" group\n", sep=""));
	if (runDiffExpr == TRUE) {
		cat("Calculate differentially expressed genes\n");
		if (DEGs_method == "eBayes") {
			dif <- eBayes_DiffGene(eset, pData);
		}
		if (DEGs_method == "RankProd") {
			dif <- RankProd_DiffGene(eset, pData);
		}
		cat("Save differentially expressed genes as txt format\n");
		filename <- paste(diffExprPath, "/all_genes_list (", groupName, ").txt", sep="");
		write.table(dif, file=filename, quote=FALSE, sep="\t");
		if (!is.null(adjPval_cutoff)) {
			dif2 <- dif[(abs(dif$logFC) >= log2(FC_cutoff)) & (dif$adj.P.Val <= adjPval_cutoff),];
		}
		else {
			dif2 <- dif[(abs(dif$logFC) >= log2(FC_cutoff)) & (dif$P.Value <= Pval_cutoff),];
		}
		dif2_up <- dif2[dif2$logFC > 0,];
		dif2_down <- dif2[dif2$logFC < 0,];
		filename <- paste(diffExprPath, "/diff_genes_list (", groupName, ").txt", sep="");
		write.table(dif2, file=filename, quote=FALSE, sep="\t");
		symbol_list <- as.data.frame(rownames(dif2));
		symbol_list_up <- as.data.frame(rownames(dif2_up));
		symbol_list_down <- as.data.frame(rownames(dif2_down));
		filename <- paste(diffExprPath, "/diff_symbol (", groupName, ").txt", sep="");
		write.table(symbol_list, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");
		filename <- paste(diffExprPath, "/diff_symbol_up (", groupName, ").txt", sep="");
		write.table(symbol_list_up, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");
		filename <- paste(diffExprPath, "/diff_symbol_down (", groupName, ").txt", sep="");
		write.table(symbol_list_down, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");
		cat("Drawing log(FC) scatterplot\n");
		filename <- paste(diffExprPath, "/log(FC) scatterplot (", groupName, ").png", sep="");
		plot_logFC(dif, FC_cutoff=FC_cutoff, savefile="png", filename=filename);
		FC_list <- dif$logFC;
		if (!is.null(adjPval_cutoff)) {
			P_list <- -log10(dif$adj.P.Val);
			P_cutoff <- adjPval_cutoff;
			ylab <- "-log10(adj.P.Val)";
		}
		else {
			P_list <- -log10(dif$P.Value);
			P_cutoff <- Pval_cutoff;
			ylab <- "-log10(P.Value)";
		}
		P_list <- P_list[!is.infinite(P_list)];
		cat("Drawing volcano plot\n");
		filename <- paste(diffExprPath, "/Volcano plot (", groupName, ").png", sep="");
		plot_volcano(FC_list, P_list, FC_cutoff=FC_cutoff, P_cutoff=P_cutoff, ylab=ylab, savefile="png", filename=filename);
		dif_count <- length(rownames(dif2));
		if (dif_count > 0) {
			dif_list <- cbind(rownames(dif2), dif2);
			colnames(dif_list)[1] <- "Symbol";
			rownames(dif_list) <- c(1:length(rownames(dif_list)));
			dif_list_up <- dif_list[dif_list$logFC>0,];
			dif_list_down <- dif_list[dif_list$logFC<0,];
			up_gene_count <- length(rownames(dif_list_up));
			down_gene_count <- length(rownames(dif_list_down));
			cat("Drawing differentially expressed genes count\n");
			filename <- paste(diffExprPath, "/DiffGene count (", groupName, ").png", sep="");
			plot_dif_count(up_gene_count, down_gene_count, savefile="png", filename=filename);
			cat("Save differentially expressed genes as html format\n");
			filename <- paste(diffExprPath, "/Differentially Expressed Genes (", groupName, ").html", sep="");
			save_as_html(dif_list, dif_list_up, dif_list_down, "DiffGene", filename);
		}
		else {
			return;
		}
	}
	if (runGO == TRUE & runGOBP == TRUE & dif_count > 0) {
		cat("Starting GO Biological Process analysis\n");
		filename <- paste(codePath, "/GOBPID_level4.txt", sep="");
		GOBPID_level4 <- read.table(file=filename, header=TRUE, sep="\t");
		GOBP_all <- enrichment_analysis(rownames(dif2), "GOBP", 4);
		GOBP_up <- enrichment_analysis(rownames(dif2_up), "GOBP", 4);
		GOBP_down <- enrichment_analysis(rownames(dif2_down), "GOBP", 4);
		if ((!is.null(GOBP_all))) {
			cat("Save GO Biological Process as txt format\n");
			filename <- paste(GOPath, "/GO Biological Process (", groupName, ").txt", sep="");
			write.table(GOBP_all, file=filename, row.names=FALSE, quote=FALSE, sep="\t");
		}
		cat("Save GO Biological Process as html format\n");
		filename <- paste(GOPath, "/GO Biological Process (", groupName, ").html", sep="");
		save_as_html(GOBP_all, GOBP_up, GOBP_down, "GOBP", filename);
		if ((!is.null(GOBP_up)) & (!is.null(GOBP_down))) {
			cat("Drawing enriched GO Biological Process\n");
			filename <- paste(GOPath, "/GO Biological Process (", groupName, ").png", sep="");
			plot_enrichment(up_gene_count, down_gene_count, GOBP_up, GOBP_down, filename=filename);
		}
	}
	if (runGO == TRUE & runGOCC == TRUE & dif_count > 0) {
		cat("Starting GO Cellular Component analysis\n");
		GOCC_all <- enrichment_analysis(rownames(dif2), "GOCC");
		GOCC_up <- enrichment_analysis(rownames(dif2_up), "GOCC");
		GOCC_down <- enrichment_analysis(rownames(dif2_down), "GOCC");
		if ((!is.null(GOCC_all))) {
			cat("Save GO Cellular Component as txt format\n");
			filename <- paste(GOPath, "/GO Cellular Component (", groupName, ").txt", sep="");
			write.table(GOCC_all, file=filename, row.names=FALSE, quote=FALSE, sep="\t");
		}
		cat("Save GO Cellular Component as html format\n");
		filename <- paste(GOPath, "/GO Cellular Component (", groupName, ").html", sep="");
		save_as_html(GOCC_all, GOCC_up, GOCC_down, "GOCC", filename);
		if ((!is.null(GOCC_up)) & (!is.null(GOCC_down))) {
			cat("Drawing enriched GO Cellular Component\n");
			filename <- paste(GOPath, "/GO Cellular Component (", groupName, ").png", sep="");
			plot_enrichment(up_gene_count, down_gene_count, GOCC_up, GOCC_down, filename=filename);
		}
	}
	if (runGO == TRUE & runGOMF == TRUE & dif_count > 0) {
		cat("Starting GO Molecular Function analysis\n");
		GOMF_all <- enrichment_analysis(rownames(dif2), "GOMF");
		GOMF_up <- enrichment_analysis(rownames(dif2_up), "GOMF");
		GOMF_down <- enrichment_analysis(rownames(dif2_down), "GOMF");
		if ((!is.null(GOMF_all))) {
			cat("Save GO Molecular Function as txt format\n");
			filename <- paste(GOPath, "/GO Molecular Function (", groupName, ").txt", sep="");
			write.table(GOMF_all, file=filename, row.names=FALSE, quote=FALSE, sep="\t");
		}
		cat("Save GO Molecular Function as html format\n");
		filename <- paste(GOPath, "/GO Molecular Function (", groupName, ").html", sep="");
		save_as_html(GOMF_all, GOMF_up, GOMF_down, "GOMF", filename);
		if ((!is.null(GOMF_up)) & (!is.null(GOMF_down))) {
			cat("Drawing enriched GO Molecular Function\n");
			filename <- paste(GOPath, "/GO Molecular Function (", groupName, ").png", sep="");
			plot_enrichment(up_gene_count, down_gene_count, GOMF_up, GOMF_down, filename=filename);
		}
	}
	if (runKEGG == TRUE & dif_count > 0) {
		cat("Starting KEGG Pathway analysis\n");
		KEGG_all <- enrichment_analysis(rownames(dif2), "KEGG");
		KEGG_up <- enrichment_analysis(rownames(dif2_up), "KEGG");
		KEGG_down <- enrichment_analysis(rownames(dif2_down), "KEGG");
		if ((!is.null(KEGG_all))) {
			cat("Save KEGG Pathway as txt format\n");
			filename <- paste(KEGGPath, "/KEGG Pathway (", groupName, ").txt", sep="");
			write.table(KEGG_all, file=filename, row.names=FALSE, quote=FALSE, sep="\t");
		}
		cat("Save KEGG Pathway as html format\n");
		filename <- paste(KEGGPath, "/KEGG Pathway (", groupName, ").html", sep="");
		save_as_html(KEGG_all, KEGG_up, KEGG_down, "KEGG", filename);
		if ((!is.null(KEGG_up)) & (!is.null(KEGG_down))) {
			cat("Drawing enriched KEGG Pathway\n");
			filename <- paste(KEGGPath, "/KEGG Pathway (", groupName, ").png", sep="");
			plot_enrichment(up_gene_count, down_gene_count, KEGG_up, KEGG_down, filename=filename);
		}
		if (length(KEGG_all$KEGGID) > 20) {
			kegg_ids <- KEGG_all$KEGGID[1:20];
		}
		else {
			kegg_ids <- KEGG_all$KEGGID;
		}
		if (length(KEGG_all$KEGGID) > 0) {
			cat("Drawing gene expression profiles in KEGG Pathway\n");
			filename <- paste(KEGGPath, "/KEGG Pathway Gene Expression (", groupName, ").png", sep="");
			plot_ave_expr_val(eset, pData, dif, kegg_ids, "KEGG", gap=0, FC_cutoff, Pval_cutoff, adjPval_cutoff, groupName, 
				map_list=NULL, inner_color=NULL, outer_color=NULL, label=c("Case", "Control"), filename=filename);
		}
		if (length(KEGG_all$KEGGID) > 5) {
			kegg_ids <- KEGG_all$KEGGID[1:5];
		}
		else {
			kegg_ids <- KEGG_all$KEGGID;
		}
		kegg_dir <- paste(codePath, "KEGG_all_pathways", sep="/");
#		cat("Drawing pathview in top KEGG Pathway\n");
#		plot_pathview(as.matrix(dif)[,1], kegg_ids, species="hsa", gene.idtype="SYMBOL", groupName, kegg_dir, KEGGPath);
	}
	if (runTF == TRUE & dif_count > 0) {
		cat("Starting transcription factor analysis\n");
		filename <- paste(codePath, "/Homo_sapiens_transcription_factor_list.txt", sep="");
		TF_list <- read.table(file=filename, header=TRUE, sep="\t");
		rownames(TF_list) <- TF_list$Symbol;
		filename <- paste(codePath, "/Homo_sapiens_transcription_factor_target.txt", sep="");
		TF_target <- read.table(file=filename, header=TRUE, sep="\t");
		dif_TF <- TFs_analysis(dif, dif2, TF_list, TF_target, top_num=100, TFPath, groupName);
		if (length(rownames(dif_TF)) > 0) {
			dif_TF_list <- cbind(rownames(dif_TF), dif_TF);
			colnames(dif_TF_list)[1] <- "Symbol";
			rownames(dif_TF_list) <- c(1:length(rownames(dif_TF_list)));
			dif_TF_list_up <- dif_TF_list[dif_TF_list$logFC>0,];
			dif_TF_list_down <- dif_TF_list[dif_TF_list$logFC<0,];
			cat("Save transcription factor as html format\n");
			filename <- paste(TFPath, "/Transcription Factor (", groupName, ").html", sep="");
			save_as_html(dif_TF_list, dif_TF_list_up, dif_TF_list_down, "TF", filename);
		}
	}
	if (runTarget == TRUE & dif_count > 0) {
		cat("Starting disease target analysis\n");
		target_analysis(eset, pData, dif, dif2, targetPath, groupName, integrityTargetPath, diseaseName);
	}
}

group_compare <- function(eset, pData, runDiffExpr, runGO, runGOBP, runGOCC, runGOMF, runKEGG, runTF, runTarget, 
	compare, compareType, diffExprPath, GOPath, KEGGPath, TFPath, targetPath) {
	groupDiff <- as.list(c());
	groupGOBP <- as.list(c());
	groupGOCC <- as.list(c());
	groupGOMF <- as.list(c());
	groupKEGG <- as.list(c());
	groupTF <- as.list(c());
	groupTarget <- as.list(c());
	groupNames <- NULL;
	if (!is.na(compare)) {
		compare_id <- 0;
		for (i in c(1:length(colnames(pData)))) {
			if (colnames(pData)[i] == compare & compare != "Category") {
				compare_id <- i;
				break;
			}
		}
		if (compare_id == 0) {
			return;
		}
		if (compareType == 1) {
			groupNames <- as.character(unique(pData[,compare_id]));
		}
		if (compareType == 2) {
			groupNames <- setdiff(as.character(unique(pData[,compare_id])), "Control");
		}
		if (compareType == 3) {
			pData_1 <- pData[!is.na(pData[,compare_id]) & pData$Category!="Control",];
			groupNames <- unique(paste(as.character(pData_1$Study), as.character(pData_1[,compare_id])));
		}
	}
	if (is.null(groupNames)) {
		return;
	}
	if (runDiffExpr == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(diffExprPath, "/diff_genes_list (", groupNames[i], ").txt", sep="");
			dif2 <- read.table(file=filename, header=TRUE, sep="\t");
			groupDiff[[i]] <- rownames(dif2);
		}
		names(groupDiff) <- groupNames;
		cat("Drawing differentially expressed genes Venn diagram\n");
		filename <- paste(diffExprPath, "/DiffGene Venn Diagram.png", sep="");
		smart_venn(groupDiff, cex.num=2, cex.lab=2, filename=filename);
	}
	if (runGO == TRUE & runGOBP == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(GOPath, "/GO Biological Process (", groupNames[i], ").txt", sep="");
			GOBP_all <- read.table(file=filename, header=TRUE, sep="\t");
			groupGOBP[[i]] <- GOBP_all$GOBPID;
		}
		names(groupGOBP) <- groupNames;
		cat("Drawing GO Biological Process Venn diagram\n");
		filename <- paste(GOPath, "/GOBP Venn Diagram.png", sep="");
		smart_venn(groupGOBP, cex.num=2, cex.lab=2, filename=filename);
	}
	if (runGO == TRUE & runGOCC == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(GOPath, "/GO Cellular Component (", groupNames[i], ").txt", sep="");
			GOCC_all <- read.table(file=filename, header=TRUE, sep="\t");
			groupGOCC[[i]] <- GOCC_all$GOCCID;
		}
		names(groupGOCC) <- groupNames;
		cat("Drawing GO Cellular Component Venn diagram\n");
		filename <- paste(GOPath, "/GOCC Venn Diagram.png", sep="");
		smart_venn(groupGOCC, cex.num=2, cex.lab=2, filename=filename);
	}
	if (runGO == TRUE & runGOMF == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(GOPath, "/GO Molecular Function (", groupNames[i], ").txt", sep="");
			GOMF_all <- read.table(file=filename, header=TRUE, sep="\t");
			groupGOMF[[i]] <- GOMF_all$GOMFID;
		}
		names(groupGOMF) <- groupNames;
		cat("Drawing GO Molecular Function Venn diagram\n");
		filename <- paste(GOPath, "/GOMF Venn Diagram.png", sep="");
		smart_venn(groupGOMF, cex.num=2, cex.lab=2, filename=filename);
	}
	if (runKEGG == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(KEGGPath, "/KEGG Pathway (", groupNames[i], ").txt", sep="");
			KEGG_all <- read.table(file=filename, header=TRUE, sep="\t");
			groupKEGG[[i]] <- KEGG_all$KEGGID;
		}
		names(groupKEGG) <- groupNames;
		cat("Drawing KEGG Pathway Venn diagram\n");
		filename <- paste(KEGGPath, "/KEGG Venn Diagram.png", sep="");
		smart_venn(groupKEGG, cex.num=2, cex.lab=2, filename=filename);
	}
	if (runTF == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(TFPath, "/dif_transcription_factor (", groupNames[i], ").txt", sep="");
			dif_TF <- read.table(file=filename, header=TRUE, sep="\t");
			groupTF[[i]] <- rownames(dif_TF);
		}
		names(groupTF) <- groupNames;
		cat("Drawing differentially transcription factors Venn diagram\n");
		filename <- paste(TFPath, "/TF Venn Diagram.png", sep="");
		smart_venn(groupTF, cex.num=2, cex.lab=2, filename=filename);
	}
	if (runTarget == TRUE) {
		for (i in c(1:length(groupNames))) {
			filename <- paste(targetPath, "/dif_target_genes (", groupNames[i], ").txt", sep="");
			dif_target <- read.table(file=filename, header=TRUE, sep="\t");
			groupTarget[[i]] <- rownames(dif_target);
		}
		names(groupTarget) <- groupNames;
		cat("Drawing differentially target genes Venn diagram\n");
		filename <- paste(targetPath, "/Target Venn Diagram.png", sep="");
		smart_venn(groupTarget, cex.num=2, cex.lab=2, filename=filename);
	}
}

celFile_to_expr_matrix <- function(dataPath, studyID=NULL, resultsPath) {
	if (is.null(studyID)) {
		studyID <- dir(dataPath);
	}
	celPaths <- paste(dataPath, studyID, sep="/");
	exprPath <- paste(resultsPath, "1. Raw Expression Matrix", sep="/");
	if (!file.exists(exprPath)) {
		dir.create(exprPath);
	}
	for (i in c(1:length(studyID))) {
		celFiles <- list.celfiles(celPaths[i], full.name=TRUE);
		cat(paste("Reading CEL files in", studyID[i], "\n"));
		celData <- read.celfiles(celFiles);
		sampleNames(celData) <- gsub(".CEL$", "", sampleNames(celData));
		cat("Robust Multichip Average (RMA) algorithm normalization:\n");
		exprMatrix <- rma(celData);
		exprMatrix <- assayDataElement(exprMatrix, "exprs");
		exprMatrix <- as.data.frame(signif(exprMatrix, 6));
		exprFileName <- paste(exprPath, studyID[i], sep="/");
		write.table(exprMatrix, file=exprFileName, sep="\t", quote=FALSE);
		cat(paste("Gene expression matrix of", studyID[i], "has written to the file.\n"));
		rm(celFiles, celData, exprMatrix);
		temp <- gc();
	}
}

integrate_multi_datasets <- function(analysisType, datasetsInfo, platformPath, samplePhenotype, resultsPath, 
	isAnnotated, isMerged, isIntegrate, isRenormalized, mergeMethod="mean", repeatSymbol=1) {
	if (analysisType == 1) {
		exprMatrixPath <- paste(resultsPath, "1. Raw Expression Matrix", sep="/");
	}
	if (analysisType == 2) {
		exprMatrixPath <- dataPath;
	}
	integrateCode <- paste("\"", paste(codePath, "integrate_multi_datasets.py", sep="/"), "\"", sep="");
	integrateCode <- paste(integrateCode, " \"", exprMatrixPath, "\"", sep="");
	datasetsCode <- paste("[", datasetsInfo$Study, ",", datasetsInfo$Platform, "]", sep="");
	for (i in c(1:length(datasetsCode))) {
		integrateCode <- paste(integrateCode, datasetsCode[i]);
	}
	integrateCode <- paste(integrateCode, " \"", platformPath, "\" \"", samplePhenotype, "\" \"", resultsPath, "\"", sep="");
	integrateCode <- paste(integrateCode, isAnnotated, isMerged, isIntegrate, isRenormalized, mergeMethod, as.character(repeatSymbol));
	system(paste("python", integrateCode));
}

matching_eset_pData <- function(eset, pData) {
	cat("Matching the expression matrix and sample phenotype\n");
	pData <- pData[order(pData$SampleID),];
	if (length(colnames(eset)) < length(rownames(pData))) {
		pData <- pData[rownames(t(eset)),];
	}
	else if (length(colnames(eset)) > length(rownames(pData))) {
		stop("Expression matrix and sample phenotype do not match!");
	}
	else {
		tempData <- pData;
		max_len <- length(colnames(eset));
		isMatch <- FALSE;
		for (i in c(1:max_len)) {
			isMatch <- FALSE;
			a <- 1;
			b <- max_len;
			count <- 0;
			while (a <= b) {
				j <- as.integer((a + b) / 2);
				tempID <- as.character(tempData$SampleID[j]);
				if (colnames(eset)[i] == tempID) {
					pData[i,] <- tempData[j,];
					isMatch <- TRUE;
					break;
				}
				else if (colnames(eset)[i] > tempID) {
					a <- j + 1;
				}
				else if (colnames(eset)[i] < tempID) {
					b <- j - 1;
				}
				count <- count + 1;
				if (count > (log2(max_len) + 1)) {
					isMatch <- FALSE;
					break;
				}
			}
		}
		rm(tempData);
		if (!isMatch) {
			stop("Expression matrix and sample phenotype do not match!");
		}
	}
}

extract_datasets <- function(pData, study=NULL, species=NULL, sex=NULL, age=NULL, tissue=NULL, tissueType=NULL, stage=NULL) {
	cat("Extracting selected dataset\n");
	if (!is.null(study)) {
		pData <- pData[!is.na(pData$Study),];
		pData <- pData[pData$Study == study,];
	}
	if (!is.null(species)) {
		pData <- pData[!is.na(pData$Species),];
		pData <- pData[pData$Species == species,];
	}
	if (!is.null(sex)) {
		pData <- pData[!is.na(pData$Sex),];
		pData <- pData[pData$Sex == sex,];
	}
	if (!is.null(age)) {
		if (length(age) != 2) {
			warning(paste("Age must be an interval (for example: c(0, 60)).\nYour input is: ",
						  as.character(age), ".\nThis information will be ignored.", sep=""));
		}
		else {
			pData <- pData[!is.na(pData$Age),];
			pData <- pData[(pData$Age>=min(age))&(pData$Age<max(age)),];
		}
	}
	if (!is.null(tissue)) {
		pData <- pData[!is.na(pData$Tissue),];
		pData <- pData[pData$Tissue == tissue,];
	}
	if (!is.null(tissueType)) {
		pData <- pData[!is.na(pData$TissueType),];
		pData <- pData[pData$TissueType == tissueType,];
	}
	if (!is.null(stage)) {
		pData <- pData[!is.na(pData$Stage),];
		pData <- pData[(pData$Stage==stage)|(pData$Stage=="Normal"),];
	}
	rownames(pData) <- pData$SampleID;
	return(pData);
}

valueDistribution <- function(eset, pData, groupName="Study", label=c("Case", "Control"), color=c("#D95F02", "#1B9E77"), 
	main=NULL, xlab=NULL, ylab=NULL, imgwd=NULL, filename="Value Distribution.png") {
	mean_case <- c();
	mean_control <- c();
	group <- as.character(unique(pData[,groupName]));
	gene_num <- length(eset[,1]);
	for (i in c(1:length(group))) {
		pData_case <- pData[(pData[,groupName]==group[i])&(pData$Category=="Case"),];
		pData_control <- pData[(pData[,groupName]==group[i])&(pData$Category=="Control"),];
		eset_case <- eset[,pData_case$SampleID];
		eset_control <- eset[,pData_control$SampleID];
		if (length(eset_case[1,]) > 1) {
			mean_case <- cbind(mean_case, apply(as.matrix(eset_case), 1, mean));
		}
		else {
			mean_case <- cbind(mean_case, rep(NA, length(gene_num)));
		}
		colnames(mean_case)[i] <- group[i];
		if (length(eset_control[1,]) > 1) {
			mean_control <- cbind(mean_control, apply(as.matrix(eset_control), 1, mean));
		}
		else {
			mean_control <- cbind(mean_control, rep(NA, length(gene_num)));
		}
		colnames(mean_control)[i] <- group[i];
	}
	leg.txt <- label;
	leg.color <- color;
	if (is.null(main)) {
		main <- "Value Distribution Boxplot";
	}
	if (is.null(xlab)) {
		xlab <- group;
	}
	if (is.null(ylab)) {
		ylab <- "Values";
	}
	n <- length(group);
	if (is.null(imgwd)) {
		imgwd <- (n + 1) * 100;
	}
	par(mar=c(3, 4, 3, 1));
	png(filename=filename, width=imgwd, height=600);
	boxplot(mean_case, boxwex=0.3, outline=FALSE, col="white", border="white", names=xlab, ylab=ylab, main=main);
	boxplot(mean_case, boxwex=0.3, add=TRUE, outline=FALSE, show=FALSE, col="#D95F02", at=1:n-0.2, lwd=2);
	boxplot(mean_control, boxwex=0.3, add=TRUE, outline=FALSE, show=FALSE, col="#1B9E77", at=1:n+0.2, lwd=2);
	legend("topright", legend=leg.txt, pch=15, col=leg.color, y.intersp=1.2, xjust=0, pt.cex=2, bty="o");
	rm(pData_case, pData_control, eset_case, eset_control);
	temp <- dev.off();
}

eBayes_DiffGene <- function(eset, pData) {
	samplelist <- factor(pData$Category);
	design <- model.matrix(~-1 + samplelist);
	colnames(design) <- c("Case", "Control");
	fit <- lmFit(eset, design);
	contrast.matrix <- makeContrasts(contrasts="Case-Control", levels=design);
	fit2 <- contrasts.fit(fit, contrast.matrix);
	fit2 <- eBayes(fit2, 0.01);
	dif <- topTable(fit2, adjust="fdr", n=nrow(fit2), sort.by="B");
	dif <- dif[!is.na(dif$logFC),];
	dif <- signif(dif, 6);
	return(dif);
}

RankProd_DiffGene <- function(eset, pData) {
	eset.cl <- as.character(pData$Category);
	for (i in c(1:length(eset.cl))) {
		if (eset.cl[i]=="Case") {
			eset.cl[i] <- 1;
		}
		else if (eset.cl[i]=="Control") {
			eset.cl[i] <- 0;
		}
	}
	eset.cl <- as.numeric(eset.cl);

	eset.origin <- as.numeric(factor(pData$Study));
	
	eset.gnames <- rownames(eset);
	rownames(eset) <- c(1:length(eset.gnames));
	
	RP.out <- RP(eset, eset.cl, eset.origin, num.perm=100, logged=TRUE, na.rm=FALSE, gene.names=eset.gnames, plot=FALSE, rand=123);
	genelist <- topGene(RP.out, cutoff=1.0, method="pval", logged=TRUE, logbase=2, gene.names=eset.gnames);
	
	uplist <- as.data.frame(genelist$Table1);
	colnames(uplist) <- c("Gene.Index", "RP", "logFC", "adj.P.Val", "P.Value");
	uplist <- uplist[uplist$logFC < 1.0,];
	downlist <- as.data.frame(genelist$Table2);
	colnames(downlist) <- c("Gene.Index", "RP", "logFC", "adj.P.Val", "P.Value");
	downlist <- downlist[downlist$logFC >= 1.0,];
	
	dif <- rbind(uplist, downlist);
	dif$logFC <- round(log2(1/dif$logFC), 4);
	dif <- dif[order(dif$RP),];
	return(dif);
}

plot_logFC <- function(dif, label="All Genes", FC_cutoff=2.0, height=600, width=800, savefile="png", filename="log(FC) scatterplot.png") {
	dif <- dif[order(rownames(dif)),];
	count <- length(dif$logFC);
	xlim <- c(1,count);
	ylim <- c(min(dif$logFC), max(dif$logFC));
	if (savefile == "png") {
		png(height=height, width=width, filename=filename);
	}
	par(mar=c(5, 5, 4, 1));
	plot(0, 0, xlim=xlim, ylim=ylim, col="white", main=label, xlab="Genes", ylab="log2(FC)", cex.main=2, cex.axis=1.5, cex.lab=1.5);
	high_count <- 0;
	low_count <- 0;
	for (i in c(1:count)) {
		val <- dif$logFC[i];
		if (val >= 0) {
			points(i, val, pch=20, col="#C3524F", cex=1);
			high_count <- high_count + 1;
		}
		else if (val < 0) {
			points(i, val, pch=20, col="#4D81BE", cex=1);
			low_count <- low_count + 1;
		}
	}
	abline(h=0);
	abline(h=c(-log2(FC_cutoff),log2(FC_cutoff)), lty=2);
	legend("topright", pch=c(16,16), col=c("#C3524F", "#4D81BE"), legend=c(high_count, low_count), pt.cex=1.5, cex=1.5);
	if (savefile == "png") {
		temp <- dev.off();
	}
}

plot_dif_count <- function(up_count, down_count, height=600, width=600, savefile="png", filename="DiffGene Count.png") {
	if (savefile == "png") {
		png(height=height, width=width, filename=filename);
	}
	y_top <- max(up_count, down_count) * 1.2;
	xlab <- c("Up", "Down");
	ylab <- "Number of Differentially Expressed Genes";
	par(mar=c(3, 5, 1, 1));
	plot(0, 0, xlim=c(0, 2), ylim=c(0, y_top), ylab=ylab, xaxt="n", xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5, col="#FFFFFF");
	rect(0.25, 0, 0.75, up_count, col="#C0504D");
	rect(1.25, 0, 1.75, down_count, col="#4F81BD");
	text(0.5, up_count, labels=as.character(up_count), pos=3, cex=1.5);
	text(1.5, down_count, labels=as.character(down_count), pos=3, cex=1.5);
	axis(1, at=c(0.5, 1.5), labels=xlab, tick=FALSE, cex.axis=1.5);
	if (savefile == "png") {
		temp <- dev.off();
	}
}

plot_volcano <- function(FC_list, P_list, FC_cutoff=2.0, P_cutoff=0.05, height=800, width=800, ylab=NA, savefile="png", filename=NA) {
	if (is.na(filename)) {
		filename <- "volcano plot.png";
	}
	xlimit <- c(-max(-FC_list),max(FC_list));
	ylimit <- c(0,max(P_list));
	if (is.na(ylab)) {
		ylab <- "-log10(P)";
	}
	if (savefile == "png") {
		png(height=height, width=width, filename=filename);
	}
	par(mar=c(5, 5, 1, 1));
	plot(0, 0, col="white", xlim=xlimit, ylim=ylimit, xlab="log2(FC)", ylab=ylab, cex.main=3, cex.axis=2, cex.lab=2);
	for(i in c(1:length(FC_list))) {
		if (FC_list[i] > log2(FC_cutoff) & P_list[i] > -log10(P_cutoff)) {
			color <- "#ED1C24";
		}
		else if (FC_list[i] < -log2(FC_cutoff) & P_list[i] > -log10(P_cutoff)) {
			color <- "#22B14C";
		}
		else {
			color <- "#C5C5C5";
		}
		points(FC_list[i], P_list[i], pch=20, col=color, cex=1);
	}
	abline(h=-log10(0.05), lty=2);
	abline(v=c(-log2(FC_cutoff),log2(FC_cutoff)), lty=2);
	if (savefile == "png") {
		temp <- dev.off();
	}
}

enrichment_analysis <- function(genelist, type="GOBP", level=NA) {
	entrezIDs <- mget(genelist, org.Hs.egSYMBOL2EG, ifnotfound=NA);
	entrezIDs <- as.character(entrezIDs);
	entrezIDs <- unique(entrezIDs[entrezIDs!="NA"]);
	if (type=="GOBP" | type=="GOCC" | type=="GOMF") {
		Anno <- get("org.Hs.egGO");
		universe <- Lkeys(Anno);
		params <- new("GOHyperGParams",
					  geneIds=entrezIDs,
					  universeGeneIds=universe,
					  annotation="org.Hs.eg.db",
					  ontology=substr(type, 3, 4),
					  pvalueCutoff=0.01,
					  conditional=FALSE,
					  testDirection="over");
	}
	if (type == "KEGG") {
		Anno <- get("org.Hs.egPATH");
		universe <- Lkeys(Anno);
		params <- new("KEGGHyperGParams",
					  geneIds=entrezIDs,
					  universeGeneIds=universe,
					  annotation="org.Hs.eg.db",
					  categoryName="KEGG",
					  pvalueCutoff=0.01,
					  testDirection="over");
	}
	try(hgOver <- hyperGTest(params), silent=TRUE);
	if (!exists("hgOver")) {
		return(NULL);
	}
	enGenes <- geneIdsByCategory(hgOver);
	enGenes <- sapply(enGenes, function(.ids) {
					.sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
					.sym[is.na(.sym)] <- .ids[is.na(.sym)]
					paste(.sym, collapse=" ")
					});
	enResult <- summary(hgOver);
	if (length(rownames(enResult))==0) {
		return(NULL);
	}
	enResult$GeneSymbols <- enGenes[as.character(enResult[,1])];
	enResult$Pvalue <- signif(enResult$Pvalue, 6);
	enResult$OddsRatio <- signif(enResult$OddsRatio, 6);
	enResult$ExpCount <- signif(enResult$ExpCount, 6);
	if (type=="GOBP" & level==4 & exists("GOBPID_level4")) {
		enResult <- filter_level4_GOBP(enResult, GOBPID_level4);
	}
	return(enResult);
}

filter_level4_GOBP <- function(GOBP_list=NULL, GOBPID_level4=NULL) {
	if (is.null(GOBP_list)) {
		stop("GO biological processes list can't be NULL!");
	}
	if (is.null(GOBPID_level4)) {
		stop("Level4 GOBP list can't be NULL!");
	}
	GOBP_filter <- intersect(GOBP_list$GOBPID, GOBPID_level4$GOBPID);
	GOBP_filter <- as.data.frame(GOBP_filter);
	colnames(GOBP_filter) <- "GOBPID";
	rownames(GOBP_filter) <- GOBP_filter$GOBPID;
	rownames(GOBP_list) <- GOBP_list$GOBPID;
	return(GOBP_list[rownames(GOBP_filter),]);
}

save_as_html <- function(all_set=NULL, up_set=NULL, down_set=NULL, type="DiffGene", filename) {
	if (type == "DiffGene") {
		html_title <- "Differentially Expressed Genes";
		color <- c("#95B3D7", "#4F81BD", "#7BA0CD", "#4F81BD", "#D3DFEE");
	}
	else if (type == "GOBP") {
		html_title <- "GO Biological Process";
		color <- c("#C3D69B", "#9BBB59", "#B3CC82", "#9BBB59", "#E6EED5");
		td_width <- c("8", "9", "8", "8", "6", "6", "20", "35");
	}
	else if (type == "GOCC") {
		html_title <- "GO Cellular Component";
		color <- c("#B2A2C7", "#8064A2", "#9F8AB9", "#8064A2", "#DFD8E8");
		td_width <- c("8", "9", "8", "8", "6", "6", "20", "35");
	}
	else if (type == "GOMF") {
		html_title <- "GO Molecular Function";
		color <- c("#92CDDC", "#4BACC6", "#78C0D4", "#4BACC6", "#D2EAF1");
		td_width <- c("8", "9", "8", "8", "6", "6", "20", "35");
	}
	else if (type == "KEGG") {
		html_title <- "KEGG Pathway";
		color <- c("#FAC08F", "#F79646", "#F9B074", "#F79646", "#FDE4D0");
		td_width <- c("8", "9", "8", "8", "6", "6", "20", "35");
	}
	else if (type == "TF") {
		html_title <- "Transcription factor";
		color <- c("#D99694", "#C0504D", "#CF7B79", "#C0504D", "#EFD3D2");
		td_width <- c("13", "13", "15", "20", "13", "13", "13");
	}
	html_line <- c("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">");
	html_line <- paste(html_line, "<html>\n\t<head>\n\t<meta charset=\"utf-8\" />\n\t<meta http-equiv=\"Content-Language\" content=\"zh-cn\" />", sep="\n");
	html_line <- paste(html_line, "\t<meta http-equiv=\"Content-Type\" content=\"text/html\"; charset=\"gb2312\" />\n\t<title>", sep="\n");
	html_line <- paste(html_line, html_title, "</title>", sep=""); # Name
	html_line <- paste(html_line, "\t<script language=\"JavaScript\" type=\"text/JavaScript\">\n\tfunction toggle(targetid){", sep="\n");
	html_line <- paste(html_line, "\t\tif (document.getElementById){\n\t\t\ttarget=document.getElementById(targetid);", sep="\n");
	html_line <- paste(html_line, "\t\t\tall=document.getElementById(\"all\");\n\t\t\tup=document.getElementById(\"up\");", sep="\n");
	html_line <- paste(html_line, "\t\t\tdown=document.getElementById(\"down\");\n\t\t\tif (targetid==\"all\") {", sep="\n");
	html_line <- paste(html_line, "\t\t\t\tall.style.display=\"block\";\n\t\t\t\tup.style.display=\"none\";", sep="\n");
	html_line <- paste(html_line, "\t\t\t\tdown.style.display=\"none\";\n\t\t\t}\n\t\t\tif (targetid==\"up\") {", sep="\n");
	html_line <- paste(html_line, "\t\t\t\tall.style.display=\"none\";\n\t\t\t\tup.style.display=\"block\";", sep="\n");
	html_line <- paste(html_line, "\t\t\t\tdown.style.display=\"none\";\n\t\t\t}\n\t\t\tif (targetid==\"down\") {", sep="\n");
	html_line <- paste(html_line, "\t\t\t\tall.style.display=\"none\";\n\t\t\t\tup.style.display=\"none\";", sep="\n");
	html_line <- paste(html_line, "\t\t\t\tdown.style.display=\"block\";\n\t\t\t}\n\t\t}\n\t}\n\t</script>", sep="\n");
	html_line <- paste(html_line, "\t<style type=\"text/css\">\n\tbody {\n\t\tpadding: 0;\n\t\tmargin: 0;\n\t\ttext-align: center;", sep="\n");
	html_line <- paste(html_line, "\t\tbackground-color: #F0F0F0;\n\t\tfont-family: Times New Roman, Calibri;\n\t}\n\t#menu {", sep="\n");
	html_line <- paste(html_line, "\t\tpadding-top: 50px;\n\t\theight: 50px;\n\t\twidth: 900px;\n\t}\n\t#menu ul {\n\t\tfloat: left;", sep="\n");
	html_line <- paste(html_line, "\t\tpadding: 0;\n\t\tmargin: 0;\n\t\theight: 50px;\n\t\tline-height: 50px;\n\t}\n\t#menu ul li {", sep="\n");
	html_line <- paste(html_line, "\t\tfloat: left;\n\t\twidth: 120px;\n\t\tmargin-left: 80px;\n\t\tmargin-right: 80px;", sep="\n");
	html_line <- paste(html_line, "\t\tfont-size: 30px;\n\t\tfont-weight: bold;\n\t\tcolor: #FFFFFF;\n\t\tbackground-color: ", sep="\n");
	html_line <- paste(html_line, color[1], ";", sep=""); # button off color
	html_line <- paste(html_line, "\t\tborder-radius: 10px 10px 0 0;\n\t\tdisplay: inline;\n\t}\n\t#menu ul li:hover {\n\t\tbackground-color: ", sep="\n");
	html_line <- paste(html_line, color[2], ";", sep=""); # button on color
	html_line <- paste(html_line, "\t}\n\thr {\n\t\tmargin: 0;\n\t\tborder-width: 3px;\n\t\tborder-color: ", sep="\n");
	html_line <- paste(html_line, color[2], ";", sep=""); # line color, as same as button on color
	html_line <- paste(html_line, "\t\tborder-style: solid;\n\t}\n\ttable {\n\t\tmargin-top: -10px;\n\t\tmargin-bottom: 20px;\n\t\twidth: 100%;", sep="\n");
	html_line <- paste(html_line, "\t\tline-height: 25px;\n\t\tborder-collapse: collapse;\n\t\tborder-width: 1px;\n\t\tborder-color: ", sep="\n");
	html_line <- paste(html_line, color[3], ";", sep=""); # table border color
	html_line <- paste(html_line, "\t\tborder-style: solid;\n\t}\n\tdiv#all {\n\t\twidth: 90%;\n\t\tdisplay: block;\n\t\tcolor: #000000;\n\t}", sep="\n");
	html_line <- paste(html_line, "\tdiv#up {\n\t\twidth: 90%;\n\t\tdisplay: none;\n\t\tcolor: #ED1C24;\n\t}\n\tdiv#down{\n\t\twidth: 90%;", sep="\n");
	html_line <- paste(html_line, "\t\tdisplay: none;\n\t\tcolor: #22B14C;\n\t}\n\ttr#row_top {\n\t\tcolor: #FFFFFF;\n\t\tfont-weight: bold;", sep="\n");
	html_line <- paste(html_line, "\t\tbackground-color: ", sep="\n");
	html_line <- paste(html_line, color[4], ";", sep=""); # table title background color
	html_line <- paste(html_line, "\t}\n\ttr#row1 {\n\t\tbackground-color: ", sep="\n");
	html_line <- paste(html_line, color[5], ";", sep=""); # table "row1" background color, "row2" is "#FFFFFF"
	html_line <- paste(html_line, "\t}\n\ttr#row2 {\n\t\tbackground-color: #FFFFFF;\n\t}\n\ttd {\n\t\ttext-align: center;\n\t}\n\tp {", sep="\n");
	html_line <- paste(html_line, "\t\ttext-align: right;\n\t\tpadding-right: 5px;\n\t\tfont-weight: bold;\n\t}\n\t</style>\n\t</head>", sep="\n");
	html_line <- paste(html_line, "\t<body>\n\t\t<center>\n\t\t<div id =\"menu\">\n\t\t\t<ul>\n\t\t\t\t<li onclick=\"toggle('all')\">All</li>", sep="\n");
	html_line <- paste(html_line, "\t\t\t\t<li onclick=\"toggle('up')\">Up</li>\n\t\t\t\t<li onclick=\"toggle('down')\">Down</li>", sep="\n");
	html_line <- paste(html_line, "\t\t\t</ul>\n\t\t</div>\n\t\t<hr />", sep="\n");
	# All list
	if (!is.null(all_set) & length(rownames(all_set))>0) {
		html_line <- paste(html_line, "\t\t<div id = \"all\">\n\t\t<p>Total Count: ", sep="\n");
		html_line <- paste(html_line, as.character(length(rownames(all_set))), "</p>", sep="");
		html_line <- paste(html_line, "\t\t<table border=\"1px\">\n\t\t\t<tr id = \"row_top\">", sep="\n");
		for (i in c(1:length(colnames(all_set)))) {
			if (type != "DiffGene") {
				temp_line <- paste("\t\t\t\t<td width=\"", td_width[i], "%\">", colnames(all_set)[i], "</td>", sep="");
			}
			else {
				temp_line <- paste("\t\t\t\t<td>", colnames(all_set)[i], "</td>", sep="");
			}
			html_line <- paste(html_line, temp_line, sep="\n");
		}
		html_line <- paste(html_line, "\t\t\t</tr>", sep="\n");
		for (i in c(1:length(rownames(all_set)))) {
			if (i %% 2 == 1) {
				html_line <- paste(html_line, "\t\t\t<tr id = \"row1\">", sep="\n");
			}
			else {
				html_line <- paste(html_line, "\t\t\t<tr id = \"row2\">", sep="\n");
			}
			for (j in c(1:length(colnames(all_set)))) {
				temp_line <- paste("\t\t\t\t<td>", all_set[i,j], "</td>", sep="");
				html_line <- paste(html_line, temp_line, sep="\n");
			}
			html_line <- paste(html_line, "\t\t\t</tr>", sep="\n");
		}
		html_line <- paste(html_line, "\t\t</table>\n\t\t</div>", sep="\n");
	}
	else {
		html_line <- paste(html_line, "\t\t<div id = \"all\">\n\t\t</div>", sep="\n");
	}
	# Up list
	if (!is.null(up_set) & length(rownames(up_set))>0) {
		html_line <- paste(html_line, "\t\t<div id = \"up\">\n\t\t<p>Total Count: ", sep="\n");
		html_line <- paste(html_line, as.character(length(rownames(up_set))), "</p>", sep="");
		html_line <- paste(html_line, "\t\t<table border=\"1px\">\n\t\t\t<tr id = \"row_top\">", sep="\n");
		for (i in c(1:length(colnames(up_set)))) {
			if (type != "DiffGene") {
				temp_line <- paste("\t\t\t\t<td width=\"", td_width[i], "%\">", colnames(up_set)[i], "</td>", sep="");
			}
			else {
				temp_line <- paste("\t\t\t\t<td>", colnames(up_set)[i], "</td>", sep="");
			}
			html_line <- paste(html_line, temp_line, sep="\n");
		}
		html_line <- paste(html_line, "\t\t\t</tr>", sep="\n");
		for (i in c(1:length(rownames(up_set)))) {
			if (i %% 2 == 1) {
				html_line <- paste(html_line, "\t\t\t<tr id = \"row1\">", sep="\n");
			}
			else {
				html_line <- paste(html_line, "\t\t\t<tr id = \"row2\">", sep="\n");
			}
			for (j in c(1:length(colnames(up_set)))) {
				temp_line <- paste("\t\t\t\t<td>", up_set[i,j], "</td>", sep="");
				html_line <- paste(html_line, temp_line, sep="\n");
			}
			html_line <- paste(html_line, "\t\t\t</tr>", sep="\n");
		}
		html_line <- paste(html_line, "\t\t</table>\n\t\t</div>", sep="\n");
	}
	else {
		html_line <- paste(html_line, "\t\t<div id = \"up\">\n\t\t</div>", sep="\n");
	}
	# Down list
	if (!is.null(down_set) & length(rownames(down_set))>0) {
		html_line <- paste(html_line, "\t\t<div id = \"down\">\n\t\t<p>Total Count: ", sep="\n");
		html_line <- paste(html_line, as.character(length(rownames(down_set))), "</p>", sep="");
		html_line <- paste(html_line, "\t\t<table border=\"1px\">\n\t\t\t<tr id = \"row_top\">", sep="\n");
		for (i in c(1:length(colnames(down_set)))) {
			if (type != "DiffGene") {
				temp_line <- paste("\t\t\t\t<td width=\"", td_width[i], "%\">", colnames(down_set)[i], "</td>", sep="");
			}
			else {
				temp_line <- paste("\t\t\t\t<td>", colnames(down_set)[i], "</td>", sep="");
			}
			html_line <- paste(html_line, temp_line, sep="\n");
		}
		html_line <- paste(html_line, "\t\t\t</tr>", sep="\n");
		for (i in c(1:length(rownames(down_set)))) {
			if (i %% 2 == 1) {
				html_line <- paste(html_line, "\t\t\t<tr id = \"row1\">", sep="\n");
			}
			else {
				html_line <- paste(html_line, "\t\t\t<tr id = \"row2\">", sep="\n");
			}
			for (j in c(1:length(colnames(down_set)))) {
				temp_line <- paste("\t\t\t\t<td>", down_set[i,j], "</td>", sep="");
				html_line <- paste(html_line, temp_line, sep="\n");
			}
			html_line <- paste(html_line, "\t\t\t</tr>", sep="\n");
		}
		html_line <- paste(html_line, "\t\t</table>\n\t\t</div>", sep="\n");
	}
	else {
		html_line <- paste(html_line, "\t\t<div id = \"down\">\n\t\t</div>", sep="\n");
	}
	html_line <- paste(html_line, "\t\t</center>\n\t</body>\n</html>", sep="\n");
	writeLines(html_line, filename);
}

plot_oval_ratio <- function(a=2, b=1, midpoint=c(0,0), ratio=0.5, direction="left", color=c("#F79646", "#4BACC6")) {
	t <- seq(0, 2*pi, pi/360);
	x <- a * cos(t) + midpoint[1];
	y <- b * sin(t) + midpoint[2];
	coordinate <- rbind(x,y);
	polygon(x, y, border="#C5C5C5");
	precision <- 1e-4;
	x1 <- 0;
	x2 <- pi;
	while (abs(x2-x1)>precision) {
		x <- (x1 + x2) / 2.0;
		if (ratio == 0.5) {
			d <- 0.0;
		}
		else if (ratio < 0.5) {
			d <- -a*b / sqrt(b^2+a^2*tan(x)^2);
		}
		else if (ratio > 0.5) {
			d <- a*b / sqrt(b^2+a^2*tan(x)^2);
		}
		y <- d*sqrt(a^2-d^2) + a^2*asin(d/a) + (0.5-ratio)*a^2*pi;
		if (y == 0) {
			break;
		}
		else if (y > 0) {
			x1 <- x;
		}
		else if (y < 0) {
			x2 <- x;
		}
	}
	d <- d + midpoint[1];
	x1 <- coordinate[,coordinate["x",]<=d]["x",];
	y1 <- coordinate[,coordinate["x",]<=d]["y",];
	x2 <- coordinate[,coordinate["x",]>=d]["x",];
	y2 <- coordinate[,coordinate["x",]>=d]["y",];
	if (direction == "left") {
		polygon(x1, y1, col=color[1], border="#C5C5C5");
		polygon(x2, y2, col=color[2], border="#C5C5C5");
	}
	else if (direction == "right") {
		polygon(2*midpoint[1]-x1, y1, col=color[1], border="#C5C5C5");
		polygon(2*midpoint[1]-x2, y2, col=color[2], border="#C5C5C5");
	}
}

plot_enrichment <- function(up_genes=0, down_genes=0, En_up=NULL, En_down=NULL, term_count=10, divide=20, cex=1.5, filename=NA) {
	if (is.na(filename)) {
		filename <- "Enrichment Result.png";
	}
	if (length(rownames(En_up)) < term_count) {
		up_count <- length(rownames(En_up));
	}
	else {
		up_count <- term_count;
	}
	if (length(rownames(En_down)) < term_count) {
		down_count <- length(rownames(En_down));
	}
	else {
		down_count <- term_count;
	}
	max_up_len <- 0;
	for (i in c(1:up_count)) {
		if (max_up_len < nchar(as.character(En_up$Term[i]))) {
			max_up_len <- nchar(as.character(En_up$Term[i]));
		}
	}
	max_down_len <- 0;
	for (i in c(1:down_count)) {
		if (max_down_len < nchar(as.character(En_down$Term[i]))) {
			max_down_len <- nchar(as.character(En_down$Term[i]));
		}
	}
	height <- 300 + 25 * (up_count + down_count);
	width <- 17 * divide + max(8.5*(70-2*divide), 8.5*max_up_len) + max(8.5*(70-2*divide), 8.5*max_down_len) + 60;
	png(filename=filename, height=height, width=width);
	par(mar=c(0, 0, 0, 0), mgp=c(0, 0, 0));
	x_lim_1 <- -(8.5 * divide + max(8.5*(70-2*divide), 8.5*max_up_len)) / 9;
	x_lim_2 <- (8.5 * divide + max(8.5*(70-2*divide), 8.5*max_down_len)) / 9;
	y_lim_1 <- -(150 + 25 * down_count) / 9;
	y_lim_2 <- (150 + 25 * up_count) / 9;
	plot(0, 0, xlim=c(x_lim_1, x_lim_2), ylim=c(y_lim_1, y_lim_2), xlab="", ylab="", axes=FALSE, asp=1, col="white");
	t <- seq(0, 2*pi, pi/180);
	x <- 25 * cos(t);
	y <- 8 * sin(t);
	polygon(x, y, col="#DBE5F1", border="#4F81BD");
	rect(-5, -4, -15, 4, col="#F2DCDB", border="#C05040");
	rect(5, -4, 15, 4, col="#EBF1DD", border="#9BBB59");
	text(-10, 1.5, labels="Up", col="#953734", cex=cex);
	text(-10, -1.5, labels=as.character(up_genes), col="#953734", cex=cex);
	text(10, 1.5, labels="Down", col="#76923C", cex=cex);
	text(10, -1.5, labels=as.character(down_genes), col="#76923C", cex=cex);

	y_bottom <- 12;
	term_width <- 3;
	max_bar_len <- 30;
	
	y_top <- y_bottom + term_width * (up_count-1);
	y_pos <- seq(y_top, y_bottom, -term_width);
	log_Pvalues <- -log10(En_up$Pvalue[1:up_count]);
	bar_lens <- log_Pvalues / max(log_Pvalues) * max_bar_len;
	polygon(c(-15, -15, -divide, -divide), c(4, -4, y_bottom-1.5, y_top+1.5), col="#F2DCDB80", border="#C05040");
	polygon(c(-5, -15, -divide, -divide+70), c(4, 4, y_top+1.5, y_top+1.5), col="#F2DCDB80", border="#C05040");
	polygon(c(-5, -5, -divide+70, -divide+70), c(4, -4, y_bottom-1.5, y_top+1.5), col="#F2DCDB80", border="#C05040");
	for (i in c(1:up_count)) {
		text(-divide, y_pos[i], labels=En_up$Term[i], pos=2, col="#953734", cex=cex);
	}
	rect(-divide, y_bottom-1.5, -divide+70, y_top+1.5, col="#F2DCDB");
	for (i in c(1:up_count)) {
		rect(-divide, y_pos[i]-1, -divide+bar_lens[i], y_pos[i]+1, col="#548DD4", border="#C5C5C5");
		text(-divide+bar_lens[i]/2, y_pos[i], labels=as.character(round(log_Pvalues[i], 2)), col="#FFFFFF", cex=cex);
	}
	rect(-divide, y_bottom-1.5, -divide+70, y_top+1.5, border="#C05040");
	for (i in seq(y_top+1.5, y_bottom-1.5, -term_width)) {
		lines(c(-divide-0.3, -divide), c(i, i), col="#C05040");
	}
	for (i in c(1:up_count)) {
		ratio <- En_up$Count[i]/En_up$Size[i];
		plot_oval_ratio(10, 1, c(-divide+52, y_pos[i]), ratio, "left");
		value <- paste(as.character(En_up$Count[i]), " (", as.character(round(ratio*100,0)), "%)", sep="");
		text(-divide+42, y_pos[i], labels=value, pos=2, col="#E36C09", cex=cex);
		text(-divide+62, y_pos[i], labels=as.character(En_up$Size[i]), pos=4, col="#31859B", cex=cex);
	}
	rect(-divide, y_top+2, -divide+30, y_top+2.5, col="#548DD4", border=NA);
	text(-divide+15, y_top+4, labels="-log10(Pvalue)", col="#548DD4", cex=cex);
	rect(-divide+31, y_top+2, -divide+41, y_top+2.5, col="#F79646", border=NA);
	text(-divide+36, y_top+4, labels="Count(%)", col="#F79646", cex=cex);
	rect(-divide+63, y_top+2, -divide+69, y_top+2.5, col="#4BACC6", border=NA);
	text(-divide+66, y_top+4, labels="Size", col="#4BACC6", cex=cex);

	text(0, 0, labels="Genes", col="#366092", cex=cex);

	y_top <- y_bottom + term_width * (down_count-1);
	y_pos <- seq(-y_bottom, -y_top, -term_width);
	log_Pvalues <- -log10(En_down$Pvalue[1:down_count]);
	bar_lens <- log_Pvalues / max(log_Pvalues) * max_bar_len;
	polygon(c(15, 15, divide, divide), c(-4, 4, -y_bottom+1.5, -y_top-1.5), col="#EBF1DD80", border="#9BBB59");
	polygon(c(5, 15, divide, divide-70), c(-4, -4, -y_top-1.5, -y_top-1.5), col="#EBF1DD80", border="#9BBB59");
	polygon(c(5, 5, divide-70, divide-70), c(-4, 4, -y_bottom+1.5, -y_top-1.5), col="#EBF1DD80", border="#9BBB59");
	for (i in c(1:down_count)) {
		text(divide, y_pos[i], labels=En_down$Term[i], pos=4, col="#76923C", cex=cex);
	}
	rect(divide, -y_bottom+1.5, divide-70, -y_top-1.5, col="#EBF1DD");
	for (i in c(1:down_count)) {
		rect(divide, y_pos[i]-1, divide-bar_lens[i], y_pos[i]+1, col="#548DD4", border="#C5C5C5");
		text(divide-bar_lens[i]/2, y_pos[i], labels=as.character(round(log_Pvalues[i], 2)), col="#FFFFFF", cex=cex);
	}
	rect(divide, -y_bottom+1.5, divide-70, -y_top-1.5, border="#9BBB59");
	for (i in seq(-y_top-1.5, -y_bottom+1.5, term_width)) {
		lines(c(divide+0.3, divide), c(i, i), col="#9BBB59");
	}
	for (i in c(1:down_count)) {
		ratio <- En_down$Count[i]/En_down$Size[i];
		plot_oval_ratio(10, 1, c(divide-52, y_pos[i]), ratio, "right");
		value <- paste(as.character(En_down$Count[i]), " (", as.character(round(ratio*100,0)), "%)", sep="");
		text(divide-42, y_pos[i], labels=value, pos=4, col="#E36C09", cex=cex);
		text(divide-62, y_pos[i], labels=as.character(En_down$Size[i]), pos=2, col="#31859B", cex=cex);
	}
	rect(divide, -y_top-2, divide-30, -y_top-2.5, col="#548DD4", border=NA);
	text(divide-15, -y_top-4, labels="-log10(Pvalue)", col="#548DD4", cex=cex);
	rect(divide-31, -y_top-2, divide-41, -y_top-2.5, col="#F79646", border=NA);
	text(divide-36, -y_top-4, labels="Count(%)", col="#F79646", cex=cex);
	rect(divide-63, -y_top-2, divide-69, -y_top-2.5, col="#4BACC6", border=NA);
	text(divide-66, -y_top-4, labels="Size", col="#4BACC6", cex=cex);
	temp <- dev.off();
}

removeNA <- function(sampleData) {
	filterData <- sampleData;
	for(i in c(1:length(sampleData[,1]))) {
		NAvalue = TRUE;
		for(j in c(1:length(sampleData[1,]))) {
			if(!is.na(sampleData[i,j])) {
				NAvalue = FALSE;
				break;
			}
		}
		if(NAvalue == TRUE) {
			filterData <- filterData[rownames(filterData) != rownames(sampleData)[i],];
		}
	}
	return(filterData);
}

get_genes_from_enrichmentID <- function(en_id, type, map_list=NULL) {
	if (type == "GO") {
		all_genes <- toTable(org.Hs.egGO);
		select <- all_genes[all_genes$go_id==en_id,];
	}
	if (type == "KEGG") {
		all_genes <- toTable(org.Hs.egPATH);
		select <- all_genes[all_genes$path_id==en_id,];
	}
	if (!is.null(map_list)) {
		select <- map_list[map_list$path_id==en_id,];
	}
	select_genes <- mget(select$gene_id, org.Hs.egSYMBOL, ifnotfound=NA);
	select_genes <- unlist(select_genes);
	return(select_genes);
}

enrichmentID_to_name <- function(en_ids, type) {
	if (type == "GO") {
		en_names <- as.list(Term(GOTERM));
	}
	if (type == "KEGG") {
		en_names <- as.list(KEGGPATHID2NAME);
	}
	select_en_names <- unlist(en_names[en_ids]);
	return(select_en_names);
}

plot_ave_expr_val <- function(eset, pData, dif, en_ids, type, gap=0, FC_cutoff=1.5, Pval_cutoff=0.05, adjPval_cutoff=NULL, groupName=NULL, 
	map_list=NULL, inner_color=NULL, outer_color=NULL, label=c("Case", "Control"), filename="Average Expression Intensity.png") {
	pData_case <- pData[pData$Category=="Case",];
	pData_control <- pData[pData$Category=="Control",];
	eset_case <- eset[,colnames(t(pData_case))];
	eset_control <- eset[,colnames(t(pData_control))];
	eset_mean <- data.frame(Case=apply(eset_case, 1, mean), Control=apply(eset_control, 1, mean));
	rownames(eset_mean) <- rownames(eset);
	rm(eset_case, eset_control);
	for (i in c(1:length(rownames(dif)))) {
		if (is.null(adjPval_cutoff)) {
			if ((abs(dif$logFC[i]) >= log2(FC_cutoff)) & (dif$P.Value[i] <= Pval_cutoff)) {
				dif$isDiffExpr[i] <- 1;
			}
			else {
				dif$isDiffExpr[i] <- 0;
			}
		}
		else {
			if ((abs(dif$logFC[i]) >= log2(FC_cutoff)) & (dif$adj.P.Val[i] <= adjPval_cutoff)) {
				dif$isDiffExpr[i] <- 1;
			}
			else {
				dif$isDiffExpr[i] <- 0;
			}
		}
	}
	en_num <- length(en_ids);
	border_color <- "#C5C5C5";
	gene_counts <- c();
	for (i in c(1:en_num)) {
		select_genes <- get_genes_from_enrichmentID(en_ids[i], type, map_list);
		select_genes <- intersect(select_genes, rownames(eset_mean));
		gene_counts[i] <- length(select_genes);
	}
	en_filter <- gene_counts <= 500;
	en_ids <- en_ids[en_filter];
	gene_counts <- gene_counts[en_filter];
	en_num <- length(en_ids);
	if (is.null(inner_color)) {
		inner_color <- rainbow(en_num);
	}
	if (is.null(outer_color)) {
		outer_color <- rainbow(en_num, alpha=0.5);
	}
	total_gene_counts <- sum(gene_counts);
	en_names <- enrichmentID_to_name(en_ids, type);
	anno_label <- c(label, "Up-regulated Gene", "Down-regulated Gene");
	max_label_len <- 0;
	for (i in c(1:length(anno_label))) {
		if (max_label_len < string_width(anno_label[i], cex=2)) {
			max_label_len <- string_width(anno_label[i], cex=2);
		}
	}
	max_name_len <- 0;
	for (i in c(1:en_num)) {
		if (max_name_len < string_width(en_names[i], cex=2)) {
			max_name_len <- string_width(en_names[i], cex=2);
		}
	}
	xlim <- c(-12, 14 + max(max_label_len, max_name_len) / 40);
	ylim <- c(-12, 12);
	outer_circle_R = 10;
	label_circle_R = c(7.6, 8.4);
	middle_circle_R = 5;
	inner_circle_R = 4;
	png(filename=filename, width=1040+max(max_label_len, max_name_len), height=960);
	par(mar=c(0, 0, 0, 0));
	plot(0, 0, xlim=xlim, ylim=ylim, xlab="", ylab="", axes=FALSE, xaxs="i", yaxs="i", col="white");
	for (i in c(1:en_num)) {
		select_genes <- get_genes_from_enrichmentID(en_ids[i], type, map_list);
		select_genes <- intersect(select_genes, rownames(eset_mean));
		select_eset_mean <- eset_mean[select_genes,];
		select_dif <- dif[select_genes,];
		interval <- 2 * pi * c(sum(gene_counts[0:(i-1)]), sum(gene_counts[0:i])) / total_gene_counts;
		t <- seq(interval[1]+gap, interval[2]-gap, (interval[2]-interval[1]-2*gap) / gene_counts[i]);
		s <- sort(t, decreasing=TRUE);
		x <- c(inner_circle_R*sin(t), outer_circle_R*sin(s));
		y <- c(inner_circle_R*cos(t), outer_circle_R*cos(s));
		polygon(x, y, col=outer_color[i], border=border_color);
		x <- c(inner_circle_R*sin(t), middle_circle_R*sin(s));
		y <- c(inner_circle_R*cos(t), middle_circle_R*cos(s));
		polygon(x, y, col=inner_color[i], border=border_color);
		for (j in c(1:gene_counts[i])) {
			case_val = outer_circle_R + (eset_mean$Case[j] / 10.0);
			x1 <- c(outer_circle_R*sin(t[j]), outer_circle_R*sin(t[j+1]), case_val*sin(t[j+1]), case_val*sin(t[j]));
			y1 <- c(outer_circle_R*cos(t[j]), outer_circle_R*cos(t[j+1]), case_val*cos(t[j+1]), case_val*cos(t[j]));
			polygon(x1, y1, col="#D95F02", border=NA);
			control_val = outer_circle_R - (eset_mean$Control[j] / 10.0);
			x2 <- c(outer_circle_R*sin(t[j]), outer_circle_R*sin(t[j+1]), control_val*sin(t[j+1]), control_val*sin(t[j]));
			y2 <- c(outer_circle_R*cos(t[j]), outer_circle_R*cos(t[j+1]), control_val*cos(t[j+1]), control_val*cos(t[j]));
			polygon(x2, y2, col="#1B9E77", border=NA);
			if (select_dif$isDiffExpr[j] == 1) {
				x3 <- c(label_circle_R[1]*sin(t[j]), label_circle_R[1]*sin(t[j+1]), label_circle_R[2]*sin(t[j+1]), label_circle_R[2]*sin(t[j]));
				y3 <- c(label_circle_R[1]*cos(t[j]), label_circle_R[1]*cos(t[j+1]), label_circle_R[2]*cos(t[j+1]), label_circle_R[2]*cos(t[j]));
				if(dif$logFC[j] > 0) {
					polygon(x3, y3, col = "#FF0000", border="#FF0000");
				}
				if(dif$logFC[j] < 0) {
					polygon(x3, y3, col = "#00FF00", border="#00FF00");
				}
			}
		}
	}
	text(0, 0, groupName, cex=3);
	rect(12, 11.4, 13.9 + max_label_len / 40, 7.4, lwd=2);
	lines(c(12.5, 13.5), c(10.75, 10.75), col="#D95F02", lwd=4);
	lines(c(12.5, 13.5), c(9.85, 9.85), col="#1B9E77", lwd=4);
	lines(c(12.5, 13.5), c(8.95, 8.95), col="#FF0000", lwd=4);
	lines(c(12.5, 13.5), c(8.05, 8.05), col="#00FF00", lwd=4);
	for (i in c(1:length(anno_label))) {
		text(13.7, 11.65-0.9*i, anno_label[i], cex=2, pos=4);
	}
	rect(12, 7.1, 13.4 + max_name_len / 40, 6.7 - 0.9 * en_num, lwd=2);
	for (i in c(1:en_num)) {
		rect(12.5, 7.6-0.9*i, 13, 7.1-0.9*i, col=inner_color[i], border=outer_color[i]);
		text(13.2, 7.35-0.9*i, en_names[i], cex=2, pos=4);
	}
	temp <- dev.off();
}

plot_pathview <- function(logFC, kegg_ids, species="hsa", gene.idtype="SYMBOL", groupName, kegg_dir, filepath) {
	cur_dir <- getwd();
	setwd(filepath);
	for (i in c(1:length(kegg_ids))) {
		pv.out <- pathview(gene.data=logFC, pathway.id=kegg_ids[i], kegg.dir=kegg_dir, out.suffix=groupName, 
			species=species, gene.idtype=gene.idtype, kegg.native=T, same.layer=F, cex=0.2);
	}
	setwd(cur_dir);
}

TFs_analysis <- function(dif, dif2, TF_list, TF_target, top_num=5, TFPath, groupName) {
	all_TF_name <- data.frame(Symbol=intersect(rownames(dif), rownames(TF_list)));
	rownames(all_TF_name) <- all_TF_name$Symbol;
	all_TF <- dif[rownames(all_TF_name),];
	cat("Save transcription factor as txt format\n");
	filename <- paste(TFPath, "/all_transcription_factor (", groupName, ").txt", sep="");
	write.table(all_TF, filename, quote=FALSE, sep="\t");
	dif_TF_name <- data.frame(Symbol=intersect(rownames(dif2), rownames(TF_list)));
	rownames(dif_TF_name) <- dif_TF_name$Symbol;
	dif_TF <- dif2[rownames(dif_TF_name),];
	filename <- paste(TFPath, "/dif_transcription_factor (", groupName, ").txt", sep="");
	write.table(dif_TF, filename, quote=FALSE, sep="\t");
	dif_TF <- dif_TF[order(abs(dif_TF$logFC), decreasing=TRUE),];
	dif_num <- length(rownames(dif_TF));
	if (dif_num >= top_num) {
		top_TF <- rownames(dif_TF)[1:top_num];
	}
	else {
		top_TF <- rownames(dif_TF);
	}
	if (length(top_TF) > 0) {
		cat("Drawing top differential expressed transcription factors\n");
		for (i in c(1:length(top_TF))) {
			filename <- paste(TFPath, "/", top_TF[i], "_target_regulation (", groupName, ").png", sep="");
			plot_TF_target(dif, top_TF[i], TF_target, filename);
		}
	}
	select_TF_list <- TF_list[rownames(dif_TF),];
	dif_TF_list <- cbind(select_TF_list[,3:5], dif_TF$logFC, dif_TF$P.Value, dif_TF$adj.P.Val);
	colnames(dif_TF_list) <- c("Family", "Species", "DNA-binding domain", "logFC", "P.Value", "adj.P.Val");
	dif_TF_list <- dif_TF_list[order(dif_TF_list$Family),];
	return(dif_TF_list);
}

FC_color <- function(FC) {
	if (FC >= 2) {
		FC <- 2;
	}
	if (FC <= -2) {
		FC <- -2;
	}
	if (FC == 0) {
		r <- as.hexmode(0);
		g <- as.hexmode(0);
	}
	if (FC > 0) {
		r <- as.hexmode(as.integer(FC*127.5));
		g <- as.hexmode(0);
	}
	if (FC < 0) {
		r <- as.hexmode(0);
		g <- as.hexmode(as.integer(-FC*127.5));
	}
	if (r < 16) {
		r <- paste("0", as.character(r), sep="");
	}
	if (g < 16) {
		g <- paste("0", as.character(g), sep="");
	}
	color <- paste("#", as.character(r), as.character(g), "00", sep="");
	color <- toupper(color);
	return(color);
}

plot_TF_target <- function(dif, TF, TF_target, filename) {
	select_TF_target <- TF_target[TF_target$Name==TF,];
	if (length(select_TF_target$Target) == 0) {
		return(0);
	}
	select_TF_target <- select_TF_target[!duplicated(select_TF_target$Target),];
	rownames(select_TF_target) <- select_TF_target$Target;
	select_TF_target <- cbind(select_TF_target[,1:3], dif[rownames(select_TF_target),]);
	select_TF_target <- select_TF_target[!is.na(select_TF_target$logFC),];
	if (length(select_TF_target$Target) == 0) {
		return(0);
	}
	select_TF_target_A <- select_TF_target[select_TF_target$Type=="Activation",];
	select_TF_target_A <- select_TF_target_A[!is.na(select_TF_target_A$logFC),];
	select_TF_target_R <- select_TF_target[select_TF_target$Type=="Repression",];
	select_TF_target_R <- select_TF_target_R[!is.na(select_TF_target_R$logFC),];
	select_TF_target_U <- select_TF_target[select_TF_target$Type=="Unknown",];
	select_TF_target_U <- select_TF_target_U[!is.na(select_TF_target_U$logFC),];
	max_len <- max(length(select_TF_target_A[,1]), length(select_TF_target_R[,1]), length(select_TF_target_U[,1]));
	if (max_len < 6) {
		max_len <- 6;
	}
	ylim <- c(-2*max_len-10, 3);
	png(filename=filename, width=640, height=130+20*max_len);
	par(mar=c(0, 0, 0, 0), xaxt="n", yaxt="n");
	plot(0, 0, xlim=c(-24, 40), ylim=ylim, xaxs="i", yaxs="i", bty="n", col="#FFFFFF", asp=1);
	rect(-20.2, -0.2, -5, 0.2, col="#FF0000", border=NA);
	rect(-20.2, 0.2, -19.8, -8, col="#FF0000", border=NA);
	polygon(c(-20.4, -19.6, -20), c(-7.6, -7.6, -8.6), col="#FF0000", border=NA);
	rect(-0.2, -2, 0.2, -8, col="#0000FF", border=NA);
	rect(-0.6, -8, 0.6, -8.4, col="#0000FF", border=NA);
	rect(5, -0.2, 20.2, 0.2, col="#AAAAAA", border=NA);
	rect(19.8, 0.2, 20.2, -8, col="#AAAAAA", border=NA);
	plot_circle(20, -8, 0.4, col="#AAAAAA", border=NA);
	rect(-5, -2, 5, 2, border="#C5C5C5", col=FC_color(dif[TF,"logFC"]), lwd=3);
	text(0, 0, TF, col="#FFFFFF", font=2);
	if (length(select_TF_target_A$Target) > 0) {
		for (i in c(1:length(select_TF_target_A$Target))) {
			rect(-23, -7-2*i, -17, -7-2*(i+1), col=FC_color(select_TF_target_A$logFC[i]), border="#C5C5C5");
			text(-17, -8-2*i, select_TF_target_A$Target[i], pos=4);
		}
	}
	if (length(select_TF_target_R$Target) > 0) {
		for (i in c(1:length(select_TF_target_R$Target))) {
			rect(-3, -7-2*i, 3, -7-2*(i+1), col=FC_color(select_TF_target_R$logFC[i]), border="#C5C5C5");
			text(3, -8-2*i, select_TF_target_R$Target[i], pos=4);
		}
	}
	if (length(select_TF_target_U$Target) > 0) {
		for (i in c(1:length(select_TF_target_U$Target))) {
			rect(17, -7-2*i, 23, -7-2*(i+1), col=FC_color(select_TF_target_U$logFC[i]), border="#C5C5C5");
			text(23, -8-2*i, select_TF_target_U$Target[i], pos=4);
		}
	}
	text(27.5, 0, "Regulation Type", font=2, pos=4);
	rect(28, -1.8, 32, -2.2, col="#FF0000", border=NA);
	polygon(c(31.8, 31.8, 32.6), c(-1.6, -2.4, -2), col="#FF0000", border=NA);
	rect(28, -3.8, 32, -4.2, col="#0000FF", border=NA);
	rect(32, -3.4, 32.4, -4.6, col="#0000FF", border=NA)
	rect(28, -5.8, 32, -6.2, col="#AAAAAA", border=NA);
	plot_circle(32, -6, 0.4, col="#AAAAAA", border=NA);
	text(32.5, -2, "Activation", pos=4);
	text(32.5, -4, "Repression", pos=4);
	text(32.5, -6, "Unknown", pos=4);
	text(31.5, -8.5, "logFC", font=2, pos=4);
	for (i in c(0:19)) {
		rect(32, -10-0.5*i, 34, -10-0.5*(i+1), col=FC_color(2-4*i/19), border=NA);
	}
	text(34, -10, " 2.0", pos=4);
	text(34, -12.5, " 1.0", pos=4);
	text(34, -15, " 0.0", pos=4);
	text(34, -17.5, "-1.0", pos=4);
	text(34, -20, "-2.0", pos=4);
	temp <- dev.off();
}

target_analysis <- function(eset, pData, dif, dif2, targetPath, groupName, integrityTargetPath, diseaseName) {
	diseaseClass_1 <- dir(integrityTargetPath);
	for (i in c(1:length(diseaseClass_1))) {
		diseaseClass_1_Path <- paste(integrityTargetPath, diseaseClass_1[i], sep="/");
		diseaseClass_2 <- dir(diseaseClass_1_Path);
		for (j in c(1:length(diseaseClass_2))) {
			if (diseaseName == diseaseClass_2[j]) {
				diseaseClass_2_Path <- paste(diseaseClass_1_Path, diseaseClass_2[j], sep="/");
				break;
			}
		}
	}
	filename <- paste(diseaseClass_2_Path, "all_targets_symbol.txt", sep="/");
	disease_target <- read.table(file=filename, stringsAsFactors=FALSE, header=FALSE, sep="\t");
	colnames(disease_target) <- "Symbol";
	rownames(disease_target) <- disease_target$Symbol;
	all_target_name <- data.frame(Symbol=intersect(rownames(eset), rownames(disease_target)));
	rownames(all_target_name) <- all_target_name$Symbol;
	all_target <- eset[rownames(all_target_name),];
	all_target_num <- length(rownames(all_target));
	all_target_genes <- dif[rownames(all_target_name),];
	cat(paste("Save ", diseaseName, " target genes as txt format\n", sep=""));
	filename <- paste(targetPath, "/all_target_genes (", groupName, ").txt", sep="");
	write.table(all_target_genes, file=filename, quote=FALSE, sep="\t");
	pData_case <- pData[pData$Category=="Case",];
	pData_control <- pData[pData$Category=="Control",];
	all_target_case <- all_target[,colnames(t(pData_case))];
	case_num <- length(rownames(pData_case));
	all_target_control <- all_target[,colnames(t(pData_control))];
	control_mean <- data.frame(Mean=apply(all_target_control, 1, mean));
	all_target_case_logFC <- all_target_case - control_mean$Mean;
	filename <- paste(targetPath, "/all_target_heatmap (", groupName, ").png", sep="");
	pheatmap(all_target_case_logFC, show_colnames=FALSE, cellwidth=400/case_num, cellheight=900/all_target_num, filename=filename);
	dif_target_name <- data.frame(Symbol=intersect(rownames(disease_target), rownames(dif2)));
	rownames(dif_target_name) <- dif_target_name$Symbol;
	dif_target <- dif2[rownames(dif_target_name),];
	filename <- paste(targetPath, "/dif_target_genes (", groupName, ").txt", sep="");
	write.table(dif_target, file=filename, quote=FALSE, sep="\t");
	dif_target_num <- length(rownames(dif_target));
	if (dif_target_num > 0) {
		dif_target_case <- all_target_case[rownames(dif_target),];
		dif_target_control <- all_target_control[rownames(dif_target),];
		cat("Drawing differential expressed target genes\n");
		for (i in c(1:dif_target_num)) {
			x <- dif_target_case[i,];
			y <- dif_target_control[i,];
			symbol <- rownames(dif_target)[i];
			label <- paste("FDR = ", signif(dif_target$P.Value[i], 4), sep="");
			filename <- paste(targetPath, "/", symbol, "_barplot (", groupName, ").png", sep="");
			plot_bar_compare(x, y, main=symbol, label=label, filename=filename);
		}
	}
}

plot_select_logFC <- function(dif, color=c("#000000", "#F0F0F0"), FC_cutoff=2.0, filename=filename) {
	len <- length(rownames(dif));
	ylim <- c(min(dif$logFC), max(dif$logFC));
	png(filename=filename, height=300, width=600);
	par(mar=c(0, 0, 0, 0), xaxt="n", yaxt="n", xaxs="i");
	plot(0, 0, xlim=c(-20, len), ylim=ylim, col="#FFFFFF", bty="n");
	rect(0, ylim[1], len, ylim[2], col=color[2], border=NA);
	for (i in c(1:len)) {
		if (dif$logFC[i] > log2(FC_cutoff)) {
			rect(i-1, 0, i, log2(FC_cutoff), col="#C5C5C5", border=NA);
			rect(i-1, log2(FC_cutoff), i, dif$logFC[i], col="#FF0000", border=NA);
		}
		else if (dif$logFC[i] < -log2(FC_cutoff)) {
			rect(i-1, 0, i, -log2(FC_cutoff), col="#C5C5C5", border=NA);
			rect(i-1, -log2(FC_cutoff), i, dif$logFC[i], col="#00FF00", border=NA);
		}
		else {
			rect(i-1, 0, i, dif$logFC[i], col="#C5C5C5", border=NA);
		}
	}
	abline(h=0);
	abline(h=c(-log2(FC_cutoff), log2(FC_cutoff)), lty=2);
	rect(-20, ylim[1], 0, ylim[2], col=color[1], border=NA);
	rect(-20, ylim[1], len, ylim[2], border=color[1], lwd=2);
	temp <- dev.off();
}

plot_bar_compare <- function(x, y, main=NA, label=NA, xlab=c("Case", "Control"), ylab="Expression Value", filename) {
	x_mean <- mean(as.numeric(x));
	x_sd <- sd(as.numeric(x));
	y_mean <- mean(as.numeric(y));
	y_sd <- sd(as.numeric(y));
	ylim <- c(0, max(x_mean + x_sd, y_mean + y_sd) * 1.25);
	png(filename=filename, height=600, width=600);
	par(mar=c(3, 5, 4, 1));
	plot(0, 0, xlim=c(0, 2), ylim=ylim, main=main, xlab="", ylab=ylab, xaxt="n", xaxs="i", yaxs="i", 
		cex.main=3, cex.lab=2, cex.axis=2, col="#FFFFFF");
	rect(0.25, 0, 0.75, x_mean, col="#F2DCDB", border="#C0504D", lwd=3);
	rect(1.25, 0, 1.75, y_mean, col="#DBE5F1", border="#4F81BD", lwd=3);
	lines(c(0.5, 0.5), c(x_mean, x_mean + x_sd), col="#C0504D", lwd=3);
	lines(c(1.5, 1.5), c(y_mean, y_mean + y_sd), col="#4F81BD", lwd=3);
	lines(c(0.35, 0.65), c(x_mean + x_sd, x_mean + x_sd), col="#C0504D", lwd=3);
	lines(c(1.35, 1.65), c(y_mean + y_sd, y_mean + y_sd), col="#4F81BD", lwd=3);
	lines(c(0.5, 0.5), c(x_mean + x_sd + ylim[2]*0.02, ylim[2]*0.85), col="#000000", lwd=3);
	lines(c(1.5, 1.5), c(y_mean + y_sd + ylim[2]*0.02, ylim[2]*0.85), col="#000000", lwd=3);
	lines(c(0.5, 1.5), c(ylim[2]*0.85, ylim[2]*0.85), col="#000000", lwd=3);
	text(1, ylim[2]*0.85, label, cex=2, pos=3);
	axis(1, at=c(-1, 3), lwd=1);
	axis(1, at=c(0.5, 1.5), labels=xlab, tick=FALSE, cex.axis=2);
	temp <- dev.off();
}

common_2 <- function(L1, L2) {
	return(intersect(L1, L2));
}

common_3 <- function(L1, L2, L3) {
	L0 <- intersect(L1, L2);
	return(intersect(L0, L3));
}

common_4 <- function(L1, L2, L3, L4) {
	L0 <- intersect(L1, L2);
	L0 <- intersect(L0, L3);
	return(intersect(L0, L4));
}

common_5 <- function(L1, L2, L3, L4, L5) {
	L0 <- intersect(L1, L2);
	L0 <- intersect(L0, L3);
	L0 <- intersect(L0, L4);
	return(intersect(L0, L5));
}

plot_venn2 <- function(L1, L2, label=c("A", "B"), color=c("#C0504D", "#4F81BD"), cex.num=3, cex.lab=3, filename) {
	if (length(c(L1, L2)) == 0) {
		return;
	}
	n12 <- length(common_2(L1, L2));
	duplicate <- common_2(L1, L2);
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	n1 <- length(L1);
	n2 <- length(L2);
	png(filename=filename, height=402, width=546);
	par(mar=c(0, 0, 0, 0));
	plot(c(0, 0), xlim=c(-9.1, 9.1), ylim=c(-6.1, 7.3), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", col="#FFFFFF");
	plot_circle(-3, 0, 6, col=paste(color[1], "80", sep=""), border=NA, lwd=2);
	plot_circle(3, 0, 6, col=paste(color[2], "80", sep=""), border=NA, lwd=2);
	plot_circle(-3, 0, 6, col=NA, border="#C5C5C5", lwd=2);
	plot_circle(3, 0, 6, col=NA, border="#C5C5C5", lwd=2);
	text(-5.5, 0, labels=as.character(n1), cex=cex.num);
	text(0, 0, labels=as.character(n12), cex=cex.num);
	text(5.5, 0, labels=as.character(n2), cex=cex.num);
	text(-4.5, 6.7, labels=label[1], col=color[1], cex=cex.lab);
	text(4.5, 6.7, labels=label[2], col=color[2], cex=cex.lab);
	temp <- dev.off();
}

plot_venn3 <- function(L1, L2, L3, label=c("A", "B", "C"), color=c("#C0504D", "#4F81BD", "#9BBB59"), cex.num=3, cex.lab=3, filename) {
	if (length(c(L1, L2, L3)) == 0) {
		return;
	}
	n123 <- length(common_3(L1, L2, L3));
	duplicate <- common_3(L1, L2, L3);
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	n12 <- length(common_2(L1, L2));
	n13 <- length(common_2(L1, L3));
	n23 <- length(common_2(L2, L3));
	duplicate <- unique(c(common_2(L1, L2), common_2(L1, L3), common_2(L2, L3)));
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	n1 <- length(L1);
	n2 <- length(L2);
	n3 <- length(L3);
	png(filename=filename, height=591, width=546);
	par(mar=c(0, 0, 0, 0));
	plot(c(0, 0), xlim=c(-9.1, 9.1), ylim=c(-12.5, 7.2), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", col="#FFFFFF");
	plot_circle(-3, 0, 6, col="#C0504D80", border=NA, lwd=2);
	plot_circle(3, 0, 6, col="#4F81BD80", border=NA, lwd=2);
	plot_circle(0, -5.2, 6, col="#9BBB5980", border=NA, lwd=2);
	plot_circle(-3, 0, 6, col=NA, border="#C5C5C5", lwd=2);
	plot_circle(3, 0, 6, col=NA, border="#C5C5C5", lwd=2);
	plot_circle(0, -5.2, 6, col=NA, border="#C5C5C5", lwd=2);
	text(-5.5, 2, labels=as.character(n1), cex=cex.num);
	text(5.5, 2, labels=as.character(n2), cex=cex.num);
	text(0, -8, labels=as.character(n3), cex=cex.num);
	text(0, 3, labels=as.character(n12), cex=cex.num);
	text(-4, -4, labels=as.character(n13), cex=cex.num);
	text(4, -4, labels=as.character(n23), cex=cex.num);
	text(0, -2, labels=as.character(n123), cex=cex.num);
	text(-4.5, 6.6, labels=label[1], col="#C0504D", cex=cex.lab);
	text(4.5, 6.6, labels=label[2], col="#4F81BD", cex=cex.lab);
	text(0, -11.8, labels=label[3], col="#9BBB59", cex=cex.lab);
	temp <- dev.off();
}

plot_venn4 <- function(L1, L2, L3, L4, label=c("A", "B", "C", "D"), 
	color=c("#C0504D", "#4F81BD", "#9BBB59", "#F79646"), cex.num=3, cex.lab=3, filename) {
	if (length(c(L1, L2, L3, L4)) == 0) {
		return;
	}
	n1234 <- length(common_4(L1, L2, L3, L4));
	duplicate <- common_4(L1, L2, L3, L4);
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	n123 <- length(common_3(L1, L2, L3));
	n124 <- length(common_3(L1, L2, L4));
	n134 <- length(common_3(L1, L3, L4));
	n234 <- length(common_3(L2, L3, L4));
	duplicate <- unique(c(common_3(L1, L2, L3), common_3(L1, L2, L4), common_3(L1, L3, L4), common_3(L2, L3, L4)));
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	n12 <- length(common_2(L1, L2));
	n13 <- length(common_2(L1, L3));
	n14 <- length(common_2(L1, L4));
	n23 <- length(common_2(L2, L3));
	n24 <- length(common_2(L2, L4));
	n34 <- length(common_2(L3, L4));
	duplicate <- unique(c(common_2(L1, L2), common_2(L1, L3), common_2(L1, L4), common_2(L2, L3), common_2(L2, L4), common_2(L3, L4)));
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	n1 <- length(L1);
	n2 <- length(L2);
	n3 <- length(L3);
	n4 <- length(L4);
	png(filename=filename, height=510, width=720);
	par(mar=c(0, 0, 0, 0));
	plot(c(0, 0), xlim=c(-12, 12), ylim=c(-10, 7), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", col="#FFFFFF");
	plot_oval(-4, -3, 9, 5, pi/5, col="#C0504D80", border=NA, lwd=2);
	plot_oval(0, -1, 9, 4, pi/5, col="#4F81BD80", border=NA, lwd=2);
	plot_oval(0, -1, 9, 4, -pi/5, col="#9BBB5980", border=NA, lwd=2);
	plot_oval(4, -3, 9, 5, -pi/5, col="#F7964680", border=NA, lwd=2);
	plot_oval(-4, -3, 9, 5, pi/5, col=NA, border="#C5C5C5", lwd=2);
	plot_oval(0, -1, 9, 4, pi/5, col=NA, border="#C5C5C5", lwd=2);
	plot_oval(0, -1, 9, 4, -pi/5, col=NA, border="#C5C5C5", lwd=2);
	plot_oval(4, -3, 9, 5, -pi/5, col=NA, border="#C5C5C5", lwd=2);
	text(-9, 0, labels=as.character(n1), cex=cex.num);
	text(-4, 4, labels=as.character(n2), cex=cex.num);
	text(4, 4, labels=as.character(n3), cex=cex.num);
	text(9, 0, labels=as.character(n4), cex=cex.num);
	text(-5.5, 2, labels=as.character(n12), cex=cex.num);
	text(-5.5, -4.5, labels=as.character(n13), cex=cex.num);
	text(0, -8, labels=as.character(n14), cex=cex.num);
	text(0, 2, labels=as.character(n23), cex=cex.num);
	text(5.5, -4.5, labels=as.character(n24), cex=cex.num);
	text(5.5, 2, labels=as.character(n34), cex=cex.num);
	text(-3, -0.5, labels=as.character(n123), cex=cex.num);
	text(2.5, -5.5, labels=as.character(n124), cex=cex.num);
	text(-2.5, -5.5, labels=as.character(n134), cex=cex.num);
	text(3, -0.5, labels=as.character(n234), cex=cex.num);
	text(0, -3, labels=as.character(n1234), cex=cex.num);
	text(-10, 4.5, labels=label[1], col="#C0504D", cex=cex.lab);
	text(-5, 6, labels=label[2], col="#4F81BD", cex=cex.lab);
	text(5, 6, labels=label[3], col="#9BBB59", cex=cex.lab);
	text(10, 4.5, labels=label[4], col="#F79646", cex=cex.lab);
	temp <- dev.off();
}

plot_venn5 <- function(L1, L2, L3, L4, L5, label=c("A", "B", "C", "D", "E"), 
	color=c("#C0504D", "#4F81BD", "#9BBB59", "#F79646", "#4BACC6"), cex.num=3, cex.lab=3, filename) {
	if (length(c(L1, L2, L3, L4, L5)) == 0) {
		return;
	}
	n12345 <- length(common_5(L1, L2, L3, L4, L5));
	duplicate <- common_5(L1, L2, L3, L4, L5);
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	L5 <- setdiff(L5, duplicate);
	n1234 <- length(common_4(L1, L2, L3, L4));
	n1235 <- length(common_4(L1, L2, L3, L5));
	n1245 <- length(common_4(L1, L2, L4, L5));
	n1345 <- length(common_4(L1, L3, L4, L5));
	n2345 <- length(common_4(L2, L3, L4, L5));
	duplicate <- unique(c(common_4(L1, L2, L3, L4), common_4(L1, L2, L3, L5), common_4(L1, L2, L4, L5), 
						  common_4(L1, L3, L4, L5), common_4(L2, L3, L4, L5)));
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	L5 <- setdiff(L5, duplicate);
	n123 <- length(common_3(L1, L2, L3));
	n124 <- length(common_3(L1, L2, L4));
	n125 <- length(common_3(L1, L2, L5));
	n134 <- length(common_3(L1, L3, L4));
	n135 <- length(common_3(L1, L3, L5));
	n145 <- length(common_3(L1, L4, L5));
	n234 <- length(common_3(L2, L3, L4));
	n235 <- length(common_3(L2, L3, L5));
	n245 <- length(common_3(L2, L4, L5));
	n345 <- length(common_3(L3, L4, L5));
	duplicate <- unique(c(common_3(L1, L2, L3), common_3(L1, L2, L4), common_3(L1, L2, L5), common_3(L1, L3, L4), common_3(L1, L3, L5), 
						  common_3(L1, L4, L5), common_3(L2, L3, L4), common_3(L2, L3, L5), common_3(L2, L4, L5), common_3(L3, L4, L5)));
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	L5 <- setdiff(L5, duplicate);
	n12 <- length(common_2(L1, L2));
	n13 <- length(common_2(L1, L3));
	n14 <- length(common_2(L1, L4));
	n15 <- length(common_2(L1, L5));
	n23 <- length(common_2(L2, L3));
	n24 <- length(common_2(L2, L4));
	n25 <- length(common_2(L2, L5));
	n34 <- length(common_2(L3, L4));
	n35 <- length(common_2(L3, L5));
	n45 <- length(common_2(L4, L5));
	duplicate <- unique(c(common_2(L1, L2), common_2(L1, L3), common_2(L1, L4), common_2(L1, L5), common_2(L2, L3),
						  common_2(L2, L4), common_2(L2, L5), common_2(L3, L4), common_2(L3, L5), common_2(L4, L5)));
	L1 <- setdiff(L1, duplicate);
	L2 <- setdiff(L2, duplicate);
	L3 <- setdiff(L3, duplicate);
	L4 <- setdiff(L4, duplicate);
	L5 <- setdiff(L5, duplicate);
	n1 <- length(L1);
	n2 <- length(L2);
	n3 <- length(L3);
	n4 <- length(L4);
	n5 <- length(L5);
	sidelen <- 28 / (cos(0.3*pi) + cos(0.1*pi));
	png(filename=filename, height=1020, width=1050);
	par(mar=c(0, 0, 0, 0));
	plot(c(0, 0), xlim=c(-17.6, 17.4), ylim=c(-17.5, 16.5), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", col="#FFFFFF");
	plot_venn_shape(0, 14, 6*pi/180, col="#C0504D80", border=NA, lwd=2);
	plot_venn_shape(sidelen*sin(0.3*pi), 14-sidelen*cos(0.3*pi), 78*pi/180, col="#4F81BD80", border=NA, lwd=2);
	plot_venn_shape(sidelen/2, -14, 150*pi/180, col="#9BBB5980", border=NA, lwd=2);
	plot_venn_shape(-sidelen/2, -14, 222*pi/180, col="#F7964680", border=NA, lwd=2);
	plot_venn_shape(-sidelen*sin(0.3*pi), 14-sidelen*cos(0.3*pi), 294*pi/180, col="#4BACC680", border=NA, lwd=2);
	plot_venn_shape(0, 14, 6*pi/180, col=NA, border="#C5C5C5", lwd=2);
	plot_venn_shape(sidelen*sin(0.3*pi), 14-sidelen*cos(0.3*pi), 78*pi/180, col=NA, border="#C5C5C5", lwd=2);
	plot_venn_shape(sidelen/2, -14, 150*pi/180, col=NA, border="#C5C5C5", lwd=2);
	plot_venn_shape(-sidelen/2, -14, 222*pi/180, col=NA, border="#C5C5C5", lwd=2);
	plot_venn_shape(-sidelen*sin(0.3*pi), 14-sidelen*cos(0.3*pi), 294*pi/180, col=NA, border="#C5C5C5", lwd=2);
	text(-0.5, 12, labels=as.character(n1), cex=cex.num);
	text(12, 2, labels=as.character(n2), cex=cex.num);
	text(7, -12, labels=as.character(n3), cex=cex.num);
	text(-8, -12, labels=as.character(n4), cex=cex.num);
	text(-12, 3, labels=as.character(n5), cex=cex.num);
	text(4, 7.5, labels=as.character(n12), cex=cex.num);
	text(-2.5, 8, labels=as.character(n13), cex=cex.num);
	text(-3.5, -10.5, labels=as.character(n14), cex=cex.num);
	text(-7.2, 5, labels=as.character(n15), cex=cex.num);
	text(9.5, -2.5, labels=as.character(n23), cex=cex.num);
	text(8, 4, labels=as.character(n24), cex=cex.num);
	text(-10, -1, labels=as.character(n25), cex=cex.num);
	text(2, -11, labels=as.character(n34), cex=cex.num);
	text(7, -8, labels=as.character(n35), cex=cex.num);
	text(-8.5, -6.2, labels=as.character(n45), cex=cex.num);
	text(0.8, 5.5, labels=as.character(n123), cex=cex.num);
	text(5.2, 5, labels=as.character(n124), cex=cex.num);
	text(-6.2, 1.5, labels=as.character(n125), cex=cex.num);
	text(-0.5, -9.8, labels=as.character(n134), cex=cex.num);
	text(-4.5, 5.5, labels=as.character(n135), cex=cex.num);
	text(-4.8, -6.5, labels=as.character(n145), cex=cex.num);
	text(6.8, 0, labels=as.character(n234), cex=cex.num);
	text(7.8, -4.5, labels=as.character(n235), cex=cex.num);
	text(-8, -3.6, labels=as.character(n245), cex=cex.num);
	text(3.4, -7.5, labels=as.character(n345), cex=cex.num);
	text(3.5, 2.5, labels=as.character(n1234), cex=cex.num);
	text(-2.8, 3.2, labels=as.character(n1235), cex=cex.num);
	text(-5.2, -2.5, labels=as.character(n1245), cex=cex.num);
	text(-0.5, -6.8, labels=as.character(n1345), cex=cex.num);
	text(4.8, -3.8, labels=as.character(n2345), cex=cex.num);
	text(0, -1.5, labels=as.character(n12345), cex=cex.num);
	text(0.1, 14.8, labels=label[1], col="#C0504D", cex=cex.lab, srt=354);
	text(15.5, 3.5, labels=label[2], col="#4F81BD", cex=cex.lab, srt=282);
	text(9.6, -14.7, labels=label[3], col="#9BBB59", cex=cex.lab, srt=210);
	text(-9.6, -14.6, labels=label[4], col="#F79646", cex=cex.lab, srt=138);
	text(-15.5, 3.7, labels=label[5], col="#4BACC6", cex=cex.lab, srt=66);
	temp <- dev.off();
}

smart_venn <- function(x, color=c("#C0504D", "#4F81BD", "#9BBB59", "#F79646", "#4BACC6"), cex.num=3, cex.lab=3, filename) {
	if (length(x) == 2) {
		plot_venn2(x[[1]], x[[2]], names(x), color, cex.num, cex.lab, filename);
	}
	if (length(x) == 3) {
		plot_venn3(x[[1]], x[[2]], x[[3]], names(x), color, cex.num, cex.lab, filename);
	}
	if (length(x) == 4) {
		plot_venn4(x[[1]], x[[2]], x[[3]], x[[4]], names(x), color, cex.num, cex.lab, filename);
	}
	if (length(x) == 5) {
		plot_venn5(x[[1]], x[[2]], x[[3]], x[[4]], x[[5]], names(x), color, cex.num, cex.lab, filename);
	}
}

plot_circle <- function(x=0, y=0, r=1, col=NA, border="#000000", lwd=1) {
	t <- seq(0, 2*pi, pi/180);
	x1 <- x + r*cos(t);
	y1 <- y + r*sin(t);
	polygon(x1, y1, col=col, border=border, lwd=lwd);
}

plot_oval <- function(x=0, y=0, a=2, b=1, φ=0, col=NA, border="#000000", lwd=1) {
	t <- seq(0, 2*pi, pi/180);
	x1 <- a*cos(t);
	y1 <- b*sin(t);
	x2 <- x + x1*cos(φ) + y1*sin(φ);
	y2 <- y + y1*cos(φ) - x1*sin(φ);
	polygon(x2, y2, col=col, border=border, lwd=lwd);
}

plot_venn_shape <- function(x=0, y=0, p=0, col=NA, border="#000000", lwd=1) {
	t <- seq(-atan(0.9), pi+atan(0.9), pi/180);
	x1 <- 7.5*cos(t);
	y1 <- 9.0*sin(t) - 9.0;
	t <- seq(pi, 2*pi, pi/180);
	x1 <- c(x1, 3.75*cos(t));
	y1 <- c(y1, 4.5*sin(t) - 22.5);
	x2 <- x + x1*cos(p) + y1*sin(p);
	y2 <- y + y1*cos(p) - x1*sin(p);
	polygon(x2, y2, col=col, border=border, lwd=lwd);
}

string_width <- function(s, cex=1) {
	character_width <- data.frame(
		Character=c("`", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0", "-", "=", "~", "!", "@", "#", "$", "%", "^", 
			"&", "*", "(", ")", "_", "+", "q", "w", "e", "r", "t", "y", "u", "i", "o", "p", "[", "]", "\\", "Q", 
			"W", "E", "R", "T", "Y", "U", "I", "O", "P", "{", "}", "|", "a", "s", "d", "f", "g", "h", "j", "k",
			"l", ";", "'", "A", "S", "D", "F", "G", "H", "J", "K", "L", ":", "\"", "z", "x", "c", "v", "b", "n", 
			"m", ",", ".", "/", "Z", "X", "C", "V", "B", "N", "M", "<", ">", "?", " "),
		Width=c(8, 11, 19, 18, 19, 18, 18, 18, 18, 18, 18, 12, 19, 19, 5, 35, 21, 18, 29, 17, 
			23, 13, 9, 9, 22, 19, 17, 25, 18, 11, 10, 17, 15, 5, 18, 17, 9, 9, 12, 26, 
			36, 21, 24, 21, 23, 22, 5, 26, 21, 12, 12, 4, 18, 16, 18, 12, 18, 15, 9, 17,  
			5, 6, 5, 25, 22, 23, 19, 26, 22, 16, 22, 18, 6, 11, 18, 16, 16, 19, 17, 15, 
			27, 6, 6, 12, 22, 23, 24, 24, 21, 22, 25, 19, 19, 18, 6));
	character_width <- character_width[order(character_width$Character),];
	inter_width <- 4;
	s <- unlist(strsplit(s, ""));
	s_width <- inter_width * (length(s) - 1);
	character_num <- length(rownames(character_width));
	for (i in c(1:length(s))) {
		a <- 1;
		b <- character_num;
		count <- 0;
		while (a <= b) {
			j <- as.integer((a + b) / 2);
			cur_character <- as.character(character_width$Character[j]);
			if (s[i] == cur_character) {
				s_width <- s_width + character_width$Width[j];
				break;
			}
			if (s[i] > cur_character) {
				a <- j + 1;
			}
			if (s[i] < cur_character) {
				b <- j - 1;
			}
			count <- count + 1;
			if (count > log2(character_num)) {
				break;
			}
		}
	}
	s_width <- s_width * cex / 3;
	s_width <- as.integer(round(s_width, 0));
	return(s_width);
}

gradient_color <- function(x, val_list, col_list) {
	if (length(val_list) != length(col_list)) {
		cat("ERROR: value list and color list do not math!");
	}
	if (is.na(x)) {
		return("#FFFFFF");
	}
	num <- length(val_list);
	if (x <= val_list[1]) {
		return(col_list[1]);
	}
	if (x >= val_list[num]) {
		return(col_list[num]);
	}
	x_col <- c(0, 0, 0);
	for (i in c(1:num)) {
		if (x > val_list[i] & x < val_list[i+1]) {
			col_1 <- col2rgb(col_list[i], alpha=FALSE);
			col_2 <- col2rgb(col_list[i+1], alpha=FALSE);
			val_diff <- val_list[i+1] - val_list[i];
			for (j in c(1:3)) {
				a <- (col_2[j] - col_1[j]) / val_diff;
				b <- (val_list[i+1]*col_1[j] - val_list[i]*col_2[j]) / val_diff;
				x_col[j] <- as.integer(a * x + b);
			}
			return(rgb(x_col[1], x_col[2], x_col[3], maxColorValue=255));
		}
	}
}

plot_P_val <- function(Pval, x, y, color=c("#FFFF00", "#F0F0F0", "#FFFFFF"), border=NA) {
	if (is.na(Pval)) {
		plot_circle(x, y, r=0.1, col=color[3], border=border, lwd=1);
	}
	else if (Pval > 0.05) {
		a <- -2 / 19;
		b <- 39 / 190;
		r <- a * Pval + b;
		plot_circle(x, y, r, col=color[2], border=border, lwd=1);
	}
	else if (Pval > 0.001 & Pval <= 0.05) {
		a <- -200 / 49;
		b <- 99 / 245;
		r <- a * Pval + b;
		plot_circle(x, y, r, col=color[1], border=border, lwd=1);
	}
	else if (Pval <= 0.001) {
		plot_circle(x, y, r=0.4, col=color[1], border=border, lwd=1);
	}
}

create_GSEA_data <- function(eset, pData, filename) {
	genes_num <- length(rownames(eset));
	sample_num <- length(rownames(pData));
	label <- data.frame(NAME=rownames(eset), DESCRIPTION=rownames(eset));
	rownames(label) <- rownames(eset);
	eset <- cbind(label, eset);
	
	filename_gct <- paste(filename, ".gct", sep="");
	cat("#1.2\n", file=filename_gct);
	cat(paste(genes_num, "\t", sample_num, "\n", sep=""), file=filename_gct, append=TRUE);
	cat(colnames(eset), sep="\t", file=filename_gct, append=TRUE);
	cat("\n", file=filename_gct, append=TRUE);
	write.table(eset, file=filename_gct, row.names=F, col.names=F, quote=F, sep="\t", append=TRUE);
	
	filename_cls <- paste(filename, ".cls", sep="");
	cat(paste(sample_num, " 2 1\n", sep=""), file=filename_cls);
	if (as.character(pData$Category[1]) == "Case") {
		cat("# Case Control\n", file=filename_cls, append=TRUE);
	}
	else {
		cat("# Control Case\n", file=filename_cls, append=TRUE);
	}
	cat(as.character(pData$Category), sep="\t", file=filename_cls, append=TRUE);
	cat("\n", file=filename_cls, append=TRUE);
}

plot_GSEA <- function(GSEA, group_title="Group", group_name=NA, group_color=NA, sig_lab="FDR", filename=NA) {
	group_num <- (length(colnames(GSEA))-3) / 2;
	pathway_num <- length(rownames(GSEA));
	if (is.na(group_name[1])) {
		group_name <- paste("G", 1:group_num, sep="");
	}
	if (is.na(group_color[1])) {
		group_color <- rainbow(group_num);
	}
	if (is.na(filename)) {
		filename <- "GSEA Resuls.pdf";
	}
	GSEA <- GSEA[order(GSEA$Name, decreasing=TRUE),];
	GSEA <- GSEA[order(GSEA$Class, decreasing=TRUE),];
	index <- seq(4, length(colnames(GSEA)), 2);
	val_list <- c(min(GSEA[,4:length(colnames(GSEA))], na.rm=TRUE), 0, max(GSEA[,4:length(colnames(GSEA))], na.rm=TRUE));
	col_list <- c("#4F81BD", "#FFFFFF", "#C0504D");
	xlim <- c(-24, group_num * 3.5 + 14);
	ylim <- c(-0.5, pathway_num + 2);
	pdf_height <- 14;
	pdf_width <- round(pdf_height * (xlim[2] - xlim[1]) / (ylim[2] - ylim[1]), 2);
	pdf(file=filename, height=pdf_height, width=pdf_width);
	par(mar=c(0, 0, 0, 0), bty="n");
	plot(0, 0, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xaxs="i", yaxs="i", col="#FFFFFF");
	for (i in c(1:pathway_num)) {
		text(0, i-0.5, GSEA$Name[i], cex=0.8, pos=2);
	}
	for (j in c(1:group_num)) {
		rect(3.5*(j-1), pathway_num+0.5, 3.5*j, pathway_num+1.5, col=group_color[j], border="#C5C5C5");
	}
	for (i in c(1:pathway_num)) {
		for (j in c(1:group_num)) {
			color <- gradient_color(GSEA[i,index[j]], val_list, col_list);
			rect(3.5*(j-1), i-1, 3.5*j, i, col=color, border="#C5C5C5");
			plot_P_val(GSEA[i,(index[j]+1)], 3.5*(j-0.5), i-0.5);
		}
	}
	anno_y <- pathway_num - 0.5;
	text(3.5*group_num, anno_y, group_title, font=2, cex=0.8, pos=4);
	for (j in c(1:group_num)) {
		rect(3.5*(group_num+0.2), anno_y-1.5*j-0.5, 3.5*(group_num+0.6), anno_y-1.5*j+0.5, col=group_color[j], border="#C5C5C5");
		text(3.5*(group_num+0.5), anno_y-1.5*j, group_name[j], cex=0.8, pos=4);
	}
	NES_y <- anno_y - 1.5 * group_num - 2;
	text(3.5*group_num, NES_y, "Normalized Enrichment", font=2, cex=0.8, pos=4);
	text(3.5*group_num, NES_y-1, "Score (NES)", font=2, cex=0.8, pos=4);
	color_y <- seq(NES_y-17, NES_y-2, 0.3);
	color_value <- seq(val_list[1], val_list[3], (val_list[3]-val_list[1])/49);
	for (i in c(1:length(color_value))) {
		color <- gradient_color(color_value[i], val_list, col_list);
		rect(3.5*(group_num+0.2), color_y[i], 3.5*(group_num+0.8), color_y[i+1], col=color, border=NA);
	}
	for (i in seq(as.integer(val_list[1]), as.integer(val_list[3]), 0.5)) {
		label_y <- 15*(i-val_list[1])/(val_list[3]-val_list[1])+color_y[1]-0.15;
		text(3.5*(group_num+0.7), label_y, as.character(i), cex=0.8, pos=4)
	}
	text(3.5*(group_num+0.7), color_y[1]-0.15, as.character(val_list[1]), cex=0.8, pos=4);
	text(3.5*(group_num+0.7), color_y[length(color_y)]-0.15, as.character(val_list[3]), cex=0.8, pos=4);
	Pval_y <- NES_y - 19;
	text(3.5*group_num, Pval_y, "Significance", font=2, cex=0.8, pos=4);
	plot_circle(3.5*(group_num+0.5), Pval_y-1.5, r=0.8, col="#FFFF00", border=NA);
	plot_circle(3.5*(group_num+0.5), Pval_y-5.5, r=0.2, col="#C5C5C5", border=NA);
	polygon_x <- c(3.5*(group_num+0.26), 3.5*(group_num+0.38), 3.5*(group_num+0.62), 3.5*(group_num+0.74));
	polygon_y <- c(Pval_y-1.5, Pval_y-3.5, Pval_y-3.5, Pval_y-1.5);
	polygon(polygon_x, polygon_y, col="#FFFF00", border=NA);
	polygon_x <- c(3.5*(group_num+0.38), 3.5*(group_num+0.44), 3.5*(group_num+0.56), 3.5*(group_num+0.62));
	polygon_y <- c(Pval_y-3.5, Pval_y-5.5, Pval_y-5.5, Pval_y-3.5);
	polygon(polygon_x, polygon_y, col="#C5C5C5", border=NA);
	text(3.5*(group_num+0.7), Pval_y-1.5, paste(sig_lab, " < 0.001", sep=""), cex=0.8, pos=4);
	text(3.5*(group_num+0.7), Pval_y-3.5, paste(sig_lab, " = 0.050", sep=""), cex=0.8, pos=4);
	text(3.5*(group_num+0.7), Pval_y-5.5, paste(sig_lab, " > 0.999", sep=""), cex=0.8, pos=4);
	dev.off();
}

eset_classifier <- function(eset, break_num) {
	total_num <- length(colnames(eset));
	mean_1 <- apply(eset[,1:break_num], 1, mean);
	mean_2 <- apply(eset[,(break_num+1):total_num], 1, mean);
	eset$score <- mean_1 - mean_2;
	min_score <- min(abs(eset$score));
	subset_1 <- eset[eset$score<min_score, 1:break_num];
	subset_2 <- eset[eset$score<min_score, (break_num+1):total_num];
	subset_3 <- eset[eset$score>=min_score, 1:break_num];
	subset_4 <- eset[eset$score>=min_score, (break_num+1):total_num];
	subset_1 <- subset_1[sample(rownames(subset_1)), sample(colnames(subset_1))];
	subset_4 <- subset_4[sample(rownames(subset_4)), sample(colnames(subset_4))];
	subset_2 <- subset_2[rownames(subset_1), colnames(subset_4)];
	subset_3 <- subset_3[rownames(subset_4), colnames(subset_1)];
	new_eset <- rbind(cbind(subset_1, subset_2), cbind(subset_3, subset_4));
	return(new_eset);
}