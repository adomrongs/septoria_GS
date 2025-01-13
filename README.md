# Septoria Genomic Selection

Remember that the column "Leaf" in phenotype does not represent different leaves (1,2 3) but different plants (1,2,3), but all leaves are scond-leaves
Additionally, Rep is not a rep is actually trial.


Script organization:
	1_PCA: Perform PCA analysis and plot PC1 vs PC2 for both septoria and wheat. Also Genomic Relationship Matrix heatmaps
	
	2_EDA:  Exploratory Data analysis for the phenotype that includes the mixes
	
	3_prepGWAS_wheat: Prepare genotype_wheat, map_wheat, k_wheat, blues_wheat to perform GWAS on wheat
		Rdata created: 3_wheat_GWAS.Rdata
		
	4_runGWAS_wheat: Perform GWAS on wheat using BLINK and MLM with 0-4 PC
		cluster: runGWAS_wheat.slurm
		results: outputs/GWAS_wheat
		
	6_preCV: Prepare everything needed for the CV using the mixes: k_mixes, I_mixes, k_wheat, I_wheat, adjusted_phenotype, map_wheat, genotype_wheat
		Rdata created: 4_CV_data.Rdata
		
	7_runCV_all: run the CV with the three differnt strategies and GRM for both septoria and wheat\
		cluster: run_cv_all.slurm
		results: data/modified_data/cv/*
		
	8_plotCV: plot boxplot for each strategy, model and mix as well as the histogram of the gwas hit per iteration
	
	9_prepTest: prepare data to make the predictions: genotype_all, k_all, cleaned_septoria_phenotype_1, cleaned_septoria_phenotype_2
		Rdata created: 5_predictions.Rdata
		
	10_modelAIC: compare the AIC using two phenotypes (all leaves vs leaf2) and 2 different modell (with or without BRep to extract the blues)
	FInally, we selected the phenotype with all leaves combined and using the model including the Brep column.
	
	11_predictions: run the predictions of the 3 traits with the 4 different models commented before on the 19 test strains.
		cluster: run_predictions.slurm 
		results: data/modified_data/predictions/
	
	12_prediction_results: Calculate the correlation between the blues of the 19 test strains and their calculated blups. 
	
	13_prep_CV_septoria: prapre the information to run the CV on the 100 train strains: genotype, map_septoria,phenotype, kinship, blues_all
		Rdata: 6_CV_septoria.Rdata
		
	14_septoriaCV: perform 30x5-fold CV on the train isolates
		cluster: run_cv_sep.slurm
		
	15_plotCVsep: plot the results from the previous analysis for the 3 traits.
	
	16_postGWAS_wheat: plot boxplots, manhattan and qqplots + identify genes and GO associated to the significant hits (2 on PLACL)
	
	17_postGWAS_sep: same as before bnuit with septoria hits. 
	
	functions_septoria_GS: script containing all the built-in function used in this project
	

