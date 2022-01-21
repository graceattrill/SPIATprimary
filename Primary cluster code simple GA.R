#This is a simple code that returns numerical values that can be analysed (Mixing score & Kcross)
#It does not return any plots, clusters, communities or heatmaps

library(data.table)
library(stringr)
#spatstat for Kcross
library(spatstat)
library(Rphenograph)
library(SPIAT)
setwd("/Volumes/Gracie/SPIAToutput")

#22551 - no B cells
#18694 - no P3, P5, P7
#20790 - no NK cells
#21461 - no P1, NK cells, LCs
#22856 - no NK cells, LCs
#27774 - no P5

#17245 - no CD8s
#16887 - no CD8s
#22244 - no P3
#32660 - no P5
#26733 - no LCs
#25541 - no P5

combi.files<-list.files(path="/Volumes/Gracie/SPIAToutput/Combifiles", pattern="_combitable.csv", full.names = T, recursive = T, include.dirs = T)
#Set immune_phenotypes vector as all phenotypes except for Melanoma - simplifies later analysis
immune_phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8", "B cell", "NK cell", "Langerhans cell")

#define columns for SPIAT to read from
dye_columns_interest <- c("P1",
                          "P2", 
                          "P3", 
                          "P4", 
                          "P5", 
                          "P6", 
                          "P7", 
                          "P8", 
                          "Melanoma", 
                          "DAPI Positive Classification",
                          "B cell",
                          "NK cell",
                          "Langerhans cell",
                          "HLAABCpMelanoma",
                          "HLAABCnMelanoma")

#Setting dye_columns_interest & intensity_columns_interest as the same columns allows manually HALO-assigned phenotypes to be used
intensity_columns_interest <- c("P1",
                                "P2", 
                                "P3", 
                                "P4", 
                                "P5", 
                                "P6", 
                                "P7", 
                                "P8", 
                                "Melanoma", 
                                "DAPI Positive Classification",
                                "B cell",
                                "NK cell",
                                "Langerhans cell",
                                "HLAABCpMelanoma",
                                "HLAABCnMelanoma")

markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "Melanoma", "DAPI", 
             "B cell", "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma")

#All HLA-ABCp & n melanoma are also Melanoma - need to set this up for SPIAT to read correctly
phenotypes_of_interest<-c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "Melanoma","B cell", "NK cell", "Langerhans cell","Melanoma,HLAABCpMelanoma", "Melanoma,HLAABCnMelanoma")
colour_vector<-c("chartreuse", "cyan", "blue", "darkorange", "deeppink", "chocolate4", "gold", "seagreen", "red", "tan4", "violet", "azure2","yellow", "tan1")

for (i in combi.files) {
  newdat<-str_replace(i, pattern = ".object_results_combitable.csv", replacement = "SPIAT")
  
  #upload file to R
  formatted_image <- format_image_to_sce(format="HALO",
  path=i,
  markers=markers,
  dye_columns_interest=dye_columns_interest,
  intensity_columns_interest=intensity_columns_interest)
  
  #Check phenotypes have been read correctly:
  print_column(formatted_image, column= "Phenotype")
    
  #Remove any double positive B cell/NK cell or LC/NK cell and check phenotypes again
  formatted_image<-select_phenotypes(formatted_image, keep=FALSE, phenotypes = c("B cell,NK cell", "NK cell,Langerhans cell"))
  print_column(formatted_image, column = "Phenotype")
  
  #Calculate mixing score - Melanoma to immune
  mixscoretab<-data.frame(phenotype=character(),mixingscore=numeric())
  for(k in immune_phenotypes){
    mixscorenum<-compute_mixing_score(formatted_image,reference_marker="Melanoma",target_marker=k,radius=20)
    looptab<-data.frame(phenotype=paste0(k, "Mix"), score=mixscorenum)
    mixscoretab<-rbind(mixscoretab,looptab)
  }
  if(i == "/Volumes/Gracie/SPIAToutput/Combifiles/0320777 A1 22614 T cell.tif_38189_job3540.object_results_combitable.csv"){
    all.mixing.scores<-mixscoretab
  } else {
    all.mixing.scores<-cbind(all.mixing.scores,mixscoretab[,2])
  }
  
#Kcross - this calculated Kcross values for all immune phenotypes to melanoma, creates file with all Kcross values but doesn't keep Kcross graphs
Ktab<-data.frame(phenotype=character(),Kscore=numeric())
for (l in immune_phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(l, "Melanoma"), column = "Phenotype")
  pdf(paste0(l, "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(l, "Kcross"), Kscore = Kcrossn)
  Ktab<-rbind(Ktab, Klooptab)
}

if(i == "/Volumes/Gracie/SPIAToutput/Combifiles/0320777 A1 22614 T cell.tif_38189_job3540.object_results_combitable.csv"){
  all.Kcross<-Ktab
} else{
  all.Kcross<-cbind(all.Kcross,Ktab[,2])
}

}

colnames(all.Kcross) = c("Phenotype", c(combi.files))
colnames(all.Kcross)<-str_replace_all(colnames(all.Kcross), pattern = ".object_results_combitable.csv", replacement = "")
colnames(all.Kcross)<-str_replace_all(colnames(all.Kcross), pattern = "\\/Volumes\\/Gracie\\/SPIAToutput\\/Combifiles\\/", replacement = "")
column_names<-c(colnames(all.Kcross))
colnames(all.mixing.scores) = column_names
alldat<-rbind(all.Kcross, all.mixing.scores)
fwrite(alldat, file = "SPIATdat.csv", quote=F)