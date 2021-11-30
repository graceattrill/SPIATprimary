library(stringr)
library(data.table)

#Move T cell files to T cell folder & Immune cell files to Immune folder
T.files<-list.files(path="/Users/gracie/Desktop/R Projects/SPIATprimary/Files/T cellxyimmuneHALO", pattern=".csv", full.names=T, recursive=T, include.dirs=T)
I.files<-list.files(path="/Users/gracie/Desktop/R Projects/SPIATprimary/Files/ImmuneHALO", pattern=".csv", full.names=T, recursive=T, include.dirs=T)

#These 3 vectors are so here so that the columns match up between the two tables in the rbind
T.colsp<-c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "CD8+ T cells", "T cells")
T.colsm<-c("CD39 (Opal 520) Positive Classification", "CD39 (Opal 520) Positive Cytoplasm Classification", "CD39 (Opal 520) Cytoplasm Intensity", "CD39 (Opal 520) Cell Intensity", 
          "PD-1 (Opal 540) Positive Classification", "PD-1 (Opal 540) Positive Cytoplasm Classification", "PD-1 (Opal 540) Cytoplasm Intensity", "PD-1 (Opal 540) Cell Intensity",
          "CD3 (Opal 570) Positive Classification", "CD3 (Opal 570) Positive Cytoplasm Classification", "CD3 (Opal 570) Cytoplasm Intensity", "CD3 (Opal 570) Cell Intensity",
          "CD103 (Opal 620) Positive Classification", "CD103 (Opal 620) Positive Cytoplasm Classification", "CD103 (Opal 620) Cytoplasm Intensity", "CD103 (Opal 620) Cell Intensity",
          "CD8 (Opal 650) Positive Classification", "CD8 (Opal 650) Positive Cytoplasm Classification", "CD8 (Opal 650) Cytoplasm Intensity", "CD8 (Opal 650) Cell Intensity")
I.cols<-c("B cell", "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma",
          "CD20 (Opal 520) Positive Classification",	"CD20 (Opal 520) Positive Cytoplasm Classification",	"CD20 (Opal 520) Cytoplasm Intensity",	"CD20 (Opal 520) Cell Intensity",
          "CD1a (Opal 540) Positive Classification",	"CD1a (Opal 540) Positive Cytoplasm Classification",	"CD1a (Opal 540) Cytoplasm Intensity",	"CD1a (Opal 540) Cell Intensity",
          "CD56 (Opal 570) Positive Classification",	"CD56 (Opal 570) Positive Cytoplasm Classification",	"CD56 (Opal 570) Cytoplasm Intensity",	"CD56 (Opal 570) Cell Intensity",
          "HLA-ABC (Opal 620) Positive Classification",	"HLA-ABC (Opal 620) Positive Cytoplasm Classification",	"HLA-ABC (Opal 620) Cytoplasm Intensity",	"HLA-ABC (Opal 620) Cell Intensity",
          "Langerin (Opal 650) Positive Classification",	"Langerin (Opal 650) Positive Cytoplasm Classification",	"Langerin (Opal 650) Cytoplasm Intensity",	"Langerin (Opal 650) Cell Intensity")

##T cell formatting
for(i in T.files){
  print(paste0("Formatting", i))
  newdat<-str_replace(i,pattern=".csv", replacement="_combitable.csv")
  
  print("replace")
  tableT<-fread(i,data.table=F)
  #get rid of Xmin Xmax Ymin Ymax
  tableT<- subset(tableT, select = -c(7:10))
  #rename XY columns
  colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered XMin",replacement="XMin")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered XMax",replacement="XMax")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered YMin",replacement="YMin")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered YMax",replacement="YMax")
  #rename P1-8 columns
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103\\+PD-1\\+",replacement="P1")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103\\+PD-1-",replacement="P2")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103-PD-1\\+",replacement="P3")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103-PD-1-",replacement="P4")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103\\+PD-1\\+",replacement="P5")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103\\+PD-1-",replacement="P6")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103-PD-1\\+",replacement="P7")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103-PD-1-",replacement="P8")
  #add columns from immune with 0 in all rows
  out1<-matrix(0,nrow=(nrow(tableT)),ncol=25)
  colnames(out1)<-c(I.cols)
  tableT<-cbind(tableT,out1)
  print(head(tableT))
  #Replace um2  
  colnames(tableT)<-str_replace_all(colnames(tableT),pattern="\\(µm²\\)",replacement="")
  colnames(tableT)<-str_replace_all(colnames(tableT),pattern="\\(µm\\)",replacement="")
  #Set max T cell object Id +1 as minimum object Id for Immune panel
  objectidn<-max(tableT$`Object Id`)+1
  
##Immune cell formatting
    for (k in I.files){
    print(paste0("Formatting", k))
    tableIp<-fread(k,data.table=F)
    #get rid of XMin XMax YMin YMax, (& Analysis Region & Layer - T & I exported from different HALO versions)
    tableIp<- subset(tableIp, select = -c(2:3,9:12))
    #rename XY columns, HLA-ABC (& title - different HALO versions)
    colnames(tableIp)<-str_replace(colnames(tableIp),pattern="Immune Registered XMin",replacement="XMin")
    colnames(tableIp)<-str_replace(colnames(tableIp),pattern="Immune Registered XMax",replacement="XMax")
    colnames(tableIp)<-str_replace(colnames(tableIp),pattern="Immune Registered YMin",replacement="YMin")
    colnames(tableIp)<-str_replace(colnames(tableIp),pattern="Immune Registered YMax",replacement="YMax")
    colnames(tableIp)<-str_replace(colnames(tableIp), pattern="HLA-ABC\\+ Melanoma", replacement="HLAABCpMelanoma")
    colnames(tableIp)<-str_replace(colnames(tableIp), pattern="HLA-ABC- Melanoma", replacement="HLAABCnMelanoma")
    colnames(tableIp)<-str_replace(colnames(tableIp),pattern="Image Location",replacement="Image File Name")
    
  #add xy columns to tableI
  tableI<-tableIp[,1:6]
  #add T cell phenotype columns to tableI
  out2<-matrix(0,nrow=(nrow(tableI)),ncol=10)
  colnames(out2)<-c(T.colsp)
  tableI<-cbind(tableI,out2)
  #Add melanoma column to tableI
  tableI$Melanoma<-tableIp[,"Melanoma"]
  #add T cell marker columns
  out3<-matrix(0,nrow=(nrow(tableI)),ncol=20)
  colnames(out3)<-c(T.colsm)
  tableI<-cbind(tableI,out3)
  #add SOX10,DAPI & cell parameters
  tableI<-cbind(tableI,tableIp[,33:46])
  #add rest of Immune columns to tableI
  tableI<-cbind(tableI, tableIp[,c(7:9, 11:32)])
  #Add object ID number to immune cells
  tableI$`Object Id`<-tableI$`Object Id` + objectidn
  
  #Replace (um2)  
  colnames(tableI)<-str_replace_all(colnames(tableI),pattern="\\(µm²\\)",replacement="")
  colnames(tableI)<-str_replace_all(colnames(tableI),pattern="\\(µm\\)",replacement="")
  
  combitable<-rbind(tableT,tableI)
  fwrite(combitable, file = newdat,quote=F,sep=",")
    }
  }

