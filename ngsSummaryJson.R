#!/usr/bin/env Rscript
# Author: OLP/01/2022
# Usage: Getting metadata for ngs file
# Parameters: Run Name ;  i.e. Rscript ngs.R UT-A01290-220528
# last modified on 22/09/12

library("easypackages");libraries("grid","rbin","data.table","progress","XML","xml2","seqinr","data.table")
sup <- suppressPackageStartupMessages
sup(library(lubridate))
sup(library(tidyverse))
options("width"=300)
args = commandArgs(trailingOnly=T);# args[1]<-c('UT-A01290-230119')
date<-ymd(substr(args[1], 11, 16));
runPath <-paste('/Volumes/NGS/Analysis/covidseq',args[1],sep="/")
runPath
df1 = read.csv(paste(runPath,"covidseq_output/Logs_Intermediates/SampleSheetValidation/SampleSheet_Intermediate.csv",sep = "/"))
dim(df1)
df1<-tail(df1,-13);df1<-df1[,c(1,3)];
colnames(df1) <- c("Sample_Accession","Sample_Type")
M<-filter(df1, Sample_Type == "PatientSample")
M<-M["Sample_Accession"]
cat(" Detected ", length(M$Sample_Accession), " patient's samples")

# Getting Pangolin lineage 
pango<-read.csv(paste(runPath,"CecretPangolin/pangolin/lineage_report.csv",sep ="/"))[,c("taxon","lineage")]
var <-paste('',args[1],sep='-'); head(var)
pango<-pango %>% mutate_at("taxon", str_replace, var, "")
M<-merge(M,pango, by.x="Sample_Accession", by.y= "taxon", all.x=T)
M$lineage[is.na(M$lineage)]= "Not able to be sequenced"
head(M, n = 120)

# From LIMS
#COVID/daily_metadata/lims/
cat(crayon::bold("Please wait, serching UPHL-LIMS","\n"))
cat(crayon::bold("1. From NGS_Covid","\n"))
NGS_covid<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/COVID/daily_metadata/lims/",
                                 pattern=glob2rx("NGS*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/COVID/daily_metadata/lims/",
                                                                                          pattern=glob2rx("NGS*.csv"))$info$mtime),])
NGS_covid <- fread(paste('/Volumes/IDGenomics_NAS/COVID/daily_metadata/lims/',NGS_covid,sep=""), 
                   select = c("submitterId","sampleNumber","lastName","firstName","DOB","collectionDate")) 
head(NGS_covid); 
NGS_covid1<-select(NGS_covid,c("sampleNumber","lastName","firstName","DOB","collectionDate"))
NGS_covid1$sampleNumber<-as.character(NGS_covid1$sampleNumber); 
NGS_covid2<-select(NGS_covid,c("submitterId","lastName","firstName","DOB","collectionDate"))
colnames(NGS_covid2)[1]<-("sampleNumber")
NGS_covid<-bind_rows(NGS_covid1,NGS_covid2)
head(NGS_covid)
NGS_covid$DOB<-as.character(NGS_covid$DOB)
NGS_covid$collectionDate<-as.character(NGS_covid$collectionDate)

cat(crayon::bold("2. From All_nocovid_NoNGS","\n"))
All_nocovid<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/COVID/daily_metadata/lims/",
                                   pattern=glob2rx("All*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/COVID/daily_metadata/lims/", 
                                                                                            pattern=glob2rx("All*.csv"))$info$mtime),])
head(All_nocovid)
All_nocovid <- fread(paste('/Volumes/IDGenomics_NAS/COVID/daily_metadata/lims/',All_nocovid,sep=""), 
                     select = c("sampleNumber","lastName","firstName","DOB","collectionDate")) 

All_nocovid$sampleNumber<-as.character(All_nocovid$sampleNumber)
All_nocovid$sampleNumber<-as.character(All_nocovid$sampleNumber)
LIMS<-bind_rows(NGS_covid,All_nocovid)
head(LIMS)     
Matrix<-merge(M,LIMS, by.x = "Sample_Accession", by.y = "sampleNumber", all.x = T)
Matrix1<-merge(Matrix$Sample_Accession[which(!is.na(Matrix$collectionDate))],Matrix, by.x = "x", by.y ="Sample_Accession", all.x = T)
colnames(Matrix1) <- c("sampleNumber","Test result","lastName","firstName","DOB","collectionDate")
head(Matrix1)
# Searching epitrax

EpiTrax_Export<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/COVID/daily_metadata/epitrax/", 
                                      pattern=glob2rx("export_1551*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/COVID/daily_metadata/epitrax/",
                                                                                                       pattern=glob2rx("export_1551*.csv"))$info$mtime),])
cat(crayon::bold("Please wait, searching EpiTrax","\n")); head(EpiTrax_Export)
M_epitrax <- fread(paste('/Volumes/IDGenomics_NAS/COVID/daily_metadata/epitrax/',EpiTrax_Export,sep=""), 
                   select = c("lab_accession_no","person_last_name","person_first_name", "patient_birth_date","lab_collection_date"))   

Mx<-merge(Matrix$Sample_Accession[which(is.na(Matrix$collectionDate))],M, by.x ="x", by.y ="Sample_Accession", all.x = T)
Matrix2<-merge(Mx,M_epitrax, by.x = "x", by.y ="lab_accession_no", all.x = T)

colnames(Matrix2) <- c("sampleNumber","Test result","lastName","firstName","DOB","collectionDate") 
Matrix2$collectionDate<-as.Date(Matrix2$collectionDate)
Matrix2$DOB<-as.Date(Matrix2$DOB)
Matrix2$collectionDate<-as.character(Matrix2$collectionDate)
Matrix2$DOB<-as.character(Matrix2$DOB)
Matrix2<-unique(Matrix2);   head(Matrix2)

ngs_file<-bind_rows(Matrix1, Matrix2)
ngs_file<-cbind(ngs_file,'Name of facility performing the test'='UPHL', 'ordering facility' ='UHPL',
                'Name or description of the test ordered'="SARS-CoV-2 Whole Genome Sequencing",
                'lab test date'=ymd(date),'lab test status' ="preliminary")
ngs_file<-select(ngs_file,"lastName",	"firstName",	"DOB",	"Name of facility performing the test","Name or description of the test ordered",
                 "collectionDate","sampleNumber","lab test date","lab test status",	"Test result",	"ordering facility")
ngs_file[is.na(ngs_file)] <- ""
head(ngs_file)

# Storing results
# Aqui va le marico JSON file for daily_metadata uumm wait for summary

d=paste(runPath,'Rsumcovidseq',sep="/")
dir.create(file.path(d,"" ), recursive = T)
LL<-paste(substring(Sys.time(), 1,10),substring(Sys.time(), 12,13),substring(Sys.time(), 15,16),substring(Sys.time(), 18,19),sep = "-")
dx<- paste(paste('ngs',args[1],sep="_"),LL,".csv",sep="")
dx<-paste(d,dx, sep ='/');write.csv(ngs_file,dx,row.names=F)
cat(crayon::bold("ngs file is complete and can be found at",d,"\n"))



# Summary
df1 = read.csv(paste(runPath,"covidseq_output/Logs_Intermediates/SampleSheetValidation/SampleSheet_Intermediate.csv",sep = "/"))
dim(df1)
head(df1,40)
df2 = read.csv(paste(runPath,"CecretPangolin/pangolin/lineage_report.csv",sep ="/"))
df3 = read.table(file = paste(runPath,"covidseq_output/Logs_Intermediates/VirusDetection/kmer/Summary.tsv",sep ="/"), sep = '\t', header = T); 
df4 = read.csv(paste(runPath,"covidseq_output/Logs_Intermediates/FastqGeneration/Reports/Demultiplex_Stats.csv",sep ="/"))
df4<-select(df4, c("SampleID","X..Reads","Mean.Quality.Score..PF."))
var <-paste('',args[1],sep='-'); head(var)
df4<-df4 %>% mutate_at("SampleID", str_replace, var, ""); df4<-df4 %>% rename(Num_reads = 2, Mean_Quality_Score = 3)
dim(df4)
#=============================================
# Info parameters about the run:
Sequencer <-substr(args[1], 4, 4); dat<-substr(args[1], 11, 16);
if (Sequencer =="A"){ TypeOfSequencer<-"NovaSeq";Pt<-"A01290"} else {TypeOfSequencer<-"NextSeq";Pt<-"NB551133"};
instrument_type<-df1[4,2];chemistry<-df1[7,2];index_adapters<-df1[6,2];wetlab_assay<-df1[5,2];pangolin_version<-df2[1,10]
software = "Illumina DRAGEN COVID pipeline (R.U.O.) version 1.0.1";# <===Where can it be found
pathOutput <-paste("/Volumes/NGS/Output/",Pt,"/",dat,"*/RunParameters.xml",sep='')

fileName<-Sys.glob(pathOutput,dirmark = TRUE);
ReadType<-xmlToDataFrame(nodes=getNodeSet(xmlParse(read_xml(fileName)),"//ReadType"))
ReadLength<-xmlToDataFrame(nodes=getNodeSet(xmlParse(read_xml(fileName)),"//Read1NumberOfCycles"))
df1<-tail(df1,-13);
head(df1)
df1<-df1[,c(1,3,11)];
colnames(df1) <- c("Sample_Accession","Sample_Type","lanes");
Nlanes<-length(unique(df1[3])$lanes)
df1<-filter(df1, Sample_Type == "PatientSample"); n<-length(df1$Sample_Type)
write.csv("")
n
df1<-merge(df1,df4, by.x = "Sample_Accession", by.y = "SampleID", all.x = T)

#=============================================
# Run parameters
RunPar <- data.frame("","");colnames(RunPar)<-c("Parameters","Info");
RunPar[1,1]<-"Wetlab assay";RunPar[2,1]<-"Instrument type";RunPar[3,1]<-"Software";RunPar[4,1]<-"Pangolin version";RunPar[5,1]<-"Read type";
RunPar[6,1]<-"Read length";RunPar[6,1]<-"Number of samples";RunPar[7,1]<-"Number of lanes";RunPar[8,1]<-"Chemistry";RunPar[9,1]<-"index adapters";
RunPar[10,1]<-"Date of the run"
RunPar[1,2]<-wetlab_assay;RunPar[2,2]<-instrument_type;RunPar[3,2]<-software;RunPar[4,2]<-pangolin_version;
RunPar[5,2]<-ReadType;RunPar[6,2]<-n;RunPar[7,2]<-Nlanes;RunPar[8,2]<-chemistry;RunPar[9,2]<-index_adapters;RunPar[10,2]<-as.character(date)
RunPar
#=============================================
# Random Number Generator for Sample_IDs
RNG<-paste(floor(runif(n, min=99, max=1000)),floor(runif(n, min=99, max=1000)),sep ="")
# Checking for duplicates in RNG
while(length(RNG)!= length(unique(RNG)))
{ RNG<-paste(floor(runif(n, min=99, max=1000)),floor(runif(n, min=99, max=1000)),sep ="");}
IDD<-paste("UT-UPHL",substr(args[1], 11, 16),sep="-")
Sample_ID<-paste(IDD,RNG,sep="")
summary1<-cbind(df1,Sample_ID)
dim(summary1)
#============================================
#join tables: SampleSheet & Combined Lineage Report
var <-paste('',args[1],sep='-');
df2 <-df2 %>% mutate_at("taxon", str_replace, var, "")
df2<-df2[,c("taxon","lineage","scorpio_call")]
df3<-df3 %>% mutate_at("Sample", str_replace, var, "")
df3<-df3[,c("Sample","nTargetsDetected.SARS.CoV2")]
colnames(df3)[2] <- 'sc2_amplicons';#
summary1<-merge(summary1,df2, by.x="Sample_Accession", by.y= "taxon", all.x=TRUE)
summary1<-merge(summary1,df3, by.x="Sample_Accession", by.y= "Sample", all.x=TRUE)

#==========================================
# Nucleotides count: (actg & n)
num_actgx<-vector();num_nx<-vector(); pass_fail<-vector()
runPath1 = paste(runPath,"covidseq_output/Sample_Analysis/",sep ="/")
num_actgx = 0; num_nx = 0;
n
for (i in 1:n){
  sample=summary1$Sample_Accession[i]
  amplicons=summary1$sc2_amplicons[i]
  if(amplicons >=50)
  {
    Lyapunovv<- Sys.glob(paste(runPath1,sample,"*/*.fasta",sep = ""),dirmark = TRUE)
    a =  str_count(read.fasta(Lyapunovv, as.string = T, forceDNAtolower = T), pattern = "a|c|t|g");num_actgx[i] = a
    b =  str_count(read.fasta(Lyapunovv, as.string = T, forceDNAtolower = T))-a; num_nx[i] = b; pass_fail[i] ="pass"; c = "pass"
  }
  else {# AmpliconThereshold <50,
    a = b =0; num_actgx[i] = 0; num_nx[i] = 0; pass_fail[i] ="fail"; c = "fail"
  }
  cat(sample," \t",a,"\t",b,"\t",c,"\n")
}


df_new <- as.data.frame(cbind(as.character(num_actgx), as.character(num_nx),pass_fail)); 
df_new<-df_new %>% rename(num_actg = 1, num_n = 2)
head(df_new)
dim(summary1)
dim(df_new)
summary<-cbind(summary1,df_new);

head(summary)
dim(summary)

summary<-select(summary,"Sample_Accession","Sample_ID","Sample_Type","lanes","Num_reads",
                "Mean_Quality_Score","sc2_amplicons","lineage","scorpio_call","num_actg","num_n","pass_fail" )

#========================================
# storing files
d=paste(runPath,'Rsumcovidseq',sep="/");dir.create(file.path(d,"" ), recursive = TRUE); p1<-Sys.time()
dx<- paste(d,'/Rcovidseq_summary-',paste0(Sys.Date(),'-',hour(Sys.time()),'-',minute(Sys.time()),'-',second(Sys.time())),".csv",sep=""); 
write.csv(summary,dx,row.names=FALSE)
dx<- paste(d,'/RunParameters-',paste0(Sys.Date(),'-',hour(Sys.time()),'-',minute(Sys.time()),'-',second(Sys.time())),".csv",sep=""); 
write.csv(RunPar,dx,row.names=FALSE)
cat("Summary & run rarameters files are complete and can be found at: ",d,"\n")

#=======================================
# calling the python script for the UPHL_metadata file (json file)
# ======================================
path2python<-paste("python3 /home/Bioinformatics/Dripping_Rock/bin/fasta_to_json.py",args[1],sep = " ")

cat(path2python,"\n")
#system('python3 /home/Bioinformatics/Dripping_Rock/bin/fasta_to_json.py UT-A01290-230118')
#cat("system('python3 /home/Bioinformatics/Dripping_Rock/bin/fasta_to_json.py UT-A01290-230118'","\n")
system(path2python)




