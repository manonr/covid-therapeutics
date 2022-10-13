

# main analysis
# this script takes all the translated protein sequences
# for each combo of interest, the aa frequencies are compared across treated an untrested pateints
# prints one csv per combo examined
# and one summary csv with all significant combinations

# M. Ragonnet
# 09/03/2022


library(ape)
library(phangorn)
library(bioseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)



runDate <- "April12"
ndays <- 5                 ########## how many days after treatment do you count sequences as post treatment?

subdir <- paste0(ndays,"dayCutoff")

refDir <- "C:/Users/mr909/OneDrive/Projects/ncov/PHE/therapeutics/reference"
refDir <- "C:/Users/manon.ragonnet/Documents/Projects/therapeutics/reference"

rootDir <- paste0("C:/Users/mr909/OneDrive/Projects/ncov/PHE/therapeutics/", runDate)
rootDir <- paste0("C:/Users/manon.ragonnet/Documents/Projects/therapeutics/", runDate)
setwd(rootDir)

wDir <- paste0(rootDir, "/",subdir)
dir.create(subdir)
setwd(wDir)

drugCombos <- read.csv("../../drug_gene_combos.csv", stringsAsFactors = F)

lineages <- read.csv("../therapeutics_with_seq_lineages.csv", stringsAsFactors = F)
lineages <- rbind(lineages[lineages$prepost=="post" & lineages$date_difference>ndays,], lineages[lineages$prepost=="pre",])

aa_files <- dir(rootDir)[grep("aa", dir(rootDir))]
variants <- unique(lineages$variant)

summary_output <- data.frame(gene="spike",variant="variant", treatment="treatment", pos=1, aminoacid="A", npost_uniqueseq=0, 
                             npre_uniqueseq=0, proppost=0, proppre=0, p=1,npost_uniquePatient=0, npre_uniquePatient=0, p2=1, stringsAsFactors = FALSE)
summary_line <- 0

for (z in 1:length(drugCombos[,1])){
  
  genename <- drugCombos$gene[z]
  drug <- drugCombos$drug[z]
  
  aafile <- aa_files[grep(paste0(genename, "\\."), aa_files)]
  
  aa_seq <- bioseq::read_fasta(paste0(rootDir, "/",aafile), type="AA")
  
  names(aa_seq) <- gsub("\r", "", names(aa_seq))
  aa_seq_uniqueID <- aa_seq[match(unique(names(aa_seq)),names(aa_seq))]
  
  
  
  for (var in variants){
    
    refgenome <- bioseq::read_fasta(paste(refDir,var,aafile, sep="/"), type="AA") 
    ref_tibble <- tibble(label = names( refgenome), sequence =  refgenome )
    
    print(paste(drug, genename, var))
    
    output <- data.frame(gene="spike",variant="variant", treatment="treatment", pos=1, aminoacid="A", npost_uniqueseq=0, 
                         npre_uniqueseq=0, proppost=0, proppre=0, p=1,npost_uniquePatient=0, npre_uniquePatient=0, p2=1, stringsAsFactors = FALSE)
    line <- 0
    
    patients_var <- unique(lineages$central_sample_id[lineages$intervention==drug & lineages$variant==var])
    drug_aa_seq <- aa_seq_uniqueID[unique(unlist(lapply(patients_var, function(x) {grep(x, names(aa_seq_uniqueID))})))]
    
    fra_data <- tibble(label = names( drug_aa_seq ), sequence =  drug_aa_seq )
    fra_data$treatment <- "unknown"
    
    fra_data$treatment[match(lineages$central_sample_id[lineages$prepost=="post"], fra_data$label)] <- "post"
    fra_data$treatment[match(lineages$central_sample_id[lineages$prepost=="pre"], fra_data$label)] <- "pre"
    
    fra_data$uniqueID <- lineages$uniq_ID[match(fra_data$label, lineages$central_sample_id)]
    
    print(table( fra_data$uniqueID,fra_data$treatment))
    
    if(nrow(fra_data)>1){
      
      for (i in 1:as.numeric(nchar(fra_data$sequence[1]))){
        
        # print(i)
        
        tab <- t(table(fra_data$treatment, unlist(lapply(fra_data$sequence, function(x) {strsplit(x, "")[[1]][i]}))))
        
        # remove rows that are X or ~
        zig <- which(rownames(tab)=="~")
        if(length(zig)>0){
          tab <- tab[-zig,]
        }
        Xs <- which(rownames(tab)=="X")
        if(length(Xs)>0){
          tab <- tab[-Xs,]
        }
        stars <- which(rownames(tab)=="*")
        if(length(stars)>0){
          tab <- tab[-stars,]
        }
        
        #print( tab)
        
        if (length(tab)>2 ){
          if(sum(tab[,1])>0){
            
          
            
            fish <- fisher.test(tab,simulate.p.value=TRUE)
            temptab <- cbind(prop.table(tab), Total = rowSums(prop.table(tab)))
            #maxAA <- names(which(temptab[,3]==max(temptab[,3])))
            
            refAA <- unlist(lapply(ref_tibble$sequence, function(x) {strsplit(x, "")[[1]][i]}))
            
            
            
            if(fish$p.value<1){ ############### you can change p value here
              #print(i)
              
              
              dftemp <- data.frame(uniqueID=fra_data$uniqueID, treat=fra_data$treatment,residue=unlist(lapply(fra_data$sequence, function(x) {strsplit(x, "")[[1]][i]})))
              
              for (j in 1:length(tab[,1])){
                line <- line+1
                AAchange <- paste0(refAA, i, rownames(tab)[j])
                
                # count the number of unique patients with each mutation (as well as numer of sequences)
                unique_pre <- length(unique(dftemp$uniqueID[dftemp$residue==rownames(tab)[j] & dftemp$treat=="pre"]))
                unique_post <- length(unique(dftemp$uniqueID[dftemp$residue==rownames(tab)[j] & dftemp$treat=="post"]))
                
                outputLine <- c(genename,var, drug, i, AAchange,tab[j,1],tab[j,2],
                                round(tab[j,1]/sum(tab[,1]),4), round(tab[j,2]/sum(tab[,2]),4),
                                round(fish$p.value,6), unique_post, unique_pre, "notyet")
                output[line,] <- outputLine
                
              }
              
              lineEnd <- line
              lineStart <- lineEnd-length(tab[,1])+1
              
              tab2 <- matrix(as.numeric(c(output$npost_uniquePatient[lineStart:lineEnd], output$npre_uniquePatient[lineStart:lineEnd])),
                             nrow = length(tab[,1]), ncol=2,byrow = F)
              fish2 <- fisher.test(tab2)
              output$p2[lineStart:lineEnd] <- round(fish2$p.value,6)
              
              if(fish2$p.value<0.01){ ############### you can change p value here
                
                summary_output <- rbind(summary_output, output[lineStart:lineEnd,])
              }
            }
          }
        }
      }
      write.csv(output, paste0(genename,"_", drug,"_",var, "_AA_changes_pre-post-treatment_", runDate,subdir, ".csv"), row.names = FALSE)
    }
  }
}

summary_output <- summary_output[2:length(summary_output[,1]),]

write.csv(summary_output, paste0("significant_AA_changes_pre-post-treatment_", runDate,"_",subdir,".csv"), row.names = FALSE)  
