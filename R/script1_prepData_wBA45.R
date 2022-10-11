
# therapeutics prep data
# takes the csv with treatment info and links to lineages
# generates list of cogids to grab from climb (then use get_climb_fasta)
# creates a csv of each id + treatment info + lineage etc


# M. Ragonnet
# 02/03/2022

# files needed for this script
### 2 csv with bluteq ids pre and post treatment
### CLIMB file /cephfs/covid/bham/results/msa/latest/alignments/cog_202X-XX-XX_all_metadata.csv with lineage assignments

runDate <- "April12" # the run date will be changed every time the script run


library(data.table)
library(gtools)



data_folder <- "C:/Users/manon.ragonnet/Documents/Projects/latest_data/"

setwd (paste0("C:/Users/manon.ragonnet/Documents/Projects/therapeutics/", runDate))


csv_pre <- read.csv(dir()[grep("pretreat", dir())])
csv_post <- read.csv(dir()[grep("posttreat", dir())])
csv_pre$prepost <- "pre"
csv_post$prepost <- "post"
csv_temp <- smartbind(csv_pre, csv_post)

therapeutics_seq <- unique(csv_temp$cog_uk_id)


lineages_csv <- paste0(data_folder, dir(data_folder)[grep("all_metadata", dir(data_folder))]) # get this file /cephfs/covid/bham/results/msa/latest/alignments/cog_202X-XX-XX_all_metadata.csv
lineages <- fread(lineages_csv)

lineages2 <- lineages[match(therapeutics_seq, lineages$central_sample_id)[!is.na(match(therapeutics_seq, lineages$central_sample_id))],]
myCols <- c("central_sample_id",  "lineage")
lineages <- lineages2[,..myCols]
lineages$variant <- "other"
lineages$seqExists <- T

# omicrons
lineages$variant[grep("BA.1", lineages$lineage)] <- "BA.1"
lineages$variant[grep("BA.2", lineages$lineage)]<- "BA.2"
lineages$variant[grep("BA.4", lineages$lineage)]<- "BA.4"
lineages$variant[grep("BA.5", lineages$lineage)]<- "BA.5"
#lineages$variant[grep("B.1.1.529", lineages$lineage)]<- "omicron"

# deltas
lineages$variant[grep("B.1.617.2", lineages$lineage)]<- "delta"
lineages$variant[grep("AY", lineages$lineage)]<- "delta"

table(lineages$lineage, lineages$variant,useNA = "always")


output <- merge(lineages, csv_temp, by.x="central_sample_id",
                by.y="cog_uk_id", all=T)
output <- output[!is.na(output$seqExists),]; output$seqExists <- NULL

write.csv(output,"therapeutics_with_seq_lineages.csv", row.names = F)
write.table(therapeutics_seq, "therapeutics_seqNames.txt", row.names = F, col.names = F, quote = F)
