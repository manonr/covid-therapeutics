


################################
# plotting for ms
##################################




runDate <- "April12"
ndays <- 10                 ########## how many days after treatment do you count sequences as post treatment?
pcutoff <- 0.001

subdir <- paste0(ndays,"dayCutoff")

rootDir <- paste0("C:/Users/mr909/OneDrive/Projects/ncov/PHE/therapeutics/", runDate)



rootDir <- paste0("C:/Users/manon.ragonnet/Documents/Projects/therapeutics/", runDate)

wDir <- paste0(rootDir, "/",subdir)
setwd(wDir)

lineages <- read.csv("../therapeutics_with_seq_lineages.csv", stringsAsFactors = F)
lineages <- rbind(lineages[lineages$prepost=="post" & lineages$date_difference>ndays,], lineages[lineages$prepost=="pre",])


# functions
DRPlot <- function(output, genelength, myPlotTitle=NULL){
  variant <- unique(output$variant)
  drug <- unique(output$treatment)
  gene <- unique(output$gene)
  
  # get counts
  npre <- length(unique(lineages$uniq_ID[lineages$intervention==drug &lineages$prepost=="pre" & lineages$variant==variant]))
  npost <- length(unique(lineages$uniq_ID[lineages$intervention==drug &lineages$prepost=="post" & lineages$variant==variant]))
  total <- npre+npost
  
  output$p2[output$p2<0.0000000001] <- 0.0000000001
  output$logp <- -log(as.numeric(output$p2),base = 10)
  
  output2 <- unique(output[,c("pos", "logp")])
  
  # colours https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=11
  
  lowp <- output[output$logp>-log(as.numeric(pcutoff),base = 10),]
  
  if(!is.null(myPlotTitle)){
    plotTitle <- paste0(myPlotTitle)#, "\n","total ", total, 
                       # " (", npre, " pre/ ", npost, " post)") #drug, "\n",variant, "\n", gene , "\n
  } else {
    plotTitle <- paste("total ", total, 
                       " (", npre, " pre/ ", npost, " post)") #drug, "\n",variant, "\n", gene , "\n
  }
  plot(output2$pos, output2$logp, pch=16,ylim=c(0, max(output2$logp)+2), xlim=c(0,genelength+200),
       xlab=paste("Residue in",gene), ylab="-10 log(p)", main=plotTitle)
  points(lowp$pos, lowp$logp, col="#E40046", pch=16)
  
  abline(h = -log(as.numeric(0.001),base = 10), col="grey")
  text(genelength+100, -log(as.numeric(0.001),base = 10)+0.1, labels="p<0.001")
  abline(h = -log(as.numeric(0.0001),base = 10), col="grey")
  text(genelength+100, -log(as.numeric(0.0001),base = 10)+0.1, labels="p<0.0001")
  abline(h = -log(as.numeric(0.00001),base = 10), col="grey")
  text(genelength+100, -log(as.numeric(0.00001),base = 10)+0.1, labels="p<0.00001")
  
  
  
}


annotatePos <- function(output, 
                        textgapx= rep(0, length= unique(output$pos[output$p2<0.01])),
                        textgapy= rep(0.2, length= unique(output$pos[output$p2<0.01])))
{
  output$p2[output$p2<0.0000000001] <- 0.0000000001
  output$logp <- -log(as.numeric(output$p2),base = 10)
  lowp <- output[output$logp>-log(as.numeric(pcutoff),base = 10),]
  tg <- 0
  
  for (lowpos in unique(lowp$pos)){
    tg <- tg+1
    posdf <- lowp[lowp$pos==lowpos,]
    aa1 <- unique(substr(posdf$aminoacid,0,1))   # edit this
    tmp <- posdf[posdf$npost_uniquePatient>0,]
    aa2a <- unique(substr(tmp$aminoacid, nchar(tmp$aminoacid[1]), nchar(tmp$aminoacid)+1))
    aa2b <- paste(aa2a[-which(aa2a==aa1)], collapse="/")
    aatxt <- paste0(aa1,lowpos,aa2b)
    text(as.numeric(lowpos)+textgapx[tg], posdf$logp[1]+textgapy[tg], aatxt, col="#E40046",cex = 0.8)
  }
}

csv <- read.csv("../../drug_interaction_sites.csv", stringsAsFactors = F)



source("C:/Users/manon.ragonnet/Documents/Rcode/figLabel.R")

pdf("Fig1_3plots.pdf", height=8, width=4)
par(mfrow=c(3,1))
par(mar=c(5.1,4.1,4.1,2.1))

## delta/cas, spike
gene <- "S"; drug <- "Casirivimab and imdevimab"; variant <- "delta"
outputs <- dir()[grep(paste(gene, drug,variant, sep="_"), dir())]
output <- read.csv(outputs[1]) 


DRPlot(output, 1300, myPlotTitle = "Delta variant, treated with casirivimab and imdevimab")
fig_label("A", region="figure", pos="topleft", cex=2)
points(csv$site[csv$drug=="casirivimab"], rep(11, length(csv$site[csv$drug=="casirivimab"])),  col="#1D57A5", pch=17)
points(csv$site[csv$drug=="imdevimab"], rep(11.8, length(csv$site[csv$drug=="imdevimab"])),  col="#8A1B61", pch=17)
text(1400, 11, labels="casirivimab",  col="#1D57A5")
text(1400,11.8, labels="imdevimab",  col="#8A1B61")
unique(output$pos[output$p2<pcutoff])
annotatePos(output, textgapx=c(-50,-20,60,70),
            textgapy=c(0.5,-0.5, 0.5,0.5))


## BA.1, spike
gene <- "S"; drug <- "Sotrovimab"; variant <- "BA.1"
outputs <- dir()[grep(paste(gene, drug, variant, sep="_"), dir())]
output <- read.csv(outputs) 
output$aminoacid <- gsub("Q493Q", "R493Q", output$aminoacid)
output$aminoacid <- gsub("Q493R", "R493R", output$aminoacid)
output$aminoacid <- gsub("Q493L", "R493L", output$aminoacid)
DRPlot(output, 1300, myPlotTitle = "BA.1 variant, treated with Sotrovimab")
fig_label("B", region="figure", pos="topleft", cex=2)
points(csv$site[csv$drug=="Sotrovimab"], rep(11, length(csv$site[csv$drug=="Sotrovimab"])),  col="#582C83", pch=17)
text(1400, 11, labels="sotrovimab",  col="#582C83")
unique(output$pos[output$p2<pcutoff])
annotatePos(output, textgapx=c(-90,-110,75,45),
            textgapy=c(-0.5,0.5, 0.4, 0.5,0.5))


## BA.2, spike
gene <- "S"; drug <- "Sotrovimab"; variant <- "BA.2"
outputs <- dir()[grep(paste(gene, drug, variant, sep="_"), dir())]
output <- read.csv(outputs) 
DRPlot(output, 1300, myPlotTitle = "BA.2 variant, treated with Sotrovimab")
fig_label("C", region="figure", pos="topleft", cex=2)
points(csv$site[csv$drug=="Sotrovimab"], rep(12, length(csv$site[csv$drug=="Sotrovimab"])),  col="#582C83", pch=17)
text(1400, 12, labels="sotrovimab",  col="#582C83")

unique(output$pos[output$p2<pcutoff])
annotatePos(output, textgapx=c(0.2),
            textgapy=c(0.2))

dev.off()
