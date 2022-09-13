library(devtools)
library(fitdistrplus)
library(depmixS4)
# install_github("rr1859/R.4Cker")
# library(R.4Cker)
# proxy won't allow install -> use dump/source
# all files associated with program in Documents>R>R-3.6.0>R.4Cker-master>R

# only need to dump one time, re-source if change .R file
dump("createR4CkerObject", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/createR4CkerObject.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/createR4CkerObject.R")

enz_file=read.table("T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/hg19_csp6i_flanking_sites_30_unique_2.bed", stringsAsFactors=FALSE)
# had to set working directory otherwise, Error in file(file, "rt") : cannot open the connection
# setwd("T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C")

# this command takes awhile to run
my_obj_1=createR4CkerObjectFromFiles(files=c("T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/sample9_1_aligned.bedGraph"),bait_chr="chr14",bait_coord=23398794,bait_name="PRMT5",primary_enz="GTAC",samples=c("PRMT5_9"),conditions="PRMT5",replicates=1,species="hs",output_dir="T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/4Cker",enz_file=enz_file)

dump("transAnalysis", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/transAnalysis.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/transAnalysis.R")

dump("buildAdaptiveWindowsTrans", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/buildAdaptiveWindowsTrans.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/buildAdaptiveWindowsTrans.R")

dump("getWindowCounts", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/getWindowCounts.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/getWindowCounts.R")

dump("removePCR", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/removePCR.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/removePCR.R")

dump("startingValues", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/startingValues.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/startingValues.R")

dump("parameterEstimation", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/parameterEstimation.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/parameterEstimation.R")

dump("differentialAnalysis", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/differentialAnalysis.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/differentialAnalysis.R")

dump("merge_windows", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/merge_windows.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/merge_windows.R")

dump("normalizeCounts", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/normalizeCounts.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/normalizeCounts.R")

dump("validateParameters", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/validateParameters.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/validateParameters.R")

dump("viterbi3State", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/viterbi3State.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/viterbi3State.R")

dump("generateSyntheticSamples", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/generateSyntheticSamples.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/generateSyntheticSamples.R")

dump("buildAdaptiveWindowsCis", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/buildAdaptiveWindowsCis.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/buildAdaptiveWindowsCis.R")

dump("cisAnalysis", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/cisAnalysis.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/cisAnalysis.R")

dump("nearBaitAnalysis", file="C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/nearBaitAnalysis.R")
source("C:/Users/brin43/Documents/R/R-3.6.0/R.4Cker-master/R/nearBaitAnalysis.R")

transAnalysis(my_obj_1,k=50)
#######################
#format my data following transAnalysis
trans_data_circos_high <- read.table("T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/4Cker/PRMT5_9_trans_highinter.bed")
trans_data_circos_low <- read.table("T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/4Cker/PRMT5_9_trans_lowinter.bed")
names(trans_data_circos_high)[1]<-"chromosome.1"
names(trans_data_circos_high)[2]<-"chromStart.1"
names(trans_data_circos_high)[3]<-"chromEnd.1"
names(trans_data_circos_low)[1]<-"chromosome.1"
names(trans_data_circos_low)[2]<-"chromStart.1"
names(trans_data_circos_low)[3]<-"chromEnd.1"
trans_data_circos_high$chromosome<-"chr14"
trans_data_circos_high$chromStart<-23398790
trans_data_circos_high$chromEnd<-23398797
trans_data_circos_low$chromosome<-"chr14"
trans_data_circos_low$chromStart<-23398790
trans_data_circos_low$chromEnd<-23398797
trans_data_circos_high_rev<-trans_data_circos_high[,c("chromosome.1","chromStart.1","chromEnd.1","chromosome","chromStart","chromEnd")]
#trans_data_circos_3<-trans_data_circos_2[order(trans_data_circos_2[,1],as.numeric(trans_data_circos_2[,2])),];
trans_data_circos_low_rev<-trans_data_circos_low[,c("chromosome.1","chromStart.1","chromEnd.1","chromosome","chromStart","chromEnd")]
#trans_data_circos_3<-trans_data_circos_2[order(trans_data_circos_2[,1],as.numeric(trans_data_circos_2[,2])),];
#######################
# make circos plot- high interactions
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude<-NULL;
cyto.info<-UCSC.HG19.Human.CytoBandIdeogram;
#tracks.inside<-3;
tracks.inside<-2;
tracks.outside<-0;
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside);
out.file<-"T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/CircosPlot_high.pdf";
pdf(file=out.file,height=8,width=8,compress=TRUE);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data);
name.col<-4;
side<-"in";
track.num<-1;
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num,side);
track.num<-2;
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num,side);
track.num<-3;
RCircos.Link.Plot(trans_data_circos_high,track.num,by.chromosome=TRUE);
dev.off()
#######################
# remove problematic chr22 interactions
trans_data_circos_low_2 <-rbind(trans_data_circos_low[1:62,],trans_data_circos_low[64:lengths(trans_data_circos_low[1]),])
#######################
# make circos plot - low interactions
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude<-NULL;
cyto.info<-UCSC.HG19.Human.CytoBandIdeogram;
#tracks.inside<-3;
tracks.inside<-2;
tracks.outside<-0;
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside);
out.file<-"T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/CircosPlot_low.pdf";
pdf(file=out.file,height=8,width=8,compress=TRUE);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data);
name.col<-4;
side<-"in";
track.num<-1;
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num,side);
track.num<-2;
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num,side);
track.num<-3;
RCircos.Link.Plot(trans_data_circos_low_2,track.num,by.chromosome=TRUE);
dev.off()
#######################
#######################
#format my data from my_obj
easy<-my_obj_1@data_trans
trans_data_circos<-easy[[1]]
trans_data_circos<-trans_data_circos[,1:3]
# remove some rows without coordiantes
trans_data_circos <-rbind(trans_data_circos[1:1298,],trans_data_circos[1303:lengths(trans_data_circos[1]),])
# remove chrM and chrUn
trans_data_circos <-rbind(trans_data_circos[1:3275,],trans_data_circos[3293:lengths(trans_data_circos[1]),])
#attempt to remove rows from chr1 that are "outside of chromosome length"
names(trans_data_circos)[1]<-"chromosome.1"
names(trans_data_circos)[2]<-"chromStart.1"
names(trans_data_circos)[3]<-"chromEnd.1"
trans_data_circos$chromosome<-"chr14"
trans_data_circos$chromStart<-23398790
trans_data_circos$chromEnd<-23398797
trans_data_circos_2<-trans_data_circos[,c("chromosome.1","chromStart.1","chromEnd.1","chromosome","chromStart","chromEnd")]
trans_data_circos_3<-trans_data_circos_2[order(trans_data_circos_2[,1],as.numeric(trans_data_circos_2[,2])),];
#remove rows that caused error for length of chr1
trans_data_circos_3<-trans_data_circos_3[123:lengths(trans_data_circos_3[1]),]
#remove rows that caused error for length of chr10
trans_data_circos_4 <-rbind(trans_data_circos_3[1:147,],trans_data_circos_3[254:lengths(trans_data_circos_3[1]),])
#remove rows that caused error for length of chr11
trans_data_circos_5 <-rbind(trans_data_circos_4[1:184,],trans_data_circos_4[305:lengths(trans_data_circos_4[1]),])
#remove rows that caused error for length of chr12
trans_data_circos_6 <-rbind(trans_data_circos_5[1:240,],trans_data_circos_5[394:lengths(trans_data_circos_5[1]),])
#remove rows that caused error for length of chr13
trans_data_circos_7 <-rbind(trans_data_circos_6[1:279,],trans_data_circos_6[361:lengths(trans_data_circos_6[1]),])
#remove rows that caused error for length of chr15
trans_data_circos_8 <-rbind(trans_data_circos_7[1:287,],trans_data_circos_7[401:lengths(trans_data_circos_7[1]),])
#remove rows that caused error for length of chr16
trans_data_circos_9 <-rbind(trans_data_circos_8[1:290,],trans_data_circos_8[313:lengths(trans_data_circos_8[1]),])
#remove rows that caused error for length of chr17
trans_data_circos_10 <-rbind(trans_data_circos_9[1:419,],trans_data_circos_9[446:lengths(trans_data_circos_9[1]),])
#remove rows that caused error for length of chr18
trans_data_circos_11 <-rbind(trans_data_circos_10[1:555,],trans_data_circos_10[571:lengths(trans_data_circos_10[1]),])
#remove rows that caused error for length of the rest of the chromosomes
trans_data_circos_12 <-rbind(trans_data_circos_11[1:631,],
                             trans_data_circos_11[649:742,],
                             trans_data_circos_11[847:981,],
                             trans_data_circos_11[1004:1093,],
                             trans_data_circos_11[1097:1238,],
                             trans_data_circos_11[1328:1439,],
                             trans_data_circos_11[1549:1639,],
                             trans_data_circos_11[1753:1845,],
                             trans_data_circos_11[1973:2032,],
                             trans_data_circos_11[2133:2185,],
                             trans_data_circos_11[2311:2349,],
                             trans_data_circos_11[2454:2517,],
                             trans_data_circos_11[2639:2685,],
                             trans_data_circos_11[2688:lengths(trans_data_circos_11[1]),])
########################
#subset of data for testing
#trans_data_circos_3<-trans_data_circos_2[1:168,]
#trans_data_circos_3 <-rbind(trans_data_circos[1:167,],trans_data_circos[179,])
#trans_data_circos_3 <-trans_data_circos[168:172,] #error
#trans_data_circos_4<-trans_data_circos_4[order(as.numeric(trans_data_circos_3[,2])),];
#trans_data_circos_4<-trans_data_circos_3[123:269,]
#trans_data_circos_5<-trans_data_circos_4[304:360,]
trans_data_circos_9<-trans_data_circos_8[312:441,]
########################
# make circos plot
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude<-NULL;
cyto.info<-UCSC.HG19.Human.CytoBandIdeogram;
#tracks.inside<-3;
tracks.inside<-2;
tracks.outside<-0;
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside);
out.file<-"T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/CircosPlot7.pdf";
pdf(file=out.file,height=8,width=8,compress=TRUE);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data);
name.col<-4;
side<-"in";
track.num<-1;
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num,side);
track.num<-2;
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num,side);
track.num<-3;
RCircos.Link.Plot(trans_data_circos_12,track.num,by.chromosome=TRUE);
dev.off()
