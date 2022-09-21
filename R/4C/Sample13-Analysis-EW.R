T_4C_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/4C/Figs"

library(todor)
library(devtools)
library(fitdistrplus)
#install.packages("depmixS4")
library(depmixS4)
#devtools::install_github("rr1859/R.4Cker")
library(R.4Cker)

enz_file <- read.table("~/network/T/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/hg19_csp6i_flanking_sites_30_unique_2.bed", stringsAsFactors=FALSE)
# had to set working directory otherwise, Error in file(file, "rt") : cannot open the connection
# setwd("T:/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C")

# this command takes awhile to run
dir.create("data")
my_obj_1 = createR4CkerObjectFromFiles(
  files = c(
    "~/network/T/Labs/EHL/Rosa/Lindsey/3_PUBLICATIONS/4C/sample13_1_aligned.bedGraph"
  ),
  bait_chr = "chr14",
  bait_coord = 23398794,
  bait_name = "PRMT5",
  primary_enz = "GTAC",
  samples = c("PRMT5_13"),
  conditions = "PRMT5",
  replicates = 1,
  species = "hs",
  output_dir = "data",
  enz_file = enz_file
)

#Near the bait analysis: For 4bp cutter experiments k=3-5 and for 6bp cutter experiments a k=10 is recommended
#DpnII is a 4bp cutter is a Csp6I 2bp cutter
###Only need to run once
#nb_results=nearBaitAnalysis(my_obj_1,k=5) 

#Trans analysis.  k=20 for 6bp experiments and a k=50+ for 4bp cutter experiments is recommended
###Only need to run once
#transAnalysis(my_obj_1,k=50) #previously 20 for Lindsey's sample9 analysis w/Csp6I & EcoRI

#######################
#format my data following transAnalysis
trans_data_circos_high <- read.table("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/4C/Sample13-1stcutter_DpnII-2ndcutter_Csp6I-Analysis_EW/4C_sample13_1st.DpnII_2nd.Csp6I_rstudio-export/PRMT5_13_trans_highinter.bed")
trans_data_circos_low <- read.table("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/4C/Sample13-1stcutter_DpnII-2ndcutter_Csp6I-Analysis_EW/4C_sample13_1st.DpnII_2nd.Csp6I_rstudio-export/PRMT5_13_trans_lowinter.bed")
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
#######################
#remove chrUn
##unique(trans_data_circos_high$chromosome.1)
trans_data_circos_high<- trans_data_circos_high |> filter(!str_detect(chromosome.1, "chrUn"))

##HighInter Filter: Remove Problematic regions chromosomal regions that dont map to UCSC.HG19.Human.CytoBandIdeogram 
data(UCSC.HG19.Human.CytoBandIdeogram)
UCSC.HG19.Human.CytoBandIdeogram |>
  group_by(Chromosome)|>
  summarise(MaxChrEnd = max(ChromEnd)) #|>filter(Chromosome =="chr3") 
#chr16
trans_data_circos_high2<- trans_data_circos_high[!(trans_data_circos_high$chromosome.1 == "chr16" & trans_data_circos_high$chromEnd.1 > 90354754),] 
##chr17
trans_data_circos_high2<- trans_data_circos_high2[!(trans_data_circos_high2$chromosome.1 == "chr17" & trans_data_circos_high2$chromEnd.1 > 81195210),] 
##chr18
trans_data_circos_high2<- trans_data_circos_high2[!(trans_data_circos_high2$chromosome.1 == "chr18" & trans_data_circos_high2$chromEnd.1 > 78077248),] 
##chr19,20,22,3
trans_data_circos_high2<- trans_data_circos_high2[!(trans_data_circos_high2$chromosome.1 == "chr19" & trans_data_circos_high2$chromEnd.1 > 59128983),] 
trans_data_circos_high2<- trans_data_circos_high2[!(trans_data_circos_high2$chromosome.1 == "chr20" & trans_data_circos_high2$chromEnd.1 > 63025520),] 
trans_data_circos_high2<- trans_data_circos_high2[!(trans_data_circos_high2$chromosome.1 == "chr22" & trans_data_circos_high2$chromEnd.1 > 51304566),] 
trans_data_circos_high2<- trans_data_circos_high2[!(trans_data_circos_high2$chromosome.1 == "chr3" & trans_data_circos_high2$chromEnd.1 > 198022430),] 

# make circos plot- high interactions
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude<-NULL;
cyto.info<-UCSC.HG19.Human.CytoBandIdeogram;
#tracks.inside<-3;
tracks.inside<-2;
tracks.outside<-0;
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside);
out.file <- fs::path(T_4C_Figs, "circos_plot_S13_highInt.pdf");
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
RCircos.Link.Plot(trans_data_circos_high2,track.num,by.chromosome=TRUE);
dev.off()

#Low Interactions
#Remove chrUn
#unique(trans_data_circos_low$chromosome.1)
trans_data_circos_low<- trans_data_circos_low |> filter(!str_detect(chromosome.1, "chrUn"))

# make circos plot - low interactions
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude<-NULL;
cyto.info<-UCSC.HG19.Human.CytoBandIdeogram;
#tracks.inside<-3;
tracks.inside<-2;
tracks.outside<-0;
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside);
out.file<-fs::path(T_4C_Figs, "circos_plot_S13_lowInt.pdf");
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
RCircos.Link.Plot(trans_data_circos_low,track.num,by.chromosome=TRUE);
dev.off()
#######################
