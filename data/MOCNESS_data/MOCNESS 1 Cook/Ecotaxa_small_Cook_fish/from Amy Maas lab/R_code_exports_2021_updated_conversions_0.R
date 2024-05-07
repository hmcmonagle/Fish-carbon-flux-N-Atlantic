
#Updated 26 July 2023 by HG to double check equations and produce figures

#Ready set go


packages <- c("tidyr","dplyr","readxl","stringr","data.table","ggplot2","LaCroixColoR","readxl", 
              "ggpubr", "gridExtra", "RColorBrewer","colorspace","OneR")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)



library(readr)
#ecotaxa_export_5410_20230615_1932 <- read_delim("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/ecotaxa_export_5410_20230615_1932.tsv", 
                                                    #   "\t", escape_double = FALSE, trim_ws = TRUE)

# Helena wrote the following code and commented out the line above, to get files to import correctly on this machine
# (adjust as needed for other file paths on your local machine)
ecotaxa_export_5410_20230615_1932 <- read.csv("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/EXPORTS_NA_fish.csv", header = TRUE)


ecotaxa_export<-ecotaxa_export_5410_20230615_1932
sub1<-ecotaxa_export[,c(1,7,8,15,18,33:34,101,134,37,72,21)]
names(sub1)=c("Label","Min_depth","Max_depth","Taxa","area","major","minor",
              "Tow_Vol","Sub_part", "feret","esd","gray_mode")
sub1$cruise<-as.factor(rep("jc0214", length(sub1$Label)))

sub1$tow<-as.factor(with(sub1, ifelse(Label %like%"^20",paste("m",(substr(Label,3,4)), sep=""),
                                      ifelse(Label %like% "m47", "m47","m48"))))
sub1$Label<-as.character(sub1$Label)
sub1$net<-as.factor(with(sub1, ifelse(Label %like%"^20",substr(Label,5,6),
                                    ifelse(Label %like% "jc0214", sapply(strsplit(Label, "_"), function(x) x[3]), "NA"))))
nets<-c("n2","n3","n4","n5","n6","n7","n8","n9","n10","n10","n9")
levels(sub1$net)<-nets

sub1$tow_net<-as.factor(paste(sub1$tow,sub1$net, sep="_"))


levels(sub1$tow_net)[levels(sub1$tow_net)=="m12_n2"] <- "m11_n10"
sub1$tow<-as.factor(as.character(substr(sub1$tow_net,1,3)))
sub1$net<-as.factor(as.character(sub("^[^_]*_([^_]*).*","\\1",sub1$tow_net)))
#sub1$net<-levels(c("n2","n3","n4","n5","n6","n7","n8","n9","n10"))
summary(sub1)

sub1$fraction<-as.factor(with(sub1, ifelse(Label %like%"^20",substr(Label,8,11),
                                           ifelse(Label %like% "jc0214", sapply(strsplit(Label, "_"), function(x) x[4]), "NA"))))



sub1$area_mm2<-sub1$"area"*0.000112 #0.00002809 for 4800
sub1$major_mm<-sub1$"major"*.010583 #.0053 for 4800
sub1$minor_mm<-as.numeric(sub1$"minor"*.010583) #.0053 for 4800
sub1$feret_mm<-as.numeric(sub1$"feret"*0.010583)
sub1$esd_mm<-as.numeric(sub1$esd*0.010583)
sub1$vol<-(4/3)*pi*((sub1$minor_mm*0.5)^2)*(sub1$major_mm/2)
sub1$hdif<-(sub1$Max_depth-sub1$Min_depth)
sub1$split<-(1/(sub1$Sub_part))
sub1$L_D<-as.factor(ifelse(sub1$Taxa %in%"not-living","not-living","living"))


#Other environmental factors can be added to this spreadsheet and indexed using these two lines
#Moc_meta_index <- read_excel("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Moc_meta_index.xlsx")

# Helena wrote:
Moc_meta_index <- read_xlsx("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/Moc_meta_index.xlsx")


sub1$temp<-as.numeric(as.character(Moc_meta_index$Temp[match(sub1$tow_net,Moc_meta_index$tow_net)]))
sub1$D_N<-as.factor(as.character(Moc_meta_index$DayNight[match(sub1$tow_net,Moc_meta_index$tow_net)]))

sub1$tow_net_bin<-as.factor(paste(sub1$tow_net,sub1$fraction, sep="_"))

#scan_split_index <- read_excel("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/split_fix_index.xlsx")

# Helena wrote: 
scan_split_index <- read_xlsx("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/split_fix_index.xlsx")


sub1$tow_net_bin<-as.factor(paste(sub1$tow_net,sub1$fraction, sep="_"))
sub1$sub_part_correct<-as.factor(as.character(scan_split_index$Acq_sub[match(sub1$tow_net_bin,scan_split_index$tow_net_fraction)]))
sub1$split_correct<-as.numeric(as.character(scan_split_index$Split[match(sub1$tow_net_bin,scan_split_index$tow_net_fraction)]))


sub1$DW<-with(sub1, ifelse(Taxa %in% "Calanoida", 0.0145*vol,
                                  ifelse(Taxa %in% "Ostracoda", 0.0517*vol,
                                         ifelse(Taxa %in% "Thecosomata",0.0647*vol,
                                                ifelse(Taxa %in% "Amphipoda",0.0285*vol,
                                                       ifelse(Taxa %in% "Euphausiacea",0.0341*vol,
                                                             0.0145*vol))))))

sub1$O2_umol<-with(sub1, ifelse(Taxa %in% "Copepoda", ((exp(-0.399+(0.801*log(DW)))+0.069*(temp))/22.4),
                                ifelse(Taxa %in% "Chaetognatha", ((exp(-0.173+(0.805*log(DW)))+0.068*(temp))/22.4),
                                       ifelse(Taxa %in% "Amphipoda", ((exp(0.407+(0.743*log(DW)))+0.037*(temp))/22.4),
                                              ifelse(Taxa %in% "Euphausiacea",((exp(0.392+(0.753*log(DW)))+0.046*(temp))/22.4),
                                                     ifelse(Taxa %in% "Mollusca",((exp(-0.56+(0.82*log(DW)))+0.046*(temp))/22.4),
                                                            ((exp(-0.399+(0.801*log(DW)))+0.069*(temp))/22.4)))))))

sub1$CO2<-with(sub1, ifelse(Taxa %in% "Copepoda", O2_umol*0.87,
                            ifelse(Taxa %in% "Chaetognatha", O2_umol*1.35,
                                   ifelse(Taxa %in% "Amphipoda", O2_umol*1.35,
                                          ifelse(Taxa %in% "Euphausiacea",O2_umol*1.35,
                                                 ifelse(Taxa %in% "Mollusca",O2_umol*0.94,
                                                        O2_umol*0.87))))))

sub2<-sub1[,c(1,13:15,17,30:31,2:3,24,8,12,18:23,4,26,28,32:34)]
colnames(sub2)


# write.csv(sub2,file=paste
          # ("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Exports 2021 long.csv", sep=""),row.names=F)

# Helena wrote: 
write.csv(sub2,file=paste("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/Exports 2021 long.csv", sep=""), row.names=FALSE)

head(sub2) #Currently removing all not-living individuals
sub3<-sub2%>%filter(L_D!="not-living")%>%droplevels


### Table for Adrian ################
summary(sub3)
sub3$BS<-as.factor(ifelse(sub3$esd_mm<5,"Small","Big"))
sub4<-sub3%>%filter(fraction=="0200")%>%droplevels

summary(sub4)

sum1<-sub3%>%group_by(tow,net,fraction,BS)%>%
  summarize(count=n(),depth_min=mean(Min_depth),depth_max=mean(Max_depth),hdif=mean(hdif),split=median(split_correct),
            freq=count/split,tv=mean(Tow_Vol), Density_m3=freq/tv,
            NBV_m3=(sum(vol)/split/tv),BM_m3=(sum(DW)/split/tv),Abundance_m2=(freq/tv*hdif),NBV_m2=NBV_m3*hdif,
            BM_m2=BM_m3*hdif) 
sum2<-sum1%>%group_by(tow,net,fraction,BS)%>%
  summarize(count=sum(count),depth_min=mean(depth_min),depth_max=mean(depth_max),Density_m3=sum(Density_m3),Abundance_m2=sum(Abundance_m2),
            Biomass_m2=sum(BM_m2),Biovol_m2=sum(NBV_m2))

# write.csv(sum2,file=paste
        #   ("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Big_Small_Zoop_Table_Fractions_4July23.csv", sep=""),row.names=F)

# Helena wrote:
write.csv(sum2,file=paste
          ("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/Big_Small_Zoop_Table_Fractions_4July23.csv", sep=""),row.names=FALSE)
          

#Currently set to drop anything outside 1-5mm feret and groups we are not interested in
#M_filt<-sub3%>%filter(esd_mm>=1&esd_mm<=5)
#summary(M_filt)
#M_filt<-sub3

#B<-as.character(c(seq(0.25, 36.25, by=0.25)))
#M_filt$bin<-cut(M_filt$esd_mm, breaks=c(seq(0.25, 36.5, by=0.25)), labels=B)
#M_filt$bin<-factor(M_filt$bin)
#M_filt$Nbin<-as.numeric(as.character(M_filt$bin))
#M_filt$net<-as.factor(M_filt$net)

sub3$log_esd<-log10(sub3$esd_mm)
summary(sub3)
B<-as.character(c(seq(-0.7, 1.9, by=0.1)))
M_filt<-sub3
M_filt$bin<-cut(M_filt$log_esd, breaks=c(seq(-0.7, 2, by=0.1)), labels=B)
M_filt$bin<-factor(M_filt$bin)
M_filt$Nbin<-as.numeric(as.character(M_filt$bin))
M_filt$net<-as.factor(M_filt$net)

levels(M_filt$fraction)

summary<-M_filt%>%group_by(tow,net,fraction,bin)%>%
  summarize(count=n(),depth=(mean(Min_depth)+mean(Max_depth))/2,hdif=median(hdif),split=median(split_correct),
            freq=count/split,tv=mean(Tow_Vol), Density_m3=freq/tv, bin2=mean(Nbin),
            NBV_m3=(sum(vol)/split/tv),BM_m3=(sum(DW)/split/tv), oxy_m3=(sum(O2_umol))/tv/split,
            CO2_m3=(sum(CO2)/tv/split),Abundance_m2=(freq/tv*hdif),NBV_m2=NBV_m3*hdif,
            BM_m2=BM_m3*hdif, oxy_m2=oxy_m3*hdif,CO2_m2=CO2_m3*hdif,gray_mode=mean(gray_mode)) 

#write.csv(summary,file=paste
    #      ("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Exports_2021_logbins.csv", sep=""),row.names=F)

# Helena wrote:
write.csv(summary,file=paste
          ("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/Exports_2021_logbins.csv", sep=""),row.names=FALSE)


### WATERFALLS #######
## Fully scanned pair is m47 (day) and m48 (night)
WF<-summary%>%filter(tow=="m48") %>% droplevels
summary(WF)


WF$bin<-as.factor(WF$bin)
WF_sub<-WF %>% group_by(net,bin) %>% summarize(Abundance_m2=sum(Abundance_m2))
WF_sub$binN<-as.numeric(as.character(WF_sub$bin))
WF_sub$net<-factor(WF_sub$net, levels=c("n2","n3","n4","n5","n6","n7","n8","n9","n10"))

# FILLS IN THE DEPTH INTERVALS FOR THE SPECIFIC MOCNESS

net_labs<-c("0-50 m","50-100 m","100-150 m","150-200m","200-300 m", "300-400 m",
            "400-500 m","500-750 m", "750-1000 m")



#NAME OF THE MOCNESS
Title<-"EXPORTS Atlantic 2021 (Night)"


WFplot<-ggplot(data=WF_sub, aes(x=binN, y=Abundance_m2,color=net))+
  geom_point(size=2)+
  scale_color_manual(values=(lacroix_palette("PeachPear", type="continuous", n=9)), labels=rev(net_labs))+
  labs(x=expression("logESD Class"~(mm)), y=expression("Abundance"~(particles~m^-2)),
       title=paste(Title), color="")+
  guides(color=guide_legend(reverse=T))+
  coord_cartesian(clip='off')+
  scale_y_log10(limits=c(0.0001,100000))+
  #you can limit your scale based on what biomass you effectively sample, but you should look to see all data first
  
  theme(plot.title=element_text(face="bold",hjust=0.5, size=14),
        legend.position='bottom')+
  theme(legend.text=element_text(size=11))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=14),axis.text.y= element_text(size=14),
        axis.title=element_text(size=14),strip.text=element_text(size=14, face="bold"))

WFplot
# set your path
# ggsave(WFplot,filename=paste(Title,"Waterfall.png", sep=" "), path=paste("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Waterfalls/", sep=""), width=6, height=5, units="in" )

# Helena commented out because we don't need to save this plot

## This code gives you linear regressions and R2 of the correlation for each line.
library(data.table)
## set the bins to exclude the nets with poor capture (particularly on the low end)
#WF_sub_trim<-WF_sub%>%filter(binN>0.1&binN<10)
dt<-data.table(WF_sub, key="net")
fits<-lapply(unique(dt$net),function(z){
  summary(lm(log(Abundance_m2)~(binN), data=dt[J(z),], y=T))
})

fits
## this saves the stats to a file. Set your path.
# capture.output(fits,file=paste("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Waterfalls/",Title,"_reg_stats.txt"))

# Helena commented out because we don't need to save this, and stopped editing code here



### Heatmaps #########################################################################


library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(pals)

path<-paste("C:/Users/localadmin.BIOS/Dropbox (BIOS)/ZoopGroup_LAB/EXPORTS/Cruise 2021/Zooscan/Heatmaps/")

#SET YOUR PATH HERE, PICKING THE DAY OF A PAIR
day2<-summary%>%filter(tow=="m47")%>%droplevels
#SET YOUR PATH HERE, PICKING THE NIGHT OF A PAIR
night2<-summary%>%filter(tow=="m48")%>%droplevels


J<-"Exports 2021 N Atlantic" #SET THE NAME OF THE PAIR
#SET THE DEPTH INTERVALS FOR THE PAIR
d_labs=c("1000","750","500","400","300","200","150","100","50","0")


day2$bin<-as.factor(day2$bin)
night2$bin<-as.factor(night2$bin)
dn.hm<-merge(day2,night2, by.x=c("net","bin","fraction"),by.y=c("net","bin","fraction"),all.x=TRUE,all.y=TRUE)
head(dn.hm)
summary(dn.hm)
dn.hm[is.na(dn.hm)]<-0
dn.hm2<-dn.hm%>%group_by(net,bin)%>%summarize(BVm3_day=sum(NBV_m3.x), BVm3_night=sum(NBV_m3.y),
                                              BVm2_day=sum(NBV_m2.x), BVm2_night=sum(NBV_m2.y),
                                              BMm2_day=sum(BM_m2.x), BMm2_night=sum(BM_m2.y),
                                              Oxm2_day=sum(oxy_m2.x),Oxm2_night=sum(oxy_m2.y),
                                              CO2m2_day=sum(CO2_m2.x),CO2m2_night=sum(CO2_m2.y),
                                              gray_mode_day=sum(gray_mode.x),gray_mode_night=sum(gray_mode.y))

#write.csv(dn.hm2,file=paste(path,J,"_DayNight_dnhm2.csv",sep=""),row.names=FALSE)

dn.hm3<-dn.hm2%>%mutate(BV_m3=(BVm3_day-BVm3_night), BV_m2=(BVm2_day-BVm2_night),BM_m2=(BMm2_day-BMm2_night),
                        Ox=(Oxm2_day-Oxm2_night),CO2=(CO2m2_day-CO2m2_night),gray=(gray_mode_day-gray_mode_night)) %>%select(,c(1:2,13:20))  


dn.hm3<-as.data.frame(dn.hm3)
#dn.hm3$binnum<-as.numeric(as.character(dn.hm3$bin))
#dn.hm4<-filter(dn.hm3,binnum>0.01&binnum<100)

dn.hm3$net<-factor(dn.hm3$net, levels=c("n2","n3","n4","n5","n6","n7","n8","n9","n10"))
summary(dn.hm3$net)
levels(dn.hm3$net)<-c("n2","n3","n4","n5","n6","n7","n8","n9","n10","n11")

summary(dn.hm3)
write.csv(dn.hm3,file=paste(path,J,"_DayNight_dnhm3.csv",sep=""),row.names=FALSE)

BV.heatmap<-ggplot(data=dn.hm3, mapping=aes(x=bin, y=net, fill=BV_m2, color=""))+
  geom_tile()+
  ylab(label="Net Depth Range(m)")+
  scale_fill_continuous_divergingx(na.value="gray45",limits=c(-5700,5700),palette = 'RdBu',
                                   rev=TRUE, mid =0, l3 = 0, p1=0.4, p3 = .4, p4 = .5)+
  scale_x_discrete(breaks=c("-1","-0.5", "1.11022302462516e-16", "0.5","1", "1.5","2"), 
                   labels=c("-1","-0.5", "0", "0.5","1", "1.5","2"))+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12, vjust=2.5),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste(J,'(Day-Night) Plankton Biovolume Shift'))+
  labs(fill=expression("Biovolume"~(mm^3/m^2)),
       x=expression("logESD Size Class (mm)"),
       color=expression("<-5700"~mm^3/m^2))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
BV.heatmap
ggsave(BV.heatmap,filename=paste(J,"BV_Shift.png", sep="_"),path=path, width=8, height=5, units="in", dpi=300)

Ox.heatmap<-ggplot(data=dn.hm3, mapping=aes(x=bin, y=net, fill=Ox, color=""))+
  geom_tile()+
  ylab(label="Net Depth Range(m)")+
  scale_fill_continuous_divergingx(na.value="gray45",limits=c(-525,525),palette = 'RdBu',
                                   rev=TRUE, mid =0, l3 = 0, p1=0.4, p3 = .4, p4 = .5)+
  scale_x_discrete(breaks=c("-1","-0.5", "1.11022302462516e-16", "0.5","1", "1.5","2"), 
                   labels=c("-1","-0.5", "0", "0.5","1", "1.5","2"))+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12, vjust=2.5),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste(J,'(Day-Night) Plankton Respiratory Shift'))+
  labs(fill=expression("O2"~(mu*mol~m^-2*h^-1)),
       x=expression("logESD Size Class (mm)"),
       color=expression("<525"~mu*mol~m^-2*h^-1))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
Ox.heatmap
ggsave(Ox.heatmap,filename=paste(J,"Ox_Shift.png", sep="_"),path=path, width=8, height=5, units="in", dpi=300)


Gray.day.heatmap<-ggplot(data=dn.hm3, mapping=aes(x=bin, y=net, fill=gray_mode_day))+
  geom_tile()+
  ylab(label="Net Depth Range(m)")+
  scale_fill_gradientn(colours=rev(ocean.ice(300)), guide = "colourbar")+
  scale_x_discrete(breaks=c("-1","-0.5", "1.11022302462516e-16", "0.5","1", "1.5","2"), 
                   labels=c("-1","-0.5", "0", "0.5","1", "1.5","2"))+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12, vjust=2.5),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste(J,'Median Gray, Day'))+
  labs(fill=expression("Gray shade"),
       x=expression("logESD (mm)"))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
Gray.day.heatmap
ggsave(Gray.day.heatmap,filename=paste(J,"Gray day.png", sep="_"),path=path, width=8, height=5, units="in", dpi=300)

Gray.night.heatmap<-ggplot(data=dn.hm3, mapping=aes(x=bin, y=net, fill=gray_mode_night))+
  geom_tile()+
  ylab(label="Net Depth Range(m)")+
  scale_fill_gradientn(colours=rev(ocean.ice(300)), guide = "colourbar")+
  scale_x_discrete(breaks=c("-1","-0.5", "1.11022302462516e-16", "0.5","1", "1.5","2"), 
                   labels=c("-1","-0.5", "0", "0.5","1", "1.5","2"))+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12, vjust=2.5),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste(J,'Median Gray, Night'))+
  labs(fill=expression("Gray shade"),
       x=expression("logESD (mm)"))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
Gray.night.heatmap
ggsave(Gray.night.heatmap,filename=paste(J,"Gray_night.png", sep="_"),path=path, width=8, height=5, units="in", dpi=300)


p1<-grid.arrange(arrangeGrob(Gray.day.heatmap+theme(legend.position="none"),Gray.night.heatmap, ncol=2),
                 nrow=2,heights=c(10, 1))

ggsave(p1,filename=paste("DN_Opacity.png", sep="_"),path=path, width=12, height=7, units="in", dpi=300)


Gray.dif.heatmap<-ggplot(data=dn.hm3, mapping=aes(x=bin, y=net, fill=gray))+
  geom_tile()+
  ylab(label="Net Depth Range(m)")+
  scale_fill_continuous_divergingx(na.value="gray45",limits=c(-600,600),palette = 'RdBu',
                                   rev=TRUE, mid =0, l3 = 0, p1=0.4, p3 = .4, p4 = .5)+
  scale_x_discrete(breaks=c("-1","-0.5", "1.11022302462516e-16", "0.5","1", "1.5","2"), 
                   labels=c("-1","-0.5", "0", "0.5","1", "1.5","2"))+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12, vjust=2.5),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold",size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste(J,'Median Gray Shift, Day-Night Difference'))+
  labs(fill=expression("Gray shade"),
       x=expression("logESD (mm)"))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
Gray.dif.heatmap
ggsave(Gray.dif.heatmap,filename=paste(J,"Gray_dif.png", sep="_"),path=path, width=8, height=5, units="in", dpi=300)


