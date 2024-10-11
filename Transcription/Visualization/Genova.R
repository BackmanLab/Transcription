#########################################################################################
## Genova ##
## This code was used to generate the pileup plots in
## "Extra-nucleolar Pol I regulates gene transcription through chromatin domain maintenance"
#########################################################################################

rm(list = ls())
gc()

# import packages
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)
library(forcats)
library(GENOVA)

## Broadly, colours of heatmaps can be changed by setting the GENOVA.colour.palette option.
## https://github.com/robinweide/GENOVA/issues/298
options("GENOVA.colour.palette" = "whitered")

## directories and variables 
dir.1 <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/juicer_analysis/"
dir.2 <- "/projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/"
dir.3 <- "/projects/p32171/juicer/work/112123_HiC/juicer_analysis/"

list.files(paste0(dir.1), all.files=F, include.dirs = FALSE)

hic ="inter_30.hic"
exps <- c("1hr_ActD",  "Pol2_6hrs_Aux","Pol1_6hrs_Aux", "WT_HCT116_CTRL")

## -- introduce data -- ##
## contact probability resolution - 5kb
wt.5kb <- load_contacts(signal_path = paste0(dir.2,exps[4],"/mega/aligned/",hic),
                        sample_name = "WT",
                        resolution = 5e3,
                        balancing = 'KR', # this is the default
                        colour = "black")

pol1.5kb <- load_contacts(signal_path = paste0(dir.1,exps[3],"/mega/aligned/",hic),
                          sample_name = "Pol1",
                          resolution = 5e3,
                          balancing = 'KR', # this is the default
                          colour = "black")

actd.5kb <- load_contacts(signal_path = paste0(dir.3,exps[1],"/mega/aligned/",hic),
                          sample_name = "ActD",
                          resolution = 5e3,
                          balancing = 'KR', # this is the default
                          colour = "black")

pol2.5kb <- load_contacts(signal_path = paste0(dir.3,exps[2],"/mega/aligned/",hic),
                          sample_name = "Pol2",
                          resolution = 5e3,
                          balancing = 'KR', # this is the default
                          colour = "black")


gc()

##-------------------------------------------------------## Loop files 
wt.loops = data.frame(read.delim('/projects/p32171/HiC2/opt/juicer/work/112123_HiC/loop_analysis/WT_HCT116_CTRL/mega/hiccups_results/merged_loops.bedpe', h = F, skip = 1))
pol1.loops = data.frame(read.delim('/projects/p32171/HiC2/opt/juicer/work/112123_HiC/loop_analysis/Pol1_6hrs_Aux/mega/hiccups_results/merged_loops.bedpe', h = F, skip = 1))

actd.loops = read.delim('/projects/p32171/HiC2/opt/juicer/work/112123_HiC/loop_analysis/1hr_ActD/mega/hiccups_results/merged_loops.bedpe', h = F, skip = 1)
pol2.loops = read.delim('/projects/p32171/HiC2/opt/juicer/work/112123_HiC/loop_analysis/Pol2_6hrs_Aux/mega/hiccups_results/merged_loops.bedpe', h = F, skip = 1)

## Diff loops
wt.dl.loops = read.delim('/projects/p32171/HiC2/opt/juicer/work/112123_HiC/loop_analysis/Pol2_6hrs_Aux/mega/hiccups_diff/pol2vswt/differential_loops2.bedpe', h = F, skip = 1)
pol2.dl.loops = read.delim("/projects/p32171/HiC2/opt/juicer/work/112123_HiC/loop_analysis/Pol2_6hrs_Aux/mega/hiccups_diff/pol2vswt/differential_loops1.bedpe", h = F, skip = 1)


cf_DFinf2NA <- function(x)
{
  for (i in 1:ncol(x)){
    x[,i][is.infinite(x[,i])] = NA
  }
  x <- na.omit(x)
  x<- x[-c(1),]
  
  return(x)
}

wt.loops <- cf_DFinf2NA(wt.loops)
pol1.loops <- cf_DFinf2NA(pol1.loops)
pol2.loops <- cf_DFinf2NA(pol2.loops)
actd.loops <- cf_DFinf2NA(actd.loops)
wt.dl.loops <- cf_DFinf2NA(wt.dl.loops)
pol2.dl.loops <- cf_DFinf2NA(pol2.dl.loops)

##-------------------------------------------------------## TAD files
wt.domains= read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_domains/arrowhead_domains/WT_HCT116_CTRL/mega/inter_30_contact_domains/5000_blocks.bedpe', h = F, skip = 1)
pol1.domains = read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_domains/arrowhead_domains/Pol1_6hrs_Aux/mega/inter_30_contact_domains/5000_blocks.bedpe', h = F, skip = 1)
actd.domains = read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_domains/arrowhead_domains/1hr_ActD/mega/inter_30_contact_domains/5000_blocks.bedpe', h = F, skip = 1)
pol2.domains = read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_domains/arrowhead_domains/Pol2_6hrs_Aux/mega/inter_30_contact_domains/5000_blocks.bedpe', h = F, skip = 1)

cf_DFinf2NA <- function(x)
{
  for (i in 1:ncol(x)){
    x[,i][is.infinite(x[,i])] = NA
  }
  x <- na.omit(x)
  x<- x[-c(1),]
  
  return(x)
}

wt.domains<- cf_DFinf2NA(wt.domains)
pol1.domains<- cf_DFinf2NA(pol1.domains)
pol2.domains <- cf_DFinf2NA(pol2.domains)
actd.domains <- cf_DFinf2NA(actd.domains)

##-------------------------------------------------------## Color scales

library(rcartocolor)
display_carto_all( colorblind_friendly = TRUE) ## get colorblind carto colors
my_colors = carto_pal(7, "Earth") 
my_colors = carto_pal(7, "ArmyRose") ## list carto colors
my_colors

## control: "#A16928" treatment: "#2887A1" ## Earth
## control: "#798234" treatment: "#D46780" ## ArmyRose

##-------------------------------------------------------## APA plotting

synced <- sync_indices(list("WT" = wt.5kb,'Pol1' = pol1.5kb))
synced.2 <- sync_indices(list("WT" = wt.5kb,'ActD' = actd.5kb))
synced.3 <- sync_indices(list("WT" = wt.5kb,'Pol2' = pol2.5kb))
synced.4 <- sync_indices(list("Pol1" = pol1.5kb,'Pol2' = pol2.5kb))

## Quantify loops | loop pileup
APA.p1<- APA(list("WT" = synced[[1]],'Pol1' = synced[[2]]),dist_thres = c(100e3, Inf),bedpe = wt.loops)
qapa.p1 <- quantify(APA.p1, size = 3)

APA.actd<- APA(list("WT" = synced.2[[1]],'ActD' = synced.2[[2]]),dist_thres = c(100e3, Inf),bedpe = wt.loops)
qapa.actd <- quantify(APA.actd, size = 3)

APA.p2<- APA(list("WT" = synced.3[[1]],'Pol2' = synced.3[[2]]),dist_thres = c(100e3, Inf),bedpe = wt.loops)
qapa.p2 <- quantify(APA.p2, size = 3)

df <- rbind(qapa.p1$per_sample, qapa.actd$per_sample, qapa.p2$per_sample)
write.csv(as.data.frame(df), "/home/lmc0633/Transcription Paper/loop_pileup_scores.csv")

gc()

## Colors for ATA part
v_colors = c("white","#ede5cf","#e0c2a2","#d39c83","#c1766f","#a65461","#813753","#541f3f", "black")

## filename: apa.pol1loops | apa.wtloops
visualise(APA, colour_lim = c(0,20),metric = "diff", focus = 1)+
  ggplot2::scale_fill_gradientn( colours = c(v_colors),limits= c(0,20), na.value="black")+
  ggplot2::scale_colour_gradientn(limits = c(-1.5,1.5),aesthetics = "altfill",colours = c(colors),na.value="#2887A1",guide = ggplot2::guide_colourbar(available_aes = "altfill"))


##-------------------------------------------------------## ATA plotting 

## aggregate TAD analysis
ATA <- ATA(list("WT" = synced.3[[1]],'Pol2' = synced.3[[2]]), bed = pol2.domains) ## Arrowhead TADs

## Colors for ATA part
v_colors = c("white","#ede5cf","#e0c2a2","#d39c83","#c1766f","#a65461","#813753","#541f3f", "black")

## filename: ata.pol1tads | ata.wt.tads
visualise(ATA, colour_lim = c(0,20),metric = "diff", focus = 1)+
  ggplot2::scale_fill_gradientn( colours = c(v_colors),limits= c(0,20), na.value="black")+
  ggplot2::scale_colour_gradientn(limits = c(-5,5),aesthetics = "altfill",colours = c(colors),na.value="#A16928",guide = ggplot2::guide_colourbar(available_aes = "altfill"))

##-------------------------------------------------------## Compartments  

## -- introduce data -- ##
## contact probability resolution - 5kb
wt.100kb <- load_contacts(signal_path = paste0(dir.2,exps[4],"/mega/aligned/",hic),
                          sample_name = "WT",
                          resolution = 100e3,
                          balancing = 'KR', # this is the default
                          colour = "black")

pol1.100kb <- load_contacts(signal_path = paste0(dir.1,exps[3],"/mega/aligned/",hic),
                            sample_name = "Pol1",
                            resolution = 100e3,
                            balancing = 'KR', # this is the default
                            colour = "black")

actd.100kb <- load_contacts(signal_path = paste0(dir.2,exps[1],"/mega/aligned/",hic),
                            sample_name = "ActD",
                            resolution = 100e3,
                            balancing = 'KR', # this is the default
                            colour = "black")

pol2.100kb <- load_contacts(signal_path = paste0(dir.2,exps[2],"/mega/aligned/",hic),
                            sample_name = "Pol2",
                            resolution = 100e3,
                            balancing = 'KR', # this is the default
                            colour = "black")

gc()

synced <- sync_indices(list("WT" = wt.5kb,'Pol1' = pol1.5kb))
synced.2 <- sync_indices(list("WT" = wt.5kb,'ActD' = actd.5kb))
synced.3 <- sync_indices(list("WT" = wt.5kb,'Pol2' = pol2.5kb))

##-------------------------------------------------------## Compartment score plots

## Get eigens
cs = compartment_score(synced)
cs.2 = compartment_score(synced.2)
cs.3 = compartment_score(synced.3)

## Saddle plot calculation
attr(cs.3, 'signed') <- T # Fake signed CS
saddle.plot = saddle(synced,
                     CS_discovery = cs,
                     bins = 50)

p <- visualise(saddle.plot)

## Colors for difference plot
colors <- colorRampPalette(c("#2887A1", "white","#A16928"))(12)

## Colors for ATA part
v_colors = c("white","#ede5cf","#e0c2a2","#d39c83","#c1766f","#a65461","#813753","#541f3f", "black")

## Remap for plotting with different color scale
p$layers[[1]]$mapping <- ggplot2::aes(fill = log2(obsexp))
lims <- range(log2(p$data$obsexp))
values <- scales::rescale(c(lims[1], 0, lims[2]))

## Plot saddle ## pol1.comp.saddle
p + ggplot2::scale_fill_gradientn(
  colours = c(v_colors),
  values  = values,
  limits = lims,
  name = "Log2(FC)"
) +ggplot2::scale_colour_gradientn(limits = c(-1,1),
                                   aesthetics = "altfill",
                                   colours = c(colors),na.value="#A16928",
                                   guide = ggplot2::guide_colourbar(available_aes = "altfill"))

##-------------------------------------------------------## Compartments Strength plots

library(PupillometryR)
library(ggrain)

CSS <- quantify(saddle.plot)
compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength')

plot.df <- compared %>% 
  gather(cond, value, 2:3) 
plot.df <- plot.df[plot.df$value < 10 & plot.df$value > 0,]

## Rain plot of compartment strength ## actd.comp.strength
plot.df %>%
  ggplot( aes(1, y=value, fill=cond, color = cond ))+
  geom_rain(alpha = .5, rain.side = 'l',trim = F,
            boxplot.args = list(color = "black", outlier.shape = NA),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = .095, width = 0.105), width = 0.075))  +
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_fill_manual(values = c("#2887A1", "#A16928" )) +
  scale_color_manual(values = c("#2887A1", "#A16928"))+   theme( plot.title = element_text(size=11)) +
  ggtitle("") + ylab("Compartment Strength") + xlab("")


## Plot differntial dot plot
with(compared, plot(WT, Pol2, xlim = c(0,5), ylim = c(0,5), pch = 20))
abline(a = 0, b = 1, lty = 5)
