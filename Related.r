###Download Related package with the tar.gz from 
##https://github.com/timothyfrasier/related

library(related) #to analyze the data
library(tidyverse) #to handle data
library(ggplot2) #to plot
library(stringr) #to break facet labels
library(patchwork) #to generate ggplot2 compositions
library(scico) #to generate color-blind friendly colors

setwd("/media/wilson/MIDGARD/Documents/working/Sibship Larvae/R/")

file.ca1 = "../Colony/LarvaeSibShipCa1/OffspringGenotype.txt"
file.ca2 = "../Colony/LarvaeSibShipCa2/OffspringGenotype.txt"
file.cin = "../Colony/LarvaeSibShipCi/OffspringGenotype.txt"

pop.ca1 <- read.table( file=file.ca1, header = T, sep = c(" "), skipNul = T)
pop.ca1 <- pop.ca1[,colnames(pop.ca1)!="X"]

pop.ca2 <- read.table( file=file.ca2, header = T, sep = c(" "), skipNul = T)
pop.ca2 <- pop.ca2[,colnames(pop.ca2)!="X"]

pop.cin <- read.table( file=file.cin, header = T, sep = c(" "), skipNul = T)
pop.cin <- pop.cin[,colnames(pop.cin)!="X"]

tail(pop.ca1)
tail(pop.ca2)
tail(pop.cin)

ca1 <- related::coancestry(pop.ca1, trioml = 1, wang = 1, lynchli = 1, lynchrd = 1, ritland = 1, quellergt = 1, dyadml = 1)

ca2 <- related::coancestry(pop.ca2, trioml = 1, wang = 1, lynchli = 1, lynchrd = 1, ritland = 1, quellergt = 1, dyadml = 1)

cin <- related::coancestry(pop.cin, trioml = 1, wang = 1, lynchli = 1, lynchrd = 1, ritland = 1, quellergt = 1, dyadml = 1)

r.ca1 <- cbind(ca1$relatedness, site="Ca1")
r.ca2 <- cbind(ca2$relatedness, site="Ca2")
r.cin <- cbind(cin$relatedness, site="Cin")

rel <- rbind(r.ca1, r.ca2,r.cin )
rel <- rel[,-c(1,4)]

plot.labels <- c("Milligan (2003)", "Li et al. (1993)", "Lynch & Ritland (1999)", 
                 "Queler & Goodnight (1989)", "Ritland (1996)", "Wang (2007)", "Wang (2002)")

ordered.labels <- c("Queler & Goodnight (1989)", "Li et al. (1993)", "Ritland (1996)",  
                 "Lynch & Ritland (1999)", "Wang (2002)", "Milligan (2003)","Wang (2007)")

qrel <- rel %>%
  pivot_longer(., cols = c("trioml", "wang", "lynchli",
                           "lynchrd", "quellergt", "dyadml","ritland"),
               names_to = "estimator", values_to = "r") %>%
  mutate(estimator=as.factor(estimator))

levels(qrel$estimator) <- plot.labels
qrel <- qrel %>% mutate(estimator = factor(estimator, levels = ordered.labels))
levels(qrel$estimator) <- stringr::str_wrap(ordered.labels, 20)

qrel.bak <- qrel
qrel <- qrel %>% 
  group_by(estimator) %>% mutate(r = sqrt((r/max(r))^2) )
qrel
## BloxPlot ####
bxplot <- ggplot(qrel, aes(x=site, y=r, group=site, fill=site))+
  geom_boxplot(alpha=0.8)+
  #geom_jitter(aes(fill=site, color=site), alpha=0.5)+
  facet_wrap(~estimator, nrow = 1, strip.position = "top")+
  labs(y="relative r-value")+
  scale_fill_scico_d(begin = 0.2, end = 0.6, palette = "batlow")+
  scale_color_scico_d(begin = 0.2, end = 0.6, palette = "batlow")+
  theme(
    strip.switch.pad.wrap = unit(0.5, "mm"),
    strip.text.y = element_text(angle = 0),
    legend.position = "none",
    axis.title.x = element_blank()
  )
bxplot
# DensPlot wraped by site~estimators ####
#qrel <- qrel %>% mutate(estimator = factor(estimator, levels = ordered.labels))
densplot <- ggplot(qrel.bak, aes(x=r, group=estimator, color=site, fill=site))+
  geom_density(size=0.4)+
  facet_grid(estimator~site)+
labs(y="kernel density distribution",
     x= "relatedness (r)")+
#  scale_fill_scico_d(palette = "batlow", begin = 0.2,  end = 0.75)+
#  scale_color_scico_d(palette = "batlow", begin = 0.2, end = 0.75)+
  scale_fill_scico_d(begin = 0.2, end = 0.6, palette = "batlow")+
  scale_color_scico_d(begin = 0.2, end = 0.6, palette = "batlow")+
  scale_x_continuous(limits=c(-0.5, 0.5))+
#  labs(title = "A")+
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "none"
    )
densplot  
#DensePlot All overlayed ####
densplotflip <- 
  ggplot(qrel.bak, aes(x=r, y=..density.., group=site, color=site, fill=site))+
  geom_density(aes(linetype=site), alpha=0)+
  facet_grid(~estimator)+
  scale_fill_scico_d(palette = "batlow", begin = 0.2,  end = 0.75, alpha=0)+
  scale_color_scico_d(palette = "batlow", begin = 0.2, end = 0.75)+
  scale_linetype_manual(values = c("dotted", "longdash", "solid"))+
  labs(y="kernel density distribution", x="relatedness (r)")+
  scale_x_continuous(limits=c(-0.5, 1))+
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "right"
  )
densplotflip
theme_common <- theme(
  panel.border = element_rect(color = "#CCCCCC", fill="transparent"),
)

## Composition ####
densplot+theme_common
bxplot+theme_common
layout <- c(
  area(t = 1, l = 1, b = 3, r = 100),
  area(t = 4, l = 1, b = 4, r = 101)
)
composition <- I(densplot+theme_common) / I(bxplot+theme_common) + plot_layout(design = layout)
composition

#HeatMap
dheat <- qrel.bak %>% 
  mutate(ind1.id = as.factor(ind1.id), ind2.id = as.factor(ind2.id)) %>%
  mutate(ind1.id = factor(ind1.id, levels = unique(ind1.id)),
         ind2.id = factor(ind2.id, levels= unique(ind2.id)))

#to create labels with interval for axes
build_discrite_labels <- function(vec=NULL, n.breaks=5){
  if(is_tibble(vec) | is.data.frame(vec)) vec <- unlist(vec)
  vec <- factor(vec, levels=unique(vec))
  vec.levels <- levels(vec)
  doreplace <- rep(c(TRUE, rep(FALSE, n.breaks-1)), ceiling(length(vec.levels)/n.breaks))
  out.labels <- vec.levels
  out.labels[!doreplace] <- ""
  out.labels <- out.labels[I(length(heatlbs)+1):length(doreplace)*-1]
  names(out.labels) <- vec.levels
  return(out.labels)
}
#To filter data for each site
ggheatmap <- function(df=NULL, statistic=NULL, locality=NULL, cust.min=round(min(df$r),2)){
  if(is.null(df)) stop("Did you forget to provide data?")
  if(!is.null(statistic)) df <- df %>% filter(estimator == statistic)
  if(!is.null(locality)) df <- df %>% filter(site==locality)
  
  outplot <- df %>% ggplot(aes(x=ind1.id, y=ind2.id, fill=r,group=site))+
    geom_tile(color=rgb(1,1,1,0.5))+
    scale_fill_stepsn(n.breaks=5, 
                      colours=scico(5,
                                    direction = -1,
                                    palette = "batlow"), 
                      limits=c(cust.min,1))+
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size=8),
      panel.background = element_blank(),
      panel.grid = element_line(color="#DDDDDD")
    )
  outplot <- outplot +
    scale_x_discrete(labels=build_discrite_labels(outplot$data$ind1.id,4))+
    scale_y_discrete(labels=build_discrite_labels(outplot$data$ind2.id,4))
  return(outplot)
}
heatplot_common_theme <-theme(
  axis.text = element_text(color="black"),
  legend.position = c(0.9,0.33),
  legend.background = element_blank(),
  legend.title = element_text(size=12, face="italic", hjust = 0.25),
  legend.key.size = unit(5, "mm"),
  legend.text.align = 0,
  legend.spacing.x = unit(1,"mm")
)
plot1 <- ggheatmap(dheat, statistic = "Queler & Goodnight\n(1989)", locality = "Ca1", cust.min = -0.3)+heatplot_common_theme+labs(title = "A")
plot2 <- ggheatmap(dheat, statistic = "Queler & Goodnight\n(1989)", locality = "Ca2", cust.min = -0.3)+heatplot_common_theme+labs(title = "B")
plot3 <- ggheatmap(dheat, statistic = "Queler & Goodnight\n(1989)", locality = "Cin", cust.min = -0.3)+heatplot_common_theme+labs(title = "C")
plot4<-ggplot(filter(dheat, estimator=="Queler & Goodnight\n(1989)"))+
  geom_boxplot(aes(x=site, y=r, fill=site))+
  labs(y="pairwise relatedness (r)")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        panel.grid = element_blank()
        )+
  scale_fill_scico_d(palette = "batlow", begin=0.2, end=0.6)+
  labs(title = "D")

plot1 + plot2 + plot3 + plot4
