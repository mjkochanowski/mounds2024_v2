library(ggplot2)
library(vegan)
library(car)
library(ggpubr)
library(dplyr)

#################################################################

### Colony forming units calculations ###

cfu=read.csv("cfu5.csv", sep=',')
cfu$substrate=factor(cfu$substrate, levels=c('forest litter', 'mound surface','mound interior'), labels = c('forest litter', 'mound surface','mound interior'))
cfu$season=factor(cfu$season, levels=c('summer','autumn'), labels = c('August', 'October'))

leveneTest(cfu ~ season*substrate, data = cfu)
aov.cfu <- aov(cfu ~ season*substrate, data = cfu)
summary(aov.cfu)
tukey.cfu=TukeyHSD(aov.cfu, data=cfu)

#################################################################

### taxonomical diveristy ###
sp=read.csv('sp5.csv')

spm = sp[,-(1:4)]
sp_meta=sp[,1:4]
rownames(spm) <- sp[,1]
sppr=specnumber(spm)
spp=cbind(sp_meta,sppr)

leveneTest(sppr ~ season*substrate, data = sp)
aov.taxa <- aov(sppr~substrate*season, data=sp)
summary(aov.taxa)
TukeyHSD(aov.taxa)

#################################################################
### shannon diversity ###

shannondiv <- diversity(spm)
df=cbind(sp_meta,shannondiv)
head(df)

leveneTest(shannondiv ~ season*substrate, data = df)
aov.shan <- aov(shannondiv~substrate*season, data=df)
summary(aov.shan)
TukeyHSD(aov.shan)

#################################################################
### PERMANOVA ###

kopce=read.csv('sp5.csv', header = T, row.names = 1)
kopce.matrix=kopce[,-(1:3)]
kopce.meta=kopce[,(1:3)]

set.seed(27041998)
dist=vegdist(kopce.matrix, method = "bray")
bd=betadisper(dist, kopce.meta$substrate)
anova(bd)
perma=adonis(dist~substrate*season+substrate*site, data=kopce.meta)

### pariwaise post hoc ###
library(pairwiseAdonis)

pair.substrat<-pairwise.adonis(dist,factors=kopce.meta$substrate , p.adjust.m = 'bonferroni')
pair.site<-pairwise.adonis(dist,factors=kopce.meta$site , p.adjust.m = 'bonferroni')
pair.substrat.season<-pairwise.adonis(dist,factors=factor(paste(kopce.meta$substrate,kopce.meta$season)) , p.adjust.m = 'bonferroni')

#################################################################
###NMDS###

kop <- read.csv(file = "sp5.csv", header = TRUE, row.names = 1)
kop.sp=kop[,-(1:3)]
kop.meta=kop[,1:3]

set.seed(123)

sp_distmat <- vegdist(kop.sp, method = "bray")

kop.envfit <- envfit(kop.mds, kop.meta, permutations = 999) # this fits environmental vectors
kop.spp.fit <- envfit(kop.mds, kop.sp, permutations = 9999) # this fits species vectors

site.scrs <- as.data.frame(scores(kop.mds, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, season = kop.meta$season) #add grouping variable "season" to dataframe
site.scrs <- cbind(site.scrs, substrate = kop.meta$substrate) #add grouping variable of cluster grouping to dataframe
site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot

spp.scrs <- as.data.frame(scores(kop.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = kop.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<0.05) #subset data to show species significant at 0.05

#######################################################

###indicator species

library(indicspecies)

pc = read.csv("sp5.csv", header= TRUE, sep="\t")
abund = pc[,5:ncol(pc)]
type = pc$substrate

inv.IV = multipatt(abund, type, func = "IndVal", control = how(nperm=999))
summary(inv.IV, alpha=0.05)
inv.rg = multipatt(abund, type, func = "r", control = how(nperm=999))
summary(inv.rg, alpha=0.05)

presence_absence_df <- abund
presence_absence_df[presence_absence_df != 0] <- 1

inv.IV = multipatt(presence_absence_df, type, func = "IndVal", control = how(nperm=999))
summary(inv.IV, alpha=0.05)
inv.rg = multipatt(presence_absence_df, type, func = "r", control = how(nperm=999))
summary(inv.rg, alpha=0.05)

##############################################

### Figures ###

### Figure 3 - CFUs ###
cfu=read.csv("cfu5.csv", sep=',')
cfu$substrate=factor(cfu$substrate, levels=c('forest litter', 'mound surface','mound interior'), labels = c('forest litter', 'mound surface','mound interior'))
cfu$season=factor(cfu$season, levels=c('summer','autumn'), labels = c('August', 'October'))

cfu_plot1=ggplot(cfu, aes(x=season,y=cfu, fill=season))+
  geom_boxplot()+
  theme_classic()+
  labs(x='', fill='season', y='CFUs per sample')+
  scale_fill_manual(values=c("#ffc20a", "#0c7bdc"))+
  ylim(0,160)+
  theme(axis.text.x=element_text(size=12),
		axis.text.y=element_text(size=12),
		axis.title = element_text(size = 12),
        legend.position = 'none',
        plot.margin = margin(0.75,1.5,0.75,0.75, 'cm'))+
  geom_text(aes(x = 1.5, y = 150, label = "p = 0.05", fontface=2), position = position_dodge(width = 0.8), vjust = -0.5, size = 4, inherit.aes = FALSE)+
  geom_text(aes(x = 0.75, y = 160, label = "A"), size = 6, fontface =2,  inherit.aes = FALSE)+
  geom_segment(aes(x = 1, xend = 2, y = 150, yend = 150), inherit.aes = FALSE)+
  geom_segment(aes(x = 1, xend = 1, y = 150, yend = 147.5), inherit.aes = FALSE)+
  geom_segment(aes(x = 2, xend = 2, y = 150, yend = 147.5), inherit.aes = FALSE)

stat_test <- compare_means(cfu ~ season, data = cfu, 
                           group.by = "substrate",
                           method = "anova",
                           p.adjust.method = "holm")
stat_test$p.adj=c(0.01,0.96,0.90)


text_annotations <- stat_test %>%
  mutate(y.position = 150,  # Adjust this value as needed
         label = paste0("p = ", p.adj)) %>%
  select(substrate, y.position, label)
text_annotations$fontface=c(2,1,1)
print(text_annotations)

line_annotations <- stat_test %>%
  mutate(y.position = 150,  # Adjust this value as needed
         x.start = as.numeric(substrate) - 0.2, 
         x.end = as.numeric(substrate) + 0.2) %>%
  select(substrate, y.position, x.start, x.end)

tick_annotations <- line_annotations %>%
  mutate(y.tick = y.position - 2.5)

cfu_plot2=ggplot(cfu, aes(x=substrate, y=cfu, fill=season)) +
  geom_boxplot()+
  theme_classic() +
  labs(x='', fill='month', y='CFUs per sample') +
  scale_fill_manual(values=c("#ffc20a", "#0c7bdc")) +
  ylim(0, 160) +
  theme(axis.text.x=element_text(size=12),
		axis.text.y=element_text(size=12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.margin = margin(0.75,0.75,0.75,0.75, 'cm'))+
  geom_text(data = text_annotations, aes(x = substrate, y = y.position, label = label, fontface=fontface), size = 4, 
            position = position_dodge(width = 0.8), vjust = -0.5, inherit.aes = FALSE)+
  geom_text(aes(x = 0.55, y = 160, label = "B"),
            size = 6, fontface =2,  inherit.aes = FALSE)+
  geom_segment(data = line_annotations, aes(x = x.start, xend = x.end, y = y.position, yend = y.position), 
               inherit.aes = FALSE)+
  geom_segment(data = tick_annotations, aes(x = x.start, xend = x.start, y = y.position, yend = y.tick), 
               inherit.aes = FALSE) +
  geom_segment(data = tick_annotations, aes(x = x.end, xend = x.end, y = y.position, yend = y.tick), 
               inherit.aes = FALSE)

cfu_combined <- ggarrange(cfu_plot1,	
                            cfu_plot2,
                            nrow = 1,
                            ncol = 2,
                            widths=c(1.2,3))

###################################################################

### Figure 4 - number of taxa per sample ###

sp=read.csv('sp5.csv')

spm = sp[,-(1:4)]
sp_meta=sp[,1:4]
rownames(spm) <- sp[,1]
sppr=specnumber(spm)
spp=cbind(sp_meta,sppr)

spp$substrate=factor(sp$substrate, levels=c('forest litter', 'mound surface','mound interior'), labels = c('forest litter', 'mound surface','mound interior'))
spp$season=factor(sp$season, levels=c('summer','autumn'), labels=c('August', 'October'))

sp_plot1=ggplot(spp, aes(x=season,y=sppr, fill=season))+
  geom_boxplot()+
  theme_classic()+
  labs(x='', fill='season', y='taxa per sample')+
  scale_fill_manual(values=c("#ffc20a", "#0c7bdc"))+
  ylim(0,12.5)+
  theme(axis.text.x=element_text(size=12),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        plot.margin = margin(0.75,1.5,0.75,0.75, 'cm'))+
  geom_text(aes(x = 1.5, y = 11.6, label = "p = 0.02", fontface=2),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4, inherit.aes = FALSE)+
  geom_text(aes(x = 0.75, y = 12.5, label = "A"),
            size = 6, fontface =2,  inherit.aes = FALSE)+
  geom_segment(aes(x = 1, xend = 2, y = 11.6, yend = 11.6), inherit.aes = FALSE)+
  geom_segment(aes(x = 1, xend = 1, y = 11.6, yend = 11.4), inherit.aes = FALSE)+
  geom_segment(aes(x = 2, xend = 2, y = 11.6, yend = 11.4), inherit.aes = FALSE)
  
stat_test <- compare_means(sppr ~ season, data = spp, 
                           group.by = "substrate",
                           method = "anova",
                           p.adjust.method = "bonferroni")

stat_test$p.adj=c(0.06,0.37,0.98)

text_annotations <- stat_test %>%
  mutate(y.position = 11.6,  # Adjust this value as needed
         label = paste0("p = ", p.adj)) %>%
  select(substrate, y.position, label)

text_annotations$fontface=c(2,1,1)

line_annotations <- stat_test %>%
  mutate(y.position = 11.6,  # Adjust this value as needed
         x.start = as.numeric(substrate) - 0.2, 
         x.end = as.numeric(substrate) + 0.2) %>%
  select(substrate, y.position, x.start, x.end)

tick_annotations <- line_annotations %>%
  mutate(y.tick = y.position - 0.2)

sp_plot2=ggplot(spp, aes(x=substrate, y=sppr, fill=season)) +
  geom_boxplot() +
  theme_classic() +
  labs(x='', fill='month', y='taxa per sample') +
  scale_fill_manual(values=c("#ffc20a", "#0c7bdc")) +
  ylim(0, 12.5) +
  theme(axis.text.x=element_text(size=12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.margin = margin(0.75,0.75,0.75,0.75, 'cm'))+
  geom_text(data = text_annotations, aes(x = substrate, y = y.position, label = label, fontface = fontface), 
            position = position_dodge(width = 0.8), vjust = -0.5,  inherit.aes = FALSE)+
  geom_text(aes(x = 0.55, y = 12.5, label = "B"),
            size = 6, fontface =2,  inherit.aes = FALSE)+
  geom_segment(data = line_annotations, aes(x = x.start, xend = x.end, y = y.position, yend = y.position), 
               inherit.aes = FALSE)+
  geom_segment(data = tick_annotations, aes(x = x.start, xend = x.start, y = y.position, yend = y.tick), 
               inherit.aes = FALSE) +
  geom_segment(data = tick_annotations, aes(x = x.end, xend = x.end, y = y.position, yend = y.tick), 
               inherit.aes = FALSE)
sp_combined <- ggarrange(sp_plot1,
                            sp_plot2,
                            nrow = 1,
                            ncol = 2,
                            widths=c(1.2,3)))

### Figure 5 NMDS ###

kop <- read.csv(file = "sp5.csv", header = TRUE, row.names = 1)
kop.sp=kop[,-(1:3)]
kop.meta=kop[,1:3]

set.seed(123)

sp_distmat <- vegdist(kop.sp, method = "bray")

kop.envfit <- envfit(kop.mds, kop.meta, permutations = 999) # this fits environmental vectors
kop.spp.fit <- envfit(kop.mds, kop.sp, permutations = 9999) # this fits species vectors

site.scrs <- as.data.frame(scores(kop.mds, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, season = kop.meta$season) #add grouping variable "season" to dataframe
site.scrs <- cbind(site.scrs, substrate = kop.meta$substrate) #add grouping variable of cluster grouping to dataframe
site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot

spp.scrs <- as.data.frame(scores(kop.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = kop.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<0.05) #subset data to show species significant at 0.05


substrate_colors <- c("forest litter" = "#A6A6A6", 
                      "mound surface" = "#F4B183", 
                      "mound interior" = "#843C0C")


nmds.plot.kop <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, colour = factor(substrate), shape = factor(season)), size = 2)+
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "substrate", shape = "season")+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))+ 
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + 
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), direction = "both", size = 4, max.overlaps = Inf)+ 
  labs(title = "") +
  scale_color_manual(values = substrate_colors)
  

### function for ellipsess ### 

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
site.scrs$substrate=factor(site.scrs$substrate)
df_ell.kop <- data.frame() 
for(g in levels(site.scrs$substrate)){
  df_ell.kop <- rbind(df_ell.kop, cbind(as.data.frame(with(site.scrs [site.scrs$substrate==g,], veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,substrate=g))
}

# data for labelling the ellipse
NMDS.mean.kop=aggregate(site.scrs[ ,c("NMDS1", "NMDS2")], 
                         list(group = site.scrs$substrate), mean)

nmds.plot.kop +
  geom_path(data = df_ell.kop, aes(x = NMDS1, y = NMDS2, group = substrate, colour= substrate))

#################################################################

### Figure 6 - abundance and prevalance ###

### figure prevalance ###

prev=read.csv('prevalence.csv', sep='\t')
prev$substrate=factor(prev$substrate, levels = c('forest litter', 'mound surface','mound interior'))
prev$taxa=factor(prev$taxa, levels=c('U. isabellina','U. curvata','A. coerulea','A. cylindrospora clade','P. verticillata/humilis clade','M. strictus','U. vinacea','U. angularis', 'E. lignicola'), labels = c('U. isabellina','U. curvata','A. coerulea',"'A. cylindrospora'","'P. verticillata/humilis'",'M. strictus','U. vinacea','U. angularis', 'E. lignicola'))


substrate_colors <- c("forest litter" = "#A6A6A6", 
                      "mound surface" = "#F4B183", 
                      "mound interior" = "#843C0C")

p = ggplot(prev, aes(x = taxa, y = prev, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6)+
  scale_fill_manual(values = substrate_colors)+
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0, 0))+
  theme_classic()+
  labs(x='', y='percentage of samples for whitch taxon was recorded', fill='')+
  theme(axis.text.x = element_text(size=11,  face = "italic", angle = 45, hjust = 1),
        legend.text = element_text(size = 11),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size = 12),
        panel.grid.major.y = element_line(color = "gray"))

###abundace###

abu=read.csv('abundance.csv', sep='\t')
abu$substrate=factor(abu$substrate, levels = c('mound interior', 'mound surface','forest litter'))
abu$taxa=factor(abu$taxa, levels=c('P. verticillata/humilis clade','U. vinacea','M. strictus','E. lignicola','U. angularis','A. coerulea','A. cylindrospora clade','U. curvata','U. isabellina'), labels = c("'P. verticillata/humilis'",'U. vinacea','M. strictus','E. lignicola','U. angularis','A. coerulea',"'A. cylindrospora'",'U. curvata','U. isabellina'))

substrate_colors <- c("forest litter" = "#A6A6A6", 
                      "mound surface" = "#F4B183", 
                      "mound interior" = "#843C0C")

a = ggplot(abu, aes(x = abu, y = taxa, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = substrate_colors)+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,250,50),
                     expand=c(0,0))+
  labs(x='cfu', y='', fill='')+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=12),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size=12, face = "italic"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        panel.grid.major.x = element_line(color = "gray"))

fig6 <- ggarrange(p,	
                  a,
                  nrow = 1
                  ncol = 2,
                  widths=c(1,1),
                  common.legend = TRUE,
                  legend = "bottom"))
