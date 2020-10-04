library(ggplot2)
library(compiler)
enableJIT(3)
library(ggplot2)
library("Cairo")

data =read.table("align.tab")

data = data[which(data$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10" )),]
data = data[which(data$V3 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10" )),]
data$V1 = factor(data$V1, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10" ))
data$V3 = factor(data$V3, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10" ))

g = ggplot(data=data, aes(x=V4, y=V2))+geom_point(size=0.5, aes(color=V5))+facet_grid(V1~V3, scales="free", space="free" )+ theme_grey(base_size = 60) +
    labs(x="Arabidopsis thaliana", y="Brassica rapa")+
    theme(axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
    axis.text.y = element_text( colour = "black"),
    legend.position='none',
    axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )

CairoPNG(file="align.png",width = 6000, height = 3500)
print(g)
dev.off()

