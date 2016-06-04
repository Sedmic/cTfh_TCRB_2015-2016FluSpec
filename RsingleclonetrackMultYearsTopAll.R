TCRdata <- read.csv("top2overtime999x.csv")
TCRdata <- TCRdata[c(1, 3, 4, 5, 6, 7)]

TCRmerge <- melt(TCRdata, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0, 2.9), breaks=c(0, 0.5, 1, 1.5, 2, 2.5)) +
  scale_color_manual(values=c("#FF0000", "#00FF00", "#0033FF", "#000000")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=25, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=20, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=25, face="bold")) +
  ggtitle("Top CXCR5- Clones of 999") +
  labs(x="Day", y="Clonotypic frequency (%)")
