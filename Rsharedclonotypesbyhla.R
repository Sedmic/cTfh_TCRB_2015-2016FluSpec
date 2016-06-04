TCRDatahla <- read.csv("Shared0401hi.csv")

TCRdata0401 <- TCRDatahla[c(1, 2, 4, 5, 6, 7, 8)]
TCRdata0401$Subject <- as.factor(TCRdata0401$Subject)
TCRmerge <- melt(TCRdata0401, id=c("CDR3.AminoAcid", "Subject"))
colnames(TCRmerge)[3] <- "Day"
colnames(TCRmerge)[4] <- "Frequency"


ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=interaction(CDR3.AminoAcid, Subject), color=Subject, linetype=CDR3.AminoAcid, shape=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=5) +
  scale_y_continuous(limits=c(0, 0.2), breaks=c(0, 0.05, 0.1, 0.15, 0.2)) +
  scale_color_manual(values=c("#FF0000", "#0033FF")) +
  scale_linetype_manual(values=c("solid","dashed")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Shared HLA-DR0401 ICOS+CD38+ Clonotypes") +
  labs(x="Day", y="Clonotypic frequency (%)")