TCROverlapMultYears <- read.csv("../R Excel DataFrames/SampleOverlapMultipleYears.csv")
library(ggplot2)

TCROverlapDatahihi <- subset(TCROverlapMultYears, Cell=="ICOS+CD38+", select=c(Subject, Days.post.vaccination, Overlap, log.p.value, LastPoint))
TCROverlapDatahihi$Subject <- as.factor(TCROverlapDatahihi$Subject)
ggplot(data=TCROverlapDatahihi, aes(Days.post.vaccination, Overlap, color=Subject, size=log.p.value)) +
  geom_point() +
  geom_line(size=1) +
  scale_y_log10(limits=c(0.002, 1.00)) +
  annotation_logticks(sides="l") +
  theme_bw() +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(0, 50, 100, 150, 200, 250, 300), range=c(0, 8)) +
  guides(color=FALSE) +
  theme(legend.key = element_blank()) +
  geom_text(aes(label=LastPoint), vjust=-2, size=6, position=position_jitter(width=0.5, height=0.5)) +
  ggtitle("TCRseq overlap of ICOS+CD38+ over time") + 
  labs(x="Days post-vaccination", y="Sample overlap score")

TCROverlapDatalolo <- subset(TCROverlapMultYears, Cell=="ICOS-CD38-", select=c(Subject, Days.post.vaccination, Overlap, log.p.value, LastPoint))
TCROverlapDatalolo$Subject <- as.factor(TCROverlapDatalolo$Subject)
ggplot(data=TCROverlapDatalolo, aes(Days.post.vaccination, Overlap, color=Subject, size=log.p.value)) +
  geom_point() +
  geom_line(size=1) +
  scale_y_log10(limits=c(0.002, 1.00)) +
  annotation_logticks(sides="l") +
  theme_bw() +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(0, 50, 100, 150, 200, 250, 300), range=c(8)) +
  guides(color=FALSE) +
  theme(legend.key = element_blank()) +
  geom_text(aes(label=LastPoint), hjust=0.8, vjust=2.5, size=6) +
  ggtitle("TCRseq overlap of ICOS-CD38- over time") + 
  labs(x="Days post-vaccination", y="Sample overlap score")

TCROverlapDatax <- subset(TCROverlapMultYears, Cell=="CXCR5-", select=c(Subject, Days.post.vaccination, Overlap, log.p.value, LastPoint))
TCROverlapDatax$Subject <- as.factor(TCROverlapDatax$Subject)
ggplot(data=TCROverlapDatax, aes(Days.post.vaccination, Overlap, color=Subject, size=log.p.value)) +
  geom_point() +
  geom_line(size=1) +
  scale_y_log10(limits=c(0.002, 1.00)) +
  annotation_logticks(sides="l") +
  theme_bw() +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(0, 50, 100, 150, 200, 250, 300), range=c(8)) +
  guides(color=FALSE) +
  theme(legend.key = element_blank()) +
  geom_text(aes(label=LastPoint), hjust=0.5, vjust=5, size=6) +
  ggtitle("TCRseq overlap of CXCR5- over time") + 
  labs(x="Days post-vaccination", y="Sample overlap score")


TCROverlapDatahihijac <- subset(TCROverlapMultYears, Subset=="ICOS+CD38+", select=c(Subject, Days.post.vaccination, Jaccard.Index, log.p.value, LastPoint))
TCROverlapDatahihijac$Subject <- as.factor(TCROverlapDatahihijac$Subject)
ggplot(data=TCROverlapDatahihijac, aes(Days.post.vaccination, Jaccard.Index, color=Subject, size=log.p.value)) +
  geom_point() +
  geom_line(size=1) +
  scale_x_continuous() +
  scale_y_log10(limits=c(0.0005, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  annotation_logticks(sides="l") +
  theme_bw() +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(0, 50, 100, 150, 200, 250, 300), range=c(0, 8)) +
  guides(color=FALSE) +
  theme(legend.key = element_blank()) +
  geom_text(aes(label=LastPoint), vjust=-4, size=6, position=position_jitter(width=0.5, height=0.5)) +
  ggtitle("TCRseq overlap of ICOS+CD38+ over time") + 
  labs(x="Days post-vaccination", y="Jaccard Index")

TCROverlapDatalolojac <- subset(TCROverlapMultYears, Subset=="ICOS-CD38-", select=c(Subject, Days.post.vaccination, Jaccard.Index, log.p.value, LastPoint))
TCROverlapDatalolojac$Subject <- as.factor(TCROverlapDatalolojac$Subject)
ggplot(data=TCROverlapDatalolojac, aes(Days.post.vaccination, Jaccard.Index, color=Subject, size=log.p.value)) +
  geom_point() +
  geom_line(size=1) +
  scale_x_continuous() +
  scale_y_log10(limits=c(0.0005, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  annotation_logticks(sides="l") +
  theme_bw() +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(0, 50, 100, 150, 200, 250, 300), range=c(7, 8)) +
  guides(color=FALSE) +
  theme(legend.key = element_blank()) +
  geom_text(aes(label=LastPoint), vjust=-4, size=6, position=position_jitter(width=0.5, height=0.5)) +
  ggtitle("TCRseq overlap of ICOS-CD38- over time") + 
  labs(x="Days post-vaccination", y="Jaccard Index")

TCROverlapDataxjac <- subset(TCROverlapMultYears, Subset=="CXCR5-", select=c(Subject, Days.post.vaccination, Jaccard.Index, log.p.value, LastPoint))
TCROverlapDataxjac$Subject <- as.factor(TCROverlapDataxjac$Subject)
ggplot(data=TCROverlapDataxjac, aes(Days.post.vaccination, Jaccard.Index, color=Subject, size=log.p.value)) +
  geom_point() +
  geom_line(size=1) +
  scale_x_continuous() +
  scale_y_log10(limits=c(0.0005, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  annotation_logticks(sides="l") +
  theme_bw() +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(0, 50, 100, 150, 200, 250, 300), range=c(7, 8)) +
  guides(color=FALSE) +
  theme(legend.key = element_blank()) +
  geom_text(aes(label=LastPoint), vjust=-4, size=6, position=position_jitter(width=0.5, height=0.5)) +
  ggtitle("TCRseq overlap of CXCR5- over time") + 
  labs(x="Days post-vaccination", y="Jaccard Index")

TCROverlapData999 <- subset(TCROverlapMultYears, Subject==999, select=c(Subset, Day.arbitrary, Overlap, log.p.value, LastPoint))
TCROverlapData999$Subset <- as.factor(TCROverlapData999$Subset)
ggplot(data=TCROverlapData999, aes(Day.arbitrary, Overlap, fill=Subset, color=Subset, size=log.p.value)) +
  geom_point(pch=21, color="black") +
  geom_line(size=1) +
  scale_x_continuous(limits=c(0,120), breaks=c(0, 37, 60, 97, 120), 
                     labels=c("Yr1 Day7", "Yr2 Day0", "Yr2 Day7", "Yr3 Day0", "Yr3 Day7")) +
  scale_y_log10(limits=c(0.001, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  theme_bw() +
  scale_fill_manual(values=c("grey", "green4", "darkorange2")) +
  scale_color_manual(values=c("grey", "green4", "darkorange2")) +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(50, 100, 150, 200, 250, 300), range=c(0.5, 6)) +
  theme(legend.key=element_blank())+
#  theme(legend.key = element_blank(), legend.title=element_text(size=2, face="bold"), 
#       legend.text=element_text(size=3, face="plain"), axis.title=element_text(size=5, face="plain"), 
#        axis.text=element_text(size=3, face="plain"  , margin=unit(0.5, "cm")), 
#        plot.title=element_text(size=5, face="plain"), axis.ticks.length=unit(0.5, "cm"), 
#        legend.key.size=unit(20, "mm")) +
  ggtitle("TCRseq overlap of 999 over time") + 
  labs(y="OVERLAP WITH Yr1 D7") +
  theme(axis.title.x=element_blank())
ggsave(filename="../R Images/RTCRoverlap999overtime2.png",height=5,width=7,units="in")

TCROverlapData999jac <- subset(TCROverlapMultYears, Subject==999, select=c(Subset, Day.arbitrary, Jaccard.Index, log.p.value, LastPoint))
TCROverlapData999jac$Subset <- as.factor(TCROverlapData999jac$Subset)
ggplot(data=TCROverlapData999jac, aes(Day.arbitrary, Jaccard.Index, fill=Subset, color=Subset, size=log.p.value)) +
  geom_point(pch=21, color="black") +
  geom_line(size=1) +
  scale_x_continuous(limits=c(0,120), breaks=c(0, 37, 60, 97, 120), labels=c("Yrl Day7", "Yr2 Day0", "Yr2 Day7", "Yr3 Day0", "Yr3 Day7")) +
  scale_y_log10(limits=c(0.001, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  theme_bw() +
  scale_fill_manual(values=c("black", "gray", "red")) +
  scale_color_manual(values=c("black", "gray", "red")) +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(50, 100, 150, 200, 250, 300), range=c(2, 16)) +
  theme(legend.key = element_blank(), legend.title=element_text(size=30, face="bold"), legend.text=element_text(size=30, face="bold"), axis.title=element_text(size=50, face="bold"), axis.text=element_text(size=30, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=50, face="bold"), axis.ticks.length=unit(0.1, "cm"), legend.key.size=unit(20, "mm")) +
  ggtitle("TCRseq jaccard of 999 over time") + 
  labs(y="Jaccard Index with Year1 Day7") +
  theme(axis.title.x=element_blank())

TCROverlapData101 <- subset(TCROverlapMultYears, Subject==101, select=c(Subset, Day.arbitrary, Overlap, log.p.value, LastPoint))
TCROverlapData101$Subset <- as.factor(TCROverlapData101$Subset)
ggplot(data=TCROverlapData101, aes(Day.arbitrary, Overlap, fill=Subset, color=Subset, size=log.p.value)) +
  geom_point(pch=21, color="black") +
  geom_line(size=1) +
  scale_x_continuous(limits=c(0,120), breaks=c(0, 37, 60, 97, 120), labels=c("Yrl Day7", "Yr2 Day0", "Yr2 Day7", "Yr3 Day0", "Yr3 Day7")) +
  scale_y_log10(limits=c(0.001, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  theme_bw() +
  scale_fill_manual(values=c("black", "gray", "red")) +
  scale_color_manual(values=c("black", "gray", "red")) +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(50, 100, 150, 200, 250, 300), range=c(2, 16)) +
  theme(legend.key = element_blank(), legend.title=element_text(size=30, face="bold"), legend.text=element_text(size=30, face="bold"), axis.title=element_text(size=50, face="bold"), axis.text=element_text(size=30, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=50, face="bold"), axis.ticks.length=unit(0.1, "cm"), legend.key.size=unit(20, "mm")) +
  ggtitle("TCRseq overlap of 101 over time") + 
  labs(y="Sample overlap score with Year1 Day7") +
  theme(axis.title.x=element_blank())

TCROverlapData101jac <- subset(TCROverlapMultYears, Subject==101, select=c(Subset, Day.arbitrary, Jaccard.Index, log.p.value, LastPoint))
TCROverlapData101jac$Subset <- as.factor(TCROverlapData101jac$Subset)
ggplot(data=TCROverlapData101jac, aes(Day.arbitrary, Jaccard.Index, fill=Subset, color=Subset, size=log.p.value)) +
  geom_point(pch=21, color="black") +
  geom_line(size=1) +
  scale_x_continuous(limits=c(0,120), breaks=c(0, 37, 60, 97, 120), labels=c("Yrl Day7", "Yr2 Day0", "Yr2 Day7", "Yr3 Day0", "Yr3 Day7")) +
  scale_y_log10(limits=c(0.001, 1.00), breaks=c(0.001, 0.01, 0.10, 1.00)) +
  theme_bw() +
  scale_fill_manual(values=c("black", "gray", "red")) +
  scale_color_manual(values=c("black", "gray", "red")) +
  scale_size_continuous(name="-log10(pvalue)", breaks=c(50, 100, 150, 200, 250, 300), range=c(2, 16)) +
  theme(legend.key = element_blank(), legend.title=element_text(size=30, face="bold"), legend.text=element_text(size=30, face="bold"), axis.title=element_text(size=50, face="bold"), axis.text=element_text(size=30, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=50, face="bold"), axis.ticks.length=unit(0.1, "cm"), legend.key.size=unit(20, "mm")) +
  ggtitle("TCRseq overlap of 101 over time") + 
  labs(y="Jaccard Index with Year1 Day7") +
  theme(axis.title.x=element_blank())
