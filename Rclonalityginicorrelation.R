TCRStatsAll <- read.csv("SampleStats1516.csv")
TCRStats <- TCRStatsAll[c(1,3,4,6,7)]
library(ggplot2)

TCRStatsD7Hi <- subset(TCRStats, Subset=="ICOS+CD38+" & Day==7, select=c(Subject, Clonality, Gini))
cor1 <- cor(TCRStatsD7Hi$Clonality, TCRStatsD7Hi$Gini, method="pearson")
cor1 <- signif(cor1, digits=3)
ols <- lm(Gini~Clonality, data=TCRStatsD7Hi)
cor.test(TCRStatsD7Hi$Clonality, TCRStatsD7Hi$Gini, paired=TRUE)

ggplot(TCRStatsD7Hi, aes(Clonality, Gini)) +
  geom_point(pch=21, fill="orange", size=7) +
  scale_x_continuous(limits=c(0, 0.050)) +
  scale_y_continuous(limits=c(0, 0.4)) +
  geom_abline(intercept= ols$coefficients[1], slope= ols$coefficients[2], color="black", size=1) +
  geom_text(aes(0.04, 0.1, label=paste("r=", cor1), size=6)) +
  geom_text(aes(0.04, 0.075, label=paste("p=", 0.00000006502), size=6)) +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=20, face="bold"), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("ICOS+CD38+ cTfh at Day 7") +
  labs(x="Clonality at day 7", y="Gini index at day 7")


TCRStatsD7Hi <- subset(TCRStats, Subset=="ICOS+CD38+" & Day==7, select=c(Subject, Clonality, Gini))
colnames(TCRStatsD7Hi)[2] <- "Clonality D7"
colnames(TCRStatsD7Hi)[3] <- "Gini D7"
TCRStatsD0Hi <- subset(TCRStats, Subset=="ICOS+CD38+" & Day==0, select=c(Subject, Clonality, Gini))
colnames(TCRStatsD0Hi)[2] <- "Clonality D0"
colnames(TCRStatsD0Hi)[3] <- "Gini D0"
TCRStatsHi <- merge(TCRStatsD0Hi, TCRStatsD7Hi, by=c("Subject"))
ClonalityHi <- TCRStatsHi[c(1, 2, 4)]
GiniHi <- TCRStatsHi[c(1, 3, 5)]

ClonalityHi$DeltaClonality <- ((ClonalityHi$`Clonality D7`)/(ClonalityHi$`Clonality D0`))
GiniHi$DeltaGini <- ((GiniHi$`Gini D7`)/(GiniHi$`Gini D0`))
DeltaHi <- merge(ClonalityHi, GiniHi, by=c("Subject"))
DeltaHi2 <- DeltaHi[c(1, 4, 7)]
cor1 <- cor(DeltaHi2$DeltaClonality, DeltaHi2$DeltaGini, method="pearson")
cor1 <- signif(cor1, digits=3)
ols <- lm(DeltaGini~DeltaClonality, data=DeltaHi2)
cor.test(DeltaHi2$DeltaClonality, DeltaHi2$DeltaGini, paired=TRUE)

ggplot(DeltaHi2, aes(DeltaClonality, DeltaGini)) +
  geom_point(pch=21, fill="orange", size=7) +
  scale_x_continuous(limits=c(0, 10)) +
  scale_y_continuous(limits=c(0, 4)) +
  geom_abline(intercept= ols$coefficients[1], slope= ols$coefficients[2], color="black", size=1) +
  geom_text(aes(7.5, 1, label=paste("r=", cor1), size=6)) +
  geom_text(aes(7.5, 0.75, label=paste("p=", 0.00000451), size=6)) +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=20, face="bold"), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("ICOS+CD38+ cTfh") +
  labs(x="DeltaClonality", y="DeltaGini index")
