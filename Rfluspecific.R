TCRfluspecific <- read.csv("FluspecificHiToLoCumulativeFreq101999.csv")
TCRfluspecific$Subject <- as.factor(TCRfluspecific$Subject)
library(ggplot2)
library(scales)

ggplot(data=TCRfluspecific, aes(x=Day_Arbitrary, y=Frequency, group=interaction(Subject,Cell), color=Cell, shape=Subject)) +
  geom_line(size=1) +
  geom_point(size=5) +
  scale_x_continuous(labels=c("Day 0", "Day 7", "Day 0", "Day 7", "Day 0", "Day 7")) +
  scale_y_log10(limits=c(0.001, 100), breaks=c(0.001, 0.01, 0.1, 1, 10, 100), labels=comma) +
  scale_color_manual(values=c("#999999", "#339900", "#FF3300")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Flu-specific Clones in Year 3") +
  labs(x="Day", y="Cumulative frequency (%)")
