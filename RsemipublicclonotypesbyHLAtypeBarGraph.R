TCRdataHLA <- read.csv("SemipublicclonotypesbyHLAtypeALL.csv")

HLAhi <- TCRdataHLA[c(1, 2)]

ggplot(data=HLAhi, aes(x=HLA.type, y=ICOShiCD38hi)) +
  geom_bar(stat="identity")
