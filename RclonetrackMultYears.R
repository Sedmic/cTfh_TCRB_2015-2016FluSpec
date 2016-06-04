TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(scales)
library(reshape2)

xduplicatedAAtf <- duplicated(x$CDR3.AminoAcid)
xrowduplicated <- which(xduplicatedAAtf, xduplicatedAAtf=="TRUE")
xduplicatedAA <- x[xrowduplicated, "CDR3.AminoAcid"]
xdupall <- x[is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
xaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=xdupall, FUN=sum)
xaggdataorder <- xaggdata[order(xaggdata$Frequency), ]
xnotdup <- x[!is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
xmerge <- rbind(xnotdup, xaggdata)
xmergeorder <- xmerge[order(xmerge$Frequency), ]
xmergetop5 <- xmergeorder[7438:7442, ]
colnames(xmergetop5)[2] <- "Day0"

yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)
yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
yaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ydupall, FUN=sum)
yaggdataorder <- yaggdata[order(yaggdata$Frequency), ]
ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
ymerge <- rbind(ynotdup, yaggdata)
ymergeorder <- ymerge[order(ymerge$Frequency), ]
ymergetop5 <- ymergeorder[8757:8761, ]
colnames(ymergetop5)[2] <- "Day7"

xAAtop5 <- xmergetop5[[1]]
xfreqd0 <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, xAAtop5), ]
xmergetop5d07 <- merge(xmergetop5, xfreqd0, all=TRUE)
xmergetop5d07[is.na(xmergetop5d07)] <- 0
colnames(xmergetop5d07)[3] <- "Day7"
yAAtop5 <- ymergetop5[[1]]
yfreqd0 <- xmergeorder[is.element(xmergeorder$CDR3.AminoAcid, yAAtop5), ]
ymergetop5d07 <- merge(ymergetop5, yfreqd0, all=TRUE)
ymergetop5d07[is.na(ymergetop5d07)] <- 0
colnames(ymergetop5d07)[3] <- "Day0"
ymergetop5d07 <- ymergetop5d07[c(1,3,2)]

xyalltop <- rbind(xmergetop5d07, ymergetop5d07)
xyallduptf <- duplicated(xyalltop$CDR3.AminoAcid)
xyallrowdup <- which(xyallduptf, xyallduptf=="TRUE")
xyalltop <- xyalltop[c(1,2,3,4,5,7,8,10), ]

TCRmerge <- melt(xyalltop, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

c("#FF99FF", "#FF0000", "#FF6600", "#FFFF00", "#00FF00", "#339900", "#00FFCC", "#0033FF", "#9966CC", "#000000")
ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0, 2.9), breaks=c(0, 0.5, 1, 1.5, 2, 2.5)) +
  scale_color_manual(values=c("#FF99FF", "#FF0000", "#FF6600", "#00FF00", "#00FFCC", "#0033FF", "#9966CC", "#000000")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=25, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=20, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=25, face="bold")) +
  ggtitle("Top 5 CXCR5- Clones of 101") +
  labs(x="Day", y="Clonotypic frequency (%)")


a4hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a4lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a4x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
a5x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

x <- a4x
y <- a5x