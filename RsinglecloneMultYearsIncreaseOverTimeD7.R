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

yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)
yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
yaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ydupall, FUN=sum)
yaggdataorder <- yaggdata[order(yaggdata$Frequency), ]
ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
ymerge <- rbind(ynotdup, yaggdata)
ymergeorder <- ymerge[order(ymerge$Frequency), ]

zduplicatedAAtf <- duplicated(z$CDR3.AminoAcid)
zrowduplicated <- which(zduplicatedAAtf, zduplicatedAAtf=="TRUE")
zduplicatedAA <- z[zrowduplicated, "CDR3.AminoAcid"]
zdupall <- z[is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
zaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=zdupall, FUN=sum)
zaggdataorder <- zaggdata[order(zaggdata$Frequency), ]
znotdup <- z[!is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
zmerge <- rbind(znotdup, zaggdata)
zmergeorder <- zmerge[order(zmerge$Frequency), ]

xyAA <- rbind(xmergeorder, ymergeorder)
xyduplicatedAAtf <- duplicated(xyAA$CDR3.AminoAcid)
xyrowduplicated <- which(xyduplicatedAAtf, xyduplicatedAAtf=="TRUE")
xyduplicatedAA <- xyAA[xyrowduplicated, "CDR3.AminoAcid"]
yshared <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, xyduplicatedAA), ]
ysharedorder <- yshared[order(yshared$CDR3.AminoAcid), ]
colnames(ysharedorder)[2] <- "Frequency.y"
xshared <- xmergeorder[is.element(xmergeorder$CDR3.AminoAcid, xyduplicatedAA), ]
xsharedorder <- xshared[order(xshared$CDR3.AminoAcid), ]
colnames(xsharedorder)[2] <- "Frequency.x"
xyshared <- merge(xshared, yshared, by=c("CDR3.AminoAcid"))
xysharedinc <- xyshared[xyshared$Frequency.x > xyshared$Frequency.y, ]

xysharedAA <- xysharedinc[[1]]
zshared <- zmergeorder[is.element(zmergeorder$CDR3.AminoAcid, xysharedAA), ]
colnames(zshared)[2] <- "Frequency.z"
zsharedAA <- zshared[[1]]
xyzshared <- xyshared[is.element(xyshared$CDR3.AminoAcid, zsharedAA), ]
xyzshared <- merge(xyzshared, zshared, by=c("CDR3.AminoAcid"))
xyzsharedinc <- xyzshared[xyzshared$Frequency.y > xyzshared$Frequency.z, ]
colnames(xyzsharedinc)[2] <- "Y3D7"
colnames(xyzsharedinc)[3] <- "Y2D7"
colnames(xyzsharedinc)[4] <- "Y1D7"
xyzsharedinc <- xyzsharedinc[c(1, 4, 3, 2)]

TCRmerge <- melt(xyzsharedinc, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1, position=position_dodge(width=0.15)) +
  geom_point(size=2, position=position_dodge(width=0.15)) +
  scale_y_log10(limits=c(0.001, 10), breaks=c(0.001, 0.01, 0.1, 1, 10), label=comma) +
  #geom_line(data=TCRmergeMed, aes(x=Day, y=Median), color="red", size=2) +
  scale_color_manual(values=c("#FF0000", "#0033FF")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold"), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Increasing ICOS+CD38+ Clones of 101") +
  labs(x="Day", y="Clonotypic frequency (%)")

c("#FF99FF", "#FF0000", "#FF6600", "#FFFF00", "#00FF00", "#339900", "#00FFCC", "#0033FF", "#9966CC", "#000000")
write.csv(xyzsharedinc, file="RsinglecloneMultYearsIncreaseOverTime101hi.csv")

a1hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a1lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a3hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a3lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))

b1hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b1lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b3hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b3lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b5hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b5lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b1x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b3x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b5x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

x <- a5hi
y <- a3hi
z <- a1hi