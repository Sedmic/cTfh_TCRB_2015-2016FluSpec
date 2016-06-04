TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(scales)

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

xysharedAA <- xyshared[[1]]
zshared <- zmergeorder[is.element(zmergeorder$CDR3.AminoAcid, xysharedAA), ]
colnames(zshared)[2] <- "Frequency.z"
zsharedAA <- zshared[[1]]
xyzshared <- xyshared[is.element(xyshared$CDR3.AminoAcid, zsharedAA), ]
xyzshared <- merge(xyzshared, zshared, by=c("CDR3.AminoAcid"))
colnames(xyzshared)[2] <- "Y3D7"
colnames(xyzshared)[3] <- "Y2D7"
colnames(xyzshared)[4] <- "Y1D7"

vduplicatedAAtf <- duplicated(v$CDR3.AminoAcid)
vrowduplicated <- which(vduplicatedAAtf, vduplicatedAAtf=="TRUE")
vduplicatedAA <- v[vrowduplicated, "CDR3.AminoAcid"]
vdupall <- v[is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
vaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=vdupall, FUN=sum)
vaggdataorder <- vaggdata[order(vaggdata$Frequency), ]
vnotdup <- v[!is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
vmerge <- rbind(vnotdup, vaggdata)
vmergeorder <- vmerge[order(vmerge$Frequency), ]

wduplicatedAAtf <- duplicated(w$CDR3.AminoAcid)
wrowduplicated <- which(wduplicatedAAtf, wduplicatedAAtf=="TRUE")
wduplicatedAA <- w[wrowduplicated, "CDR3.AminoAcid"]
wdupall <- w[is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
waggdata <- aggregate(Frequency~CDR3.AminoAcid, data=wdupall, FUN=sum)
waggdataorder <- waggdata[order(waggdata$Frequency), ]
wnotdup <- w[!is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
wmerge <- rbind(wnotdup, waggdata)
wmergeorder <- wmerge[order(wmerge$Frequency), ]

xyzsharedAA <- xyzshared[[1]]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, xyzsharedAA), ]
wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, xyzsharedAA), ]
colnames(vshared)[2] <- "Y3D0"
colnames(wshared)[2] <- "Y2D0"
vxyzshared <- merge(xyzshared, vshared, all=TRUE)
vxyzshared[is.na(vxyzshared)] <- 0
vwxyzshared <- merge(vxyzshared, wshared, all=TRUE)
vwxyzshared[is.na(vwxyzshared)] <- 0
#xyzshared[5] <- 0
#colnames(xyzshared)[5] <- "Y3D0"
#xyzshared[6] <- 0
#colnames(xyzshared)[6] <- "Y2D0"
#xyzshared <- xyzshared[c(1, 4, 6, 3, 5, 2)]
#xyzsharedtop <- xyzshared[c(123,40,93,70,30,57,98,96,126,17,116,87), ]
vwxyzshared <- vwxyzshared[c(1, 4, 6, 3, 5, 2)]
vwxyzsharedtop <- vwxyzshared[c(127,68,176,4,14,108,12,118,101,114,131), ]

TCRmerge <- melt(vwxyzsharedtop, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_log10(limits=c(0.001, 10), breaks=c(0.001, 0.01, 0.1, 1, 10), label=comma) +
  scale_color_manual(values=c("#FF99FF", "#990000", "#FF0000", "#FF3300", "#FF9900", "#FFFF00", "#00FF00", "#339900", "#00FFCC", "#0033FF", "#9966CC", "#9900CC", "#663300", "#999999", "#000000")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold"), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Top Shared CXCR5- Clones of 999") +
  labs(x="Day", y="Clonotypic frequency (%)")

a1hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a1lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a2hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a2lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a3hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a3lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a4hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a4lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))

b1hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b1lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b2hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b2lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b3hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b3lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b4hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b4lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b5hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b5lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b1x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b2x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b3x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b4x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b5x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

x <- b5x
y <- b3x
z <- b1x
v <- b4x
w <- b2x