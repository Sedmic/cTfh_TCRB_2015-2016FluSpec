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
xmergetop <- xmergeorder[7110:7112, ]
colnames(xmergetop)[2] <- "Y3D7"

yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)
yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
yaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ydupall, FUN=sum)
yaggdataorder <- yaggdata[order(yaggdata$Frequency), ]
ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
ymerge <- rbind(ynotdup, yaggdata)
ymergeorder <- ymerge[order(ymerge$Frequency), ]
ymergetop <- ymergeorder[11450:11452, ]
colnames(ymergetop)[2] <- "Y2D7"

zduplicatedAAtf <- duplicated(z$CDR3.AminoAcid)
zrowduplicated <- which(zduplicatedAAtf, zduplicatedAAtf=="TRUE")
zduplicatedAA <- z[zrowduplicated, "CDR3.AminoAcid"]
zdupall <- z[is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
zaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=zdupall, FUN=sum)
zaggdataorder <- zaggdata[order(zaggdata$Frequency), ]
znotdup <- z[!is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
zmerge <- rbind(znotdup, zaggdata)
zmergeorder <- zmerge[order(zmerge$Frequency), ]
zmergetop <- zmergeorder[10960:10962, ]
colnames(zmergetop)[2] <- "Y1D7"

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

xAAtop <- xmergetop[[1]]
yshared <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, xAAtop), ]
zshared <- zmergeorder[is.element(zmergeorder$CDR3.AminoAcid, xAAtop), ]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, xAAtop), ]
wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, xAAtop), ]
colnames(yshared)[2] <- "Y2D7"
colnames(zshared)[2] <- "Y1D7"
colnames(vshared)[2] <- "Y3D0"
colnames(wshared)[2] <- "Y2D0"
xshared77 <- merge(xmergetop, yshared, all=TRUE)
xshared77[is.na(xshared77)] <- 0
xshared777 <- merge(xshared77, zshared, all=TRUE)
xshared777[is.na(xshared777)] <- 0
xshared7707 <- merge(xshared777, vshared, all=TRUE)
xshared7707[is.na(xshared7707)] <- 0
xshared70707 <- merge(xshared7707, wshared, all=TRUE)
xshared70707[is.na(xshared70707)] <- 0
xshared70707 <- xshared70707[c(1, 4, 6, 3, 5, 2)]

yAAtop <- ymergetop[[1]]
xshared <- xmergeorder[is.element(xmergeorder$CDR3.AminoAcid, yAAtop), ]
zshared <- zmergeorder[is.element(zmergeorder$CDR3.AminoAcid, yAAtop), ]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, yAAtop), ]
wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, yAAtop), ]
colnames(xshared)[2] <- "Y3D7"
colnames(zshared)[2] <- "Y1D7"
colnames(vshared)[2] <- "Y3D0"
colnames(wshared)[2] <- "Y2D0"
yshared77 <- merge(ymergetop, xshared, all=TRUE)
yshared77[is.na(yshared77)] <- 0
yshared777 <- merge(yshared77, zshared, all=TRUE)
yshared777[is.na(yshared777)] <- 0
yshared7707 <- merge(yshared777, vshared, all=TRUE)
yshared7707[is.na(yshared7707)] <- 0
yshared70707 <- merge(yshared7707, wshared, all=TRUE)
yshared70707[is.na(yshared70707)] <- 0
yshared70707 <- yshared70707[c(1, 4, 6, 2, 5, 3)]

zAAtop <- zmergetop[[1]]
xshared <- xmergeorder[is.element(xmergeorder$CDR3.AminoAcid, zAAtop), ]
yshared <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, zAAtop), ]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, zAAtop), ]
wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, zAAtop), ]
colnames(xshared)[2] <- "Y3D7"
colnames(yshared)[2] <- "Y2D7"
colnames(vshared)[2] <- "Y3D0"
colnames(wshared)[2] <- "Y2D0"
zshared77 <- merge(zmergetop, xshared, all=TRUE)
zshared77[is.na(zshared77)] <- 0
zshared777 <- merge(zshared77, yshared, all=TRUE)
zshared777[is.na(zshared777)] <- 0
zshared7707 <- merge(zshared777, vshared, all=TRUE)
zshared7707[is.na(zshared7707)] <- 0
zshared70707 <- merge(zshared7707, wshared, all=TRUE)
zshared70707[is.na(zshared70707)] <- 0
zshared70707 <- zshared70707[c(1, 2, 6, 4, 5, 3)]

xyzshared <- rbind(xshared70707, yshared70707, zshared70707)
xyzduplicatedTF <- duplicated(xyzshared$CDR3.AminoAcid)
xyzrowduplicated <- which(xyzduplicatedTF, xyzduplicatedTF=="TRUE")
xyzshared <- xyzshared[c(1,2,3,5,6), ]
TCRmerge <- melt(xyzshared, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0, 2.9), breaks=c(0, 0.5, 1, 1.5, 2, 2.5)) +
  scale_color_manual(values=c("#FF0000", "#00FF00", "#00FFCC", "#9966CC", "#000000")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Top ICOS-CD38- Clones of 999") +
  labs(x="Day", y="Clonotypic frequency (%)")


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

x <- b5lo
y <- b3lo
z <- b1lo
v <- b4lo
w <- b2lo