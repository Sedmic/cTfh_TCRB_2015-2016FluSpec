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

uduplicatedAAtf <- duplicated(u$CDR3.AminoAcid)
urowduplicated <- which(uduplicatedAAtf, uduplicatedAAtf=="TRUE")
uduplicatedAA <- u[urowduplicated, "CDR3.AminoAcid"]
udupall <- u[is.element(u$CDR3.AminoAcid, uduplicatedAA), ]
uaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=udupall, FUN=sum)
uaggdataorder <- uaggdata[order(uaggdata$Frequency), ]
unotdup <- u[!is.element(u$CDR3.AminoAcid, uduplicatedAA), ]
umerge <- rbind(unotdup, uaggdata)
umergeorder <- umerge[order(umerge$Frequency), ]

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

aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)
arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)
aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
amerge <- rbind(anotdup, aaggdata)
amergeorder <- amerge[order(amerge$Frequency), ]

bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)
browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)
baggdataorder <- baggdata[order(baggdata$Frequency), ]
bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
bmerge <- rbind(bnotdup, baggdata)
bmergeorder <- bmerge[order(bmerge$Frequency), ]

cduplicatedAAtf <- duplicated(c$CDR3.AminoAcid)
crowduplicated <- which(cduplicatedAAtf, cduplicatedAAtf=="TRUE")
cduplicatedAA <- c[crowduplicated, "CDR3.AminoAcid"]
cdupall <- c[is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
caggdata <- aggregate(Frequency~CDR3.AminoAcid, data=cdupall, FUN=sum)
caggdataorder <- caggdata[order(caggdata$Frequency), ]
cnotdup <- c[!is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
cmerge <- rbind(cnotdup, caggdata)
cmergeorder <- cmerge[order(cmerge$Frequency), ]


xyAAshared <- xyshared[[1]]
ushared <- umergeorder[is.element(umergeorder$CDR3.AminoAcid, xyAAshared), ]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, xyAAshared), ]
wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, xyAAshared), ]
ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, xyAAshared), ]
bshared <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, xyAAshared), ]
cshared <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, xyAAshared), ]
colnames(ushared)[2] <- "ICOS+CD38+"
ushared$Day <- "0"
colnames(vshared)[2] <- "ICOS-CD38-"
vshared$Day <- "0"
colnames(wshared)[2] <- "CXCR5-"
wshared$Day <- "0"
colnames(ashared)[2] <- "ICOS+CD38+"
ashared$Day <- "7"
colnames(bshared)[2] <- "ICOS-CD38-"
bshared$Day <- "7"
colnames(cshared)[2] <- "CXCR5-"
cshared$Day <- "7"
xysharedu <- merge(xyshared, ushared, all=TRUE)
xysharedu[is.na(xysharedu)] <- 0
xyshareduv <- merge(xysharedu, vshared, all=TRUE)
xyshareduv[is.na(xyshareduv)] <- 0
xyshareduvw <- merge(xyshareduv, wshared, all=TRUE)
xyshareduvw[is.na(xyshareduvw)] <- 0
xyshareduvw <- xyshareduvw[c(1, 6, 4, 5, 7)]

xyshareda <- merge(xyshared, ashared, all=TRUE)
xyshareda[is.na(xyshareda)] <- 0
xyshareda$Day <- "7"
xysharedab <- merge(xyshareda, bshared, all=TRUE)
xysharedab[is.na(xysharedab)] <- 0
xysharedab$Day <- "7"
xysharedabc <- merge(xysharedab, cshared, all=TRUE)
xysharedabc[is.na(xysharedabc)] <- 0
xysharedabc$Day <- "7"
xysharedabc <- xysharedabc[c(1, 2, 5, 6, 7)]

xysharedabcuvw <- rbind(xyshareduvw, xysharedabc)
write.csv(xysharedabcuvw, file="101fluspecificHiToLo.csv")

#TCRmerge <- melt(xysharedabcuvw, id=c("CDR3.AminoAcid", "Day"))
colnames(TCRmerge)[3] <- "Cell"
colnames(TCRmerge)[4] <- "Frequency"


#ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=interaction(CDR3.AminoAcid, Cell), color=Cell)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0, 0.3), breaks=c(0, 0.1, 0.2, 0.3)) +
  scale_color_manual(values=c("red", "gray", "black")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Flu-specific 101 Clones in Year 3") +
  labs(x="Day", y="Clonotypic frequency (%)")

TCRfluspecific <- read.csv("101fluspecificHiToLoCumulativeFreq.csv")
TCRfluspecific$Day <- as.factor(TCRfluspecific$Day)
ggplot(data=TCRfluspecific, aes(x=Day, y=Frequency, group=Cell, color=Cell)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_log10(limits=c(0.001, 100), breaks=c(0.001, 0.01, 0.1, 1, 10, 100), labels=comma) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title=element_text(size=20, face="bold"), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Flu-specific 101 Clones in Year 3") +
  labs(x="Day", y="Cumulative frequency (%)")

sum(xyshareduvw$`ICOS+CD38+`)
sum(xyshareduvw$`ICOS-CD38-`)
sum(xyshareduvw$`CXCR5-`)
sum(xysharedabc$`ICOS+CD38+`)
sum(xysharedabc$`ICOS-CD38-`)
sum(xysharedabc$`CXCR5-`)


a1hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a3hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a4hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a4lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a4x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
a5x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

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

x <- a1hi
y <- a3hi

u <- a4hi
v <- a4lo
w <- a4x
a <- a5hi
b <- a5lo
c <- a5x