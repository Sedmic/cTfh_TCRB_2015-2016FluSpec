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
xnotdup <- xnotdup[c(1,2)]
xmerge <- rbind(xnotdup, xaggdata)
xmergeorder <- xmerge[order(xmerge$Frequency), ]

yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)
yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
yaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ydupall, FUN=sum)
yaggdataorder <- yaggdata[order(yaggdata$Frequency), ]
ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
ynotdup <- ynotdup[c(1,2)]
ymerge <- rbind(ynotdup, yaggdata)
ymergeorder <- ymerge[order(ymerge$Frequency), ]

zduplicatedAAtf <- duplicated(z$CDR3.AminoAcid)
zrowduplicated <- which(zduplicatedAAtf, zduplicatedAAtf=="TRUE")
zduplicatedAA <- z[zrowduplicated, "CDR3.AminoAcid"]
zdupall <- z[is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
zaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=zdupall, FUN=sum)
zaggdataorder <- zaggdata[order(zaggdata$Frequency), ]
znotdup <- z[!is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
znotdup <- znotdup[c(1,2)]
zmerge <- rbind(znotdup, zaggdata)
zmergeorder <- zmerge[order(zmerge$Frequency), ]

vduplicatedAAtf <- duplicated(v$CDR3.AminoAcid)
vrowduplicated <- which(vduplicatedAAtf, vduplicatedAAtf=="TRUE")
vduplicatedAA <- v[vrowduplicated, "CDR3.AminoAcid"]
vdupall <- v[is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
vaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=vdupall, FUN=sum)
vaggdataorder <- vaggdata[order(vaggdata$Frequency), ]
vnotdup <- v[!is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
vnotdup <- vnotdup[c(1,2)]
vmerge <- rbind(vnotdup, vaggdata)
vmergeorder <- vmerge[order(vmerge$Frequency), ]

wduplicatedAAtf <- duplicated(w$CDR3.AminoAcid)
wrowduplicated <- which(wduplicatedAAtf, wduplicatedAAtf=="TRUE")
wduplicatedAA <- w[wrowduplicated, "CDR3.AminoAcid"]
wdupall <- w[is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
waggdata <- aggregate(Frequency~CDR3.AminoAcid, data=wdupall, FUN=sum)
waggdataorder <- waggdata[order(waggdata$Frequency), ]
wnotdup <- w[!is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
wnotdup <- wnotdup[c(1,2)]
wmerge <- rbind(wnotdup, waggdata)
wmergeorder <- wmerge[order(wmerge$Frequency), ]

aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)
arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)
aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
anotdup <- anotdup[c(1,2)]
amerge <- rbind(anotdup, aaggdata)
amergeorder <- amerge[order(amerge$Frequency), ]

bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)
browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)
baggdataorder <- baggdata[order(baggdata$Frequency), ]
bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
bnotdup <- bnotdup[c(1,2)]
bmerge <- rbind(bnotdup, baggdata)
bmergeorder <- bmerge[order(bmerge$Frequency), ]

axshared <- rbind(xmergeorder, amergeorder)
axduplicatedTF <- duplicated(axshared$CDR3.AminoAcid)
axrowduplicated <- which(axduplicatedTF, axduplicatedTF=="TRUE")
axduplicatedAA <- axshared[axrowduplicated, "CDR3.AminoAcid"]
ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, axduplicatedAA), ]
asharedorder <- ashared[order(ashared$CDR3.AminoAcid), ]
colnames(asharedorder)[2] <- "Frequency.a"
xshared <- xmergeorder[is.element(xmergeorder$CDR3.AminoAcid, axduplicatedAA), ]
xsharedorder <- xshared[order(xshared$CDR3.AminoAcid), ]
colnames(xsharedorder)[2] <- "Frequency.x"
axshared <- merge(xshared, ashared, by=c("CDR3.AminoAcid"))
colnames(axshared)[2] <- "Frequency.x"
colnames(axshared)[3] <- "Frequency.a"

axAA <- axshared[[1]]
yshared <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, axAA), ]
zshared <- zmergeorder[is.element(zmergeorder$CDR3.AminoAcid, axAA), ]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, axAA), ]
wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, axAA), ]
colnames(yshared)[2] <- "Y2D7"
colnames(zshared)[2] <- "Y1D7"
colnames(vshared)[2] <- "Y3D0"
colnames(wshared)[2] <- "Y2D0"
axsharedy <- merge(axshared, yshared, all=TRUE)
axsharedy[is.na(axsharedy)] <- 0
axsharedyz <- merge(axsharedy, zshared, all=TRUE)
axsharedyz[is.na(axsharedyz)] <- 0
axsharedyzv <- merge(axsharedyz, vshared, all=TRUE)
axsharedyzv[is.na(axsharedyzv)] <- 0
axsharedyzvw <- merge(axsharedyzv, wshared, all=TRUE)
axsharedyzvw[is.na(axsharedyzvw)] <- 0
colnames(axsharedyzvw)[2] <- "Y3D7"
colnames(axsharedyzvw)[3] <- "a"
axsharedyzvw <- axsharedyzvw[c(1, 5, 7, 4, 6, 2)]
TCRmerge <- melt(axsharedyzvw, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

#axAA <- axshared[[1]]
yshared <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, axAA), ]
vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, axAA), ]
colnames(yshared)[2] <- "Y1D7"
colnames(vshared)[2] <- "Y2D0"
axsharedy <- merge(axshared, yshared, all=TRUE)
axsharedy[is.na(axsharedy)] <- 0
axsharedyv <- merge(axsharedy, vshared, all=TRUE)
axsharedyv[is.na(axsharedyv)] <- 0
colnames(axsharedyv)[2] <- "Y2D7"
colnames(axsharedyv)[3] <- "a"
axsharedyv <- axsharedyv[c(1, 4, 5, 2)]
TCRmerge <- melt(axsharedyv, id=c("CDR3.AminoAcid"))
colnames(TCRmerge)[2] <- "Day"
colnames(TCRmerge)[3] <- "Frequency"

ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0, 2.9), breaks=c(0, 0.5, 1, 1.5, 2, 2.5)) +
  scale_color_manual(values=c("#FF0000", "#0033FF")) +
  theme_bw() +
  theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold", margin=unit(0.5, "cm")), plot.title=element_text(size=20, face="bold")) +
  ggtitle("Shared 0401 ICOS+CD38+ in 999") +
  labs(x="Day", y="Clonotypic frequency (%)")
write.csv(axsharedyzvw, file="Shared0401hiIn999.csv")

a1hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a2hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a3hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a1lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a2lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a3lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

b1hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b1lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b1x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b6x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day=="6wk", select=c(CDR3.AminoAcid, Frequency, Count))

c3hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c3lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c3x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c6x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day=="6wk", select=c(CDR3.AminoAcid, Frequency, Count))

d3hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d3lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d3x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

e3hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e3lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e3x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

f4hi <- subset(TCRdataMultYearsIn, Subject==102 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
f5hi <- subset(TCRdataMultYearsIn, Subject==102 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
f4lo <- subset(TCRdataMultYearsIn, Subject==102 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
f5lo <- subset(TCRdataMultYearsIn, Subject==102 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

h4hi <- subset(TCRdataMultYearsIn, Subject==113 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
h5hi <- subset(TCRdataMultYearsIn, Subject==113 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
h4lo <- subset(TCRdataMultYearsIn, Subject==113 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
h5lo <- subset(TCRdataMultYearsIn, Subject==113 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))


x <- f5hi
y <- f3hi
v <- f4hi

a <- d5hi
b <- e5hi