TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))

xduplicatedAAtf <- duplicated(x$CDR3.AminoAcid)
xrowduplicated <- which(xduplicatedAAtf, xduplicatedAAtf=="TRUE")
xduplicatedAA <- x[xrowduplicated, "CDR3.AminoAcid"]
xdupall <- x[is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
xaggdata1 <- aggregate(Frequency~CDR3.AminoAcid, data=xdupall, FUN=sum)
xaggdata2 <- aggregate(Count~CDR3.AminoAcid, data=xdupall, FUN=sum)
xaggdatamerge <- merge(xaggdata1, xaggdata2, by=c("CDR3.AminoAcid"))
xaggdataorder <- xaggdatamerge[order(xaggdatamerge$Frequency), ]
xnotdup <- x[!is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
xmerge <- rbind(xnotdup, xaggdatamerge)
xmergeorder <- xmerge[order(xmerge$Frequency), ]

yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)
yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
yaggdata1 <- aggregate(Frequency~CDR3.AminoAcid, data=ydupall, FUN=sum)
yaggdata2 <- aggregate(Count~CDR3.AminoAcid, data=ydupall, FUN=sum)
yaggdatamerge <- merge(yaggdata1, yaggdata2, by=c("CDR3.AminoAcid"))
yaggdataorder <- yaggdatamerge[order(yaggdatamerge$Frequency), ]
ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
ymerge <- rbind(ynotdup, yaggdatamerge)
ymergeorder <- ymerge[order(ymerge$Frequency), ]

xyAA <- rbind(xmergeorder, ymergeorder)
xyduplicatedAAtf <- duplicated(xyAA$CDR3.AminoAcid)
xyrowduplicated <- which(xyduplicatedAAtf, xyduplicatedAAtf=="TRUE")
xyduplicatedAA <- xyAA[xyrowduplicated, "CDR3.AminoAcid"]
yshared <- ymergeorder[is.element(ymergeorder$CDR3.AminoAcid, xyduplicatedAA), ]
ysharedorder <- yshared[order(yshared$CDR3.AminoAcid), ]
colnames(ysharedorder)[2] <- "Frequency.y"
colnames(ysharedorder)[3] <- "Count.y"
xshared <- xmergeorder[is.element(xmergeorder$CDR3.AminoAcid, xyduplicatedAA), ]
xsharedorder <- xshared[order(xshared$CDR3.AminoAcid), ]
colnames(xsharedorder)[2] <- "Frequency.x"
colnames(xsharedorder)[3] <- "Count.x"
xyshared <- merge(xshared, yshared, by=c("CDR3.AminoAcid"))
xyshared$Group <- "Shared"


shared <- nrow(xyshared)
shared7 <- sum(xyshared$Frequency.x)
shared0 <- sum(xyshared$Frequency.x)

a1hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a2hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a3hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a1lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a2lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a3lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

b1hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b1lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b1x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

c3hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c3lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c3x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

d3hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d3lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d3x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

e3hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e3lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e3x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

x <- b5x
y <- b3x
