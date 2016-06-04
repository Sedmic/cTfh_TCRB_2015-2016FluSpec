TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdata1415In <- subset(TCRdataMultYears, Year==1415 & SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(entropy)
library(ineq)

aa.overlap <- function(x) {
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
  
  xcount <- unlist(xmergeorder$Count)
  xfreqs <- ((xcount)/(sum(xcount)))
  
  entropyx <- entropy.empirical(xfreqs, unit="log2")
  entropyx <- signif(entropyx, digits=6)
  entropynormx <- ((entropyx)/(log2(nrow(xmergeorder))))
  clonalityx <- (1 - (entropynormx))
  clonalityx <- signif(clonalityx, digits=6)
  ginix <- ineq(xfreqs, type="Gini")
  ginix <- signif(ginix, digits=6)
  
  n <- xcount
  m <- (n*(n-1))
  mn <- sum(m)
  N <- sum(n)
  Dx <- ((mn)/(N*(N-1)))
  simpsonsx <- (1 - (Dx))
  simpsonsx <- signif(simpsonsx, digits=6)
  
  return(list(entropyx, clonalityx, ginix, simpsonsx))
}

a2hi <- subset(TCRdata1415In, Subject==101 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a3hi <- subset(TCRdata1415In, Subject==101 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a2lo <- subset(TCRdata1415In, Subject==101 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a3lo <- subset(TCRdata1415In, Subject==101 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

b2hi <- subset(TCRdata1415In, Subject==999 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3hi <- subset(TCRdata1415In, Subject==999 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2lo <- subset(TCRdata1415In, Subject==999 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3lo <- subset(TCRdata1415In, Subject==999 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b2x <- subset(TCRdata1415In, Subject==999 & Year==1415 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b3x <- subset(TCRdata1415In, Subject==999 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

c2hi <- subset(TCRdata1415In, Subject==108 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c3hi <- subset(TCRdata1415In, Subject==108 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c2lo <- subset(TCRdata1415In, Subject==108 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c3lo <- subset(TCRdata1415In, Subject==108 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c2x <- subset(TCRdata1415In, Subject==108 & Year==1415 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c3x <- subset(TCRdata1415In, Subject==108 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

d2hi <- subset(TCRdata1415In, Subject==106 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d3hi <- subset(TCRdata1415In, Subject==106 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d2lo <- subset(TCRdata1415In, Subject==106 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d3lo <- subset(TCRdata1415In, Subject==106 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d2x <- subset(TCRdata1415In, Subject==106 & Year==1415 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d3x <- subset(TCRdata1415In, Subject==106 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

e2hi <- subset(TCRdata1415In, Subject==117 & Year==1415 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e3hi <- subset(TCRdata1415In, Subject==117 & Year==1415 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e2lo <- subset(TCRdata1415In, Subject==117 & Year==1415 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e3lo <- subset(TCRdata1415In, Subject==117 & Year==1415 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e2x <- subset(TCRdata1415In, Subject==117 & Year==1415 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e3x <- subset(TCRdata1415In, Subject==117 & Year==1415 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

aa.overlap(a2hi)