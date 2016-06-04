TCRdata1516 <- read.csv("../Samples_1516Data_AsOf02242016.csv")
TCRdata1516In <- subset(TCRdata1516, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(entropy)
library(ineq)

#calculating clonality and gini index (for any given subset, time, and individual)
aa.overlap <- function(x) {
  xduplicatedAAtf <- duplicated(x$CDR3.AminoAcid)   #condense rows with same amino acid sequence but different nucleotide sequences
  xrowduplicated <- which(xduplicatedAAtf, xduplicatedAAtf=="TRUE")
  xduplicatedAA <- x[xrowduplicated, "CDR3.AminoAcid"]
  xdupall <- x[is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
  xaggdata1 <- aggregate(Frequency~CDR3.AminoAcid, data=xdupall, FUN=sum)   #add up frequencies for rows with same amino acid sequence
  xaggdata2 <- aggregate(Count~CDR3.AminoAcid, data=xdupall, FUN=sum)   #add up counts for rows with same amino acid sequence
  xaggdatamerge <- merge(xaggdata1, xaggdata2, by=c("CDR3.AminoAcid"))
  xaggdataorder <- xaggdatamerge[order(xaggdatamerge$Frequency), ]
  xnotdup <- x[!is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
  xmerge <- rbind(xnotdup, xaggdatamerge)
  xmergeorder <- xmerge[order(xmerge$Frequency), ]
  
  xcount <- unlist(xmergeorder$Count)
  xfreqs <- ((xcount)/(sum(xcount)))
  
  entropyx <- entropy.empirical(xfreqs, unit="log2")    #entropy score calculation
  entropyx <- signif(entropyx, digits=6)
  entropynormx <- ((entropyx)/(log2(nrow(xmergeorder))))    #normalized entropy score calculation
  clonalityx <- (1 - (entropynormx))    #clonality score calculation
  clonalityx <- signif(clonalityx, digits=6)
  ginix <- ineq(xfreqs, type="Gini")    #gini index calculation
  ginix <- signif(ginix, digits=6)
  
  return(list(clonalityx, ginix))
}

aa.overlap(a4hi)

#values were recorded in an Excel csv file until values for all individuals were calculated and recorded, then plotted separately

a4hi <- subset(TCRdata1516In, Subject==101 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5hi <- subset(TCRdata1516In, Subject==101 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4lo <- subset(TCRdata1516In, Subject==101 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5lo <- subset(TCRdata1516In, Subject==101 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
a4x <- subset(TCRdata1516In, Subject==101 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
a5x <- subset(TCRdata1516In, Subject==101 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

b4hi <- subset(TCRdata1516In, Subject==999 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5hi <- subset(TCRdata1516In, Subject==999 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4lo <- subset(TCRdata1516In, Subject==999 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5lo <- subset(TCRdata1516In, Subject==999 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
b4x <- subset(TCRdata1516In, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
b5x <- subset(TCRdata1516In, Subject==999 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

c4hi <- subset(TCRdata1516In, Subject==108 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5hi <- subset(TCRdata1516In, Subject==108 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4lo <- subset(TCRdata1516In, Subject==108 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5lo <- subset(TCRdata1516In, Subject==108 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
c4x <- subset(TCRdata1516In, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
c5x <- subset(TCRdata1516In, Subject==108 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

d4hi <- subset(TCRdata1516In, Subject==106 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5hi <- subset(TCRdata1516In, Subject==106 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4lo <- subset(TCRdata1516In, Subject==106 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5lo <- subset(TCRdata1516In, Subject==106 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
d4x <- subset(TCRdata1516In, Subject==106 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
d5x <- subset(TCRdata1516In, Subject==106 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

e4hi <- subset(TCRdata1516In, Subject==117 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5hi <- subset(TCRdata1516In, Subject==117 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4lo <- subset(TCRdata1516In, Subject==117 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5lo <- subset(TCRdata1516In, Subject==117 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
e4x <- subset(TCRdata1516In, Subject==117 & Year==1516 & Cell=="CXCR5-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
e5x <- subset(TCRdata1516In, Subject==117 & Year==1516 & Cell=="CXCR5-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

f4hi <- subset(TCRdata1516In, Subject==102 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
f5hi <- subset(TCRdata1516In, Subject==102 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
f4lo <- subset(TCRdata1516In, Subject==102 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
f5lo <- subset(TCRdata1516In, Subject==102 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

g4hi <- subset(TCRdata1516In, Subject==110 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
g5hi <- subset(TCRdata1516In, Subject==110 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
g4lo <- subset(TCRdata1516In, Subject==110 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
g5lo <- subset(TCRdata1516In, Subject==110 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

h4hi <- subset(TCRdata1516In, Subject==113 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
h5hi <- subset(TCRdata1516In, Subject==113 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
h4lo <- subset(TCRdata1516In, Subject==113 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
h5lo <- subset(TCRdata1516In, Subject==113 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

i4hi <- subset(TCRdata1516In, Subject==111 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
i5hi <- subset(TCRdata1516In, Subject==111 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
i4lo <- subset(TCRdata1516In, Subject==111 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
i5lo <- subset(TCRdata1516In, Subject==111 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

j4hi <- subset(TCRdata1516In, Subject==112 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
j5hi <- subset(TCRdata1516In, Subject==112 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
j4lo <- subset(TCRdata1516In, Subject==112 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
j5lo <- subset(TCRdata1516In, Subject==112 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

k4hi <- subset(TCRdata1516In, Subject==114 & Year==1516 & Cell=="ICOS+CD38+" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
k5hi <- subset(TCRdata1516In, Subject==114 & Year==1516 & Cell=="ICOS+CD38+" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))
k4lo <- subset(TCRdata1516In, Subject==114 & Year==1516 & Cell=="ICOS-CD38-" & Day==0, select=c(CDR3.AminoAcid, Frequency, Count))
k5lo <- subset(TCRdata1516In, Subject==114 & Year==1516 & Cell=="ICOS-CD38-" & Day==7, select=c(CDR3.AminoAcid, Frequency, Count))

matrixResults <- data.frame(entropy=0,clonality=0,gini=0,simpsons=0)
matrixResults <- rbind(matrixResults,aa.overlap(a4hi)); matrixResults <- rbind(matrixResults,aa.overlap(a5hi))
matrixResults <- rbind(matrixResults,aa.overlap(a4lo)); matrixResults <- rbind(matrixResults,aa.overlap(a5lo))
matrixResults <- rbind(matrixResults,aa.overlap(a4x)); matrixResults <- rbind(matrixResults,aa.overlap(a5x))

matrixResults <- rbind(matrixResults,aa.overlap(b4hi)); matrixResults <- rbind(matrixResults,aa.overlap(b5hi))
matrixResults <- rbind(matrixResults,aa.overlap(b4lo)); matrixResults <- rbind(matrixResults,aa.overlap(b5lo))
matrixResults <- rbind(matrixResults,aa.overlap(b4x)); matrixResults <- rbind(matrixResults,aa.overlap(b5x))

