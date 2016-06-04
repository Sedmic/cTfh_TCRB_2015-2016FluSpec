TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(scales)
library(reshape2)

#calculate cumulative frequency of clones in each subset found in previous day 7 subset (for each individual separately)
aa.subsettrack <- function(x, u, v, w, a, b, c) {
  xduplicatedAAtf <- duplicated(x$CDR3.AminoAcid)   #condense rows in x (Y1D7 or Y2D7) with same amino acid sequence
  xrowduplicated <- which(xduplicatedAAtf, xduplicatedAAtf=="TRUE")
  xduplicatedAA <- x[xrowduplicated, "CDR3.AminoAcid"]
  xdupall <- x[is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
  xaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=xdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  xaggdataorder <- xaggdata[order(xaggdata$Frequency), ]
  xnotdup <- x[!is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
  xmerge <- rbind(xnotdup, xaggdata)
  xmergeorder <- xmerge[order(xmerge$Frequency), ]    #list of unique aa sequences in x
  
  
  uduplicatedAAtf <- duplicated(u$CDR3.AminoAcid)   #condense rows in u (Y2D0 or Y3D0 ICOS+CD38+) with same aa sequence
  urowduplicated <- which(uduplicatedAAtf, uduplicatedAAtf=="TRUE")
  uduplicatedAA <- u[urowduplicated, "CDR3.AminoAcid"]
  udupall <- u[is.element(u$CDR3.AminoAcid, uduplicatedAA), ]
  uaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=udupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  uaggdataorder <- uaggdata[order(uaggdata$Frequency), ]
  unotdup <- u[!is.element(u$CDR3.AminoAcid, uduplicatedAA), ]
  umerge <- rbind(unotdup, uaggdata)
  umergeorder <- umerge[order(umerge$Frequency), ]    #list of unique aa sequences in u
  
  vduplicatedAAtf <- duplicated(v$CDR3.AminoAcid)   #condense rows in v (Y2D0 or Y3D0 ICOS-CD38-) with same aa sequence
  vrowduplicated <- which(vduplicatedAAtf, vduplicatedAAtf=="TRUE")
  vduplicatedAA <- v[vrowduplicated, "CDR3.AminoAcid"]
  vdupall <- v[is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
  vaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=vdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  vaggdataorder <- vaggdata[order(vaggdata$Frequency), ]
  vnotdup <- v[!is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
  vmerge <- rbind(vnotdup, vaggdata)
  vmergeorder <- vmerge[order(vmerge$Frequency), ]    #list of unique aa sequences in v
  
  wduplicatedAAtf <- duplicated(w$CDR3.AminoAcid)   #condense rows in w (Y2D0 or Y3D0 CXCR5-) with same aa sequence
  wrowduplicated <- which(wduplicatedAAtf, wduplicatedAAtf=="TRUE")
  wduplicatedAA <- w[wrowduplicated, "CDR3.AminoAcid"]
  wdupall <- w[is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
  waggdata <- aggregate(Frequency~CDR3.AminoAcid, data=wdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  waggdataorder <- waggdata[order(waggdata$Frequency), ]
  wnotdup <- w[!is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
  wmerge <- rbind(wnotdup, waggdata)
  wmergeorder <- wmerge[order(wmerge$Frequency), ]    #list of unique aa sequences in w
  
  aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)   #condense rows in a (Y2D7 or Y3D7 ICOS+CD38+) with same aa sequence
  arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
  aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
  adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
  anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  amerge <- rbind(anotdup, aaggdata)
  amergeorder <- amerge[order(amerge$Frequency), ]    #list of unique aa sequences in a
  
  bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)   #condense rows in b (Y2D7 or Y3D7 ICOS-CD38-) with same aa sequence
  browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
  bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
  bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  baggdataorder <- baggdata[order(baggdata$Frequency), ]
  bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  bmerge <- rbind(bnotdup, baggdata)
  bmergeorder <- bmerge[order(bmerge$Frequency), ]    #list of unique aa sequences in b
  
  cduplicatedAAtf <- duplicated(c$CDR3.AminoAcid)   #condense rows in c (Y2D7 or Y3D7 CXCR5-) with same aa sequence
  crowduplicated <- which(cduplicatedAAtf, cduplicatedAAtf=="TRUE")
  cduplicatedAA <- c[crowduplicated, "CDR3.AminoAcid"]
  cdupall <- c[is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  caggdata <- aggregate(Frequency~CDR3.AminoAcid, data=cdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  caggdataorder <- caggdata[order(caggdata$Frequency), ]
  cnotdup <- c[!is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  cmerge <- rbind(cnotdup, caggdata)
  cmergeorder <- cmerge[order(cmerge$Frequency), ]    #list of unique aa sequences in c
  
  xAA <- xmergeorder[[1]]
  ushared <- umergeorder[is.element(umergeorder$CDR3.AminoAcid, xAA), ]   #determine if any aa sequences in x are present in u
  vshared <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, xAA), ]   #determine if any aa sequences in x are present in v
  wshared <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, xAA), ]   #determine if any aa sequences in x are present in w
  ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, xAA), ]   #determine if any aa sequences in x are present in a
  bshared <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, xAA), ]   #determine if any aa sequences in x are present in b
  cshared <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, xAA), ]   #determine if any aa sequences in x are present in c
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
  
  xsharedu <- merge(xmergeorder, ushared, all=TRUE)
  xsharedu[is.na(xsharedu)] <- 0
  xshareduv <- merge(xsharedu, vshared, all=TRUE)
  xshareduv[is.na(xshareduv)] <- 0
  xshareduvw <- merge(xshareduv, wshared, all=TRUE)
  xshareduvw[is.na(xshareduvw)] <- 0
  xshareduvw <- xshareduvw[c(1, 2, 4, 5, 6)]
  
  xshareda <- merge(xmergeorder, ashared, all=TRUE)
  xshareda[is.na(xshareda)] <- 0
  xshareda$Day <- "7"
  xsharedab <- merge(xshareda, bshared, all=TRUE)
  xsharedab[is.na(xsharedab)] <- 0
  xsharedab$Day <- "7"
  xsharedabc <- merge(xsharedab, cshared, all=TRUE)
  xsharedabc[is.na(xsharedabc)] <- 0
  xsharedabc$Day <- "7"
  xsharedabc <- xsharedabc[c(1, 2, 4, 5, 6)]
  
  xsharedabcuvw <- rbind(xshareduvw, xsharedabc)
  
  hid0 <- sum(xshareduvw$`ICOS+CD38+`)    #calculate cumulative frequency of Y2D0 or Y3D0 ICOS+CD38+ clones also present in x (Y1D7 or Y2D7)
  lod0 <- sum(xshareduvw$`ICOS-CD38-`)    #calculate cumulative frequency of Y2D0 or Y3D0 ICOS-CD38- clones also present in x (Y1D7 or Y2D7)
  xd0 <- sum(xshareduvw$'CXCR5-')       #calculate cumulative frequency of Y2D0 or Y3D0 CXCR5- clones also present in x (Y1D7 or Y2D7)
  hid7 <- sum(xsharedabc$`ICOS+CD38+`)    #calculate cumulative frequency of Y2D0 or Y3D0 ICOS+CD38+ clones also present in x (Y1D7 or Y2D7)
  lod7 <- sum(xsharedabc$`ICOS-CD38-`)    #calculate cumulative frequency of Y2D0 or Y3D0 ICOS-CD38- clones also present in x (Y1D7 or Y2D7)
  xd7 <- sum(xsharedabc$'CXCR5-')       #calculate cumulative frequency of Y2D0 or Y3D0 CXCR5- clones also present in x (Y1D7 or Y2D7)
  
  return (list(hid0, lod0, xd0, hid7, lod7, xd7))
}

aa.subsettrack(b3hi, b4hi, b4lo, b4x, b5hi, b5lo, b5x)

#the values were recorded in an Excel csv and plotted once values for all individuals calculated and recorded

a1hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a1lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a1x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1314 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
a2hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a2lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a2x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
a3hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a3lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a3x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
a4hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a4lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a4x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
a5hi <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
a5lo <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
a5x <- subset(TCRdataMultYearsIn, Subject==101 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

b1hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b1lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b1x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1314 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b2hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b2lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b2x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b3hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b3lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b3x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b4hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b4lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b4x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
b5hi <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
b5lo <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
b5x <- subset(TCRdataMultYearsIn, Subject==999 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

c3hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
c3lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
c3x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
c4hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
c4lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
c4x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
c5hi <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
c5lo <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
c5x <- subset(TCRdataMultYearsIn, Subject==108 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

d3hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
d3lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
d3x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
d4hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
d4lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
d4x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
d5hi <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
d5lo <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
d5x <- subset(TCRdataMultYearsIn, Subject==106 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

e3hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
e3lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
e3x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
e4hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
e4lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
e4x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
e5hi <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
e5lo <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
e5x <- subset(TCRdataMultYearsIn, Subject==117 & Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))

