TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(scales)
library(reshape2)

#plot the frequency of all clonotypes present at all day 7 timepoints (only day 7 timepoints plotted)
aa.sharedd7 <- function(x, y, z) {
  xduplicatedAAtf <- duplicated(x$CDR3.AminoAcid)   #condense rows in x (Y3D7) with same amino acid sequence
  xrowduplicated <- which(xduplicatedAAtf, xduplicatedAAtf=="TRUE")
  xduplicatedAA <- x[xrowduplicated, "CDR3.AminoAcid"]
  xdupall <- x[is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
  xaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=xdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  xaggdataorder <- xaggdata[order(xaggdata$Frequency), ]
  xnotdup <- x[!is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
  xmerge <- rbind(xnotdup, xaggdata)
  xmergeorder <- xmerge[order(xmerge$Frequency), ]    #list of unique aa sequences in x with frequencies
  
  yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)   #condense rows in y (Y2D7) with same aa sequence
  yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
  yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
  ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
  yaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ydupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  yaggdataorder <- yaggdata[order(yaggdata$Frequency), ]
  ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
  ymerge <- rbind(ynotdup, yaggdata)
  ymergeorder <- ymerge[order(ymerge$Frequency), ]    #list of unique aa sequences in y with frequencies
  
  zduplicatedAAtf <- duplicated(z$CDR3.AminoAcid)   #condense rows in z (Y1D7) with same aa sequence
  zrowduplicated <- which(zduplicatedAAtf, zduplicatedAAtf=="TRUE")
  zduplicatedAA <- z[zrowduplicated, "CDR3.AminoAcid"]
  zdupall <- z[is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
  zaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=zdupall, FUN=sum)    #add up frequencies for rows with same aa sequence
  zaggdataorder <- zaggdata[order(zaggdata$Frequency), ]
  znotdup <- z[!is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
  zmerge <- rbind(znotdup, zaggdata)
  zmergeorder <- zmerge[order(zmerge$Frequency), ]    #list of unique aa sequences in z with frequencies
  
  xyAA <- rbind(xmergeorder, ymergeorder)   #find overlapping sequences between x and y with frequencies for each
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
  
  xysharedAA <- xyshared[[1]]   #find overlapping sequences between x, y, and z (all D7's) with frequencies for each
  zshared <- zmergeorder[is.element(zmergeorder$CDR3.AminoAcid, xysharedAA), ]
  colnames(zshared)[2] <- "Frequency.z"
  zsharedAA <- zshared[[1]]
  xyzshared <- xyshared[is.element(xyshared$CDR3.AminoAcid, zsharedAA), ]
  xyzshared <- merge(xyzshared, zshared, by=c("CDR3.AminoAcid"))
  colnames(xyzshared)[2] <- "Y3D7"
  colnames(xyzshared)[3] <- "Y2D7"
  colnames(xyzshared)[4] <- "Y1D7"
  xyzshared <- xyzshared[c(1, 4, 3, 2)]
  
  TCRmerge <- melt(xyzshared, id=c("CDR3.AminoAcid"))
  colnames(TCRmerge)[2] <- "Day"
  colnames(TCRmerge)[3] <- "Frequency"
  TCRmerge$Median <- NA
  
  tempDataY1D7 <- subset(TCRmerge, Day=="Y1D7")
  a <- median(tempDataY1D7$Frequency)   #calculation of Y1D7 median
  tempDataY2D7 <- subset(TCRmerge, Day=="Y2D7")
  b <- median(tempDataY2D7$Frequency)   #calculation of Y2D7 median
  tempDataY3D7 <- subset(TCRmerge, Day=="Y3D7")
  c <- median(tempDataY3D7$Frequency)   #calculation of Y3D7 median
  
  dataMedian <- as.data.frame(c(a, b, c))
  colnames(dataMedian)[1] <- "Median"
  dataMedian$CDR3.AminoAcid <- NA
  dataMedian$Day <- c("Y1D7", "Y2D7", "Y3D7")
  dataMedian$Frequency <- NA
  dataMedian <- dataMedian[c(2, 3, 4, 1)]
  
  TCRmergeMed <- rbind(TCRmerge, dataMedian)
  
  p.sharedd7 <- ggplot(data=TCRmergeMed, aes(x=Day, y=Frequency, group=CDR3.AminoAcid)) +
    geom_line(size=1, color="#FFCC66", alpha=0.4) +           #color: #FFCC66 (light orange) for ICOS+CD38+, #C0D8C0 (light green) for ICOS-CD38-, #999999 (light gray) for CXCR5-
    geom_point(size=2, color="#FFCC66", alpha=0.4) +           #color: #FFCC66 (light orange) for ICOS+CD38+, #C0D8C0 (light green) for ICOS-CD38-, #999999 (light gray) for CXCR5-
    scale_y_log10(limits=c(0.001, 10), breaks=c(0.001, 0.01, 0.1, 1, 10), label=comma) +
    geom_line(data=TCRmergeMed, aes(x=Day, y=Median), color="darkorange2", size=2) +      #color: darkorange2 (dark orange) for ICOS+CD38+, green4 (dark green) for ICOS-CD38, #333333 (dark gray) for CXCR5-
    theme_bw() +
    theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold"), plot.title=element_text(size=20, face="bold")) +
    ggtitle("Shared ICOS+CD38+ Clones of 101") +
    labs(x="Day", y="Clonotypic frequency (%)")
  
  return(p.sharedd7)
}

aa.sharedd7(a5hi, a3hi, a1hi)

#same function was used to separately graph ICOS+CD38+, ICOS-CD38-, and CXCR5-

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