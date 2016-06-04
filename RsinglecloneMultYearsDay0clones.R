TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(scales)
library(reshape2)

aa.sharedd7zerod0 <- function(x, y, z, v, w, a, b, c, d, e) {
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
  vwxyzshared <- vwxyzshared[c(1, 4, 6, 3, 5, 2)]
  
  TCRd0zero <- vwxyzshared[vwxyzshared$Y2D0==0 | vwxyzshared$Y3D0==0, ]
  
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
  
  dduplicatedAAtf <- duplicated(d$CDR3.AminoAcid)
  drowduplicated <- which(dduplicatedAAtf, dduplicatedAAtf=="TRUE")
  dduplicatedAA <- d[drowduplicated, "CDR3.AminoAcid"]
  ddupall <- d[is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  daggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ddupall, FUN=sum)
  daggdataorder <- daggdata[order(daggdata$Frequency), ]
  dnotdup <- d[!is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  dmerge <- rbind(dnotdup, daggdata)
  dmergeorder <- dmerge[order(dmerge$Frequency), ]
  
  eduplicatedAAtf <- duplicated(e$CDR3.AminoAcid)
  erowduplicated <- which(eduplicatedAAtf, eduplicatedAAtf=="TRUE")
  eduplicatedAA <- e[erowduplicated, "CDR3.AminoAcid"]
  edupall <- e[is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  eaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=edupall, FUN=sum)
  eaggdataorder <- eaggdata[order(eaggdata$Frequency), ]
  enotdup <- e[!is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  emerge <- rbind(enotdup, eaggdata)
  emergeorder <- emerge[order(emerge$Frequency), ]
  
  TCRd0zeroAA <- TCRd0zero[[1]]
  TCRd0zero <- TCRd0zero[1]
  asharedd0 <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, TCRd0zeroAA), ]
  bsharedd0 <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, TCRd0zeroAA), ]
  csharedd0 <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, TCRd0zeroAA), ]
  dsharedd0 <- dmergeorder[is.element(dmergeorder$CDR3.AminoAcid, TCRd0zeroAA), ]
  esharedd0 <- emergeorder[is.element(emergeorder$CDR3.AminoAcid, TCRd0zeroAA), ]
  colnames(asharedd0)[2] <- "Y1D7"
  colnames(bsharedd0)[2] <- "Y2D0"
  colnames(csharedd0)[2] <- "Y2D7"
  colnames(dsharedd0)[2] <- "Y3D0"
  colnames(esharedd0)[2] <- "Y3D7"
  TCRd0zeroshareda <- merge(TCRd0zero, asharedd0, all=TRUE)
  TCRd0zeroshareda[is.na(TCRd0zeroshareda)] <- 0
  TCRd0zerosharedab <- merge(TCRd0zeroshareda, bsharedd0, all=TRUE)
  TCRd0zerosharedab[is.na(TCRd0zerosharedab)] <- 0
  TCRd0zerosharedabc <- merge(TCRd0zerosharedab, csharedd0, all=TRUE)
  TCRd0zerosharedabc[is.na(TCRd0zerosharedabc)] <- 0
  TCRd0zerosharedabcd <- merge(TCRd0zerosharedabc, dsharedd0, all=TRUE)
  TCRd0zerosharedabcd[is.na(TCRd0zerosharedabcd)] <- 0
  TCRd0zerosharedabcde <- merge(TCRd0zerosharedabcd, esharedd0, all=TRUE)
  TCRd0zerosharedabcde[is.na(TCRd0zerosharedabcde)] <- 0
  TCRd0zerosharedabcde2 <- TCRd0zerosharedabcde[!(TCRd0zerosharedabcde$Y1D7==0 & TCRd0zerosharedabcde$Y2D0==0 & TCRd0zerosharedabcde$Y2D7==0 & TCRd0zerosharedabcde$Y3D0==0 & TCRd0zerosharedabcde$Y3D7==0), ]
  
  TCRmerge <- melt(TCRd0zerosharedabcde2, id=c("CDR3.AminoAcid"))
  colnames(TCRmerge)[2] <- "Day"
  colnames(TCRmerge)[3] <- "Frequency"
  #TCRmerge$Median <- NA
  
  #tempDataY1D7 <- subset(TCRmerge, Day=="Y1D7")
  #amed <- median(tempDataY1D7$Frequency)   #calculation of Y1D7 median
  #tempDataY2D0 <- subset(TCRmerge, Day=="Y2D0")
  #bmed <- median(tempDataY2D0$Frequency)   #calculation of Y2D0 median
  #tempDataY2D7 <- subset(TCRmerge, Day=="Y2D7")
  #cmed <- median(tempDataY2D7$Frequency)   #calculation of Y2D7 median
  #tempDataY3D0 <- subset(TCRmerge, Day=="Y3D0")
  #dmed <- median(tempDataY3D0$Frequency)   #calculation of Y3D0 median
  #tempDataY3D7 <- subset(TCRmerge, Day=="Y3D7")
  #emed <- median(tempDataY3D7$Frequency)   #calculation of Y3D7 median
  
  #dataMedian <- as.data.frame(c(amed, bmed, cmed, dmed, emed))
  #colnames(dataMedian)[1] <- "Median"
  #dataMedian$CDR3.AminoAcid <- NA
  #dataMedian$Day <- c("Y1D7", "Y2D0", "Y2D7", "Y3D0", "Y3D7")
  #dataMedian$Frequency <- NA
  #dataMedian <- dataMedian[c(2, 3, 4, 1)]
  
  #TCRmergeMed <- rbind(TCRmerge, dataMedian)
  
  p.sharedd7zerod0 <- ggplot(data=TCRmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid)) +
    geom_line(size=1, color="#C0D8C0", alpha=0.4) +           #color: #FFCC66 (light orange) for ICOS+CD38+, #C0D8C0 (light green) for ICOS-CD38-, #999999 (light gray) for CXCR5-
    geom_point(size=2, color="#C0D8C0", alpha=0.4) +           #color: #FFCC66 (light orange) for ICOS+CD38+, #C0D8C0 (light green) for ICOS-CD38-, #999999 (light gray) for CXCR5-
    scale_y_log10(limits=c(0.001, 10), breaks=c(0.001, 0.01, 0.1, 1, 10), label=comma) +
    #geom_line(data=TCRmergeMed, aes(x=Day, y=Median), color="darkorange2", size=2) +      #color: darkorange2 (dark orange) for ICOS+CD38+, green4 (dark green) for ICOS-CD38, #333333 (dark gray) for CXCR5-
    scale_color_manual(values=c()) +
    theme_bw() +
    theme(legend.title=element_text(size=17), legend.text=element_text(size=17), axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_blank(), axis.text=element_text(size=18, face="bold"), plot.title=element_text(size=20, face="bold")) +
    ggtitle("Day 0 CXCR5- Clones in ICOS-CD38- in 999") +
    labs(x="Day", y="Clonotypic frequency (%)")

  return(list(sum(asharedd0$Y1D7), sum(bsharedd0$Y2D0), sum(csharedd0$Y2D7), sum(dsharedd0$Y3D0), sum(esharedd0$Y3D7), p.sharedd7zerod0, TCRd0zerosharedabcde2))
}

aa.sharedd7zerod0(b5x, b3x, b1x, b4x, b2x, b1lo, b2lo, b3lo, b4lo, b5lo)


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
a <- b1hi
b <- b2hi
c <- b3hi
d <- b4hi
e <- b5hi