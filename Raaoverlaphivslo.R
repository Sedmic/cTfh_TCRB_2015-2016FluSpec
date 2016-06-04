TCRdataMultYears <- read.csv("../Samples_MultipleYearsData_AsOf02242016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(scales)


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

x <- a2hi; y <- a1hi
x <- a2lo; y <- a1lo


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

x <- b2hi; y <- b1hi
x <- b2lo; y <- b1lo


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

  ynotshared <- ymergeorder[!is.element(ymergeorder$CDR3.AminoAcid, xyduplicatedAA), ]
  ynotsharedorder <- ynotshared[order(ynotshared$CDR3.AminoAcid), ]
  colnames(ynotsharedorder)[2] <- "Frequency.y"
  ynotsharedorder$Frequency.x <- 0
  ynotsharedorder <- ynotsharedorder[c(1,3,2)]
  xnotshared <- xmergeorder[!is.element(xmergeorder$CDR3.AminoAcid, xyduplicatedAA), ]
  xnotsharedorder <- xnotshared[order(xnotshared$CDR3.AminoAcid), ]
  colnames(xnotsharedorder)[2] <- "Frequency.x"
  xnotsharedorder$Frequency.y <- 0
  xynotshared <- rbind(xnotsharedorder, ynotsharedorder)
  
  xyallhi <- rbind(xyshared, xynotshared)
  xyallhi$Cell <- "ICOS+CD38+ cTfh"
  sharedhi <- nrow(xyshared)
  totalhi <- nrow(xyallhi)
  
  xyalllo <- rbind(xyshared, xynotshared)
  xyalllo$Cell <- "ICOS-CD38- cTfh"
  sharedlo <- nrow(xyshared)
  totallo <- nrow(xyalllo)
  
  xyall <- rbind(xyallhi, xyalllo)
  
  ggplot() +
    geom_point(data=xyalllo, aes(Frequency.x, Frequency.y, fill="#C0D8C0"), size=5, pch=21, color="black", alpha=0.5) +
    geom_point(data=xyallhi, aes(Frequency.x, Frequency.y, fill="darkorange2"), size=5, pch=21, color="black", alpha=0.8) +
    scale_x_log10(limits=c(0.0001, 10), breaks=c(0.0001, 0.001, 0.010, 0.100, 1.00, 10.00), labels=comma) +
    scale_y_log10(limits=c(0.0001, 10), breaks=c(0.0001, 0.001, 0.010, 0.100, 1.00, 10.00), labels=comma) +
    theme_bw() +
    scale_fill_manual(values=c("#C0D8C0", "darkorange2")) +
    theme(legend.position="none") + 
    ggtitle("Year 2 Day 0") +
#   labs(x="Clonotypic frequency in Year3 Day7 (%)", y="Clonotypic frequency in Year1 Day7 (%)")
#    labs(x="Clonotypic frequency in Year3 Day0 (%)", y="Clonotypic frequency in Year1 Day7 (%)")
#   labs(x="Clonotypic frequency in Year2 Day7 (%)", y="Clonotypic frequency in Year1 Day7 (%)")
   labs(x="Clonotypic frequency in Year2 Day0 (%)", y="Clonotypic frequency in Year1 Day7 (%)") 

#ggsave(filename = "../R Images/R101Y3D7vY1D7.png",height=4,width=4.5,units="in")
#ggsave(filename = "../R Images/R101Y3D0vY1D7.png",height=4,width=4.5,units="in")
#ggsave(filename = "../R Images/R101Y2D7vY1D7.png",height=4,width=4.5,units="in")
#ggsave(filename = "../R Images/R101Y2D0vY1D7.png",height=4,width=4.5,units="in")

#ggsave(filename = "../R Images/R999Y3D7vY1D7.png",height=4,width=4.5,units="in")
#ggsave(filename = "../R Images/R999Y3D0vY1D7.png",height=4,width=4.5,units="in")
#ggsave(filename = "../R Images/R999Y2D7vY1D7.png",height=4,width=4.5,units="in")
#ggsave(filename = "../R Images/R999Y2D0vY1D7.png",height=4,width=4.5,units="in")
