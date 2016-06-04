TCRdata1516 <- read.csv("Samples_1516Data_AsOf02242016.csv")
TCRdata1516In <- subset(TCRdata1516, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
library(ggplot2)
library(reshape2)

x <- subset(TCRdata1516In, Subject==101 & Cell=="ICOS+CD38+" & Day==7, select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency))
v <- subset(TCRdata1516In, Subject==101 & Cell=="ICOS+CD38+" & Day==0, select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency))

xduplicatedAAtf <- duplicated(x$CDR3.AminoAcid)
xrowduplicated <- which(xduplicatedAAtf, xduplicatedAAtf=="TRUE")
xduplicatedAA <- x[xrowduplicated, "CDR3.AminoAcid"]
xdupall <- x[is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
xaggdata <- aggregate(Frequency~CDR3.AminoAcid+Day, data=xdupall, FUN=sum)
xaggdataorder <- xaggdata[order(xaggdata$Frequency), ]
xnotdup <- x[!is.element(x$CDR3.AminoAcid, xduplicatedAA), ]
xnotdup <- subset(xnotdup, select=c(CDR3.AminoAcid, Day, Frequency))
xnotduporder <- xnotdup[order(xnotdup$Frequency), ]
xmerge <- rbind(xnotduporder, xaggdataorder)
xmergeorder <- xmerge[order(xmerge$Frequency), ]
xmergetop <- xmergeorder[767:776, ]
xtopclones7 <- c("CASRELNTEAFF", "CASRALMGQTNSNQPQHF", "CASSQLSTGAYEQYF", "CASSHDSPYNSPLHF", "CASQRTSGGRRSYEQYF", "CSAPTSGSSYEQYF", "CASSWTVAEQYF", "CASSLSRASEQYF", "CASRAQKGNQPQHF", "CATLGPYSNQPQHF")

vduplicatedAAtf <- duplicated(v$CDR3.AminoAcid)
vrowduplicated <- which(vduplicatedAAtf, vduplicatedAAtf=="TRUE")
vduplicatedAA <- v[vrowduplicated, "CDR3.AminoAcid"]
vdupall <- v[is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
vaggdata <- aggregate(Frequency~CDR3.AminoAcid+Day, data=vdupall, FUN=sum)
vaggdataorder <- vaggdata[order(vaggdata$Frequency), ]
vnotdup <- v[!is.element(v$CDR3.AminoAcid, vduplicatedAA), ]
vnotdup <- subset(vnotdup, select=c(CDR3.AminoAcid, Day, Frequency))
vnotduporder <- vnotdup[order(vnotdup$Frequency), ]
vmerge <- rbind(vnotduporder, vaggdataorder)
vmergeorder <- vmerge[order(vmerge$Frequency), ]

xtopclones0 <- vmergeorder[is.element(vmergeorder$CDR3.AminoAcid, xtopclones7), ]
colnames(xmergetop)[3] <- "7"
xmergetop$Day0 <- 0
colnames(xmergetop)[4] <- "0"
xmergetop <- xmergetop[c(1, 4, 3)]
xmergetop2 <- melt(xmergetop, id=c("CDR3.AminoAcid"))
colnames(xmergetop2)[2] <- "Day"
colnames(xmergetop2)[3] <- "Frequency"

ggplot(data=xmergetop2, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, label=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("Single clone frequencies of ICOS+CD38+ of Subject 101") +
  labs(x="Day of Year 1516", y="Clonotypic frequency (%)")


y <- subset(TCRdata1516In, Subject==101 & Cell=="ICOS-CD38-" & Day==7, select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency))
w <- subset(TCRdata1516In, Subject==101 & Cell=="ICOS-CD38-" & Day==0, select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency))

yduplicatedAAtf <- duplicated(y$CDR3.AminoAcid)
yrowduplicated <- which(yduplicatedAAtf, yduplicatedAAtf=="TRUE")
yduplicatedAA <- y[yrowduplicated, "CDR3.AminoAcid"]
ydupall <- y[is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
yaggdata <- aggregate(Frequency~CDR3.AminoAcid+Day, data=ydupall, FUN=sum)
yaggdataorder <- yaggdata[order(yaggdata$Frequency), ]
ynotdup <- y[!is.element(y$CDR3.AminoAcid, yduplicatedAA), ]
ynotdup <- subset(ynotdup, select=c(CDR3.AminoAcid, Day, Frequency))
ynotduporder <- ynotdup[order(ynotdup$Frequency), ]
ymerge <- rbind(ynotduporder, yaggdataorder)
ymergeorder <- ymerge[order(ymerge$Frequency), ]
ymergetop <- ymergeorder[1520:1529, ]
ytopclones7 <- c("CASSQEFGYTF", "CASSFVQLNEQFF", "CASSSTSEGETQYF", "CASSLERGTAYGYTF", "CASSRTTGPIGYTF")

wduplicatedAAtf <- duplicated(w$CDR3.AminoAcid)
wrowduplicated <- which(wduplicatedAAtf, wduplicatedAAtf=="TRUE")
wduplicatedAA <- w[wrowduplicated, "CDR3.AminoAcid"]
wdupall <- w[is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
waggdata <- aggregate(Frequency~CDR3.AminoAcid+Day, data=wdupall, FUN=sum)
waggdataorder <- waggdata[order(waggdata$Frequency), ]
wnotdup <- w[!is.element(w$CDR3.AminoAcid, wduplicatedAA), ]
wnotdup <- subset(wnotdup, select=c(CDR3.AminoAcid, Day, Frequency))
wnotduporder <- wnotdup[order(wnotdup$Frequency), ]
wmerge <- rbind(wnotduporder, waggdataorder)
wmergeorder <- wmerge[order(wmerge$Frequency), ]

ytopclones0 <- wmergeorder[is.element(wmergeorder$CDR3.AminoAcid, ytopclones7), ]
colnames(ymergetop)[3] <- "Frequency"
ywmerge <- rbind(ymergetop, ytopclones0)

ggplot(data=ywmerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, label=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  ylim(c(0, 2.85)) +
  theme_bw() +
  ggtitle("Single clone frequencies of ICOS-CD38- of Subject 101") +
  labs(x="Day of Year 1516", y="Clonotypic frequency (%)")


z <- subset(TCRdata1516In, Subject==101 & Cell=="CXCR5-" & Day==7, select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency))
u <- subset(TCRdata1516In, Subject==101 & Cell=="CXCR5-" & Day==0, select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency))

zduplicatedAAtf <- duplicated(z$CDR3.AminoAcid)
zrowduplicated <- which(zduplicatedAAtf, zduplicatedAAtf=="TRUE")
zduplicatedAA <- z[zrowduplicated, "CDR3.AminoAcid"]
zdupall <- z[is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
zaggdata <- aggregate(Frequency~CDR3.AminoAcid+Day, data=zdupall, FUN=sum)
zaggdataorder <- zaggdata[order(zaggdata$Frequency), ]
znotdup <- z[!is.element(z$CDR3.AminoAcid, zduplicatedAA), ]
znotdup <- subset(znotdup, select=c(CDR3.AminoAcid, Day, Frequency))
znotduporder <- znotdup[order(znotdup$Frequency), ]
zmerge <- rbind(znotduporder, zaggdataorder)
zmergeorder <- zmerge[order(zmerge$Frequency), ]
zmergetop <- zmergeorder[8752:8761, ]
ztopclones7 <- c("CASSQGASNQPQHF", "CASSLGQSAFF", "CASSQGGQGNGYTF", "CASSQSGGFSYEQYF", "CASSEWGEQYF")

uduplicatedAAtf <- duplicated(u$CDR3.AminoAcid)
urowduplicated <- which(uduplicatedAAtf, uduplicatedAAtf=="TRUE")
uduplicatedAA <- u[urowduplicated, "CDR3.AminoAcid"]
udupall <- u[is.element(u$CDR3.AminoAcid, uduplicatedAA), ]
uaggdata <- aggregate(Frequency~CDR3.AminoAcid+Day, data=udupall, FUN=sum)
uaggdataorder <- uaggdata[order(uaggdata$Frequency), ]
unotdup <- u[!is.element(u$CDR3.AminoAcid, uduplicatedAA), ]
unotdup <- subset(unotdup, select=c(CDR3.AminoAcid, Day, Frequency))
unotduporder <- unotdup[order(unotdup$Frequency), ]
umerge <- rbind(unotduporder, uaggdataorder)
umergeorder <- umerge[order(umerge$Frequency), ]

ztopclones0 <- umergeorder[is.element(umergeorder$CDR3.AminoAcid, ztopclones7), ]
a <- c("CASSLGQSAFF", 0, 0)
a <- as.data.frame(a)
a <- t(a)
a <- as.data.frame(a)
colnames(a)[1] <- "CDR3.AminoAcid"
colnames(a)[2] <- "Day"
colnames(a)[3] <- "Frequency"
row.names(a) [1] <- 1
a$Frequency <- as.numeric(a$Frequency)
a$Frequency <- 0
ztopclones0 <- rbind(ztopclones0, a)
colnames(zmergetop)[3] <- "Frequency"
zumerge <- rbind(zmergetop, ztopclones0)

ggplot(data=zumerge, aes(x=Day, y=Frequency, group=CDR3.AminoAcid, label=CDR3.AminoAcid, color=CDR3.AminoAcid)) +
  geom_line(size=1) +
  geom_point(size=3) +
  ylim(c(0, 2.85)) +
  theme_bw() +
  ggtitle("Single clone frequencies of CXCR5- of Subject 101") +
  labs(x="Day of Year 1516", y="Clonotypic frequency (%)")
