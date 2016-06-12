library(ggplot2); library(gplots); library(RColorBrewer); 

dataMatrix <- read.csv(file="../Samples_1516Data_AsOf02242016.csv",stringsAsFactors = FALSE)

just108 <- subset(dataMatrix, Subject == 117) #; subsetData <- subset(subsetData, Day == "7")
just108hihi <- subset(just108, Cell == "ICOS+CD38+")   ## now subsetData has just 101, last year, ICOS+CD38+


## display counts as a histogram
hist(just108hihi$Count[which(just108hihi$Day=="0")])
hist(just108hihi$Count[which(just108hihi$Day=="7")])

## invert the histogram 
hist(1/just108hihi$Count[which(just108hihi$Day == "0")])
hist(1/just108hihi$Count[which(just108hihi$Day == "7")])

set.seed(50)
subsetData <- subset(just108hihi,Day=="7")
plotData <- data.frame(x=sample(1:nrow(subsetData)*3,nrow(subsetData),replace=F),y=sample(1:nrow(subsetData)*3,nrow(subsetData),replace=F))
subsetData <- cbind(subsetData,plotData)
S101d7 <- ggplot(data=subsetData,aes(x=subsetData$x,y=subsetData$y,size=subsetData$Frequency,fill=subsetData$Frequency)) + theme_bw() + 
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank()) +
  xlab("arbitrary X") + ylab("arbitrary Y")+ggtitle("ICOS+CD38+ at d7 for 1516-117") +
  geom_point(colour="black", pch=21) + theme(legend.key = element_blank()) + 
  scale_fill_distiller(palette="Oranges",trans="reverse",limits=c(2,0),name="Clone Freq") +
  scale_size_continuous(limits = c(0,2),range=c(0,10),guide=guide_legend(title="Clone Freq")) 
S101d7
ggsave(S101d7,filename="S117d7.jpg",width=6,height=5)

subsetData <- subset(just108hihi,Day=="0")
plotData <- data.frame(x=sample(1:nrow(subsetData)*100,nrow(subsetData),replace=F),y=sample(1:nrow(subsetData)*100,nrow(subsetData),replace=F))
subsetData <- cbind(subsetData,plotData)
S101d0 <- ggplot(data=subsetData,aes(x=subsetData$x,y=subsetData$y,size=subsetData$Frequency,fill=subsetData$Frequency)) + theme_bw() +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank()) +
  xlab("arbitrary X") + ylab("arbitrary Y")+ggtitle("ICOS+CD38+ at d0 for 1516-117") +
  geom_point(colour="black", pch=21) + theme(legend.key = element_blank()) + 
  scale_fill_distiller(palette="Oranges",trans="reverse",limits=c(2,0),name="Clone Freq") +
  scale_size_continuous(limits = c(0,2),range=c(0,10),guide=guide_legend(title="Clone Freq")) 
S101d0
ggsave(S101d0,filename="S117d0.jpg", width=6, height=5)



dataMatrix2 <- read.csv(file="../Samples_MultipleYearsData_AsOf02242016.csv",stringsAsFactors = FALSE)
just101 <- subset(dataMatrix2, Subject == 101) #; subsetData <- subset(subsetData, Day == "7")
just101hihi <- subset(just101, Cell == "ICOS+CD38+")   ## now subsetData has just 101, last year, ICOS+CD38+


###   ******************** reshape the matrix to put freq for matching CDR3 together  **********

library(reshape2)
dataMatrix <- read.csv(file="../Samples_MultipleYearsData_AsOf02242016.csv",stringsAsFactors = FALSE)

allsubsets999 <- subset(dataMatrix,Subject == 999)
subsetData <- allsubsets999[which(allsubsets999$SequenceStatus=="In"),]  # only want in-frame sequences
subsetData <- subsetData[,-c(9,10,11,12,13)] # take away columns on V, D, and J usage
subsetData <- subsetData[-which(subsetData$Day==c("6wk")),]  # get rid of tetramer sequences

mapUnique <- unique(subsetData$CDR3.Nucleotide)

subsetData$timeIndex <- 0  # initialize column by setting to zero first
subsetData$timeIndex[which(subsetData$Year==1314 & subsetData$Day=="7")] <- 8
subsetData$timeIndex[which(subsetData$Year==1415 & subsetData$Day=="0")] <- 15
subsetData$timeIndex[which(subsetData$Year==1415 & subsetData$Day=="7")] <- 17
subsetData$timeIndex[which(subsetData$Year==1516 & subsetData$Day=="0")] <- 24
subsetData$timeIndex[which(subsetData$Year==1516 & subsetData$Day=="7")] <- 26
mydata <- dcast(subsetData, CDR3.Length+CDR3.Nucleotide+CDR3.AminoAcid+timeIndex~Cell,value.var="Frequency")
#mydata <- dcast(subsetData, CDR3.Length+CDR3.Nucleotide+CDR3.AminoAcid~Cell+timeIndex,value.var="Frequency")
mydata[is.na(mydata)] <- 0; 
#colnames(mydata) <- c("CDR3Length","CDR3Nucl","CDR3AA","timeIndex","X5lo1314_7","X5lo1415_0",
#                      "X5lo1415_7","X5lo1516_0", "X5lo1516_7","Lolo1314_7","Lolo1415_0","Lolo1415_7",
#                      "Lolo1516_0","Lolo1516_7","Hihi1314_7","Hihi1415_0","Hihi1415_7","Hihi1516_0",
#                      "Hihi1516_7")
colnames(mydata) <- c("CDR3Length","CDR3Nucl","CDR3AA","timeIndex","X5lo","Lolo","Hihi")
#mydata <- mydata[order(mydata$timeIndex,-mydata$Hihi1516_7, na.last=TRUE),]
mydata <- mydata[order(mydata$timeIndex,-mydata$Hihi, na.last=TRUE),]

# now i need to add CDR3 sequences with 0 freq that did not appear in the particular timepoint
mydata8 <- mydata[which(mydata$timeIndex==8),]
mydata15 <- mydata[which(mydata$timeIndex==15),]
mydata17 <- mydata[which(mydata$timeIndex==17),]
mydata24 <- mydata[which(mydata$timeIndex==24),]
mydata26 <- mydata[which(mydata$timeIndex==26),]

# each block will now match data from columns against complete CDR3 series, for each timeIndex, NA where not matched
mydata8compl <- data.frame(CDR3Length=0,CDR3Nucl=mapUnique,CDR3AA="z",timeIndex=0,X5lo=0,Lolo=0,Hihi=0,stringsAsFactors = FALSE)
mydata8compl$timeIndex <- mydata8[match(mydata8compl$CDR3Nucl,mydata8$CDR3Nucl),4]
mydata8compl$X5lo <- mydata8[match(mydata8compl$CDR3Nucl,mydata8$CDR3Nucl),5]
mydata8compl$Lolo <- mydata8[match(mydata8compl$CDR3Nucl,mydata8$CDR3Nucl),6]
mydata8compl$Hihi <- mydata8[match(mydata8compl$CDR3Nucl,mydata8$CDR3Nucl),7]
mydata8compl$timeIndex <- 8


mydata15compl <- data.frame(CDR3Length=0,CDR3Nucl=mapUnique,CDR3AA="z",timeIndex=0,X5lo=0,Lolo=0,Hihi=0,stringsAsFactors = FALSE)
mydata15compl$timeIndex <- mydata15[match(mydata15compl$CDR3Nucl,mydata15$CDR3Nucl),4]
mydata15compl$X5lo <- mydata15[match(mydata15compl$CDR3Nucl,mydata15$CDR3Nucl),5]
mydata15compl$Lolo <- mydata15[match(mydata15compl$CDR3Nucl,mydata15$CDR3Nucl),6]
mydata15compl$Hihi <- mydata15[match(mydata15compl$CDR3Nucl,mydata15$CDR3Nucl),7]
mydata15compl$timeIndex <- 15


mydata17compl <- data.frame(CDR3Length=0,CDR3Nucl=mapUnique,CDR3AA="z",timeIndex=0,X5lo=0,Lolo=0,Hihi=0,stringsAsFactors = FALSE)
mydata17compl$timeIndex <- mydata17[match(mydata17compl$CDR3Nucl,mydata17$CDR3Nucl),4]
mydata17compl$X5lo <- mydata17[match(mydata17compl$CDR3Nucl,mydata17$CDR3Nucl),5]
mydata17compl$Lolo <- mydata17[match(mydata17compl$CDR3Nucl,mydata17$CDR3Nucl),6]
mydata17compl$Hihi <- mydata17[match(mydata17compl$CDR3Nucl,mydata17$CDR3Nucl),7]
mydata17compl$timeIndex <- 17


mydata24compl <- data.frame(CDR3Length=0,CDR3Nucl=mapUnique,CDR3AA="z",timeIndex=0,X5lo=0,Lolo=0,Hihi=0,stringsAsFactors = FALSE)
mydata24compl$timeIndex <- mydata24[match(mydata24compl$CDR3Nucl,mydata24$CDR3Nucl),4]
mydata24compl$X5lo <- mydata24[match(mydata24compl$CDR3Nucl,mydata24$CDR3Nucl),5]
mydata24compl$Lolo <- mydata24[match(mydata24compl$CDR3Nucl,mydata24$CDR3Nucl),6]
mydata24compl$Hihi <- mydata24[match(mydata24compl$CDR3Nucl,mydata24$CDR3Nucl),7]
mydata24compl$timeIndex <- 24


mydata26compl <- data.frame(CDR3Length=0,CDR3Nucl=mapUnique,CDR3AA="z",timeIndex=0,X5lo=0,Lolo=0,Hihi=0,stringsAsFactors = FALSE)
mydata26compl$timeIndex <- mydata26[match(mydata26compl$CDR3Nucl,mydata26$CDR3Nucl),4]
mydata26compl$X5lo <- mydata26[match(mydata26compl$CDR3Nucl,mydata26$CDR3Nucl),5]
mydata26compl$Lolo <- mydata26[match(mydata26compl$CDR3Nucl,mydata26$CDR3Nucl),6]
mydata26compl$Hihi <- mydata26[match(mydata26compl$CDR3Nucl,mydata26$CDR3Nucl),7]
mydata26compl$timeIndex <- 26

fullver <- rbind(mydata8compl, mydata15compl, mydata17compl, mydata24compl, mydata26compl)
unique_mydata <- mydata[!duplicated(mydata$CDR3Nucl),]

fullver$mapUnique <- match(fullver$CDR3Nucl,unique_mydata$CDR3Nucl)
fullver$CDR3Length <- mydata[match(fullver$CDR3Nucl,mydata$CDR3Nucl),1]
fullver$CDR3AA <- mydata[match(fullver$CDR3Nucl,mydata$CDR3Nucl),3]
fullver[is.na(fullver)] <- 0; 

## now fullver is the complete matrix containing complete data set for 999

write.csv(fullver, file="fullver999.csv")

smallver <- fullver[sample(519000,10000,replace=FALSE),]



###   alternative where matrix is in blocks by subset to facilitate coloring and showing of all data

blockMat <- melt(fullver,id.vars=c("CDR3Length","CDR3Nucl","CDR3AA","timeIndex","mapUnique"),value.name="Freq")
smallver <- blockMat[sample(1557060,50000,replace=FALSE),]

library(rgl)

x <- smallver$mapUnique
y <- smallver$timeIndex
z <- smallver$Freq

rgl.open(); rgl.bg(color="white")
plot3d(x=x,y=y,z=z, col=as.numeric(smallver$variable), type="p")
rgl.bbox(color="white")



x <- blockMat$mapUnique
z <- blockMat$timeIndex
y <- blockMat$Freq

rgl.open(); rgl.bg(color="white")
plot3d(x=x,y=y,z=z, col=as.numeric(blockMat$variable)+1, type="p")
rgl.bbox(color="white")
par3d(windowRect = c(3290, 140, 4051, 694 ),zoom=0.65 )

#play3d(spin3d(axis=c(0,1,0),rpm=5),duration=5)

movie3d(spin3d(axis=c(0,1,0),rpm=4), duration=15, dev = rgl.cur(), fps = 20, dir="./movie",
 #       movie = "movie", frames = movie, dir = tempdir(), 
        convert = FALSE, clean = TRUE, verbose = TRUE,
        top = TRUE, type = "gif", startTime = 0) 


rgl.close()


require(akima) ; require(rgl)
x=seq(1,10,1)
y=sample(10,10,replace=TRUE)
z=rnorm(10) + 5
s=interp(x,y,z)
dim(s$z)

surface3d(s$x,s$y,s$z)



