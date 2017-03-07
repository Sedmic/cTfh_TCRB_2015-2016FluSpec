TCRdataMultYears <- read.csv("Samples_MultipleYearsData_AsOf09162016.csv")
TCRdataMultYearsIn <- subset(TCRdataMultYears, SequenceStatus=="In", select=c(Subject, Year, Cell, Day, CDR3.AminoAcid, Frequency, Count))
Fluspecific <- read.csv("RepeatingFluSpecificClones.csv")
AIMs <- read.csv("AIMsUnique.csv")
library("gplots")
library("RColorBrewer")
RdYlBu <- colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(400)
library("ggExtra")
library("scales")


TCRdata999 <- subset(TCRdataMultYearsIn, Subject==999, select=c(Year, Cell, Day, CDR3.AminoAcid, Frequency))
ROR999 <- subset(Fluspecific, Subject==999, select=c(CDR3.AminoAcid)) #recurring oligoclonal response clonotypes
aim999 <- subset(AIMs, Subject==999, select=c(CDR3.AminoAcid)) # AIM clonotypes

#Heatmap: Subject 999 - Recurrent oligoclonal response clones over time and subset
#Enter ROR for 999 as "x" and all in-sequence data for 999 as "y"
aa.RORheatmap999 <- function(x,y) {
  a <- subset(y, Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  b <- subset(y, Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  c <- subset(y, Year==1314 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  d <- subset(y, Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  e <- subset(y, Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  f <- subset(y, Year==1415 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  g <- subset(y, Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  h <- subset(y, Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  i <- subset(y, Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  j <- subset(y, Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  k <- subset(y, Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  l <- subset(y, Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  m <- subset(y, Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  n <- subset(y, Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  o <- subset(y, Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  
  
  aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)   #condense rows in a with same amino acid sequence but different nucleotide sequences
  arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
  aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
  adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)
  aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
  anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  amerge <- rbind(anotdup, aaggdata)
  amergeorder <- amerge[order(amerge$Frequency), ]    #list of unique aa sequences in a
  
  bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)   #condense rows in b with same amino acid sequence but different nucleotide sequences
  browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
  bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
  bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)
  baggdataorder <- baggdata[order(baggdata$Frequency), ]
  bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  bmerge <- rbind(bnotdup, baggdata)
  bmergeorder <- bmerge[order(bmerge$Frequency), ]    #list of unique aa sequences in b
  
  cduplicatedAAtf <- duplicated(c$CDR3.AminoAcid)   #condense rows in c with same amino acid sequence but different nucleotide sequences
  crowduplicated <- which(cduplicatedAAtf, cduplicatedAAtf=="TRUE")
  cduplicatedAA <- c[crowduplicated, "CDR3.AminoAcid"]
  cdupall <- c[is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  caggdata <- aggregate(Frequency~CDR3.AminoAcid, data=cdupall, FUN=sum)
  caggdataorder <- caggdata[order(caggdata$Frequency), ]
  cnotdup <- c[!is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  cmerge <- rbind(cnotdup, caggdata)
  cmergeorder <- cmerge[order(cmerge$Frequency), ]    #list of unique aa sequences in c
  
  dduplicatedAAtf <- duplicated(d$CDR3.AminoAcid)   #condense rows in d with same amino acid sequence but different nucleotide sequences
  drowduplicated <- which(dduplicatedAAtf, dduplicatedAAtf=="TRUE")
  dduplicatedAA <- d[drowduplicated, "CDR3.AminoAcid"]
  ddupall <- d[is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  daggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ddupall, FUN=sum)
  daggdataorder <- daggdata[order(daggdata$Frequency), ]
  dnotdup <- d[!is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  dmerge <- rbind(dnotdup, daggdata)
  dmergeorder <- dmerge[order(dmerge$Frequency), ]    #list of unique aa sequences in d
  
  eduplicatedAAtf <- duplicated(e$CDR3.AminoAcid)   #condense rows in e with same amino acid sequence but different nucleotide sequences
  erowduplicated <- which(eduplicatedAAtf, eduplicatedAAtf=="TRUE")
  eduplicatedAA <- e[erowduplicated, "CDR3.AminoAcid"]
  edupall <- e[is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  eaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=edupall, FUN=sum)
  eaggdataorder <- eaggdata[order(eaggdata$Frequency), ]
  enotdup <- e[!is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  emerge <- rbind(enotdup, eaggdata)
  emergeorder <- emerge[order(emerge$Frequency), ]    #list of unique aa sequences in e
  
  fduplicatedAAtf <- duplicated(f$CDR3.AminoAcid)   #condense rows in f with same amino acid sequence but different nucleotide sequences
  frowduplicated <- which(fduplicatedAAtf, fduplicatedAAtf=="TRUE")
  fduplicatedAA <- f[frowduplicated, "CDR3.AminoAcid"]
  fdupall <- f[is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  faggdata <- aggregate(Frequency~CDR3.AminoAcid, data=fdupall, FUN=sum)
  faggdataorder <- faggdata[order(faggdata$Frequency), ]
  fnotdup <- f[!is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  fmerge <- rbind(fnotdup, faggdata)
  fmergeorder <- fmerge[order(fmerge$Frequency), ]    #list of unique aa sequences in f
  
  gduplicatedAAtf <- duplicated(g$CDR3.AminoAcid)   #condense rows in g with same amino acid sequence but different nucleotide sequences
  growduplicated <- which(gduplicatedAAtf, gduplicatedAAtf=="TRUE")
  gduplicatedAA <- g[growduplicated, "CDR3.AminoAcid"]
  gdupall <- g[is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=gdupall, FUN=sum)
  gaggdataorder <- gaggdata[order(gaggdata$Frequency), ]
  gnotdup <- g[!is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gmerge <- rbind(gnotdup, gaggdata)
  gmergeorder <- gmerge[order(gmerge$Frequency), ]    #list of unique aa sequences in g
  
  hduplicatedAAtf <- duplicated(h$CDR3.AminoAcid)   #condense rows in h with same amino acid sequence but different nucleotide sequences
  hrowduplicated <- which(hduplicatedAAtf, hduplicatedAAtf=="TRUE")
  hduplicatedAA <- h[hrowduplicated, "CDR3.AminoAcid"]
  hdupall <- h[is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  haggdata <- aggregate(Frequency~CDR3.AminoAcid, data=hdupall, FUN=sum)
  haggdataorder <- haggdata[order(haggdata$Frequency), ]
  hnotdup <- h[!is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  hmerge <- rbind(hnotdup, haggdata)
  hmergeorder <- hmerge[order(hmerge$Frequency), ]    #list of unique aa sequences in h
  
  iduplicatedAAtf <- duplicated(i$CDR3.AminoAcid)   #condense rows in i with same amino acid sequence but different nucleotide sequences
  irowduplicated <- which(iduplicatedAAtf, iduplicatedAAtf=="TRUE")
  iduplicatedAA <- i[irowduplicated, "CDR3.AminoAcid"]
  idupall <- i[is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  iaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=idupall, FUN=sum)
  iaggdataorder <- iaggdata[order(iaggdata$Frequency), ]
  inotdup <- i[!is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  imerge <- rbind(inotdup, iaggdata)
  imergeorder <- imerge[order(imerge$Frequency), ]    #list of unique aa sequences in i
  
  jduplicatedAAtf <- duplicated(j$CDR3.AminoAcid)   #condense rows in j with same amino acid sequence but different nucleotide sequences
  jrowduplicated <- which(jduplicatedAAtf, jduplicatedAAtf=="TRUE")
  jduplicatedAA <- j[jrowduplicated, "CDR3.AminoAcid"]
  jdupall <- j[is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=jdupall, FUN=sum)
  jaggdataorder <- jaggdata[order(jaggdata$Frequency), ]
  jnotdup <- j[!is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jmerge <- rbind(jnotdup, jaggdata)
  jmergeorder <- jmerge[order(jmerge$Frequency), ]    #list of unique aa sequences in j
  
  kduplicatedAAtf <- duplicated(k$CDR3.AminoAcid)   #condense rows in k with same amino acid sequence but different nucleotide sequences
  krowduplicated <- which(kduplicatedAAtf, kduplicatedAAtf=="TRUE")
  kduplicatedAA <- k[krowduplicated, "CDR3.AminoAcid"]
  kdupall <- k[is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=kdupall, FUN=sum)
  kaggdataorder <- kaggdata[order(kaggdata$Frequency), ]
  knotdup <- k[!is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kmerge <- rbind(knotdup, kaggdata)
  kmergeorder <- kmerge[order(kmerge$Frequency), ]    #list of unique aa sequences in k
  
  lduplicatedAAtf <- duplicated(l$CDR3.AminoAcid)   #condense rows in l with same amino acid sequence but different nucleotide sequences
  lrowduplicated <- which(lduplicatedAAtf, lduplicatedAAtf=="TRUE")
  lduplicatedAA <- l[lrowduplicated, "CDR3.AminoAcid"]
  ldupall <- l[is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  laggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ldupall, FUN=sum)
  laggdataorder <- laggdata[order(laggdata$Frequency), ]
  lnotdup <- l[!is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  lmerge <- rbind(lnotdup, laggdata)
  lmergeorder <- lmerge[order(lmerge$Frequency), ]    #list of unique aa sequences in l
  
  mduplicatedAAtf <- duplicated(m$CDR3.AminoAcid)   #condense rows in m with same amino acid sequence but different nucleotide sequences
  mrowduplicated <- which(mduplicatedAAtf, mduplicatedAAtf=="TRUE")
  mduplicatedAA <- m[mrowduplicated, "CDR3.AminoAcid"]
  mdupall <- m[is.element(m$CDR3.AminoAcid, mduplicatedAA), ]
  maggdata <- aggregate(Frequency~CDR3.AminoAcid, data=mdupall, FUN=sum)
  maggdataorder <- maggdata[order(maggdata$Frequency), ]
  mnotdup <- m[!is.element(m$CDR3.AminoAcid, mduplicatedAA), ]
  mmerge <- rbind(mnotdup, maggdata)
  mmergeorder <- mmerge[order(mmerge$Frequency), ]    #list of unique aa sequences in m
  
  nduplicatedAAtf <- duplicated(n$CDR3.AminoAcid)   #condense rows in n with same amino acid sequence but different nucleotide sequences
  nrowduplicated <- which(nduplicatedAAtf, nduplicatedAAtf=="TRUE")
  nduplicatedAA <- n[nrowduplicated, "CDR3.AminoAcid"]
  ndupall <- n[is.element(n$CDR3.AminoAcid, nduplicatedAA), ]
  naggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ndupall, FUN=sum)
  naggdataorder <- naggdata[order(naggdata$Frequency), ]
  nnotdup <- n[!is.element(n$CDR3.AminoAcid, nduplicatedAA), ]
  nmerge <- rbind(nnotdup, naggdata)
  nmergeorder <- nmerge[order(nmerge$Frequency), ]    #list of unique aa sequences in n
  
  oduplicatedAAtf <- duplicated(o$CDR3.AminoAcid)   #condense rows in o with same amino acid sequence but different nucleotide sequences
  orowduplicated <- which(oduplicatedAAtf, oduplicatedAAtf=="TRUE")
  oduplicatedAA <- o[orowduplicated, "CDR3.AminoAcid"]
  odupall <- o[is.element(o$CDR3.AminoAcid, oduplicatedAA), ]
  oaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=odupall, FUN=sum)
  oaggdataorder <- oaggdata[order(oaggdata$Frequency), ]
  onotdup <- o[!is.element(o$CDR3.AminoAcid, oduplicatedAA), ]
  omerge <- rbind(onotdup, oaggdata)
  omergeorder <- omerge[order(omerge$Frequency), ]    #list of unique aa sequences in o
  
  xAA <- x[[1]]
  ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, xAA), ] #find overlapping sequences between ROR and each subset
  bshared <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, xAA), ]
  cshared <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, xAA), ]
  dshared <- dmergeorder[is.element(dmergeorder$CDR3.AminoAcid, xAA), ]
  eshared <- emergeorder[is.element(emergeorder$CDR3.AminoAcid, xAA), ]
  fshared <- fmergeorder[is.element(fmergeorder$CDR3.AminoAcid, xAA), ]
  gshared <- gmergeorder[is.element(gmergeorder$CDR3.AminoAcid, xAA), ]
  hshared <- hmergeorder[is.element(hmergeorder$CDR3.AminoAcid, xAA), ]
  ishared <- imergeorder[is.element(imergeorder$CDR3.AminoAcid, xAA), ]
  jshared <- jmergeorder[is.element(jmergeorder$CDR3.AminoAcid, xAA), ]
  kshared <- kmergeorder[is.element(kmergeorder$CDR3.AminoAcid, xAA), ]
  lshared <- lmergeorder[is.element(lmergeorder$CDR3.AminoAcid, xAA), ]
  mshared <- mmergeorder[is.element(mmergeorder$CDR3.AminoAcid, xAA), ]
  nshared <- nmergeorder[is.element(nmergeorder$CDR3.AminoAcid, xAA), ]
  oshared <- omergeorder[is.element(omergeorder$CDR3.AminoAcid, xAA), ]
  colnames(ashared)[2] <- "Hi.Y1D7"
  colnames(bshared)[2] <- "Lo.Y1D7"
  colnames(cshared)[2] <- "X5.Y1D7"
  colnames(dshared)[2] <- "Hi.Y2D0"
  colnames(eshared)[2] <- "Lo.Y2D0"
  colnames(fshared)[2] <- "X5.Y2D0"
  colnames(gshared)[2] <- "Hi.Y2D7"
  colnames(hshared)[2] <- "Lo.Y2D7"
  colnames(ishared)[2] <- "X5.Y2D7"
  colnames(jshared)[2] <- "Hi.Y3D0"
  colnames(kshared)[2] <- "Lo.Y3D0"
  colnames(lshared)[2] <- "X5.Y3D0"
  colnames(mshared)[2] <- "Hi.Y3D7"
  colnames(nshared)[2] <- "Lo.Y3D7"
  colnames(oshared)[2] <- "X5.Y3D7"
  
  xshareda <- merge(x, ashared, all=TRUE) #merge all data into one data frame
  xshareda[is.na(xshareda)] <- 0
  xsharedab <- merge(xshareda, bshared, all=TRUE)
  xsharedab[is.na(xsharedab)] <- 0
  xsharedabc <- merge(xsharedab, cshared, all=TRUE)
  xsharedabc[is.na(xsharedabc)] <- 0
  xsharedabcd <- merge(xsharedabc, dshared, all=TRUE)
  xsharedabcd[is.na(xsharedabcd)] <- 0
  xsharedabcde <- merge(xsharedabcd, eshared, all=TRUE)
  xsharedabcde[is.na(xsharedabcde)] <- 0
  xsharedabcdef <- merge(xsharedabcde, fshared, all=TRUE)
  xsharedabcdef[is.na(xsharedabcdef)] <- 0
  xsharedabcdefg <- merge(xsharedabcdef, gshared, all=TRUE)
  xsharedabcdefg[is.na(xsharedabcdefg)] <- 0
  xsharedabcdefgh <- merge(xsharedabcdefg, hshared, all=TRUE)
  xsharedabcdefgh[is.na(xsharedabcdefgh)] <- 0
  xsharedabcdefghi <- merge(xsharedabcdefgh, ishared, all=TRUE)
  xsharedabcdefghi[is.na(xsharedabcdefghi)] <- 0
  xsharedabcdefghij <- merge(xsharedabcdefghi, jshared, all=TRUE)
  xsharedabcdefghij[is.na(xsharedabcdefghij)] <- 0
  xsharedabcdefghijk <- merge(xsharedabcdefghij, kshared, all=TRUE)
  xsharedabcdefghijk[is.na(xsharedabcdefghijk)] <- 0
  xsharedabcdefghijkl <- merge(xsharedabcdefghijk, lshared, all=TRUE)
  xsharedabcdefghijkl[is.na(xsharedabcdefghijkl)] <- 0
  xsharedabcdefghijklm <- merge(xsharedabcdefghijkl, mshared, all=TRUE)
  xsharedabcdefghijklm[is.na(xsharedabcdefghijklm)] <- 0
  xsharedabcdefghijklmn <- merge(xsharedabcdefghijklm, nshared, all=TRUE)
  xsharedabcdefghijklmn[is.na(xsharedabcdefghijklmn)] <- 0
  xsharedabcdefghijklmno <- merge(xsharedabcdefghijklmn, oshared, all=TRUE)
  xsharedabcdefghijklmno[is.na(xsharedabcdefghijklmno)] <- 0
  
  RORHeatmap <- xsharedabcdefghijklmno
  row.names(RORHeatmap) <- RORHeatmap$CDR3.AminoAcid
  RORHeatmap <- RORHeatmap[,2:16]
  RORHeatmap_matrix <- data.matrix(RORHeatmap)
  
  ROR999_heatmap <- heatmap.2(RORHeatmap_matrix, col=RdYlBu, scale="row", margins = c(8,5),labRow="", Colv=FALSE,
                              trace="none", key.title=NA, key.ylab=NA, key.xlab=NA, density.info="none", keysize=1)
}
aa.RORheatmap999(ROR999, TCRdata999)

#Heatmap: Subject 999 - AIM clones over time and subset
#Enter AIMs for 999 as "x" and all in-sequence data for 999 as "y"
aa.AIMheatmap999 <- function(x,y) {
  a <- subset(y, Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  b <- subset(y, Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  c <- subset(y, Year==1314 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  d <- subset(y, Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  e <- subset(y, Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  f <- subset(y, Year==1415 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  g <- subset(y, Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  h <- subset(y, Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  i <- subset(y, Year==1415 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  j <- subset(y, Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  k <- subset(y, Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  l <- subset(y, Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  m <- subset(y, Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  n <- subset(y, Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  o <- subset(y, Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  
  
  aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)   #condense rows in a with same amino acid sequence but different nucleotide sequences
  arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
  aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
  adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)
  aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
  anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  amerge <- rbind(anotdup, aaggdata)
  amergeorder <- amerge[order(amerge$Frequency), ]    #list of unique aa sequences in a
  
  bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)   #condense rows in b with same amino acid sequence but different nucleotide sequences
  browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
  bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
  bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)
  baggdataorder <- baggdata[order(baggdata$Frequency), ]
  bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  bmerge <- rbind(bnotdup, baggdata)
  bmergeorder <- bmerge[order(bmerge$Frequency), ]    #list of unique aa sequences in b
  
  cduplicatedAAtf <- duplicated(c$CDR3.AminoAcid)   #condense rows in c with same amino acid sequence but different nucleotide sequences
  crowduplicated <- which(cduplicatedAAtf, cduplicatedAAtf=="TRUE")
  cduplicatedAA <- c[crowduplicated, "CDR3.AminoAcid"]
  cdupall <- c[is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  caggdata <- aggregate(Frequency~CDR3.AminoAcid, data=cdupall, FUN=sum)
  caggdataorder <- caggdata[order(caggdata$Frequency), ]
  cnotdup <- c[!is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  cmerge <- rbind(cnotdup, caggdata)
  cmergeorder <- cmerge[order(cmerge$Frequency), ]    #list of unique aa sequences in c
  
  dduplicatedAAtf <- duplicated(d$CDR3.AminoAcid)   #condense rows in d with same amino acid sequence but different nucleotide sequences
  drowduplicated <- which(dduplicatedAAtf, dduplicatedAAtf=="TRUE")
  dduplicatedAA <- d[drowduplicated, "CDR3.AminoAcid"]
  ddupall <- d[is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  daggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ddupall, FUN=sum)
  daggdataorder <- daggdata[order(daggdata$Frequency), ]
  dnotdup <- d[!is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  dmerge <- rbind(dnotdup, daggdata)
  dmergeorder <- dmerge[order(dmerge$Frequency), ]    #list of unique aa sequences in d
  
  eduplicatedAAtf <- duplicated(e$CDR3.AminoAcid)   #condense rows in e with same amino acid sequence but different nucleotide sequences
  erowduplicated <- which(eduplicatedAAtf, eduplicatedAAtf=="TRUE")
  eduplicatedAA <- e[erowduplicated, "CDR3.AminoAcid"]
  edupall <- e[is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  eaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=edupall, FUN=sum)
  eaggdataorder <- eaggdata[order(eaggdata$Frequency), ]
  enotdup <- e[!is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  emerge <- rbind(enotdup, eaggdata)
  emergeorder <- emerge[order(emerge$Frequency), ]    #list of unique aa sequences in e
  
  fduplicatedAAtf <- duplicated(f$CDR3.AminoAcid)   #condense rows in f with same amino acid sequence but different nucleotide sequences
  frowduplicated <- which(fduplicatedAAtf, fduplicatedAAtf=="TRUE")
  fduplicatedAA <- f[frowduplicated, "CDR3.AminoAcid"]
  fdupall <- f[is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  faggdata <- aggregate(Frequency~CDR3.AminoAcid, data=fdupall, FUN=sum)
  faggdataorder <- faggdata[order(faggdata$Frequency), ]
  fnotdup <- f[!is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  fmerge <- rbind(fnotdup, faggdata)
  fmergeorder <- fmerge[order(fmerge$Frequency), ]    #list of unique aa sequences in f
  
  gduplicatedAAtf <- duplicated(g$CDR3.AminoAcid)   #condense rows in g with same amino acid sequence but different nucleotide sequences
  growduplicated <- which(gduplicatedAAtf, gduplicatedAAtf=="TRUE")
  gduplicatedAA <- g[growduplicated, "CDR3.AminoAcid"]
  gdupall <- g[is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=gdupall, FUN=sum)
  gaggdataorder <- gaggdata[order(gaggdata$Frequency), ]
  gnotdup <- g[!is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gmerge <- rbind(gnotdup, gaggdata)
  gmergeorder <- gmerge[order(gmerge$Frequency), ]    #list of unique aa sequences in g
  
  hduplicatedAAtf <- duplicated(h$CDR3.AminoAcid)   #condense rows in h with same amino acid sequence but different nucleotide sequences
  hrowduplicated <- which(hduplicatedAAtf, hduplicatedAAtf=="TRUE")
  hduplicatedAA <- h[hrowduplicated, "CDR3.AminoAcid"]
  hdupall <- h[is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  haggdata <- aggregate(Frequency~CDR3.AminoAcid, data=hdupall, FUN=sum)
  haggdataorder <- haggdata[order(haggdata$Frequency), ]
  hnotdup <- h[!is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  hmerge <- rbind(hnotdup, haggdata)
  hmergeorder <- hmerge[order(hmerge$Frequency), ]    #list of unique aa sequences in h
  
  iduplicatedAAtf <- duplicated(i$CDR3.AminoAcid)   #condense rows in i with same amino acid sequence but different nucleotide sequences
  irowduplicated <- which(iduplicatedAAtf, iduplicatedAAtf=="TRUE")
  iduplicatedAA <- i[irowduplicated, "CDR3.AminoAcid"]
  idupall <- i[is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  iaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=idupall, FUN=sum)
  iaggdataorder <- iaggdata[order(iaggdata$Frequency), ]
  inotdup <- i[!is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  imerge <- rbind(inotdup, iaggdata)
  imergeorder <- imerge[order(imerge$Frequency), ]    #list of unique aa sequences in i
  
  jduplicatedAAtf <- duplicated(j$CDR3.AminoAcid)   #condense rows in j with same amino acid sequence but different nucleotide sequences
  jrowduplicated <- which(jduplicatedAAtf, jduplicatedAAtf=="TRUE")
  jduplicatedAA <- j[jrowduplicated, "CDR3.AminoAcid"]
  jdupall <- j[is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=jdupall, FUN=sum)
  jaggdataorder <- jaggdata[order(jaggdata$Frequency), ]
  jnotdup <- j[!is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jmerge <- rbind(jnotdup, jaggdata)
  jmergeorder <- jmerge[order(jmerge$Frequency), ]    #list of unique aa sequences in j
  
  kduplicatedAAtf <- duplicated(k$CDR3.AminoAcid)   #condense rows in k with same amino acid sequence but different nucleotide sequences
  krowduplicated <- which(kduplicatedAAtf, kduplicatedAAtf=="TRUE")
  kduplicatedAA <- k[krowduplicated, "CDR3.AminoAcid"]
  kdupall <- k[is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=kdupall, FUN=sum)
  kaggdataorder <- kaggdata[order(kaggdata$Frequency), ]
  knotdup <- k[!is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kmerge <- rbind(knotdup, kaggdata)
  kmergeorder <- kmerge[order(kmerge$Frequency), ]    #list of unique aa sequences in k
  
  lduplicatedAAtf <- duplicated(l$CDR3.AminoAcid)   #condense rows in l with same amino acid sequence but different nucleotide sequences
  lrowduplicated <- which(lduplicatedAAtf, lduplicatedAAtf=="TRUE")
  lduplicatedAA <- l[lrowduplicated, "CDR3.AminoAcid"]
  ldupall <- l[is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  laggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ldupall, FUN=sum)
  laggdataorder <- laggdata[order(laggdata$Frequency), ]
  lnotdup <- l[!is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  lmerge <- rbind(lnotdup, laggdata)
  lmergeorder <- lmerge[order(lmerge$Frequency), ]    #list of unique aa sequences in l
  
  mduplicatedAAtf <- duplicated(m$CDR3.AminoAcid)   #condense rows in m with same amino acid sequence but different nucleotide sequences
  mrowduplicated <- which(mduplicatedAAtf, mduplicatedAAtf=="TRUE")
  mduplicatedAA <- m[mrowduplicated, "CDR3.AminoAcid"]
  mdupall <- m[is.element(m$CDR3.AminoAcid, mduplicatedAA), ]
  maggdata <- aggregate(Frequency~CDR3.AminoAcid, data=mdupall, FUN=sum)
  maggdataorder <- maggdata[order(maggdata$Frequency), ]
  mnotdup <- m[!is.element(m$CDR3.AminoAcid, mduplicatedAA), ]
  mmerge <- rbind(mnotdup, maggdata)
  mmergeorder <- mmerge[order(mmerge$Frequency), ]    #list of unique aa sequences in m
  
  nduplicatedAAtf <- duplicated(n$CDR3.AminoAcid)   #condense rows in n with same amino acid sequence but different nucleotide sequences
  nrowduplicated <- which(nduplicatedAAtf, nduplicatedAAtf=="TRUE")
  nduplicatedAA <- n[nrowduplicated, "CDR3.AminoAcid"]
  ndupall <- n[is.element(n$CDR3.AminoAcid, nduplicatedAA), ]
  naggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ndupall, FUN=sum)
  naggdataorder <- naggdata[order(naggdata$Frequency), ]
  nnotdup <- n[!is.element(n$CDR3.AminoAcid, nduplicatedAA), ]
  nmerge <- rbind(nnotdup, naggdata)
  nmergeorder <- nmerge[order(nmerge$Frequency), ]    #list of unique aa sequences in n
  
  oduplicatedAAtf <- duplicated(o$CDR3.AminoAcid)   #condense rows in o with same amino acid sequence but different nucleotide sequences
  orowduplicated <- which(oduplicatedAAtf, oduplicatedAAtf=="TRUE")
  oduplicatedAA <- o[orowduplicated, "CDR3.AminoAcid"]
  odupall <- o[is.element(o$CDR3.AminoAcid, oduplicatedAA), ]
  oaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=odupall, FUN=sum)
  oaggdataorder <- oaggdata[order(oaggdata$Frequency), ]
  onotdup <- o[!is.element(o$CDR3.AminoAcid, oduplicatedAA), ]
  omerge <- rbind(onotdup, oaggdata)
  omergeorder <- omerge[order(omerge$Frequency), ]    #list of unique aa sequences in o
  
  xAA <- x[[1]]
  ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, xAA), ] #find overlapping sequences between AIM and each subset
  bshared <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, xAA), ]
  cshared <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, xAA), ]
  dshared <- dmergeorder[is.element(dmergeorder$CDR3.AminoAcid, xAA), ]
  eshared <- emergeorder[is.element(emergeorder$CDR3.AminoAcid, xAA), ]
  fshared <- fmergeorder[is.element(fmergeorder$CDR3.AminoAcid, xAA), ]
  gshared <- gmergeorder[is.element(gmergeorder$CDR3.AminoAcid, xAA), ]
  hshared <- hmergeorder[is.element(hmergeorder$CDR3.AminoAcid, xAA), ]
  ishared <- imergeorder[is.element(imergeorder$CDR3.AminoAcid, xAA), ]
  jshared <- jmergeorder[is.element(jmergeorder$CDR3.AminoAcid, xAA), ]
  kshared <- kmergeorder[is.element(kmergeorder$CDR3.AminoAcid, xAA), ]
  lshared <- lmergeorder[is.element(lmergeorder$CDR3.AminoAcid, xAA), ]
  mshared <- mmergeorder[is.element(mmergeorder$CDR3.AminoAcid, xAA), ]
  nshared <- nmergeorder[is.element(nmergeorder$CDR3.AminoAcid, xAA), ]
  oshared <- omergeorder[is.element(omergeorder$CDR3.AminoAcid, xAA), ]
  colnames(ashared)[2] <- "Hi.Y1D7"
  colnames(bshared)[2] <- "Lo.Y1D7"
  colnames(cshared)[2] <- "X5.Y1D7"
  colnames(dshared)[2] <- "Hi.Y2D0"
  colnames(eshared)[2] <- "Lo.Y2D0"
  colnames(fshared)[2] <- "X5.Y2D0"
  colnames(gshared)[2] <- "Hi.Y2D7"
  colnames(hshared)[2] <- "Lo.Y2D7"
  colnames(ishared)[2] <- "X5.Y2D7"
  colnames(jshared)[2] <- "Hi.Y3D0"
  colnames(kshared)[2] <- "Lo.Y3D0"
  colnames(lshared)[2] <- "X5.Y3D0"
  colnames(mshared)[2] <- "Hi.Y3D7"
  colnames(nshared)[2] <- "Lo.Y3D7"
  colnames(oshared)[2] <- "X5.Y3D7"
  
  xshareda <- merge(x, ashared, all=TRUE) #merge all data into one data frame
  xshareda[is.na(xshareda)] <- 0
  xsharedab <- merge(xshareda, bshared, all=TRUE)
  xsharedab[is.na(xsharedab)] <- 0
  xsharedabc <- merge(xsharedab, cshared, all=TRUE)
  xsharedabc[is.na(xsharedabc)] <- 0
  xsharedabcd <- merge(xsharedabc, dshared, all=TRUE)
  xsharedabcd[is.na(xsharedabcd)] <- 0
  xsharedabcde <- merge(xsharedabcd, eshared, all=TRUE)
  xsharedabcde[is.na(xsharedabcde)] <- 0
  xsharedabcdef <- merge(xsharedabcde, fshared, all=TRUE)
  xsharedabcdef[is.na(xsharedabcdef)] <- 0
  xsharedabcdefg <- merge(xsharedabcdef, gshared, all=TRUE)
  xsharedabcdefg[is.na(xsharedabcdefg)] <- 0
  xsharedabcdefgh <- merge(xsharedabcdefg, hshared, all=TRUE)
  xsharedabcdefgh[is.na(xsharedabcdefgh)] <- 0
  xsharedabcdefghi <- merge(xsharedabcdefgh, ishared, all=TRUE)
  xsharedabcdefghi[is.na(xsharedabcdefghi)] <- 0
  xsharedabcdefghij <- merge(xsharedabcdefghi, jshared, all=TRUE)
  xsharedabcdefghij[is.na(xsharedabcdefghij)] <- 0
  xsharedabcdefghijk <- merge(xsharedabcdefghij, kshared, all=TRUE)
  xsharedabcdefghijk[is.na(xsharedabcdefghijk)] <- 0
  xsharedabcdefghijkl <- merge(xsharedabcdefghijk, lshared, all=TRUE)
  xsharedabcdefghijkl[is.na(xsharedabcdefghijkl)] <- 0
  xsharedabcdefghijklm <- merge(xsharedabcdefghijkl, mshared, all=TRUE)
  xsharedabcdefghijklm[is.na(xsharedabcdefghijklm)] <- 0
  xsharedabcdefghijklmn <- merge(xsharedabcdefghijklm, nshared, all=TRUE)
  xsharedabcdefghijklmn[is.na(xsharedabcdefghijklmn)] <- 0
  xsharedabcdefghijklmno <- merge(xsharedabcdefghijklmn, oshared, all=TRUE)
  xsharedabcdefghijklmno[is.na(xsharedabcdefghijklmno)] <- 0
  
  
  AIMHeatmap <- xsharedabcdefghijklmno
  AIMHeatmap2 <- AIMHeatmap[rowSums(AIMHeatmap[,-1])>0, ]
  
  row.names(AIMHeatmap2) <- AIMHeatmap2$CDR3.AminoAcid
  AIMHeatmap2 <- AIMHeatmap2[,2:16]
  AIMHeatmap2_matrix <- data.matrix(AIMHeatmap2)
  
  AIM999_heatmap <- heatmap.2(AIMHeatmap2_matrix, col=RdYlBu, scale="row", margins = c(8,5),labRow="", Colv=FALSE,
                              trace="none", key.title=NA, key.ylab=NA, key.xlab=NA, density.info="none", keysize=1)
}
aa.AIMheatmap999(aim999, TCRdata999)


TCRdata101 <- subset(TCRdataMultYearsIn, Subject==101, select=c(Year, Cell, Day, CDR3.AminoAcid, Frequency))
ROR101 <- subset(Fluspecific, Subject==101, select=c(CDR3.AminoAcid)) #recurring oligoclonal response clonotypes
aim101 <- subset(AIMs, Subject==101, select=c(CDR3.AminoAcid)) # AIM clonotypes

#Heatmap: Subject 101 - Recurrent oligoclonal response clones over time and subset
#Enter ROR for 101 as "x" and all in-sequence data for 999 as "y"
aa.RORheatmap101 <- function(x,y) {
  a <- subset(y, Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  b <- subset(y, Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  c <- subset(y, Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  d <- subset(y, Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  e <- subset(y, Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  f <- subset(y, Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  g <- subset(y, Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  h <- subset(y, Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  i <- subset(y, Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  j <- subset(y, Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  k <- subset(y, Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  l <- subset(y, Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  
  aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)   #condense rows in a with same amino acid sequence but different nucleotide sequences
  arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
  aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
  adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)
  aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
  anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  amerge <- rbind(anotdup, aaggdata)
  amergeorder <- amerge[order(amerge$Frequency), ]    #list of unique aa sequences in a
  
  bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)   #condense rows in b with same amino acid sequence but different nucleotide sequences
  browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
  bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
  bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)
  baggdataorder <- baggdata[order(baggdata$Frequency), ]
  bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  bmerge <- rbind(bnotdup, baggdata)
  bmergeorder <- bmerge[order(bmerge$Frequency), ]    #list of unique aa sequences in b
  
  cduplicatedAAtf <- duplicated(c$CDR3.AminoAcid)   #condense rows in c with same amino acid sequence but different nucleotide sequences
  crowduplicated <- which(cduplicatedAAtf, cduplicatedAAtf=="TRUE")
  cduplicatedAA <- c[crowduplicated, "CDR3.AminoAcid"]
  cdupall <- c[is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  caggdata <- aggregate(Frequency~CDR3.AminoAcid, data=cdupall, FUN=sum)
  caggdataorder <- caggdata[order(caggdata$Frequency), ]
  cnotdup <- c[!is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  cmerge <- rbind(cnotdup, caggdata)
  cmergeorder <- cmerge[order(cmerge$Frequency), ]    #list of unique aa sequences in c
  
  dduplicatedAAtf <- duplicated(d$CDR3.AminoAcid)   #condense rows in d with same amino acid sequence but different nucleotide sequences
  drowduplicated <- which(dduplicatedAAtf, dduplicatedAAtf=="TRUE")
  dduplicatedAA <- d[drowduplicated, "CDR3.AminoAcid"]
  ddupall <- d[is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  daggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ddupall, FUN=sum)
  daggdataorder <- daggdata[order(daggdata$Frequency), ]
  dnotdup <- d[!is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  dmerge <- rbind(dnotdup, daggdata)
  dmergeorder <- dmerge[order(dmerge$Frequency), ]    #list of unique aa sequences in d
  
  eduplicatedAAtf <- duplicated(e$CDR3.AminoAcid)   #condense rows in e with same amino acid sequence but different nucleotide sequences
  erowduplicated <- which(eduplicatedAAtf, eduplicatedAAtf=="TRUE")
  eduplicatedAA <- e[erowduplicated, "CDR3.AminoAcid"]
  edupall <- e[is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  eaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=edupall, FUN=sum)
  eaggdataorder <- eaggdata[order(eaggdata$Frequency), ]
  enotdup <- e[!is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  emerge <- rbind(enotdup, eaggdata)
  emergeorder <- emerge[order(emerge$Frequency), ]    #list of unique aa sequences in e
  
  fduplicatedAAtf <- duplicated(f$CDR3.AminoAcid)   #condense rows in f with same amino acid sequence but different nucleotide sequences
  frowduplicated <- which(fduplicatedAAtf, fduplicatedAAtf=="TRUE")
  fduplicatedAA <- f[frowduplicated, "CDR3.AminoAcid"]
  fdupall <- f[is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  faggdata <- aggregate(Frequency~CDR3.AminoAcid, data=fdupall, FUN=sum)
  faggdataorder <- faggdata[order(faggdata$Frequency), ]
  fnotdup <- f[!is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  fmerge <- rbind(fnotdup, faggdata)
  fmergeorder <- fmerge[order(fmerge$Frequency), ]    #list of unique aa sequences in f
  
  gduplicatedAAtf <- duplicated(g$CDR3.AminoAcid)   #condense rows in g with same amino acid sequence but different nucleotide sequences
  growduplicated <- which(gduplicatedAAtf, gduplicatedAAtf=="TRUE")
  gduplicatedAA <- g[growduplicated, "CDR3.AminoAcid"]
  gdupall <- g[is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=gdupall, FUN=sum)
  gaggdataorder <- gaggdata[order(gaggdata$Frequency), ]
  gnotdup <- g[!is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gmerge <- rbind(gnotdup, gaggdata)
  gmergeorder <- gmerge[order(gmerge$Frequency), ]    #list of unique aa sequences in g
  
  hduplicatedAAtf <- duplicated(h$CDR3.AminoAcid)   #condense rows in h with same amino acid sequence but different nucleotide sequences
  hrowduplicated <- which(hduplicatedAAtf, hduplicatedAAtf=="TRUE")
  hduplicatedAA <- h[hrowduplicated, "CDR3.AminoAcid"]
  hdupall <- h[is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  haggdata <- aggregate(Frequency~CDR3.AminoAcid, data=hdupall, FUN=sum)
  haggdataorder <- haggdata[order(haggdata$Frequency), ]
  hnotdup <- h[!is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  hmerge <- rbind(hnotdup, haggdata)
  hmergeorder <- hmerge[order(hmerge$Frequency), ]    #list of unique aa sequences in h
  
  iduplicatedAAtf <- duplicated(i$CDR3.AminoAcid)   #condense rows in i with same amino acid sequence but different nucleotide sequences
  irowduplicated <- which(iduplicatedAAtf, iduplicatedAAtf=="TRUE")
  iduplicatedAA <- i[irowduplicated, "CDR3.AminoAcid"]
  idupall <- i[is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  iaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=idupall, FUN=sum)
  iaggdataorder <- iaggdata[order(iaggdata$Frequency), ]
  inotdup <- i[!is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  imerge <- rbind(inotdup, iaggdata)
  imergeorder <- imerge[order(imerge$Frequency), ]    #list of unique aa sequences in i
  
  jduplicatedAAtf <- duplicated(j$CDR3.AminoAcid)   #condense rows in j with same amino acid sequence but different nucleotide sequences
  jrowduplicated <- which(jduplicatedAAtf, jduplicatedAAtf=="TRUE")
  jduplicatedAA <- j[jrowduplicated, "CDR3.AminoAcid"]
  jdupall <- j[is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=jdupall, FUN=sum)
  jaggdataorder <- jaggdata[order(jaggdata$Frequency), ]
  jnotdup <- j[!is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jmerge <- rbind(jnotdup, jaggdata)
  jmergeorder <- jmerge[order(jmerge$Frequency), ]    #list of unique aa sequences in j
  
  kduplicatedAAtf <- duplicated(k$CDR3.AminoAcid)   #condense rows in k with same amino acid sequence but different nucleotide sequences
  krowduplicated <- which(kduplicatedAAtf, kduplicatedAAtf=="TRUE")
  kduplicatedAA <- k[krowduplicated, "CDR3.AminoAcid"]
  kdupall <- k[is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=kdupall, FUN=sum)
  kaggdataorder <- kaggdata[order(kaggdata$Frequency), ]
  knotdup <- k[!is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kmerge <- rbind(knotdup, kaggdata)
  kmergeorder <- kmerge[order(kmerge$Frequency), ]    #list of unique aa sequences in k
  
  lduplicatedAAtf <- duplicated(l$CDR3.AminoAcid)   #condense rows in l with same amino acid sequence but different nucleotide sequences
  lrowduplicated <- which(lduplicatedAAtf, lduplicatedAAtf=="TRUE")
  lduplicatedAA <- l[lrowduplicated, "CDR3.AminoAcid"]
  ldupall <- l[is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  laggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ldupall, FUN=sum)
  laggdataorder <- laggdata[order(laggdata$Frequency), ]
  lnotdup <- l[!is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  lmerge <- rbind(lnotdup, laggdata)
  lmergeorder <- lmerge[order(lmerge$Frequency), ]    #list of unique aa sequences in l
  
  xAA <- x[[1]]
  ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, xAA), ] #find overlapping sequences between ROR and each subset
  bshared <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, xAA), ]
  cshared <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, xAA), ]
  dshared <- dmergeorder[is.element(dmergeorder$CDR3.AminoAcid, xAA), ]
  eshared <- emergeorder[is.element(emergeorder$CDR3.AminoAcid, xAA), ]
  fshared <- fmergeorder[is.element(fmergeorder$CDR3.AminoAcid, xAA), ]
  gshared <- gmergeorder[is.element(gmergeorder$CDR3.AminoAcid, xAA), ]
  hshared <- hmergeorder[is.element(hmergeorder$CDR3.AminoAcid, xAA), ]
  ishared <- imergeorder[is.element(imergeorder$CDR3.AminoAcid, xAA), ]
  jshared <- jmergeorder[is.element(jmergeorder$CDR3.AminoAcid, xAA), ]
  kshared <- kmergeorder[is.element(kmergeorder$CDR3.AminoAcid, xAA), ]
  lshared <- lmergeorder[is.element(lmergeorder$CDR3.AminoAcid, xAA), ]
  colnames(ashared)[2] <- "Hi.Y1D7"
  colnames(bshared)[2] <- "Lo.Y1D7"
  colnames(cshared)[2] <- "Hi.Y2D0"
  colnames(dshared)[2] <- "Lo.Y2D0"
  colnames(eshared)[2] <- "Hi.Y2D7"
  colnames(fshared)[2] <- "Lo.Y2D7"
  colnames(gshared)[2] <- "Hi.Y3D0"
  colnames(hshared)[2] <- "Lo.Y3D0"
  colnames(ishared)[2] <- "X5.Y3D0"
  colnames(jshared)[2] <- "Hi.Y3D7"
  colnames(kshared)[2] <- "Lo.Y3D7"
  colnames(lshared)[2] <- "X5.Y3D7"

  xshareda <- merge(x, ashared, all=TRUE) #merge all data into one data frame
  xshareda[is.na(xshareda)] <- 0
  xsharedab <- merge(xshareda, bshared, all=TRUE)
  xsharedab[is.na(xsharedab)] <- 0
  xsharedabc <- merge(xsharedab, cshared, all=TRUE)
  xsharedabc[is.na(xsharedabc)] <- 0
  xsharedabcd <- merge(xsharedabc, dshared, all=TRUE)
  xsharedabcd[is.na(xsharedabcd)] <- 0
  xsharedabcde <- merge(xsharedabcd, eshared, all=TRUE)
  xsharedabcde[is.na(xsharedabcde)] <- 0
  xsharedabcdef <- merge(xsharedabcde, fshared, all=TRUE)
  xsharedabcdef[is.na(xsharedabcdef)] <- 0
  xsharedabcdefg <- merge(xsharedabcdef, gshared, all=TRUE)
  xsharedabcdefg[is.na(xsharedabcdefg)] <- 0
  xsharedabcdefgh <- merge(xsharedabcdefg, hshared, all=TRUE)
  xsharedabcdefgh[is.na(xsharedabcdefgh)] <- 0
  xsharedabcdefghi <- merge(xsharedabcdefgh, ishared, all=TRUE)
  xsharedabcdefghi[is.na(xsharedabcdefghi)] <- 0
  xsharedabcdefghij <- merge(xsharedabcdefghi, jshared, all=TRUE)
  xsharedabcdefghij[is.na(xsharedabcdefghij)] <- 0
  xsharedabcdefghijk <- merge(xsharedabcdefghij, kshared, all=TRUE)
  xsharedabcdefghijk[is.na(xsharedabcdefghijk)] <- 0
  xsharedabcdefghijkl <- merge(xsharedabcdefghijk, lshared, all=TRUE)
  xsharedabcdefghijkl[is.na(xsharedabcdefghijkl)] <- 0

  RORHeatmap <- xsharedabcdefghijkl
  row.names(RORHeatmap) <- RORHeatmap$CDR3.AminoAcid
  RORHeatmap <- RORHeatmap[,2:13]
  RORHeatmap_matrix <- data.matrix(RORHeatmap)
  
  ROR101_heatmap <- heatmap.2(RORHeatmap_matrix, col=RdYlBu, scale="row", margins = c(8,5),labRow="", Colv=FALSE,
                              trace="none", key.title=NA, key.ylab=NA, key.xlab=NA, density.info="none", keysize=1)
}
aa.RORheatmap101(ROR101, TCRdata101)

#Heatmap: Subject 101 - AIM clones over time and subset
#Enter AIMs for 101 as "x" and all in-sequence data for 999 as "y"
aa.AIMheatmap101 <- function(x,y) {
  a <- subset(y, Year==1314 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  b <- subset(y, Year==1314 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  c <- subset(y, Year==1415 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  d <- subset(y, Year==1415 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  e <- subset(y, Year==1415 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  f <- subset(y, Year==1415 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  g <- subset(y, Year==1516 & Day==0 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  h <- subset(y, Year==1516 & Day==0 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  i <- subset(y, Year==1516 & Day==0 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  j <- subset(y, Year==1516 & Day==7 & Cell=="ICOS+CD38+", select=c(CDR3.AminoAcid, Frequency))
  k <- subset(y, Year==1516 & Day==7 & Cell=="ICOS-CD38-", select=c(CDR3.AminoAcid, Frequency))
  l <- subset(y, Year==1516 & Day==7 & Cell=="CXCR5-", select=c(CDR3.AminoAcid, Frequency))
  
  aduplicatedAAtf <- duplicated(a$CDR3.AminoAcid)   #condense rows in a with same amino acid sequence but different nucleotide sequences
  arowduplicated <- which(aduplicatedAAtf, aduplicatedAAtf=="TRUE")
  aduplicatedAA <- a[arowduplicated, "CDR3.AminoAcid"]
  adupall <- a[is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  aaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=adupall, FUN=sum)
  aaggdataorder <- aaggdata[order(aaggdata$Frequency), ]
  anotdup <- a[!is.element(a$CDR3.AminoAcid, aduplicatedAA), ]
  amerge <- rbind(anotdup, aaggdata)
  amergeorder <- amerge[order(amerge$Frequency), ]    #list of unique aa sequences in a
  
  bduplicatedAAtf <- duplicated(b$CDR3.AminoAcid)   #condense rows in b with same amino acid sequence but different nucleotide sequences
  browduplicated <- which(bduplicatedAAtf, bduplicatedAAtf=="TRUE")
  bduplicatedAA <- b[browduplicated, "CDR3.AminoAcid"]
  bdupall <- b[is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  baggdata <- aggregate(Frequency~CDR3.AminoAcid, data=bdupall, FUN=sum)
  baggdataorder <- baggdata[order(baggdata$Frequency), ]
  bnotdup <- b[!is.element(b$CDR3.AminoAcid, bduplicatedAA), ]
  bmerge <- rbind(bnotdup, baggdata)
  bmergeorder <- bmerge[order(bmerge$Frequency), ]    #list of unique aa sequences in b
  
  cduplicatedAAtf <- duplicated(c$CDR3.AminoAcid)   #condense rows in c with same amino acid sequence but different nucleotide sequences
  crowduplicated <- which(cduplicatedAAtf, cduplicatedAAtf=="TRUE")
  cduplicatedAA <- c[crowduplicated, "CDR3.AminoAcid"]
  cdupall <- c[is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  caggdata <- aggregate(Frequency~CDR3.AminoAcid, data=cdupall, FUN=sum)
  caggdataorder <- caggdata[order(caggdata$Frequency), ]
  cnotdup <- c[!is.element(c$CDR3.AminoAcid, cduplicatedAA), ]
  cmerge <- rbind(cnotdup, caggdata)
  cmergeorder <- cmerge[order(cmerge$Frequency), ]    #list of unique aa sequences in c
  
  dduplicatedAAtf <- duplicated(d$CDR3.AminoAcid)   #condense rows in d with same amino acid sequence but different nucleotide sequences
  drowduplicated <- which(dduplicatedAAtf, dduplicatedAAtf=="TRUE")
  dduplicatedAA <- d[drowduplicated, "CDR3.AminoAcid"]
  ddupall <- d[is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  daggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ddupall, FUN=sum)
  daggdataorder <- daggdata[order(daggdata$Frequency), ]
  dnotdup <- d[!is.element(d$CDR3.AminoAcid, dduplicatedAA), ]
  dmerge <- rbind(dnotdup, daggdata)
  dmergeorder <- dmerge[order(dmerge$Frequency), ]    #list of unique aa sequences in d
  
  eduplicatedAAtf <- duplicated(e$CDR3.AminoAcid)   #condense rows in e with same amino acid sequence but different nucleotide sequences
  erowduplicated <- which(eduplicatedAAtf, eduplicatedAAtf=="TRUE")
  eduplicatedAA <- e[erowduplicated, "CDR3.AminoAcid"]
  edupall <- e[is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  eaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=edupall, FUN=sum)
  eaggdataorder <- eaggdata[order(eaggdata$Frequency), ]
  enotdup <- e[!is.element(e$CDR3.AminoAcid, eduplicatedAA), ]
  emerge <- rbind(enotdup, eaggdata)
  emergeorder <- emerge[order(emerge$Frequency), ]    #list of unique aa sequences in e
  
  fduplicatedAAtf <- duplicated(f$CDR3.AminoAcid)   #condense rows in f with same amino acid sequence but different nucleotide sequences
  frowduplicated <- which(fduplicatedAAtf, fduplicatedAAtf=="TRUE")
  fduplicatedAA <- f[frowduplicated, "CDR3.AminoAcid"]
  fdupall <- f[is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  faggdata <- aggregate(Frequency~CDR3.AminoAcid, data=fdupall, FUN=sum)
  faggdataorder <- faggdata[order(faggdata$Frequency), ]
  fnotdup <- f[!is.element(f$CDR3.AminoAcid, fduplicatedAA), ]
  fmerge <- rbind(fnotdup, faggdata)
  fmergeorder <- fmerge[order(fmerge$Frequency), ]    #list of unique aa sequences in f
  
  gduplicatedAAtf <- duplicated(g$CDR3.AminoAcid)   #condense rows in g with same amino acid sequence but different nucleotide sequences
  growduplicated <- which(gduplicatedAAtf, gduplicatedAAtf=="TRUE")
  gduplicatedAA <- g[growduplicated, "CDR3.AminoAcid"]
  gdupall <- g[is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=gdupall, FUN=sum)
  gaggdataorder <- gaggdata[order(gaggdata$Frequency), ]
  gnotdup <- g[!is.element(g$CDR3.AminoAcid, gduplicatedAA), ]
  gmerge <- rbind(gnotdup, gaggdata)
  gmergeorder <- gmerge[order(gmerge$Frequency), ]    #list of unique aa sequences in g
  
  hduplicatedAAtf <- duplicated(h$CDR3.AminoAcid)   #condense rows in h with same amino acid sequence but different nucleotide sequences
  hrowduplicated <- which(hduplicatedAAtf, hduplicatedAAtf=="TRUE")
  hduplicatedAA <- h[hrowduplicated, "CDR3.AminoAcid"]
  hdupall <- h[is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  haggdata <- aggregate(Frequency~CDR3.AminoAcid, data=hdupall, FUN=sum)
  haggdataorder <- haggdata[order(haggdata$Frequency), ]
  hnotdup <- h[!is.element(h$CDR3.AminoAcid, hduplicatedAA), ]
  hmerge <- rbind(hnotdup, haggdata)
  hmergeorder <- hmerge[order(hmerge$Frequency), ]    #list of unique aa sequences in h
  
  iduplicatedAAtf <- duplicated(i$CDR3.AminoAcid)   #condense rows in i with same amino acid sequence but different nucleotide sequences
  irowduplicated <- which(iduplicatedAAtf, iduplicatedAAtf=="TRUE")
  iduplicatedAA <- i[irowduplicated, "CDR3.AminoAcid"]
  idupall <- i[is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  iaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=idupall, FUN=sum)
  iaggdataorder <- iaggdata[order(iaggdata$Frequency), ]
  inotdup <- i[!is.element(i$CDR3.AminoAcid, iduplicatedAA), ]
  imerge <- rbind(inotdup, iaggdata)
  imergeorder <- imerge[order(imerge$Frequency), ]    #list of unique aa sequences in i
  
  jduplicatedAAtf <- duplicated(j$CDR3.AminoAcid)   #condense rows in j with same amino acid sequence but different nucleotide sequences
  jrowduplicated <- which(jduplicatedAAtf, jduplicatedAAtf=="TRUE")
  jduplicatedAA <- j[jrowduplicated, "CDR3.AminoAcid"]
  jdupall <- j[is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=jdupall, FUN=sum)
  jaggdataorder <- jaggdata[order(jaggdata$Frequency), ]
  jnotdup <- j[!is.element(j$CDR3.AminoAcid, jduplicatedAA), ]
  jmerge <- rbind(jnotdup, jaggdata)
  jmergeorder <- jmerge[order(jmerge$Frequency), ]    #list of unique aa sequences in j
  
  kduplicatedAAtf <- duplicated(k$CDR3.AminoAcid)   #condense rows in k with same amino acid sequence but different nucleotide sequences
  krowduplicated <- which(kduplicatedAAtf, kduplicatedAAtf=="TRUE")
  kduplicatedAA <- k[krowduplicated, "CDR3.AminoAcid"]
  kdupall <- k[is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kaggdata <- aggregate(Frequency~CDR3.AminoAcid, data=kdupall, FUN=sum)
  kaggdataorder <- kaggdata[order(kaggdata$Frequency), ]
  knotdup <- k[!is.element(k$CDR3.AminoAcid, kduplicatedAA), ]
  kmerge <- rbind(knotdup, kaggdata)
  kmergeorder <- kmerge[order(kmerge$Frequency), ]    #list of unique aa sequences in k
  
  lduplicatedAAtf <- duplicated(l$CDR3.AminoAcid)   #condense rows in l with same amino acid sequence but different nucleotide sequences
  lrowduplicated <- which(lduplicatedAAtf, lduplicatedAAtf=="TRUE")
  lduplicatedAA <- l[lrowduplicated, "CDR3.AminoAcid"]
  ldupall <- l[is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  laggdata <- aggregate(Frequency~CDR3.AminoAcid, data=ldupall, FUN=sum)
  laggdataorder <- laggdata[order(laggdata$Frequency), ]
  lnotdup <- l[!is.element(l$CDR3.AminoAcid, lduplicatedAA), ]
  lmerge <- rbind(lnotdup, laggdata)
  lmergeorder <- lmerge[order(lmerge$Frequency), ]    #list of unique aa sequences in l
  
  xAA <- x[[1]]
  ashared <- amergeorder[is.element(amergeorder$CDR3.AminoAcid, xAA), ] #find overlapping sequences between AIM and each subset
  bshared <- bmergeorder[is.element(bmergeorder$CDR3.AminoAcid, xAA), ]
  cshared <- cmergeorder[is.element(cmergeorder$CDR3.AminoAcid, xAA), ]
  dshared <- dmergeorder[is.element(dmergeorder$CDR3.AminoAcid, xAA), ]
  eshared <- emergeorder[is.element(emergeorder$CDR3.AminoAcid, xAA), ]
  fshared <- fmergeorder[is.element(fmergeorder$CDR3.AminoAcid, xAA), ]
  gshared <- gmergeorder[is.element(gmergeorder$CDR3.AminoAcid, xAA), ]
  hshared <- hmergeorder[is.element(hmergeorder$CDR3.AminoAcid, xAA), ]
  ishared <- imergeorder[is.element(imergeorder$CDR3.AminoAcid, xAA), ]
  jshared <- jmergeorder[is.element(jmergeorder$CDR3.AminoAcid, xAA), ]
  kshared <- kmergeorder[is.element(kmergeorder$CDR3.AminoAcid, xAA), ]
  lshared <- lmergeorder[is.element(lmergeorder$CDR3.AminoAcid, xAA), ]
  colnames(ashared)[2] <- "Hi.Y1D7"
  colnames(bshared)[2] <- "Lo.Y1D7"
  colnames(cshared)[2] <- "Hi.Y2D0"
  colnames(dshared)[2] <- "Lo.Y2D0"
  colnames(eshared)[2] <- "Hi.Y2D7"
  colnames(fshared)[2] <- "Lo.Y2D7"
  colnames(gshared)[2] <- "Hi.Y3D0"
  colnames(hshared)[2] <- "Lo.Y3D0"
  colnames(ishared)[2] <- "X5.Y3D0"
  colnames(jshared)[2] <- "Hi.Y3D7"
  colnames(kshared)[2] <- "Lo.Y3D7"
  colnames(lshared)[2] <- "X5.Y3D7"
  
  xshareda <- merge(x, ashared, all=TRUE) #merge all data into one data frame
  xshareda[is.na(xshareda)] <- 0
  xsharedab <- merge(xshareda, bshared, all=TRUE)
  xsharedab[is.na(xsharedab)] <- 0
  xsharedabc <- merge(xsharedab, cshared, all=TRUE)
  xsharedabc[is.na(xsharedabc)] <- 0
  xsharedabcd <- merge(xsharedabc, dshared, all=TRUE)
  xsharedabcd[is.na(xsharedabcd)] <- 0
  xsharedabcde <- merge(xsharedabcd, eshared, all=TRUE)
  xsharedabcde[is.na(xsharedabcde)] <- 0
  xsharedabcdef <- merge(xsharedabcde, fshared, all=TRUE)
  xsharedabcdef[is.na(xsharedabcdef)] <- 0
  xsharedabcdefg <- merge(xsharedabcdef, gshared, all=TRUE)
  xsharedabcdefg[is.na(xsharedabcdefg)] <- 0
  xsharedabcdefgh <- merge(xsharedabcdefg, hshared, all=TRUE)
  xsharedabcdefgh[is.na(xsharedabcdefgh)] <- 0
  xsharedabcdefghi <- merge(xsharedabcdefgh, ishared, all=TRUE)
  xsharedabcdefghi[is.na(xsharedabcdefghi)] <- 0
  xsharedabcdefghij <- merge(xsharedabcdefghi, jshared, all=TRUE)
  xsharedabcdefghij[is.na(xsharedabcdefghij)] <- 0
  xsharedabcdefghijk <- merge(xsharedabcdefghij, kshared, all=TRUE)
  xsharedabcdefghijk[is.na(xsharedabcdefghijk)] <- 0
  xsharedabcdefghijkl <- merge(xsharedabcdefghijk, lshared, all=TRUE)
  xsharedabcdefghijkl[is.na(xsharedabcdefghijkl)] <- 0

  AIMHeatmap <- xsharedabcdefghijkl
  AIMHeatmap2 <- AIMHeatmap[rowSums(AIMHeatmap[,-1])>0, ]
  
  row.names(AIMHeatmap2) <- AIMHeatmap2$CDR3.AminoAcid
  AIMHeatmap2 <- AIMHeatmap2[,2:13]
  AIMHeatmap2_matrix <- data.matrix(AIMHeatmap2)
  
  AIM101_heatmap <- heatmap.2(AIMHeatmap2_matrix, col=RdYlBu, scale="row", margins = c(8,5),labRow="", Colv=FALSE,
                              trace="none", key.title=NA, key.ylab=NA, key.xlab=NA, density.info="none", keysize=1)
}
aa.AIMheatmap101(aim101, TCRdata101)