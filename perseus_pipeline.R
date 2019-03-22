############################ Data Cleaning: ####################################
##read two datasets
pep <- as.data.frame(peptides)
pep<- peptides[,c(1,6)]
rawdata<- as.data.frame(shashi_allgroup)#from maxquant, PD, etc.
r<- melt(rawdata)
##extract accession numbers from database
library(splitstackshape)
peps <- as.data.frame(cSplit(pep, splitCols = c("Proteins"), 
                                 sep = ";", direction = "wide", fixed=T))
peptides = as.data.frame(peps[!(names(peps)%in% c("X1"))])
#colnames(db)<- c("ProteinID","Accession", "F1", "F2")
db$ <- paste(db$F1,"_",db$F2)
db = as.data.frame(db[!(names(db)%in% c("F1","F2"))])
i <- sapply(db, is.factor)
db[i] <- lapply(db[i], as.character)
#sapply(woimp, class)
#save(shashi_2, file= change.RData')# for database data
load(file="change.RData")

final.df<- data.frame(change.RData)
final.df[is.na(final.df)]<- NA ##converting miss to NA
##to remove NA
colSums(is.na(final.df)) #getting NAs
final.2 <- as.data.frame(final.df[rowSums(is.na(final.df)) <= 5, ]) #min 5 NA's
colSums(is.na(final.2))
change.RData[change.RData==0]<-NA #converting 0's to NA
str(change.RData)
summary(change.RData)

wimp<- as.data.frame(change.RData)#with imputation
db<- as.data.frame(db)
##change to simpler names
names(wimp2)[names(wimp2) == colnames(wimp2)[10]] <- "ProteinID"
#*********************************************************#
##source: https://www.r-bloggers.com/the-splitstackshape-package-for-r/
#install.packages('splitstackshape_1.4.2.zip', lib='C:\Users\Jitendra\Documents\R\win-library\3.3',repos = NULL)
#*********************************************************#
##seperate strings into columns
library(splitstackshape)
#help("splitstackshape")
mergeS <- as.data.frame(cSplit(wimp2, splitCols = c("ProteinID"), 
                               sep = ";", direction = "wide", fixed=T))
sapply(mergeS, class)

##convert variables Protein IDs_n from factor to character; applicable 
                                                    ##for any factor conversion
i <- sapply(mergeS, is.factor)
mergeS[i] <- lapply(mergeS[i], as.character)
mergeS[is.na(mergeS)] <- " " #convert all NA to blanks/0/NA
#sapply(merge, class)

##Merge using dplyr --------- 1(B)
library(dplyr)
CSW.f2<- rawdata %>% left_join(pep,by = c("Protein IDs"="Proteins")) %>% 
  left_join(db, by=c("ProteinID_2"="ProteinID")) %>%
  left_join(db, by=c("ProteinID_3"="ProteinID")) %>%
  left_join(db, by=c("ProteinID_4"="ProteinID")) %>%
  left_join(db, by=c("ProteinID_5"="ProteinID")) %>%
  left_join(db, by=c("ProteinID_6"="ProteinID")) %>%
  left_join(db, by=c("ProteinID_7"="ProteinID")) 

##RESULT: saves environment space and could be used to join multiple dataframes using %>%
CSW.f2[is.na(CSW.f2)] <- " "

##combine Protein ID and fasta headers
library(tidyr)
CSW.final <- unite_(CSW.f2, paste0(colnames(CSW.f2)[12:18], collapse = ""
), sep= ';', colnames(CSW.f2)[12:18])
#traceback()
names(CSW.final)[names(CSW.final) == colnames(CSW.final)[12]] <- "ProteinID" ##----- 1(B)

##for fasta headers
CSW.final <- unite_(CSW.final, 
                       paste0(colnames(CSW.final)[c(14,16,18,20,22,24,26)],
                       collapse = " "), sep = ";",
                       colnames(CSW.final)[c(14,16,18,20,22,24,26)])
names(CSW.final)[names(CSW.final) == colnames(CSW.final)[14]] <- "Fasta"

#for accession numbers
CSW.final <- unite_(CSW.final,
                       paste0(colnames(CSW.final)[c(13,15,16,17,18,19,20)],
                       collapse = " "), sep = ";",
                       colnames(CSW.final)[c(13,15,16,17,18,19,20)])
names(CSW.final)[names(CSW.final) == colnames(CSW.final)[13]] <- "Accession"

CSW.final <- as.data.frame(CSW.final)
CSW.final<- CSW.final[-c(10:11)]# drop variables maj PID,fasta headers
CSW.final
gsub("^;*|(?<=;);|;*$", "", CSW.final, perl=T)

#library(xlsx)
#library(rJava)
#write.xlsx(CSW.final, file = "merge.xlsx", sheetName= NULL, 
#row.names =T, col.names = T)
write.csv(CSW.f2, file = "CSW_final.csv", col.names = T)
##**************************************************************************##

###################### Data Preprocessing: #####################################
##read new data 
dat.a<- as.data.frame(dat_final)

fasta<- dat.a[,12]
group1<- (dat.a[ ,1:3])
group2<- (dat.a[ ,4:6])
group3<- (dat.a[ ,7:9])

##log transform ---------1(A)
dat.log<- log(dat.a[,3:8],2)
##column and row Correlation & fast pearson column correlation ---------- 1(C),(D),(E) 
#a) using MKmisc - column corrrleation 
library(MKmisc)
dat.cor<- cor((dat.a[,2:10]), method = "pearson", use = "pairwise.complete.obs")
cols<- colnames(dat.a[ ,2:10])
par(mar=c(5,7,4,2)+0.1)
corPlot(dat.cor,minCor = 0.95, cex.axis = 0.7, labels = cols, lab.both.axes = T,
        las =2, title = "pearson correlation matrix")
CD <- corDist(t(dat.a[,2:10]), method = "pearson")
plot(hclust(CD))

#b) Using corrplot package - column correlation
library(corrplot)
dat.Corr<- cor((dat.a[ ,2:10]), method = "pearson", use = "pairwise.complete.obs")
corrplot(dat.Corr, method="shade", order="hclust",addrect= 5, type= c("full"), is.corr=F, 
         sig.level = 0.05)

#i <- sapply(dat.a, is.factor) ##Protein ID is factor?
#dat.a[i] <- lapply(dat.a[i], as.character)

#for fast pearson column correlation
library(WGCNA)
cor(x, y = NULL, 
    use = "all.obs", 
    method = c("pearson", "kendall", "spearman"),
    quick = 0, 
    cosine = FALSE, 
    cosineX = cosine,
    cosineY = cosine, 
    nThreads = 0, 
    verbose = 0, indent = 0)
## for row correlation transpose the matrix and use cor() function 
corFast(x, y = NULL, 
        use = "all.obs", 
        quick = 0, nThreads = 0, 
        verbose = 0, indent = 0)

cor1(x, use = "all.obs", verbose = 0, indent = 0)

##Density Kernel Estimation ----------- 1(I)
density(x, ...)
## Default S3 method:
density(dat.a, bw = "nrd0", adjust = 1,
        kernel = c("gaussian", "epanechnikov", "rectangular",
                   "triangular", "biweight",
                   "cosine", "optcosine"),
        weights = NULL, window = kernel, width,
        give.Rkern = FALSE,
        n = 512, from, to, cut = 3, na.rm = FALSE, ...)

## Perfromance Curves using package - ROCR --------------1(J)
## cloning using package - dclone -------------1(L)


#Visualizations
#Boxplot 
boxplot(dat.a[ ,1:9],col = ("yellow"), las= 2, beside=T)
#extracting outliers
#N2<- boxplot.stats(dat.a[, 1])$out
#N3<- boxplot.stats(dat.a[, 2])$out
#N4<- boxplot.stats(dat.a[, 3])$out
#SR1<- boxplot.stats(dat.a[, 4])$out
#SR3<- boxplot.stats(dat.a[, 5])$out
#SR4<- boxplot.stats(dat.a[, 6])$out
#SS1<- boxplot.stats(dat.a[, 7])$out
#SS2<- boxplot.stats(dat.a[, 8])$out
#SS3<- boxplot.stats(dat.a[, 9])$out
#outlier <- table(c(N2,N3,N4,SR1,SR3, SR4, SS1,SS2, SS3))

##Density Plots
library(ggplot2)
density.df<- stack(dat.a[,1:9])
is.factor(density.df[ ,2])
ggplot(density.df, aes(x=values)) + geom_density()
ggplot(density.df, aes(x=values)) + geom_density(aes(group=ind))
ggplot(density.df, aes(x=values)) + geom_density(aes(group=ind, colour=ind))
##for filling densities
ggplot(density.df, aes(x=values)) + geom_density(aes(group=ind, colour=ind, fill=ind),alpha=0.3)

##missing data 
library(VIM)
dat.miss <- aggr(dat.a, col =c("Yellow","Blue"), numbers=T,
sortVars=T, cex.axis=.5,gap=1,
    ylab=c("Missingdata","Pattern"), labels=names(dat.a))

##or
library(Amelia)
#par(mar = c(4, 2, 2, 2)+0.1)
rows1<- val.mis[,1]
missmap(val.mis[,3:8], legend=T, col = c("yellow", "Blue"), x.cex = 0.5, 
                                        y.cex = 0.3, y.lables = rows1 )

##A)imputation
#a) using MICE
library(mice)# for more robustness
md.pattern(CleanData)
md.pairs(CleanData)# freq of missing vals b/w pairs 
#m=1 number of multiple imputations
#maxit=2 number of iterations. 10-20 is sufficient.
CleanData.imp <- mice(CleanData, m=5, maxit=2, printFlag=TRUE) 
summary(CleanData.imp)
#************************************************************#
#Others
# b) using AMELIA
#library(Amelia)
#idvars <- c("Gene names", "Protein names")
#CleanData.am<- amelia(x= CleanData, m=2, idvars = idvars)

#c) Using missForest
#library(missForest)
#CleanData.miss <- missForest(CleanData, verbose=T, maxiter=2)
#CleanData.miss$CleanData
#CleanData.miss$OOBerror

#d) Using Hmisc*not functional
#library(Hmisc)
#CleanData.h<- with(CleanData, 'random')

#e) using mi
#library(mi)
#CleanData.mi <- mi(CleanData, n.iter = 2)
#***********************************************#

#B) Filtering
#1) Extract imputed dataframes
CleanData2<- complete(CleanData.imp, 2)
boxplot(CleanData2[, 3:18],col = ("yellow"), las= 2, beside=T)
foo<- boxplot(CleanData2[, 3:18],col = ("yellow"), las= 2, beside=T)
boxplot.stats(CleanData2[, 3])$out
CleanData2<- as.data.frame(CleanData2)

################## Visualizations2: ######################################

##Scatter plot
#no package
par(mfrow=c(2,3))
plot(dat.a$N2, dat.a$N3, pch=2,col= "red")
plot(dat.a$N2, dat.a$N4, pch=3,col= "red")
plot(dat.a$N3, dat.a$N4, pch=4,col= "red")
plot(dat.a$SR1, dat.a$SR3, pch=5,col= "blue")
plot(dat.a$SR1, dat.a$SR4, pch=6,col= "blue")
plot(dat.a$SR3, dat.a$SR4, pch=7,col= "blue")
plot(dat.a$SS1, dat.a$SS2, pch=8,col= "green")
plot(dat.a$SS1, dat.a$SS3, pch=9,col= "green")
plot(dat.a$SS2, dat.a$SS3, pch=1,col= "green")
dev.off()

#or
#using car
library(car)
dat.a[1:9]
scatterplotMatrix(dat.a[1:9])


library(RColorBrewer)
samples <- c("N2", "N3","N4", "SR1", "SR3", "SR4", "SS1", "SS2", "SS3")
sample1<- c("N2", "N3","N4")
sample2<- c("SR1", "SR3", "SR4")
sample3<- c("SS1", "SS2", "SS3")
dat.list <- list(dat.a$N2,dat.a$N3,dat.a$N4,dat.a$SR1,dat.a$SR3,dat.a$SR4,
                 dat.a$SS1, dat.a$SS2, dat.a$SS3)
dat.list1<- list(dat.a$N2,dat.a$N3,dat.a$N4)
dat.list2<- list(dat.a$SR1,dat.a$SR3,dat.a$SR4)
dat.list3<- list(dat.a$SS1, dat.a$SS2, dat.a$SS3)
makeProfilePlot(dat.list1,sample1)
makeProfilePlot(dat.list2, sample2)
makeProfilePlot(dat.list3, sample3)

################## Clustering ###################################
#Using pvclust
#complete linkage
library(pvclust)
dat.pv<-pvclust(dat.a[ ,1:9], nboot=100, method.hclust= "complete", 
                method.dist="correlation",use.cor="pairwise.complete.obs")
rnames <- dat.a[,11]
plot(dat.pv, print.pv=T, print.num=T)
table(cutree(dat.pv$hclust, k=3))
pvrect(dat.pv, alpha = 0.95)
#average linkage
dat.pv2<-pvclust(dat.a[ ,1:9], nboot=100, method.hclust= "average", 
                 method.dist="correlation",use.cor="pairwise.complete.obs")
rnames <- dat.a[,11]
plot(dat.pv2, print.pv=T, print.num=T)
#single linkage
dat.pv3<-pvclust(dat.a[ ,1:9], nboot=100, method.hclust= "single", 
                 method.dist="correlation",use.cor="pairwise.complete.obs")
rnames <- dat.a[,11]
plot(dat.pv3, print.pv=T, print.num=T)

##kmeans clustering; needed rows w/o NA

#dat.num <- as.data.frame(dat.a[rowSums(is.na(dat.a)) <= 6,])
#cannot perform bacause every row has NA
#impute using mice
#library(mice)
#md.pattern(dat.a)
#md.pairs(dat.a)# freq of missing vals b/w pairs 
#dat.imp <- mice(dat.a, m=5, maxit=2, printFlag=TRUE) 
#summary(dat.imp)
#dat.imp1<- complete(dat.imp, 1)
#dat.imp1<- as.data.frame(dat.imp1) 

set.seed(123)
dat.kmeans<- kmeans(dat.a[ ,1:9],5, 
                    nstart=25, iter.max = 500) #nstart good b/w 25-1000
attributes(dat.kmeans)
dat.kmeans$centres
dat.kmeans$totss
dat.kmeans$withinss
dat.kmeans$betweenss
dat.kmeans$iter
dat.kmeans$ifault
dat.kmeans$cluster # <- important
dat.kmeans$cluster
dat.kmeans$size
#table(dat.kmeans$cluster, dat.a$N2) <- not really needed
# Apply k-means with k=3
dat.km2 <- kmeans(dat.a[1:9], 3, nstart=25, iter.max=1000)
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(dat.a[,1:9], col=dat.km2$clust, pch=16)

# Cluster sizes
sort(table(dat.km2$clust))
dat.clust <- sort(table(dat.km2$clust))

#***********************************************************************#
#using package cluster
library(cluster)
clusplot(dat.a[ ,1:9], dat.kmeans$cluster, main='2D Cluster solution',
         color=TRUE, shade=TRUE, labels=2, lines=0, stand=F)

#using package factoextra
library(factoextra)
fviz_cluster(dat.kmeans, dat.a[ ,1:9], main = "k-means cluster plot; k=5", 
             outlier.color = "black")

################################################################################
#)PCA plot
#normal plot: without package
dat.pc <- prcomp(dat.a[ ,1:9], scale. = F, center = T)
dat.pcv<- princomp(dat.a[ ,1:9], cor = T, score = T)
summary(dat.pc)
print(dat.pc)
dat.pc$sdev
dat.pc$rotation
pccomp <- data.frame(colnames(dat.a[ ,1:9]), dat.pc$rotation)
names(pccomp) [names(pccomp)==colnames(pccomp)[1]]<- "Category"
pccomp <- 
  
  plot(dat.pc, type = "barplot", las= 2, scale.=T)
plot(dat.pc, type = "lines", las= 2, label=T)
dat.df<- data.frame(dat.pc$x[ ,1:9]) # dataframe made of first 2 PCs
range(dat.df[,"PC1"]) 
range(dat.df[,"PC2"]) 
range(dat.df[,"PC3"]) 
range(dat.df[,"PC4"]) 

#using factoextra
library(factoextra)
dat.group2 <- as.factor(colnames(dat.df[,1:9]))
fviz_pca_ind(dat.pc, addEllipses = TRUE, ellipse.level = 0.68) +
  theme_minimal()

#using ggbiplot
library(ggbiplot)
strains<- colnames(dat.a[ ,1:9])

ggbiplot(dat.pc, choices = 1:2, 
         scale = 1, obs.scale = 1 - scale, groups = strains,
         var.scale = 1, var.axes = T, ellipse = TRUE, circle = TRUE)
+ scale_color_discrete(name = '')
+ theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)
print(ggscreeplot(dat.pc))

#using ggfortify package
library(ggfortify)
category<- colnames(dat.a[ ,1:9])
autoplot(prcomp(dat.a[ ,1:9]),label= T, loadings=T, colour = category, 
         loadings.label=T)

# using pca2d package
library(pca3d)
category<- colnames(dat.a[ ,1:9])
pca3d(dat.pc, fancy = T, bg = "black", 
      axes.color = "white", show.ellipses = T, new =T)
pca2d(dat.pc, group = dat.a[ ,12])