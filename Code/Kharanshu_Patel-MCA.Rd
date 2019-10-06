---
title: "Kharanshu_Patel-MCA"
author: "Kharanshu Patel"
date: "11/6/2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# devtools::install_github('HerveAbdi/PTCA4CATA', force = TRUE)
# devtools::install_github('HerveAbdi/data4PCCAR')

library(ExPosition)
library(corrplot)
library(ggplot2)
library(InPosition)
library(dplyr)
library(PTCA4CATA)
library(data4PCCAR)
library(classInt)
library(Hmisc)
library(psych)
library(TInPosition)
library(RColorBrewer)
library(plotly)
library(heatmaply)
library(TExPosition)

```
# Method: MCA

Multiple Correspondence Analysis (MCA) is an extension of Correspondence Analysis(CA). It helps to analyze the pattern of relationship among several nominal variables. It is used to analyze a set of observations described by a set of nominal variables. Each nominal variable is comprised of several levels and each of these levels is coded as a binary variable.

MCA can also accommodate quantitative variables by recoding them as “bins.”

# Dataset : oasis_cross-sectional

The Open Access Series of Imaging Studies (OASIS) is a project aimed at making MRI data sets of the brain freely available to the scientific community.
It is a cross-sectional MRI Data in Young, Middle Aged, Nondemented and Demented Older Adults.
This set consists of a cross-sectional collection of 216 subjects (rows) aged 33 to 96. The subjects include both men and women. 
It measures 216 subjects (rows) on 4 quanitiative variables (columns). The 5 qualitative variables have been used to interpret the data. 

```{r data_set}

library(readxl)
oasis_cross_sectional <- read_excel("oasis_cross_sectional.xlsx")

head(oasis_cross_sectional)
oasis.dummy <- oasis_cross_sectional[c(6,8,9)]

# Recode Variable Educ
# Used dummy values 50 and 100 for the sake of recode

table(oasis_cross_sectional$Educ)
oasis_cross_sectional$Educ[oasis_cross_sectional$Educ == 1] <- "50" 
oasis_cross_sectional$Educ[oasis_cross_sectional$Educ == 5] <- "1" 
oasis_cross_sectional$Educ[oasis_cross_sectional$Educ == 50] <- "5" 

oasis_cross_sectional$Educ[oasis_cross_sectional$Educ == 2] <- "100" 
oasis_cross_sectional$Educ[oasis_cross_sectional$Educ == 4] <- "2" 
oasis_cross_sectional$Educ[oasis_cross_sectional$Educ == 100] <- "4" 

# Check Educ
table(oasis_cross_sectional$Educ)
```
### Recode Columns for MCA

```{r Recode Columns for MCA}

# Plot histogram

### Age
hist(oasis_cross_sectional[,3], breaks = 10, main="Age", xlab = "Values") 
qts.age <- quantile(as.numeric(unlist(oasis_cross_sectional[,3])))[2:4]

# cut Age
Age_recode <- cut(as.numeric(unlist(oasis_cross_sectional[,3])), breaks= c(min(oasis_cross_sectional[,3])-1,qts.age, max(oasis_cross_sectional[,3])+1),labels=c(1,2,3,4))

# check correlation
cor(as.numeric(Age_recode), as.numeric(unlist(oasis_cross_sectional[,3])), method = "spearman")

table(Age_recode)

### MMSE
hist(oasis.dummy[,1], breaks = 1, main="MMSE", xlab = "Values") 

# Cut MMSE
MMSE_recode <- cut(as.numeric(unlist(oasis.dummy[,1])), breaks= c(min(oasis.dummy[,1])-1,27, 29, max(oasis.dummy[,1])+1),labels=c(1,2,3))

table(MMSE_recode)                   
# check correlation
cor(as.numeric(MMSE_recode), as.numeric(unlist(oasis.dummy[,1])), method = "spearman")

### eTIV
min(oasis.dummy[,2])
max(oasis.dummy[,2])

hist(oasis.dummy[,2], breaks = 10, main = "eTIV", xlab = "values")
qts.e <- quantile(as.numeric(unlist(oasis.dummy[,2])))[2:4]

# cut eTIV
eTIV_recode <- cut(as.numeric(unlist(oasis.dummy[,2])), breaks= c(min(oasis.dummy[,2])-1,qts.e, max(oasis.dummy[,2])+1),labels=c(1,2,3,4))

# check correlation
cor(as.numeric(eTIV_recode), as.numeric(unlist(oasis.dummy[,2])), method = "spearman")

table(eTIV_recode)

### nWBV
min(oasis.dummy[,3])
max(oasis.dummy[,3])

hist(oasis.dummy[,3], breaks = 10, main = "nWBV", xlab = "values")

# cut nWBV
qts.n <- quantile(as.numeric(unlist(oasis.dummy[,3])))[2:4]

nWBV_recode <- cut(as.numeric(unlist(oasis.dummy[,3])), breaks= c(min(oasis.dummy[,3])-1,qts.n, max(oasis.dummy[,3])+1),labels=c(1,2,3,4))

# check correlation
cor(as.numeric(nWBV_recode), as.numeric(unlist(oasis.dummy[,3])), method = "spearman")

table(nWBV_recode)

```
## Bind Columns 

```{r Bind Columns}

Col.bind <- cbind(Age_recode, MMSE_recode, eTIV_recode, nWBV_recode)

colnames(Col.bind) <- c("Age", "MMSE", "eTIV", "nWBV")

head(Col.bind)
```
# Compute MCA
```{r MCA}
mca.res <- epMCA(Col.bind, make_data_nominal = TRUE, graphs = FALSE)

resMCA.sym <- epMCA(Col.bind, make_data_nominal = TRUE, graphs = FALSE, symmetric = TRUE)

resMCA.asym <- epMCA(Col.bind, make_data_nominal = TRUE, graphs = FALSE, symmetric = FALSE)

# MCA Factor Scores
mcaFj.a <- resMCA.asym$ExPosition.Data$fj
mcaFj <- resMCA.sym$ExPosition.Data$fj
mcaFi <- resMCA.sym$ExPosition.Data$fi

# MCA Constraints 
mca.constraints.sym <- minmaxHelper(mat1 = mcaFi, mat2  = mcaFj)
mca.constraints.asym <- minmaxHelper(mat1 = mcaFi, mat2  = mcaFj.a)

# Bootstrap for Confidence Interval
BootCube.oasis <- PTCA4CATA::Boot4Mean(mca.res$ExPosition.Data$fi, 
                                 design = oasis_cross_sectional$CDR,
                                 niter = 100,
                                 suppressProgressBar = TRUE)

```
# Design Colors
```{r Colors}
COLPAL <- list()

COLPAL$rows$CDR$vec <- oasis_cross_sectional$CDR
COLPAL$rows$MF$vec <- oasis_cross_sectional$`M/F`
COLPAL$rows$Educ$vec <- oasis_cross_sectional$Educ
COLPAL$rows$SES$vec <- oasis_cross_sectional$SES

COLPAL$rows$CDR$mat <- makeNominalData(as.matrix(COLPAL$rows$CDR$vec))
COLPAL$rows$MF$mat <- makeNominalData(as.matrix(COLPAL$rows$MF$vec))
COLPAL$rows$Educ$mat <- makeNominalData(as.matrix(COLPAL$rows$Educ$vec))
COLPAL$rows$SES$mat <- makeNominalData(as.matrix(COLPAL$rows$SES$vec))


# Reorder the columns in ascending order
COLPAL$rows$CDR$mat <- COLPAL$rows$CDR$mat[,order(colnames(COLPAL$rows$CDR$mat))]
COLPAL$rows$Educ$mat <- COLPAL$rows$Educ$mat[,order(colnames(COLPAL$rows$Educ$mat))]
COLPAL$rows$SES$mat <- COLPAL$rows$SES$mat[,order(colnames(COLPAL$rows$SES$mat))]
COLPAL$rows$MF$mat <- COLPAL$rows$MF$mat[,order(colnames(COLPAL$rows$MF$mat))]

COLPAL$rows$CDR$labels <- colnames(COLPAL$rows$CDR$mat)
COLPAL$rows$Educ$labels <- colnames(COLPAL$rows$Educ$mat)
COLPAL$rows$SES$labels <- colnames(COLPAL$rows$SES$mat)
COLPAL$rows$MF$labels <- unique(oasis_cross_sectional$`M/F`)

## CDR original --> 
COLPAL$rows$CDR$color_groups <- c("#037520","#daeb20","#f08300","#780000")
COLPAL$rows$Educ$color_groups <-c("#03022f","#03a6a8","#84026a","#f802ab","#d40108")
COLPAL$rows$SES$color_groups <- c("#013600","#8c8b00","#cbcc00","#ff6858","#ff0200")
COLPAL$rows$MF$color_groups <- c("#fc33e2", "#282be2")

#Next we need to create a vector of the group colors.
#We'll take one group at a time and trade the group names for the color names...

#First, copy the group names
COLPAL$rows$CDR$color_obs <- as.matrix(COLPAL$rows$CDR$vec)
COLPAL$rows$Educ$color_obs <- as.matrix(COLPAL$rows$Educ$vec)
COLPAL$rows$SES$color_obs <- as.matrix(COLPAL$rows$SES$vec)
COLPAL$rows$MF$color_obs <- as.matrix(COLPAL$rows$MF$vec)

#Then, for each group, replace the group name with the respective group color

COLPAL$rows$CDR$color_obs[which(COLPAL$rows$CDR$vec == "0")] <- 
COLPAL$rows$CDR$color_groups[1]
COLPAL$rows$CDR$color_obs[which(COLPAL$rows$CDR$vec == "0.5")] <- COLPAL$rows$CDR$color_groups[2]
COLPAL$rows$CDR$color_obs[which(COLPAL$rows$CDR$vec == "1")] <-
COLPAL$rows$CDR$color_groups[3]
COLPAL$rows$CDR$color_obs[which(COLPAL$rows$CDR$vec == "2")] <-
COLPAL$rows$CDR$color_groups[4]

COLPAL$rows$Educ$color_obs[which(COLPAL$rows$Educ$vec == "1")] <- 
COLPAL$rows$Educ$color_groups[1] #darkblue
COLPAL$rows$Educ$color_obs[which(COLPAL$rows$Educ$vec == "2")] <- 
COLPAL$rows$Educ$color_groups[2] #cyan
COLPAL$rows$Educ$color_obs[which(COLPAL$rows$Educ$vec == "3")] <- 
COLPAL$rows$Educ$color_groups[3] #magenta
COLPAL$rows$Educ$color_obs[which(COLPAL$rows$Educ$vec == "4")] <- 
COLPAL$rows$Educ$color_groups[4] #pink
COLPAL$rows$Educ$color_obs[which(COLPAL$rows$Educ$vec == "5")] <- 
COLPAL$rows$Educ$color_groups[5] #red

COLPAL$rows$SES$color_obs[which(COLPAL$rows$SES$vec == "1")] <- 
COLPAL$rows$SES$color_groups[1] #darkgreen
COLPAL$rows$SES$color_obs[which(COLPAL$rows$SES$vec == "2")] <- 
COLPAL$rows$SES$color_groups[2] #olivegreen
COLPAL$rows$SES$color_obs[which(COLPAL$rows$SES$vec == "3")] <- 
COLPAL$rows$SES$color_groups[3] #yellow
COLPAL$rows$SES$color_obs[which(COLPAL$rows$SES$vec == "4")] <- 
COLPAL$rows$SES$color_groups[4] #lightpink
COLPAL$rows$SES$color_obs[which(COLPAL$rows$SES$vec == "5")] <- 
COLPAL$rows$SES$color_groups[5] #red

COLPAL$rows$MF$color_obs[which(COLPAL$rows$MF$vec == "F")] <- 
COLPAL$rows$MF$color_groups[1]
COLPAL$rows$MF$color_obs[which(COLPAL$rows$MF$vec == "M")] <- 
COLPAL$rows$MF$color_groups[2]


col4jmca <- c("darkgreen","darkblue","darkred","orange")
col4mca <- prettyGraphsColorSelection(n.colors = nrow(mcaFj))
```
# Plot Results for MCA

### Heat Map(Not Correlation)

```{r Heat Map}

heatmaply(scale(Col.bind), k_row = 3, k_col = 2, dendrogram = FALSE, 
          main = "MCA Heat Map" )

```
### Scree Plot MCA + Inference
```{r Scree Plot MCA}

# Compute Inference
oasis.infMCA <- InPosition::epMCA.inference.battery(DATA = Col.bind,
                    DESIGN = oasis_cross_sectional$CDR,
                    graphs =  FALSE 
                    )


# Create Function for Scree Plot
PlotScreeFix <- function (ev, p.ev = NULL, max.ev = NULL, alpha = 0.05, col.ns = "#006D2C", 
          col.sig = "#54278F", title = "Explained Variance per Dimension", 
          plotKaiser = FALSE, color4Kaiser = "darkorchid4", lwd4Kaiser = 2.5) 
{
  val.tau = (100 * ev/sum(ev))
  Top.y = ceiling(max(val.tau) * 0.1) * 10
  if (!is.null(max.ev)) {
    ev = ev * (max.ev/ev[1])
  }
  par(mar = c(5, 6, 4, 4))
  plot(x = seq(1, length(val.tau)), y = val.tau, xlab = "Dimensions", 
       ylab = "Percentage of Explained Variance", main = title, 
       type = "l", col = col.ns, lwd = 1, xlim = c(1, length(val.tau)), 
       ylim = c(0, Top.y))
  points(x = seq(1, length(val.tau)), y = val.tau, pch = 16, 
         cex = 1, col = col.ns, lwd = 2.5)
  if (!is.null(p.ev)) {
    signi.vp = which(p.ev < alpha)
    ##### These are the lines I changed #####
    lines(x = seq(1, max(signi.vp)), y = val.tau[1:max(signi.vp)], 
          type = "l", col = col.sig, lwd = 1.5)
    points(x = signi.vp, y = val.tau[signi.vp], 
           pch = 16, cex = 1.5, col = col.sig, lwd = 3.5)
    #########################################
  }
  par(new = TRUE)
  par(mar = c(5, 6, 4, 4) + 0.5)
  le.max.vp = Top.y * (ev[1]/val.tau[1])
  plot(ev, ann = FALSE, axes = FALSE, type = "n", ylim = c(0, 
                                                           le.max.vp))
  if (plotKaiser) {
    abline(h = sum(ev)/length(ev), col = color4Kaiser, lwd = lwd4Kaiser)
  }
  mtext("Inertia Extracted by the Components", side = 4, line = 3)
  axis(4)
}

# Plot the Scree
PlotScreeFix(ev = oasis.infMCA$Inference.Data$components$eigs,
          p.ev = oasis.infMCA$Inference.Data$components$p.vals,
       title = 'Eigenvalues Inference',
       plotKaiser = TRUE
       )
       
```
# Symmetric MCA Factor Scores

### MCA Row Factors
```{r MCA Row Factors}

mcaMap.i <- createFactorMap(mcaFi, constraints = mca.constraints.sym,
                            col.points = COLPAL$rows$CDR$color_obs,
                            display.points = TRUE,
                            display.names = TRUE,
                            display.labels = TRUE,
                            alpha.points = 1,
                            pch = 20,
                            cex = 4,
                            col.labels = COLPAL$rows$CDR$color_obs,
                            text.cex = 2,
                            alpha.labels = 0.6,
                            col.axes = "black",
                            alpha.axes = 1,
                            col.background = "gray50",
                            width.axes = 1,
                            title = "Symmetric Row Factor Scores")

# Labels for Inertia
label4Map.oasis <- createxyLabels.gen(1,2,
lambda =mca.res$ExPosition.Data$eigs,
tau = mca.res$ExPosition.Data$t)

mca.oasis.ggMap <- mcaMap.i$zeMap + label4Map.oasis +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(mca.oasis.ggMap)

# Plot Group Means
CDR.Means <- PTCA4CATA::getMeans(mca.res$ExPosition.Data$fi, oasis_cross_sectional$CDR)

col4Means <- unique(COLPAL$rows$CDR$color_obs)

MapGroup <- PTCA4CATA::createFactorMap(CDR.Means,
                            constraints = mcaMap.i$constraints,
                            display.labels = FALSE,
                            col.points = col4Means,
                            alpha.points = 0.7,
                            cex = 6,  # size of the dot (bigger)
                            col.labels = col4Means,
                            text.cex = 4,
                            dev.new=FALSE,
                            new.plot = FALSE)

mca.oasis.withMeans <- mca.oasis.ggMap + MapGroup$zeMap_dots + MapGroup$zeMap_text + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(mca.oasis.withMeans)

# Plot Confidence Intervals using ellipses
GraphElli <- PTCA4CATA::MakeCIEllipses(BootCube.oasis$BootCube[,1:2,],
                            names.of.factors = c("Dimension 1","Dimension 2"),
                            col = col4Means,
                            p.level = .95
)

mca.oasis.withCI <-  mca.oasis.withMeans +  GraphElli
print(mca.oasis.withCI)
```

### MCA Column Factors

```{r MCA COlumn Factor Scores}

# Code Color
col4mca <- c("orange","orange","orange","orange","darkgreen","darkgreen","darkgreen",
             "darkblue","darkblue","darkblue","darkblue","darkred","darkred","darkred","darkred")


mcaMap.j <- createFactorMap(mcaFj, constraints = mca.constraints.sym,  
col.points = gplots::col2hex(col4mca),
col.labels = gplots::col2hex(col4mca) ,
display.points = TRUE, pos = 2,
display.names = TRUE,
display.labels = TRUE,
alpha.points = 0.8,
pch = 20,
cex = 4,
col.background = "#e9e8ff",
text.cex = 3,
font.face = "bold",
col.axes = "black",
alpha.axes = 0.3,
width.axes = 0.7,
title  = "Symmetric Column Factor Scores "
)

# Labels for Inertia
label4Map.oasis <- createxyLabels.gen(1,2,
lambda =mca.res$ExPosition.Data$eigs,
tau = mca.res$ExPosition.Data$t)

mca.oasis.ggMap2 <- mcaMap.j$zeMap + label4Map.oasis +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(mca.oasis.ggMap2)

Fj <- mca.res$ExPosition.Data$fj
colnames(Fj) <- paste0('Dimension ',1:ncol(Fj))
seg.1 <- as.data.frame(Fj[c(1,2,3,4),])
seg.2 <- as.data.frame(Fj[c(5,6,7),])
seg.3 <- as.data.frame(Fj[c(8,9,10,11),])
seg.4 <- as.data.frame(Fj[c(12,13,14,15),])

axis1 = 1
axis2 = 2

colnames(seg.1) <- paste0('Dimension ',1:ncol(Fj))
colnames(seg.2) <- paste0('Dimension ',1:ncol(Fj))
colnames(seg.3) <- paste0('Dimension ',1:ncol(Fj))
colnames(seg.4) <- paste0('Dimension ',1:ncol(Fj))


mca.line.map <- mca.oasis.ggMap2 +
geom_line(data = seg.1, col = col4jmca[c(4)], alpha =1) +
geom_line(data = seg.2, col = col4jmca[c(1)], alpha =1) + 
geom_line(data = seg.3, col = col4jmca[c(2)], alpha =1) +
geom_line(data = seg.4, col = col4jmca[c(3)], alpha =1)
print(mca.line.map)
```
### Symmetric IJ Map

```{r Symmetric IJ Map}

mca.symMap  <- createFactorMapIJ(mcaFi, mcaFj,
                             col.points.i = COLPAL$rows$CDR$color_obs,
                             col.labels.i = COLPAL$rows$CDR$color_obs,
                             col.points.j = col4mca,
                             col.labels.j = col4mca,
                             pch.i = 15, pch.j = 19,
                             cex.i = 1.2, cex.j = 3,
                             alpha.labels.i = 0.5,
                             alpha.labels.j = 1,
                             alpha.axes = 0.2,
                             alpha.points.i = 0.8,
                             alpha.points.j = 1,
                             col.axes = "#333434",
                             col.background = "gray 50",
                             text.cex.i = 1,
                             text.cex.j = 2.5,
                             title = "Symmetric I-J Map"
                             )

labels4MCA <- PTCA4CATA::createxyLabels.gen(1,2, lambda = mca.res$ExPosition.Data$eigs, tau = mca.res$ExPosition.Data$t)

map.IJ.sym <- mca.symMap$baseMap + mca.symMap$I_points +
  mca.symMap$J_labels + mca.symMap$J_points + labels4MCA + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(map.IJ.sym)
```

# Asymmetric MCA Factor Scores

### MCA Factor Scores : Asymmetric
```{r MCA Factor Scores : Assymetric}

mcaMap.j.a <- createFactorMap(mcaFj.a, constraints = mca.constraints.asym,  
color.points = col4mca ,
color.labels = col4mca ,
display.points = TRUE, pos = 2,
display.names = TRUE,
display.labels = TRUE,
alpha.points = 0.8,
pch = 20,
cex = 3.5,
text.cex = 2.5,
font.face = "bold",
col.axes = "red",
alpha.axes = 0.4,
col.background = "gray50",
width.axes = 0.8,
title = "Asymmetric Column Factor Scores"
)

# Labels for Inertia
label4Map.oasis <- createxyLabels.gen(1,2,
lambda =mca.res$ExPosition.Data$eigs,
tau = mca.res$ExPosition.Data$t)

mca.oasis.ggMap3 <- mcaMap.j.a$zeMap + label4Map.oasis +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(mca.oasis.ggMap3)
```
### Asymmetric IJ Map
```{r Asymmetric IJ Map}
mca.asymMap  <- createFactorMapIJ(mcaFi, mcaFj.a,
                                  display.labels = FALSE,
                                  display.names = TRUE,
                                  col.points.i = COLPAL$rows$CDR$color_obs,
                             col.labels.i = COLPAL$rows$CDR$color_obs,
                             col.points.j = col4mca,
                             col.labels.j = col4mca,
                             pch.i = 15, pch.j = 19,
                             cex.i = 1, cex.j = 3.5,
                             alpha.labels.i = 0.5,
                             alpha.labels.j = 0.8,
                             alpha.axes = 0.4,
                             col.axes = "#333434",
                             col.background = "#fff3e1",
                             text.cex.i = 1.5,
                             text.cex.j = 2,
                             title = "Asymmetric I-J Map"
                              )

map.IJ.asym <- mca.asymMap$baseMap + mca.asymMap$I_points +  
  mca.asymMap$J_points + mca.asymMap$J_labels + labels4MCA + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(map.IJ.asym)
```
# Contribution Bars

```{r MCA - Contribution Bars}

create.ctrPlot <- function(res, axis = 1,
col = NULL,
set = 'J', ...){
if (set == 'I'){
ctr.s <- res$ExPosition.Data$ci[,axis] *
sign(res$ExPosition.Data$fi[,axis])
} else {
ctr.s <- res$ExPosition.Data$cj[,axis] *
sign(res$ExPosition.Data$fj[,axis])
}
zebars <- PrettyBarPlot2(ctr.s,
threshold = 1 / length(ctr.s),
color4bar = gplots::col2hex(col), # we need hex code
main = paste0('Signed Contributions. Dimension ',axis),
ylab = 'Contributions',
ylim = c(1.2*min(ctr.s), 1.2*max(ctr.s) ), ...)
return(zebars)
} # end of function create.ctrPlot


### Dimension 1
m1.ctrJ <- create.ctrPlot(mca.res, axis = 1, font.size = 5,
col = col4mca, horizontal = FALSE)
print(m1.ctrJ)

### Dimension 2
m2.ctrJ <- create.ctrPlot(mca.res, axis = 2, font.size = 5,
col = col4mca, horizontal = FALSE)
print(m2.ctrJ)

### Dimension 3
m3.ctrJ <- create.ctrPlot(mca.res, axis = 3, font.size = 5,
col = col4mca, horizontal = FALSE)
print(m3.ctrJ)

### Dimension 4
m4.ctrJ <- create.ctrPlot(mca.res, axis = 4, font.size = 5,
col = col4mca, horizontal = FALSE)
print(m4.ctrJ)

### Dimension 5
m5.ctrJ <- create.ctrPlot(mca.res, axis = 5, font.size = 5,
col = col4mca, horizontal = FALSE)
print(m5.ctrJ)

```
