library("vegan")
library("iNEXT")
library("purrr")
library("BiodiversityR")

setwd("C:\Users\Sam\Dropbox\MSc Oceanography Work\SOES6039 - MSc Research Project\MSc_code\R.msc\.RData")

reps <-as.data.frame(read.table(filepath<-file.choose(),header=T,sep=',', row.names = 1))     ### 9reps
threereps <- as.data.frame(read.table(filepath<-file.choose(),header=T,sep=',', row.names = 1)) ## 3reps
alldata <-as.data.frame(read.table(filepath<-file.choose(),header=T,sep=',', row.names = 1)) ##alldata 


#WORMS taxonomy - database - 



#number of species

#cleaning alldata
alldata$UNKNOWN <- NULL
alldata$Unkn_21 <- NULL
alldata$XENO <- NULL
alldata$Plate.shaped <- NULL
alldata$Reticulated <- NULL
alldata$Tube.shaped <- NULL
alldata$UNDETERMINED_EPIFAUNA <- NULL
alldata$No_fauna <- NULL
alldata$`Partial.Specimen` <- NULL
alldata$ON_HardSub <- NULL
alldata$OFFbottom <- NULL
alldata$Whalebone <- NULL
alldata$TOOSMALL <- NULL
alldata$Interesting <- NULL
alldata$xeno_annotated <- NULL

alldataclean <- as.data.frame(alldata)

#turning NAs into 0s
alldataclean[is.na(alldataclean)] <- 0

write.csv(alldataclean, "alldataclean.csv" )

#subset alldata into its own DF
Area1<- subset(alldataclean, Area == "1")
Area1$Area_m2 <- cumsum(Area1$Area_m2)
Area2<- subset(alldataclean, Area == "2")
Area2$Area_m2 <- cumsum(Area2$Area_m2)
Area3<- subset(alldataclean, Area == "3")
Area3$Area_m2 <- cumsum(Area3$Area_m2)

#calc soecies accumulaton curve for area_m2
curve_all = specaccum(alldataclean[, 10:ncol(alldataclean)], method = "random", permutations = 100, xvar = area_m2)
curve_area1 = specaccum(Area1[, 10:ncol(Area1) ], method = "random", permutations = 100, xvar = area_m2)
curve_area2 = specaccum(Area2[, 10:ncol(Area2) ], method = "random", permutations = 100, xvar = area_m2)
curve_area3 = specaccum(Area3[, 10:ncol(Area3) ], method = "random", permutations = 100, xvar = area_m2)


#plot curves (with just areas)
plot(curve_area1, ci = 2, ci.type = c("polygon"), ci.lty = 0, xlab = "Area sampled (m2)", ylab = "Number of species", ylim=c(0,100), col = 2)
plot(curve_area2, add = TRUE, col = 3, ci = 2, ci.type = c("polygon"), ci.lty = 0)
plot(curve_area3, add = TRUE, col = 4, ci = 2, ci.type = c("polygon"), ci.lty = 0)
plot(curve_area1, add = TRUE, col = 1, ci = 2, ci.type = c("line"), ci.lty = 0)
plot(curve_area2, add = TRUE, col = 1, ci = 2, ci.type = c("line"), ci.lty = 0)
plot(curve_area3, add = TRUE, col = 1, ci = 2, ci.type = c("line"), ci.lty = 0)
legend('bottomright', c( 'Area1', 'Area2', 'Area3'), col=2:4, lty=1, bty='n', inset=0.025)

#SAC for all
plot(curve_all, ci = 1, ci.type = c("polygon"), ci.lty = 3, xlab = "Area sampled (m2)", ylab = "Number of species", ylim=c(0,120), col = 2)
plot(curve_all, add = TRUE, col = 1, ci = 1, ci.type = c("line"), ci.lty = 3)
legend('bottomright', c( 'All'), col=2, lty=1, bty='n', inset=0.025)

#estimating asymptote species diversity
estimatedspecies <-data.frame(matrix(nrow = 3, ncol=9))
rownames(estimatedspecies) <- c("Area1","Area2","Area3")
colnames(estimatedspecies) <- c( "Species","chao","chao.se","jack1","jack1.se","jack2","boot","boot.se","n")
repsforpool <- as.data.frame(threereps)
est1 <- repsforpool
est1 <- specpool(Area1[, 10:ncol(Area1) ])
est2 <- specpool(Area2[, 10:ncol(Area2) ])
est3 <- specpool(Area3[, 10:ncol(Area3) ])


# iNEXT SC and extrapolations
i1 <- as.data.frame(Area1[, 10:ncol(Area1) ])
i1<- t(i1)
i2 <- as.data.frame(Area2[, 10:ncol(Area2) ])
i2 <- t(i2)
i3 <- as.data.frame(Area3[, 10:ncol(Area3) ])
i3 <- t(i3)

library(data.table)

#iNEXT calculating all data together for SC etc
iall <- as.data.frame(alldataclean)
iall$filename <- NULL
iall$Dive <- NULL
iall$Lat <- NULL
iall$Long <- NULL
iall$Px_mm <- NULL
iall$Depth <- NULL
iall$Area_m2 <- NULL
iall$POC_flux <- NULL
#iall$Area <- NULL
#iall <- t(iall)
isplit <- split(iall[-1], f = iall$Area)
isplit <- lapply(isplit,function(x){
  t(x)
})

#isplit <- transpose(isplit)

iNEXT(i1, q=0, datatype="incidence_raw")
iNEXT(i2, q=0, datatype="incidence_raw")
iNEXT(i3, q=0, datatype="incidence_raw")
iNEXT(isplit, q=0, datatype="incidence_raw")
out1 <- iNEXT(i1, q= c(0,1,2), datatype="incidence_raw")
out2 <- iNEXT(i2, q= c(0,1,2), datatype="incidence_raw")
out3 <- iNEXT(i3, q= c(0,1,2), datatype="incidence_raw")
outall <- iNEXT(isplit, q=c(0,1,2), datatype="incidence_raw")
ggiNEXT(out1, type=1, se=TRUE, facet.var="none", color.var="order", grey=FALSE)
ggiNEXT(out2, type=1, se=TRUE, facet.var="none", color.var="order", grey=FALSE)
ggiNEXT(out3, type=1, se=TRUE, facet.var="none", color.var="order", grey=FALSE)
ggiNEXT(outall, type=1, se=TRUE, facet.var="site", color.var="order", grey = FALSE)
ggiNEXT(outall, type=2, facet.var="order", color.var="site")

#
RAD <- rankabundance(threereps[1,])
RAD2 <- rankabundance(threereps[2,])
RAD3 <- rankabundance(threereps[3,])
rankabunplot(RAD, scale = "abundance", col = 2, specnames = 0, type = "l")
rankabunplot(RAD2, addit = T , scale = "abundance", col = 3, specnames = 0, type = "l"  )
rankabunplot(RAD3, addit = T , scale = "abundance", col = 4, specnames = 0, type = "l"   )
legend('topright', c( 'Area1', 'Area2', 'Area3'), col=2:4, lty=1, bty='n', inset=0.025)

#all areas
RAD <-  radfit(threereps, model = "Mandelbrot")
plot(RAD, ylab = "Number of individuals", xlab = "Rank (Morphotype)", scale = "proportion")

#extract the iNEXT data
sapply(names(outall), 
       function (x) write.csv(outall[[x]], file=paste(x, "txt", sep=".") ))

#calc rarefaction curve for number of individuals
curve_alls = specaccum(alldataclean[, 10:ncol(alldataclean)], method = "rarefaction", xvar = individuals)
curve_area1s = specaccum(Area1[, 10:ncol(Area1) ], method = "rarefaction", xvar = individuals)
curve_area2s = specaccum(Area2[, 10:ncol(Area2) ], method = "rarefaction", xvar = individuals)
curve_area3s = specaccum(Area3[, 10:ncol(Area3) ], method = "rarefaction", xvar = individuals)

rareindiv1 <- curve_area1s$individuals
rareindiv2 <- curve_area2s$individuals
rareindiv3 <- curve_area3s$individuals
rareerror1 <- curve_area1s$sd
rareerror2 <- curve_area2s$sd
rareerror3 <- curve_area3s$sd
rarerich1 <- curve_area1s$richness
rarerich2 <- curve_area2s$richness
rarerich3 <- curve_area3s$richness

#rarefaction curve for species discovered vs number of individuals = species richness
plot(curve_area1s$individuals, curve_area1s$richness, xlab = "Number of Individuals", ylab = "Number of Species", ylim=c(0,100), xlim =c(0,800), col = 0)
with(curve_area1s, arrows(rareindiv1, rarerich1 - (2 * rareerror1), rareindiv1, rarerich1 + (2 * rareerror1),
                          angle = 90, code = 1, length = 0.05,col = 2))
with(curve_area1s, arrows(rareindiv3, rarerich3 - (2 * rareerror3), rareindiv3, rarerich3 + (2 * rareerror3),
                          angle = 90, code = 1, length = 0.05,col = 4))
with(curve_area1s, arrows(rareindiv2, rarerich2 - (2 * rareerror2), rareindiv2, rarerich2 + (2 * rareerror2),
                          angle = 90, code = 1, length = 0.05,col = 3))
lines(curve_area2s$individuals, curve_area2s$richness, col = 1, ylim=c(0,100), xlim =c(0,800))
lines(curve_area3s$individuals, curve_area3s$richness, col = 1, ylim=c(0,100), xlim =c(0,800))
lines(curve_area1s$individuals, curve_area1s$richness, col = 1, ylim=c(0,100), xlim =c(0,800))
legend('bottomright', c( 'Area1', 'Area2', 'Area3'), col=2:4, lty=1, bty='n', inset=0.025)




# number of individuals per transect
indiv1<- subset(alldataclean, Area == "1")
indiv2<- subset(alldataclean, Area == "2")
indiv3<- subset(alldataclean, Area == "3")

#total no of individuals
indiv1$total<- (rowSums(indiv1[, 10:ncol(indiv1)]))
indiv2$total<- (rowSums(indiv2[, 10:ncol(indiv2)]))
indiv3$total<- (rowSums(indiv3[, 10:ncol(indiv3)]))

#individuals per area
sum(indiv1$total)
sum(indiv2$total)
sum(indiv3$total)
#individuals per m2
sum(indiv1$total/(1656*3.6))
sum(indiv2$total/(1170*3.6))
sum(indiv3$total/(1228*3.6))

#number of individuals overall 
alldataclean$total<- (rowSums(alldataclean[, 10:ncol(alldataclean)]))
sum(alldataclean$total)
#number of individuals per m2 overall
sum(alldataclean$total/(4054*3.6))
alldataclean$total <- NULL


#9 reps species numbers
specnumber(reps)
#3 reps species numbers
specnumber(threereps)




#labels
Annelida <- c("Echiura_3", "Polychaete_1", "Polychaete_2", "Polychaete_3", "Polychaete_5")
Arthropoda <- c("Amphipod_4", "Aristeid_1", "Aristeid_2", "Aristeid_3", "Aristeid_4","Aristeid_5","ARTROPODA","Cirriped_2","Crustaceans", "Isopod_1", "Isopod_2","Mysidid_1", "Podocerid_1", "SpiderShr_1","Tanaidacea_1")
Bryozoa <- c("Bryozoa_1","Bryozoa_8")
Chordata <- c("FISH","Ipnops_1","Ophidiidae_4")
Cnidaria <- c("Actinia_10", "Actinia_15", "Actinia_16", "Actinia_22", "Actinia_36", "Actinia_9", "Actinia_new1", "Actinia_new2", "Actiniaria", "Alcyon_5", "Alcyon_8", "Alcyon_new1", "Alcyon_new2", "Alcyonacea", "ANTHOZOA", "Antipatharia", "Ceriantharia", "Hydroid_4", "Hydroid_epigrowth", "Pennatulacea", "Umbellula_1", "Zoantharia_1")       
Ctenophora <- c("Cteno_3","CTENOPHORA")
Echinodermata <- c("Amperima_1", "Amperima_2", "ASTEROIDEA", "Benthodytes_1", "Benthodytes_2", "Benthothur_2", "Brisingid_1", "CRINOIDEA","Deima","Deimatidae_2", "Echinoid_10","Enypniastes_1","HOL_002","HOL_006","HOL_044","Holo_new1","Holo_new2","HOLOTHUROIDEA", "Mesothuria_1", "Molpadiod_1", "Ophiuroid_1", "Ophiuroid_2", "OPHIUROIDEA", "Paelopatides_fl", "Paelopatides_wh", "Peniagone_1",  "Peniagone_lea", "Psychrop_3", "Psychrop_4", "Psychrop_6", "Synallact_1", "Velatid_3", "Velatid_new1", "Velatid_new2")
Mollusca <- c("Gastropod_1","Octopod_1","Scaphopod_1", "Teuthoid_1")
Porifera <- c("Desmo_1","Desmo_2","Desmo_4","Desmo_5","Desmo_6","Desmo_8","Desmo_9","Hexact_1","Hexact_10","Hexact_17","Hexact_3","Hexact_4","Hexact_5","Hexact_9", "PORIFERA", "Porifera_15","Porifera_23","Porifera_3","Porifera_5","Porifera_6","Porifera_new1")
Tunicata <- c("Ascidian_1","Ascidian_2","TUNICATA")
#Unknown


#phylaperarea -  turning column names into just phyla
Phylaperarea <- data.frame(matrix(nrow = 3, ncol = 10))
rownames(Phylaperarea) <- c("Area1","Area2","Area3")
colnames(Phylaperarea) <- c( "Annelida", "Arthropoda", "Bryozoa", "Chordata", "Cnidaria", "Ctenophora", "Echinodermata", "Mollusca", "Porifera", "Tunicata")
Phylaperarea$Arthropoda <- rowSums(threereps[c(Arthropoda)], na.rm = TRUE)
Phylaperarea$Annelida <- rowSums(threereps[c(Annelida)], na.rm = TRUE)
Phylaperarea$Echinodermata <- rowSums(threereps[c(Echinodermata)], na.rm = TRUE)
Phylaperarea$Bryozoa <- rowSums(threereps[c(Bryozoa)], na.rm = TRUE)
Phylaperarea$Chordata <- rowSums(threereps[c(Chordata)], na.rm = TRUE)
Phylaperarea$Cnidaria <- rowSums(threereps[c(Cnidaria)], na.rm = TRUE)
Phylaperarea$Ctenophora <- rowSums(threereps[c(Ctenophora)], na.rm = TRUE)
Phylaperarea$Mollusca <- rowSums(threereps[c(Mollusca)], na.rm = TRUE)
Phylaperarea$Porifera <- rowSums(threereps[c(Porifera)], na.rm = TRUE)
Phylaperarea$Tunicata <- rowSums(threereps[c(Tunicata)], na.rm = TRUE)
write.csv(Phylaperarea, "Phylaperarea.csv")

#allphyla
allphyla <- data.frame(matrix(nrow = 4054, ncol = 11))
rownames(allphyla) <- c(rownames(alldataclean))
colnames(allphyla) <- c( "Area", "Annelida", "Arthropoda", "Bryozoa", "Chordata", "Cnidaria", "Ctenophora", "Echinodermata", "Mollusca", "Porifera", "Tunicata")
allphyla$Arthropoda <- rowSums(alldataclean[c(Arthropoda)], na.rm = TRUE)
allphyla$Annelida <- rowSums(alldataclean[c(Annelida)], na.rm = TRUE)
allphyla$Echinodermata <- rowSums(alldataclean[c(Echinodermata)], na.rm = TRUE)
allphyla$Bryozoa <- rowSums(alldataclean[c(Bryozoa)], na.rm = TRUE)
allphyla$Chordata <- rowSums(alldataclean[c(Chordata)], na.rm = TRUE)
allphyla$Cnidaria <- rowSums(alldataclean[c(Cnidaria)], na.rm = TRUE)
allphyla$Ctenophora <- rowSums(alldataclean[c(Ctenophora)], na.rm = TRUE)
allphyla$Mollusca <- rowSums(alldataclean[c(Mollusca)], na.rm = TRUE)
allphyla$Porifera <- rowSums(alldataclean[c(Porifera)], na.rm = TRUE)
allphyla$Tunicata <- rowSums(alldataclean[c(Tunicata)], na.rm = TRUE)
allphyla$Area <- alldataclean$Area

allphylatotals <- data.frame(matrix(nrow = 3, ncol = 10))
rownames(Phylaperarea) <- c("Area1","Area2","Area3")
colnames(allphyla) <- c( "Annelida", "Arthropoda", "Bryozoa", "Chordata", "Cnidaria", "Ctenophora", "Echinodermata", "Mollusca", "Porifera", "Tunicata")
allphylatotals$Arthropoda <- colSums(allphyla[c(Arthropoda)] & allphyla[, c("Area")] = "1", na.rm = TRUE)
write.csv(allphylatotals, "allphylatotals.csv")


#phylaperrep
Phylaperrep <- data.frame(matrix(nrow = 9, ncol = 10))
rownames(Phylaperrep) <- c("A1R1","A1R2","A1R3","A2R1","A2R2","A2R3","A3R1","A3R2","A3R3")
colnames(Phylaperrep) <- c( "Annelida", "Arthropoda", "Bryozoa", "Chordata", "Cnidaria", "Ctenophora", "Echinodermata", "Mollusca", "Porifera", "Tunicata")
Phylaperrep$Arthropoda <- rowSums(reps[c(Arthropoda)], na.rm = TRUE)
Phylaperrep$Annelida <- rowSums(reps[c(Annelida)], na.rm = TRUE)
Phylaperrep$Echinodermata <- rowSums(reps[c(Echinodermata)], na.rm = TRUE)
Phylaperrep$Bryozoa <- rowSums(reps[c(Bryozoa)], na.rm = TRUE)
Phylaperrep$Chordata <- rowSums(reps[c(Chordata)], na.rm = TRUE)
Phylaperrep$Cnidaria <- rowSums(reps[c(Cnidaria)], na.rm = TRUE)
Phylaperrep$Ctenophora <- rowSums(reps[c(Ctenophora)], na.rm = TRUE)
Phylaperrep$Mollusca <- rowSums(reps[c(Mollusca)], na.rm = TRUE)
Phylaperrep$Porifera <- rowSums(reps[c(Porifera)], na.rm = TRUE)
Phylaperrep$Tunicata <- rowSums(reps[c(Tunicata)], na.rm = TRUE)
write.csv(Phylaperrep, "Phylaperrep.csv")



#analysing cnidarians
Cnidariananalysis <- data.frame(matrix(nrow = 3:length(Cnidaria), ncol =22))
rownames(Cnidariananalysis) <- c("Area1","Area2","Area3")
colnames(Cnidariananalysis) <- c(Cnidaria)
Cnidarians <- as.data.frame(Cnidariananalysis)
Cnidarians <- subset(threereps, select = c(Cnidaria))
Cnidariananalysis <- merge(threereps, Cnidariananalysis, by = colnames(Cnidaria))
Cnidariananalysis <- Cnidariananalysis[-(1:6),]
write.csv(Cnidarians, "Cnidarians.csv")
remove(Cnidariananalysis)
remove(Cnidarians)

#analysing porifera
Poriferanalysis <- data.frame(matrix(nrow = 3:length(Porifera), ncol =22))
rownames(Poriferanalysis) <- c("Area1","Area2","Area3")
colnames(Poriferanalysis) <- c(Porifera)
Poriferas <- as.data.frame(Poriferanalysis)
Poriferas <- subset(threereps, select = c(Porifera))
Poriferanalysis <- merge(threereps, Poriferanalysis, by = colnames(Porifera))
Poriferanalysis <- Poriferanalysis[-(1:6),]
write.csv(Poriferas, "Poriferas.csv")
remove(Poriferas)
remove(Poriferanalysis)

#analysing echinoderms
echinodermnalysis <- data.frame(matrix(nrow = 3:length(Echinodermata), ncol =34))
rownames(echinodermnalysis) <- c("Area1","Area2","Area3")
colnames(echinodermnalysis) <- c(Echinodermata)
echinoderms <- as.data.frame(echinodermnalysis)
echinoderms <- subset(threereps, select = c(Echinodermata))
echinodermnalysis <- merge(threereps, echinodermnalysis, by = colnames(Echinodermata))
echinodermnalysis <- echinodermnalysis[-(1:6),]
write.csv(echinoderms, "echinoderms.csv")
remove(echinoderms)
remove(echinodermnalysis)

#analysing arthropoda
arthropodanalysis <- data.frame(matrix(nrow = 3:length(Arthropoda), ncol =34))
rownames(arthropodanalysis) <- c("Area1","Area2","Area3")
colnames(arthropodanalysis) <- c(Arthropoda)
Arthropods <- as.data.frame(arthropodanalysis)
Arthropods <- subset(threereps, select = c(Arthropoda))
arthropodanalysis <- merge(threereps, arthropodanalysis, by = colnames(Arthropods))
arthropodanalysis <- arthropodanalysis[-(1:6),]
write.csv(Arthropods, "arthropods.csv")
remove(Arthropods)
remove(arthropodanalysis)

#analysis of functional types

SF <- c("Ascidian_1","Ascidian_2","TUNICATA", "Desmo_1","Desmo_2","Desmo_8","Desmo_9","Hexact_1","Hexact_10","Hexact_17","Hexact_3","Hexact_4","Hexact_5","Hexact_9", "PORIFERA", "Porifera_15","Porifera_23","Porifera_3","Porifera_5","Porifera_new1", "Actinia_10", "Actinia_15", "Actinia_16", "Actinia_22", "Actinia_36", "Actinia_9", "Actinia_new1", "Actinia_new2", "Actiniaria", "Alcyon_5", "Alcyon_8", "Alcyon_new1", "Alcyon_new2", "Alcyonacea", "ANTHOZOA", "Antipatharia", "Ceriantharia", "Hydroid_4", "Hydroid_epigrowth", "Pennatulacea", "Umbellula_1", "Zoantharia_1", "Bryozoa_1","Bryozoa_8", "Cirriped_2", "CRINOIDEA", "Echiura_3", "Polychaete_1", "Polychaete_2", "Polychaete_3", "Polychaete_5")
DF <- c( "ARTROPODA" ,"Tanaidacea_1", "SpiderShr_1", "Mysidid_1", "Amperima_1", "Amperima_2", "ASTEROIDEA", "Benthodytes_1", "Benthodytes_2", "Benthothur_2", "Brisingid_1","Deima","Deimatidae_2", "Echinoid_10","Enypniastes_1","HOL_002","HOL_006","HOL_044","Holo_new1","Holo_new2","HOLOTHUROIDEA", "Mesothuria_1", "Molpadiod_1", "Ophiuroid_1", "Ophiuroid_2", "OPHIUROIDEA", "Paelopatides_fl", "Paelopatides_wh", "Peniagone_1",  "Peniagone_lea", "Psychrop_3", "Psychrop_4", "Psychrop_6", "Synallact_1", "Velatid_3", "Velatid_new1", "Velatid_new2", "Scaphopod_1", "Gastropod_1", "Aristeid_1", "Aristeid_2", "Aristeid_3", "Aristeid_4","Aristeid_5", "Isopod_1", "Isopod_2", "Crustaceans")
PSC <- c( "FISH","Ipnops_1","Ophidiidae_4", "Cteno_3","CTENOPHORA", "Podocerid_1", "Desmo_4","Desmo_5","Desmo_6", "Porifera_6", "Octopod_1", "Teuthoid_1", "Amphipod_4") 


functionalgroups <- data.frame(matrix(nrow = 3:length(Porifera), ncol =3))
rownames(functionalgroups) <- c("Area1","Area2","Area3")
colnames(functionalgroups) <- c("SF","DF","PSC")
functionalgroups$SF <- rowSums(threereps[c(SF)], na.rm = TRUE)
functionalgroups$DF <- rowSums(threereps[c(DF)], na.rm = TRUE)
functionalgroups$PSC <- rowSums(threereps[c(PSC)], na.rm = TRUE)

write.csv(functionalgroups, "functionalgroups.csv")



#taxonomic richness by phyla
#taxonomic richness by functional group



#alternative for doing diversity tests
# Q = 1
test <- matrix (nrow = 9, ncol = 5)
colnames(test) <- c("Shannon", "Invsimpson", "Fisher", "Hulbert_rarefac", "pielous_evenness")
rownames(test) <- rownames(reps)
testresults <- data.frame(test)
testresults$Shannon <- diversity(reps, index = "shannon")

# Q = 2
testresults$Invsimpson <- diversity(reps, index = "invsimpson")

#commonly used fishers alpha value
testresults$Fisher <- fisher.alpha(reps)

#hulburt rarefaction
rarefy(threereps, min(rowSums(threereps)))

raref
#pielous eveness
H <- diversity(reps)
J <- H/log(specnumber(reps))
testresults$pielous_evenness <- J

testresults$Abund <- rowSums(reps) / (390*3.6)
write.csv(testresults, "testresults.csv")





# code for MDS plots
similarity <- metaMDS(reps^0.25, trace= FALSE, distance="bray", zerodist = "add")
areaMDS <- scores(similarity,display="sites")
plot(areaMDS,col=rep(2:4, each=3),xlab="",ylab="", xaxt='n',yaxt='n', cex.lab=1.5, cex.axis=1.5) 
legend("bottomright",c( 'Area1', 'Area2', 'Area3'), pch=1, col= 2:4)
similarity$stress
meta
#cluster tree (bray)
Bray.Curtis.Distances <- vegdist(reps, method = "bray")
clust.res<-hclust(dist.mat)
Bray.Curtis.Distance <- hclust(Bray.Curtis.Distances)
plot(clust.res)
plot(Bray.Curtis.Distance)
plot(Bray.Curtis.Distance, hang =-1, xlab = "Replicates", ylab="Bray Curtis Similarity", axes=FALSE)
axis(2, seq(1, 0.2, by=-0.1), seq(0.2, 1, by=0.1))





