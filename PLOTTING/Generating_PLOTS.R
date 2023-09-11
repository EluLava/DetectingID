library("vioplot")

load("../1KG/data/ID_EST_Datasets/WORLD_ID3_simus_COMPLETE_ID_EST.RData")
load("../1KG/data/ID_EST_Datasets/EAS_ID3_simus_COMPLETE_ID_EST.RData")
load("../1KG/data/ID_EST_Datasets/AFR_ID3_simus_COMPLETE_ID_EST.RData")
load("../PEDIGREE/data/ID_EST_Datasets/PGPED_ID3_simus_COMPLETE_ID_EST.RData")

#Load the data Fdatasets
dtaEAS = read.table("../1KG/data/F_datasets/EAS_allchr_allFs.txt", header = T)
dtaAFR = read.table("../1KG/data/F_datasets/AFR_allchr_allFs.txt", header = T)
dtaWORLD = read.table("../1KG/data/F_datasets/WORLD_allchr_allFs.txt", header = T)
dtaPOLYPED = read.table("../PEDIGREE/data/F_datasets/POLYGAMOUS_PED_allchr_allFs.txt", header = T)

#Load polypeds subsampling AFR RANGED
load("../PEDIGREE/data/ID_EST_Datasets/PGPED_AFR_RANGEDSUB.RData")
WORLD_SUB_RANGED_POLYPED_ID3 = foreachoutput
#subsample the scenario
WORLD_SUB_RANGED_POLYPED_ID3 = WORLD_SUB_RANGED_POLYPED_ID3[WORLD_SUB_RANGED_POLYPED_ID3$Simuscen == "LM_ID3_selhT_selpT_demaT",]
#Load polypeds subsampling AFR RANDOM
load("../PEDIGREE/data/ID_EST_Datasets/PGPED_WORLD_RANDOMSUB.RData")
WORLD_SUB_RANDOM_POLYPED_ID3 = foreachoutput
#subsample the scenario
WORLD_SUB_RANDOM_POLYPED_ID3 = WORLD_SUB_RANDOM_POLYPED_ID3[WORLD_SUB_RANDOM_POLYPED_ID3$Simuscen == "LM_ID3_selhT_selpT_demaT",]

load("../PEDIGREE/data/ID_EST_Datasets/PGPED_EAS_RANGEDSUB.RData")
EAS_SUB_RANGED_POLYPED_ID3 = foreachoutput
#subsample the scenario
EAS_SUB_RANGED_POLYPED_ID3 =EAS_SUB_RANGED_POLYPED_ID3[EAS_SUB_RANGED_POLYPED_ID3$Simuscen == "LM_ID3_selhT_selpT_demaT",]
#Load polypeds subsampling AFR RANDOM
load("../PEDIGREE/data/ID_EST_Datasets/PGPED_EAS_RANDOMSUB.RData")
EAS_SUB_RANDOM_POLYPED_ID3 = foreachoutput
#subsample the scenario
EAS_SUB_RANDOM_POLYPED_ID3 = EAS_SUB_RANDOM_POLYPED_ID3[EAS_SUB_RANDOM_POLYPED_ID3$Simuscen == "LM_ID3_selhT_selpT_demaT",]

load("../PEDIGREE/data/ID_EST_Datasets/PGPED_AFR_RANGEDSUB.RData")
AFR_SUB_RANGED_POLYPED_ID3 = foreachoutput
#subsample the scenario
AFR_SUB_RANGED_POLYPED_ID3 = AFR_SUB_RANGED_POLYPED_ID3[AFR_SUB_RANGED_POLYPED_ID3$Simuscen == "LM_ID3_selhT_selpT_demaT",]
#Load polypeds subsampling AFR RANDOM
load("../PEDIGREE/data/ID_EST_Datasets/PGPED_AFR_RANDOMSUB.RData")
AFR_SUB_RANDOM_POLYPED_ID3 = foreachoutput
#subsample the scenario
AFR_SUB_RANDOM_POLYPED_ID3 = AFR_SUB_RANDOM_POLYPED_ID3[AFR_SUB_RANDOM_POLYPED_ID3$Simuscen == "LM_ID3_selhT_selpT_demaT",]


#rm
rm(foreachoutput)

#########################################################
########### FIGURE 1 POLYPED LM VS LMM GRM AS ###########

#PLOT
pdf("./FIG1_LM_LMMAS_IDest_POLYPED.pdf", height = 50, width = 32)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5),0,5,0,6,0,rep(0,5)),
                    ncol = 5, nrow = 7, byrow = T), heights = c(0.2,1.0,0.1,1.0,0.1,1.5,0), widths = c(0.5,1,.3,1,0))

par(mar = c(1,2,5,1))

#### PANNEL A: POLYPED FULL SAMPLE SIZE LM ####

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  5, cex = 10)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 11,924)\ncomplete sample size", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL B: POLYPED FULL SAMPLE SIZE LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 5, cex = 10)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL C: POLYPED REDUCED SAMPLE SIZE LM ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 2,500)\nrandom sub-sampling", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL D: POLYPED REDUCED SAMPLE SIZE LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL C: POLYPED REDUCED SAMPLE SIZE LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 2,500)\nranged sub-sampling", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL D: POLYPED REDUCED SAMPLE SIZE LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

#####################################################
########### FIGURE 2 1KG LM VS LMM GRM AS ###########

#PLOT
pdf("./PLOTvIII/FIG2_LM_LMMAS_IDest_1KG.pdf", height = 50, width = 32)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5),0,5,0,6,0,rep(0,5)),
                    ncol = 5, nrow = 7, byrow = T), heights = c(0.2,1.0,0.1,1.0,0.1,1.5,0), widths = c(0.5,1,.3,1,0))

par(mar = c(1,2,5,1))

#### PANNEL A: EAS LM ####

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:11])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[5:11]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[5:11])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  5, cex = 10)
#ADD MTEXT FOR POP
mtext(text = "1KG: EAS (n = 504)", side = 2, outer = F, line =  35, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL B: EAS LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:21])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[15:21]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[15:21])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR MODEL
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  5, cex = 10)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL C: AFR LM ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:11])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[5:11]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[5:11])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR POP
mtext(text = "1KG: AFR (n = 661)", side = 2, outer = F, line =  35, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL D: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:21])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[15:21]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[15:21])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL E: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:11])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[5:11]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[5:11])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR POP
mtext(text = "1KG: WORLD (n = 2,504)", side = 2, outer = F, line =  35, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL F: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:21])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[15:21]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[15:21])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

####################################################
#### FIGURE 3 PLOTTING F.UNI.WE for comparisons ####

pdf("./PLOTvIII/FIG3_FuniWE_ALL.pdf", height = 30, width = 25)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)),
                    ncol = 5, nrow = 5, byrow = T), heights = c(0.1,1.0,0.2,1.5,0), widths = c(0.3,1,.3,1,0))

#margins
par(mar = c(1,2,5,1))

#### PANNEL A: PEDIGREE ####

#Set colors
colours = c(rep("#3366FF",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,9])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,9]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,9])),
              log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,18])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,18]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,18])),
              log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,27])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,27]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[27])),
              log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[36])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,36]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,36]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[UNI]^w ~ " PEDIGREE")), line = 0, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -.75, cex = 7, line = -3)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL B: EAS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,7])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,7]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,7])),
              log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,17])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,17]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,17])),
              log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,27])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,27]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,27])),
              log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,37])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,37]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,37]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[UNI]^w ~ " EAS")), line = 0, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -.75, cex = 7, line = -3)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL C: AFR ####

#change margins for x axis labels
par(mar = c(20,1,2,1))

#Set models names
labs = c("LM", expression(LMM[AS]), expression(LMM[GCTA]^{w}), expression(LMM[GCTA]^{u}))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,7])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,7]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,7])),
              log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,17])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,17]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,17])),
              log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,27])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,27]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,27])),
              log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,37])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,37]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,37]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[UNI]^w ~ " AFR")), line = 0, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -.75, cex = 7, line = -3)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL D: WORLD ####

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,7])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,7]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,7])),
              log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,17])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,17]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,17])),
              log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,27])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,27]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,27])),
              log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,37])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,37]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,37]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500),log10(2500)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(8000), labels = labs , srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[UNI]^w ~ " WORLD")), line = 0, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -.75, cex = 7, line = -3)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

###################################################
#### FIG S2 POLYPED Comparison of different Fs ####

#Subset only the F we are interested in
dta = dtaPOLYPED[,c(4,5,7:12)]

#Set F names
Fnames = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
           expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
           expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Open pdf
pdf("./PLOTvIII/FIGS2_F_COMP_POLYPED.pdf", height = 40, width = 41)

#Create layout
layout(mat = matrix(c(rep(0,17),0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,rep(0,17),0,0,0,9,0,10,0,11,0,12,0,13,0,14,0,15,0,rep(0,17),0,0,0,0,0,16,0,17,0,18,0,19,0,20,0,21,0,
                      rep(0,17),0,0,0,0,0,0,0,22,0,23,0,24,0,25,0,26,0,rep(0,17),0,0,0,0,0,0,0,0,0,27,0,28,0,29,0,30,0,rep(0,17),0,0,0,0,0,0,0,0,0,0,0,31,0,32,0,33,0,
                      rep(0,17),0,0,0,0,0,0,0,0,0,0,0,0,0,34,0,35,0,rep(0,17),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,0,rep(0,17)),
                    nrow = 17, ncol = 17, byrow = F), heights = c(.2,1,rep(c(0,1),7),.2), widths = c(.3,rep(c(1,0),7),1,.2))

#Loop through F
for (fmaj in 1:ncol(dta)) {
  #First plot is always to display distributions of individual F
  hist(dta[,fmaj], ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.1, 0.5), ylim = c(0,6500), col = "deepskyblue3", main = NULL, breaks = seq(-0.1, 0.5, by = 0.02))
  box(col = "black", lwd = 1)
  
  #If last F, add axis labels, else only add thick
  if(fmaj < ncol(dta)){
    #Only add the thicks in the axis
    axis(1, xlab = NULL, at=seq(-0.1, 0.5, by = 0.05), labels = NA, hadj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.5, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(0, 0.4, by = 0.2), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
  }
  
  
  #Add F name
  mtext(Fnames[fmaj], side = 3, outer = FALSE, line = 1.5, cex = 5)
  
  
  #Loop through the other F IF fmaj < 9
  if(fmaj < ncol(dta)) {
    
    for(fmin in (fmaj + 1):ncol(dta)){
      
      #Plot first F again the others
      plot(dta[,fmaj] ~ dta[,fmin], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.1, 0.5), ylim = c(-0.1, 0.5), cex = 4)
      
      #Add the equality line
      abline(0,1, lwd = 3)
      
      #Add axis 1 if fmin = ncol(dta)
      if(fmin == ncol(dta)){
        #Add the thicks
        axis(1, xlab = NULL, at=seq(-0.1, 0.5, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(1, xlab = NULL, at=seq(0, 0.4, by = 0.2), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
      } else {
        #only add the ticks
        axis(1, xlab = NULL, at=seq(-0.1, 0.5, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
        
      }
      
      #Add axis 2 if fmaj = 1
      if(fmaj == 1){
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.1, 0.5, by = 0.05), labels = NA, hadj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(2, xlab = NULL, at=seq(-0.1, 0.5, by = 0.2), hadj = 1.5, cex.axis= 5, lwd.ticks = 4.5, las = 1, tck = -0.03)
      } else {
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.1, 0.5, by = 0.05), labels = NA, hadj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
      }
    }
  }   
}

dev.off()


###############################################
#### FIG S3 EAS Comparison of different Fs ####

#Subset only the F we are interested in
dta = dtaEAS[,c(4,6,8,10,12,14,16)]

#Set F names
Fnames = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
           expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
           expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))


#Open pdf
pdf("./PLOTvIII/FIGS3_F_COMP_EAS.pdf", height = 35, width = 35)

#Create layout
layout(mat = matrix(c(rep(0,15),0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,rep(0,15),0,0,0,8,0,9,0,10,0,11,0,12,0,13,0,rep(0,15),0,0,0,0,0,14,0,15,0,16,0,17,0,18,0,
                      rep(0,15),0,0,0,0,0,0,0,19,0,20,0,21,0,22,0,rep(0,15),0,0,0,0,0,0,0,0,0,23,0,24,0,25,0,rep(0,15),0,0,0,0,0,0,0,0,0,0,0,26,0,27,0,
                      rep(0,15),0,0,0,0,0,0,0,0,0,0,0,0,0,28,0,rep(0,15)), nrow = 15, ncol = 15, byrow = F), heights = c(.2,1,rep(c(0,1),6),.2), widths = c(.2,rep(c(1,0),6),1,.2))

#Loop through F
for (fmaj in 1:ncol(dta)) {
  #First plot is always to display distributions of individual F
  hist(dta[,fmaj], ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.05, 0.35), ylim = c(0,500), col = "deepskyblue3", main = NULL, breaks = seq(-0.1, 0.4, by = 0.01))
  box(col = "black", lwd = 1)
  
  #If last F, add axis labels, else only add thick
  if(fmaj < ncol(dta)){
    #Only add the thicks in the axis
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, hadj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
  }
  
  
  #Add F name
  mtext(Fnames[fmaj], side = 3, outer = FALSE, line = 1.5, cex = 5)
  
  
  #Loop through the other F IF fmaj < 7
  if(fmaj < 7) {
    
    for(fmin in (fmaj + 1):ncol(dta)){
      
      #Plot first F again the others
      plot(dta[,fmaj] ~ dta[,fmin], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.05, 0.35), ylim = c(-0.05, 0.35), cex = 4)
      
      #Add the equality line
      abline(0,1, lwd = 3)
      
      #Add axis 1 if fmin = ncol(dta)
      if(fmin == ncol(dta)){
        #Add the thicks
        axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(1, xlab = NULL, at=seq(0, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
      } else {
        #just ad the ticks
        axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
      }
      
      #Add axis 2 if fmaj = 1
      if(fmaj == 1){
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(2, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), hadj = 1.5, cex.axis= 5, lwd.ticks = 3.5, las = 1, tck = -0.03)
      } else {
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
      }
    }
  }   
}

dev.off()


###############################################
#### FIG S4 AFR Comparison of different Fs ####

#Subset only the F we are interested in
dta = dtaAFR[,c(4,6,8,10,12,14,16)]

#Set F names
Fnames = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
           expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
           expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Open pdf
pdf("./PLOTvIII/FIGS4_F_COMP_AFR.pdf", height = 35, width = 35)

#Create layout
layout(mat = matrix(c(rep(0,15),0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,rep(0,15),0,0,0,8,0,9,0,10,0,11,0,12,0,13,0,rep(0,15),0,0,0,0,0,14,0,15,0,16,0,17,0,18,0,
                      rep(0,15),0,0,0,0,0,0,0,19,0,20,0,21,0,22,0,rep(0,15),0,0,0,0,0,0,0,0,0,23,0,24,0,25,0,rep(0,15),0,0,0,0,0,0,0,0,0,0,0,26,0,27,0,
                      rep(0,15),0,0,0,0,0,0,0,0,0,0,0,0,0,28,0,rep(0,15)), nrow = 15, ncol = 15, byrow = F), heights = c(.2,1,rep(c(0,1),6),.2), widths = c(.2,rep(c(1,0),6),1,.2))

#Loop through F
for (fmaj in 1:ncol(dta)) {
  #First plot is always to display distributions of individual F
  hist(dta[,fmaj], ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.05, 0.35), ylim = c(0,750), col = "deepskyblue3", main = NULL, breaks = seq(-0.1, 0.4, by = 0.01))
  box(col = "black", lwd = 1)
  
  #If last F, add axis labels, else only add thick
  if(fmaj < ncol(dta)){
    #Only add the thicks in the axis
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, hadj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
  }
  
  
  #Add F name
  mtext(Fnames[fmaj], side = 3, outer = FALSE, line = 1.5, cex = 5)
  
  
  #Loop through the other F IF fmaj < 7
  if(fmaj < 7) {
    
    for(fmin in (fmaj + 1):ncol(dta)){
      
      #Plot first F again the others
      plot(dta[,fmaj] ~ dta[,fmin], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.05, 0.35), ylim = c(-0.05, 0.35), cex = 4)
      
      #Add the equality line
      abline(0,1, lwd = 3)
      
      #Add axis 1 if fmin = ncol(dta)
      if(fmin == ncol(dta)){
        #Add the thicks
        axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(1, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
      } else {
        #just ad the ticks
        axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
      }
      
      #Add axis 2 if fmaj = 1
      if(fmaj == 1){
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(2, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), hadj = 1.5, cex.axis= 5, lwd.ticks = 3.5, las = 1, tck = -0.03)
      } else {
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
      }
    }
  }   
}

dev.off()


#################################################
#### FIG S5 WORLD Comparison of different Fs ####

#Subset only the F we are interested in
dta = dtaWORLD[,c(5,7,9,11,13,15,17)]

#Set F names
Fnames = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
           expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
           expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))
#Open pdf
pdf("./PLOTvIII/FIGS5_F_COMP_WORLD.pdf", height = 40, width = 41)

#Create layout
layout(mat = matrix(c(rep(0,15),0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,rep(0,15),0,0,0,8,0,9,0,10,0,11,0,12,0,13,0,rep(0,15),0,0,0,0,0,14,0,15,0,16,0,17,0,18,0,
                      rep(0,15),0,0,0,0,0,0,0,19,0,20,0,21,0,22,0,rep(0,15),0,0,0,0,0,0,0,0,0,23,0,24,0,25,0,rep(0,15),0,0,0,0,0,0,0,0,0,0,0,26,0,27,0,
                      rep(0,15),0,0,0,0,0,0,0,0,0,0,0,0,0,28,0,rep(0,15)), nrow = 15, ncol = 15, byrow = F), heights = c(.2,1,rep(c(0,1),6),.2), widths = c(.3,rep(c(1,0),6),1,.2))

#Loop through F
for (fmaj in 1:ncol(dta)) {
  #First plot is always to display distributions of individual F
  hist(dta[,fmaj], ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.2, 0.4), ylim = c(0,2250), col = "deepskyblue3", main = NULL, breaks = seq(-0.2, 0.4, by = 0.01))
  box(col = "black", lwd = 1)
  
  #If last F, add axis labels, else only add thick
  if(fmaj < ncol(dta)){
    #Only add the thicks in the axis
    axis(1, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, hadj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(0, 0.4, by = 0.2), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
  }
  
  
  #Add F name
  mtext(Fnames[fmaj], side = 3, outer = FALSE, line = 1.5, cex = 5)
  
  
  #Loop through the other F IF fmaj < 7
  if(fmaj < 7) {
    
    for(fmin in (fmaj + 1):ncol(dta)){
      
      #Plot first F again the others
      plot(dta[,fmaj] ~ dta[,fmin], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.2, 0.4), ylim = c(-0.2, 0.4), cex = 4)
      
      #Add the equality line
      abline(0,1, lwd = 3)
      
      #Add axis 1 if fmin = ncol(dta)
      if(fmin == ncol(dta)){
        #Add the thicks
        axis(1, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(1, xlab = NULL, at=seq(-0.1, 0.3, by = 0.2), padj = 1.5, lwd.ticks = 3, cex.axis = 5, las = 1, tck = -0.03)
      } else {
        #just ad the ticks
        axis(1, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
      }
      
      #Add axis 2 if fmaj = 1
      if(fmaj == 1){
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
        #Add the text
        axis(2, xlab = NULL, at=seq(-0.2, 0.4, by = 0.2), hadj = 1.5, cex.axis= 5, lwd.ticks = 3.5, las = 1, tck = -0.03)
      } else {
        #Add the thicks
        axis(2, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
      }
    }
  }   
}

dev.off()

####################################################
#### FIG S6 PLOTTING F.W VS F.C for comparisons ####

pdf("./PLOTvIII/FIGS6_Fw_VS_Fc.pdf", height = 44, width = 22)

layout(mat = matrix(c(rep(0,15),0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,rep(0,15),0,8,0,9,0,10,0,11,0,12,0,13,0,14,0,rep(0,15),0,15,0,16,0,17,0,18,0,19,0,20,0,21,0,rep(0,15)),
                    nrow = 15, ncol = 7, byrow = F), heights = c(.5,rep(c(1,.1),6),1,.7), widths = c(1,rep(c(1,.1),2),1,.2))


#Set F names
Fnames = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
           expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
           expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#### EAS

dta = dtaEAS

#Start count for different F for names
j = 1

#Loop through Fs
for (i in seq(4, 16, by = 2)) {
  
  #Plot first F again the others
  plot(dta[,i+1] ~ dta[,i], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.05, 0.35), ylim = c(-0.05, 0.35), cex = 4)
  
  #Add the equality line
  abline(0,1, lwd = 3)
  
  #Add EAS LEGEND IF first graph
  if(i == 4){
    #Add F.WORLD legend
    mtext(text = "EAS", side = 3, line = 12, outer = F, cex = 8)
  }
  
  #Add F.WORLD LEGEND IF middle graph
  if(i == 10){
    #Add F.WORLD legend
    mtext(text = expression(F[WORLD]), side = 2, at = 0, line = 30, outer = F, cex = 8)
  }
  
  #Add axis 1 if i = last = 16
  if(i == 16){
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5.5, las = 1, tck = -0.03)
    #else just add the ticks
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  }
  
  #Add axis because EAS
  
  #Add the thicks
  axis(2, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
  #Add the text
  axis(2, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), cex.axis= 5.5, lwd.ticks = 3.5, las = 1, tck = -0.03, hadj = 1.5)
  
  #Add Y axis label becasue EAS
  mtext(text = Fnames[j], side = 2, line = 15, outer = F, cex = 6)
  
  #increment count for which F
  j = j + 1
}

#### AFR

dta = dtaAFR

#Loop through Fs
for (i in seq(4, 16, by = 2)) {
  
  #Plot first F again the others
  plot(dta[,i+1] ~ dta[,i], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.15, 0.35), ylim = c(-0.15, 0.35), cex = 4)
  
  #Add the equality line
  abline(0,1, lwd = 3)
  
  #Add AFR LEGEND IF first graph
  if(i == 4){
    #Add AFR legend
    mtext(text = "AFR", side = 3, line = 12, outer = F, cex = 8)
  }
  
  #Add F.CONTINENT LEGEND IF last graph
  if(i == 16){
    #Add F.WORLD legend
    mtext(text = expression(F[CONTINENT]), side = 1, line = 25, outer = F, cex = 8)
  }
  
  #Add axis 1 if i = last = 16
  if(i == 16){
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(-0.1, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5.5, las = 1, tck = -0.03)
    #else just add the ticks
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  }
  
  #Add axis becasue EAS
  
  #Add the thicks
  axis(2, xlab = NULL, at=seq(-0.1, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
  
}

#### WORLD

dta = dtaWORLD

#Loop through Fs
for (i in seq(4, 16, by = 2)) {
  
  #Plot first F again the others
  plot(dta[,i+1] ~ dta[,i], pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.26, 0.35), ylim = c(-0.26, 0.35), cex = 4)
  
  #Add the equality line
  abline(0,1, lwd = 3)
  
  #Add WORLD LEGEND IF first graph
  if(i == 4){
    #Add WORLD legend
    mtext(text = "WORLD", side = 3, line = 12, outer = F, cex = 8)
  }
  
  #Add axis 1 if i = last = 16
  if(i == 16){
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(-0.2, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 5.5, las = 1, tck = -0.03)
    #else just add the ticks
  } else {
    #Add the thicks 
    axis(1, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  }
  
  #Add axis because EAS
  
  #Add the thicks
  axis(2, xlab = NULL, at=seq(-0.2, 0.4, by = 0.05), labels = NA, lwd.ticks = 3, las = 1, tck = -0.03)
  
}

dev.off()



########################################################
#### FIG S7 PLOTTING FUNI.WE.W VS FUNI.WE.C per pop ####

pdf("./PLOTvIII/FIGS7_F_UNI.WE_Fw_VS_Fc_perPop.pdf", width = 30, height = 33)

layout(mat = matrix(c(rep(0,11),0,1,0,2,0,3,0,4,0,5,0,rep(0,11),0,6,0,7,0,8,0,9,0,10,0,rep(0,11),0,11,0,12,0,13,0,14,0,15,0,rep(0,11),0,16,0,17,0,18,0,19,0,20,0,rep(0,11),
                      0,21,0,22,0,23,0,24,0,25,0,rep(0,11),0,26,0,0,0,0,0,0,0,0,0,rep(0,11)), nrow = 13, ncol = 11, byrow = T),
       heights = c(rep(c(0 ,1),6),.1), widths = c(.5,rep(c(1,0),5)))

#Loop through POPs
for(pop in sort(unique(dtaWORLD$pop))){
  
  #Subset the data
  data = dtaWORLD[dtaWORLD$pop == pop,]
  
  #Plot the data
  plot(data$Funi.we.w ~ data$Funi.we.c, pch = 19, col = "darkorange2", ylab = NA, xlab = NA, yaxt = 'n', xaxt = 'n', xlim = c(-0.05, 0.31), ylim = c(-0.05, 0.31), cex = 3)
  
  #Add the equality line
  abline(0,1, lwd = 3)
  
  #Add POP ID
  mtext(pop, at = 0 , line = -4, cex = 3)
  
  #IF first column, add y axis
  
  if(pop %in% c("ACB", "CHB", "GBR", "JPT", "PEL", "YRI")){
    
    #Add the thicks
    axis(2, xlab = NULL, at=seq(-0.05, 0.3, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(2, xlab = NULL, at=seq(0, 0.3, by = 0.1), lwd.ticks = 3, cex.axis = 4.5, las = 1, tck = -0.03, hadj = 1.5)
    
    #else just add the ticks
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(0, 0.3, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  }
  
  #If last column, add x axis
  if(pop %in% c("YRI", "PJL", "PUR", "STU", "TSI")){
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.05, 0.3, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
    #Add the text
    axis(1, xlab = NULL, at=seq(0, 0.3, by = 0.1), padj = 1.5, lwd.ticks = 3, cex.axis = 4.5, las = 1, tck = -0.03)
    
    #else just add the ticks
  } else {
    #Add the thicks
    axis(1, xlab = NULL, at=seq(-0.05, 0.4, by = 0.05), labels = NA, padj = 1.5, lwd.ticks = 3, las = 1, tck = -0.03)
  }
  
  #Add F.WORLD LEGEND IF pop = GBR
  if(pop == "GBR"){
    #Add F.WORLD legend
    mtext(text = expression(italic(F[UNI]^w ~ "WORLD")), side = 2, line = 15, at = -0.125, outer = F, cex = 6)
  }
  
  #Add F.CONTINENT LEGEND IF pop = PUR
  if(pop == "PUR"){
    #Add F.WORLD legend
    mtext(text = expression(italic(F[UNI]^w ~ "CONTINENT")), side = 1, line = 25, outer = F, cex = 6)
  }
}

dev.off()

######################################################################
########### FIGURE S8 POLYPED LMM: GCTA WE VS LMM: GCTA UN ###########

#PLOT
pdf("./PLOTvIII/FIGS8_LMMGCTAWE_LMMGCTAUN_IDest_POLYPED.pdf", height = 50, width = 32)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5),0,5,0,6,0,rep(0,5)),
                    ncol = 5, nrow = 7, byrow = T), heights = c(0.2,1.0,0.1,1.0,0.1,1.5,0), widths = c(0.5,1,.3,1,0))

par(mar = c(1,2,5,1))

#### PANNEL A: POLYPED FULL SAMPLE SIZE LMM GCTA WE ####

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 6, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR MODEL
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  5, cex = 10)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 11,924)\ncomplete sample size", side = 2, outer = F, line =  30, cex = 5)

#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL B: POLYPED FULL SAMPLE SIZE LMM GCTA UN ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 6, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  5, cex = 10)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL C: POLYPED REDUCED SAMPLE SIZE RANDOM LMM GCTA WE ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(6000), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 2,500)\nrandom sub-sampling", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL D: POLYPED REDUCED SAMPLE SIZE RANDOM LMM GCTA UN ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(6000), labels = NULL, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL E: POLYPED REDUCED SAMPLE SIZE LMM GCTA WE ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(10000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 2,500)\nranged sub-sampling", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL F: POLYPED REDUCED SAMPLE SIZE LMM GCTA UN ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(10000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7, line = -5)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

##################################################################
########### FIGURE S9 1KG LMM: GCTA WE VS LMM: GCTA UN ###########

#PLOT
pdf("./PLOTvIII/FIGS9_LMMGCTAWE_LMMGCTAUN_IDest_1KG.pdf", height = 50, width = 32)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5),0,5,0,6,0,rep(0,5)),
                    ncol = 5, nrow = 7, byrow = T), heights = c(0.2,1.0,0.1,1.0,0.1,1.5,0), widths = c(0.5,1,.3,1,0))

par(mar = c(1,2,5,1))

#### PANNEL A: EAS LM ####

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:31])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[25:31]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[25:31])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR MODEL
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  5, cex = 10)
#ADD MTEXT FOR POP
mtext(text = "1KG: EAS (n = 504)", side = 2, outer = F, line =  35, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL B: EAS LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:41])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[35:41]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[35:41])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR MODEL
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  5, cex = 10)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL C: AFR LM ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:31])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[25:31]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[25:31])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT FOR POP
mtext(text = "1KG: AFR (n = 661)", side = 2, outer = F, line =  35, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL D: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:41])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[35:41]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[35:41])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

#### PANNEL E: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:31])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[25:31]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[25:31])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(10000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT FOR POP
mtext(text = "1KG: WORLD (n = 2,504)", side = 2, outer = F, line =  35, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


#### PANNEL F: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:41])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[35:41]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[35:41])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(10000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 6, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()


######################################################
#### FIGURE S10 ID EST ALL pops standard scenario ####

pdf("./PLOTvIII/FIGS10_IDest_ALLpops_standard_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



##################################################
#### FIGURE S11 ID EST ALL pops h scenario #### 

pdf("./PLOTvIII/FIGS11_IDest_ALLpops_h_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaF[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

##################################################
#### FIGURE S12 ID EST ALL pops p scenario #### 

pdf("./PLOTvIII/FIGS12_IDest_ALLpops_p_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaF[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25000) ,log10(25000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(10000,1000,100,10)),log10(c(1,10,100,1000,10000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis labels
axis(2, at = c(-log10(c(10000,100)),log10(c(1,100,10000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10000,-100,0,100,10000))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25000) ,log10(25000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(10000,1000,100,10)),log10(c(1,10,100,1000,10000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis labels
axis(2, at = c(-log10(c(10000,100)),log10(c(1,100,10000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10000,-100,0,100,10000))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25000) ,log10(25000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(10000,1000,100,10)),log10(c(1,10,100,1000,10000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis labels
axis(2, at = c(-log10(c(10000,100)),log10(c(1,100,10000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10000,-100,0,100,10000))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25000) ,log10(25000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(10000,1000,100,10)),log10(c(1,10,100,1000,10000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis labels
axis(2, at = c(-log10(c(10000,100)),log10(c(1,100,10000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10000,-100,0,100,10000))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



##################################################
#### FIGURE S13 ID EST ALL pops DEMA scenario ####

pdf("./PLOTvIII/FIGS13_IDest_ALLpops_dema_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaT[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhF_selpF_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhF_selpF_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpF_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



##################################################
#### FIGURE S14 ID EST ALL pops h & p scenario ####

pdf("./PLOTvIII/FIGS14_IDest_ALLpops_hp_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaF[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(6000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000)),log10(c(1,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,0,1000,100000))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(6000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000)),log10(c(1,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,0,1000,100000))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(6000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000)),log10(c(1,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,0,1000,100000))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaF[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(6000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000)),log10(c(1,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,0,1000,100000))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



######################################################
#### FIGURE S15 ID EST ALL pops h & dema scenario ####

pdf("./PLOTvIII/FIGS15_IDest_ALLpops_hdema_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpF_demaT[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhT_selpF_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhT_selpF_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpF_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(1000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

######################################################
#### FIGURE S16 ID EST ALL pops p & dema scenario ####

pdf("./PLOTvIII/FIGS16_IDest_ALLpops_pdema_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000,10)),log10(c(1,10,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,-10,0,10,1000,100000))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000,10)),log10(c(1,10,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,-10,0,10,1000,100000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[, c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[, c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[, c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000,10)),log10(c(1,10,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,-10,0,10,1000,100000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhF_selpT_demaT[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250000) ,log10(250000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(100000,1000,10)),log10(c(1,10,1000,100000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100000,-1000,-10,0,10,1000,100000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

#set margins
par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhF_selpT_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhF_selpT_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500000) ,log10(2500000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(100000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(1000000,100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(1000000,10000,100)),log10(c(1,100,10000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000000,-10000,-100,0,100,10000,1000000))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500000) ,log10(2500000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(100000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(1000000,100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(1000000,10000,100)),log10(c(1,100,10000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000000,-10000,-100,0,100,10000,1000000))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500000) ,log10(2500000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(100000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(1000000,100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(1000000,10000,100)),log10(c(1,100,10000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000000,-10000,-100,0,100,10000,1000000))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhF_selpT_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500000) ,log10(2500000)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(100000000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis thicks
axis(2, at = c(-log10(c(1000000,100000,10000,1000,100,10)),log10(c(1,10,100,1000,10000,100000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = NA)
#Add y axis text
axis(2, at = c(-log10(c(1000000,10000,100)),log10(c(1,100,10000,1000000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000000,-10000,-100,0,100,10000,1000000))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



#################################################
#### FIGURE S17 ID EST ALL pops ALL scenario ####

pdf("./PLOTvIII/FIGS17_IDest_ALLpops_ALL_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#SET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[, c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[, c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[, c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

#set margins
par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4),rep("springgreen",3))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:14])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:14]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "EAS", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:24])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:24]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:34])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:34]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:44])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:44]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:14])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:14]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "AFR", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:24])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:24]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:34])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:34]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:44])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:44]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[AS[LDMS]])), expression(italic(F[UNI[LDMS]]^u)), expression(italic(F[UNI[LDMS]]^w)))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:14])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:14]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:24])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:24]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:24])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:34])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:34]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:44])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:44]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:44])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

#######################################################################
#### FIGURE S18 ID EST ALL pops ALL scenario filtering on MAF 0.05 ####

#Load data
load("./ID_EST/WORLD_1000_genomes_LMM/WORLDMAF0.05_ID3_simus_COMPLETE_ID_EST.RData")

#PLOT
pdf("./PLOTvIII/FIGS18_IDest_WORLD_MAF0.05_ALL_scenario.pdf", height = 14, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9)),
                    ncol = 9, nrow = 3, byrow = T), heights = c(0.2,1,0), widths = c(0.1,1,.3,1,.3,1,.3,1,0))

par(mar = c(1,2,5,1))


#### PANNEL A: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

colours = c(rep("#3366FF",3),rep("darkorange2",5))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:10])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:10]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:10])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(35) ,log10(35)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(100), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,11:17])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,11:17]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,11:17])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(35) ,log10(35)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(100), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,21:27])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,21:27]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,21:27])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(35) ,log10(35)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(100), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[31:37])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,31:37]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,31:37])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(35) ,log10(35)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(100), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()


###############################################
#### FIGURE S19 ID EST WORLD BCFTools Qual ####

#Load data
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_FROHBCFToolsQual_ID3_simus_COMPLETE_ID_EST.RData")
#because of me and names problems
ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT = ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT
rm(list = ls(pattern = "ID_EST_WORLD_LMM"))
#Load World
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_ID3_simus_COMPLETE_ID_EST.RData")

#Merge WORLD all Fs + BCFTools.Qual
ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT = cbind(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,1:7], FROHs.BCFTools.Qual = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,1],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,8:17], FROHs.BCFTools.Qual.lme.AS = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,2],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,18:27], FROHs.BCFTools.Qual.lme.AS = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,3],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,28:37], FROHs.BCFTools.Qual.lme.AS = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,4],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,38:44])

#PLOT
pdf("./PLOTvIII/FIGS19_IDest_WORLD_BCFTools.Qual_ALL_scenario.pdf", height = 14, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9)),
                    ncol = 9, nrow = 3, byrow = T), heights = c(0.2,1,0), widths = c(0.1,1,.3,1,.3,1,.3,1,0))

par(mar = c(1,2,5,1))


#### PANNEL I: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"Qual."})),expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

colours = c(rep("#3366FF",3),rep("darkorange2",5))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,5:12])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,5:12]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,5:12])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,16:23])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,16:23]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,16:23])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,27:34])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,27:34]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,27:34])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,38:45])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,38:45]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,38:45])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

#################################################################
#### FIGURE S20 ID EST: Intermediate frequencies causal loci ####

load("./ID_EST/EAS_1000_genomes_LMM/EAS_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")
load("./ID_EST/AFR_1000_genomes_LMM/AFR_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")
load("./ID_EST/POLYPED_20cM_LMM/POLYPED_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")

pdf("./PLOTvIII/FIGS20_IDest_ALLpops_IntermediateCausalLoci_ALL_scenario.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1.5,0,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: POLYPED LM ####

#SET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])),expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(5,6,8:13)])+1)*(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(5,6,8:13)]/abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(5,6,8:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE (n = 11,924)\ncomplete sample size", side = 2, outer = F, line =  25, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: POLYPED LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(14,15,17:22)])+1)*(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(14,15,17:22)]/abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(14,15,17:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: POLYPED LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[, c(23,24,26:31)])+1)*(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[, c(23,24,26:31)]/abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[, c(23,24,26:31)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: POLYPED LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(32,33,35:40)])+1)*(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(32,33,35:40)]/abs(ID_EST_PGPED_LMM_ID3_intermediateLoci_selhT_selpT_demaT[,c(32,33,35:40)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: EAS LM ####

#set margins
par(mar = c(1,2,5,1))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4))

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:11])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:11]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,5:11])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR POP
mtext(text = "1KG: EAS\n(n = 504)", side = 2, outer = F, line =  25, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:21])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:21]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,15:21])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:31])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:31]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,25:31])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:41])+1)*(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:41]/abs(ID_EST_EAS_LMM_ID3_selhT_selpT_demaT[,35:41])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: AFR LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:11])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:11]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,5:11])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

## ADD MTEXT FOR POP
mtext(text = "1KG: AFR\n(n = 661)", side = 2, outer = F, line =  25, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:21])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:21]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,15:21])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:31])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:31]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,25:31])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:41])+1)*(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:41]/abs(ID_EST_AFR_LMM_ID3_selhT_selpT_demaT[,35:41])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:11])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:11]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,5:11])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "1KG: WORLD\n(n = 2,504)", side = 2, outer = F, line =  25, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:21])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:21]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,15:21])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:31])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:31]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,25:31])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:41])+1)*(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:41]/abs(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,35:41])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(25) ,log10(25)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(60), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(10)),log10(c(1,10))), las = 1, cex.axis = 4, hadj = 1, labels = c(-10,0,10))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

#######################################################################
#### FIGURE SXX ID EST TwoPops SecondPop SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/SecondAdmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3_simus_FasbothFuni_hpdema_ID_EST.RData")

pdf("./PLOTvIII/FIGSZ_TwoPops_admixed_secondAdmixedPop.pdf", height = 30, width = 29)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)),
                    ncol = 5, nrow = 5, byrow = T), heights = c(0.3,1.5,0.2,1.5,0), widths = c(0.5,1,.6,1,0))

#change margins for x axis labels
par(mar = c(20,1,2,1))

#Set colors
colours = c(rep("#3366FF",4))
#Set lables
labs = c("LM", expression(LMM[AS]), expression(LMM[GCTA]^w), expression(LMM[GCTA]^u))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,7])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,7]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,7])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,10])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,10]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,10])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,13])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,13]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,13])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,16])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,16]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,16]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 10
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[AS])), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -.75, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

######################## FuniUN ######################## 

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,8])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,8]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,8])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,11])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,11]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,11])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,14])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,14]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,14])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,17])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,17]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,17]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[UNI]^u)), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -.75, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

######################## FuniWE ######################## 

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,9])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,9]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,9])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,12])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,12]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,12])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,15])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,15]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,15])),
              log10(abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,18])+1)*(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,18]/abs(SecondadmixedPop_SUB_POLYPED_TwoPopsadmixed_ID3[,18]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(italic(F[UNI]^w)), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -.75, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()











###########################################################################
#### FIGURE SXX ID EST TwoPops admixed RANDOM SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/WORLD_RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3_simus_Funiwe_hpdema_ID_EST.RData")

pdf("./PLOTvIII/TEstingTwoPopsFig.pdf", height = 30, width = 29)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)),
                    ncol = 5, nrow = 5, byrow = T), heights = c(0.3,1.5,0.2,1.5,0), widths = c(0.7,1,.6,1,0))

#change margins for x axis labels
par(mar = c(25,1,2,1))

#Set colors
colours = c(rep("#3366FF",4))
#Set lables
labs = c("LM", expression(LMM[AS]), expression(LMM[GCTA[WE]]), expression(LMM[GCTA[UN]]))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,7])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,7]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,7])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,8])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,8]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,8])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,9])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,9]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,9])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,10])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,10]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,10]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = "TwoPops admixed\nRANDOM SUB", line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)


##########################################################################
#### FIGURE SXX ID EST TwoPopsadmixed RANGED SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/WORLD_RANGED_SUB_POLYPED_TwoPopsadmixed_ID3_simus_Funiwe_hpdema_ID_EST.RData")

#Set colors
colours = c(rep("#3366FF",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,7])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,7]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,7])),
              log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,8])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,8]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,8])),
              log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,9])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,9]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,9])),
              log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,10])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,10]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,10]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = "TwoPops admixed\nRANGED SUB", line = 3, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

###########################################################################
#### FIGURE SXX ID EST TwoPops COALESC RANDOM SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/WORLD_RANDOM_SUB_POLYPED_TwoPopscoalesc_ID3_simus_Funiwe_hpdema_ID_EST.RData")

#Set colors
colours = c(rep("#3366FF",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,7])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,7]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,7])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,8])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,8]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,8])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,9])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,9]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,9])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,10])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,10]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID3[,10]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = "TwoPops coaslesc\nRANDOM SUB", line = 4, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

###########################################################################
#### FIGURE SXX ID EST TwoPops COALESC RANGED SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/WORLD_RANGED_SUB_POLYPED_TwoPopscoalesc_ID3_simus_Funiwe_hpdema_ID_EST.RData")

#Set colors
colours = c(rep("#3366FF",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,7])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,7]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,7])),
              log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,8])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,8]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,8])),
              log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,9])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,9]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,9])),
              log10(abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,10])+1)*(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,10]/abs(RANGED_SUB_POLYPED_TwoPopsadmixed_ID3[,10]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = "TwoPops coaslesc\nRANGED SUB", line = 5, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

############################################################################
#### FIGURE SXX ID EST TwoPops ADMIXED LASTGEN SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/WORLD_LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3_simus_FasbothFuni_hpdema_ID_EST.RData")

pdf("./PLOTvIII/TwoPops_admixed_LastGenSUB.pdf", height = 30, width = 29)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)),
                    ncol = 5, nrow = 5, byrow = T), heights = c(0.3,1.5,0.2,1.5,0), widths = c(0.7,1,.6,1,0))

#change margins for x axis labels
par(mar = c(25,1,2,1))

#Set colors
colours = c(rep("#3366FF",4))
#Set lables
labs = c("LM", expression(LMM[AS]), expression(LMM[GCTA[WE]]), expression(LMM[GCTA[UN]]))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,7])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,7]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,7])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,10])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,10]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,10])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,13])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,13]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,13])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,16])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,16]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,16]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(F[AS]), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

######################## FuniUN ######################## 

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,8])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,8]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,8])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,11])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,11]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,11])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,14])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,14]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,14])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,17])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,17]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,17]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(F[UN.UN]), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

######################## FuniWE ######################## 

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,9])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,9]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,9])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,12])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,12]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,12])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,15])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,15]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,15])),
              log10(abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,18])+1)*(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,18]/abs(LASTGEN_SUB_POLYPED_TwoPopsadmixed_ID3[,18]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(F[UN.WE]), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

#########################################################################
#### FIGURE SXX ID EST TwoPops ADMIXED b=10 SUBSAMPLING ADD DOM DEMA ####

load("./ID_EST/TwoPops_SUBSAMPLING/WORLD_RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10_simus_FasbothFuni_hpdema_ID_EST.RData")

pdf("./PLOTvIII/TwoPops_admixed_b10.pdf", height = 30, width = 29)

#layout
layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5)),
                    ncol = 5, nrow = 5, byrow = T), heights = c(0.3,1.5,0.2,1.5,0), widths = c(0.7,1,.6,1,0))

#change margins for x axis labels
par(mar = c(25,1,2,1))

#Set colors
colours = c(rep("#3366FF",4))
#Set lables
labs = c("LM", expression(LMM[AS]), expression(LMM[GCTA[WE]]), expression(LMM[GCTA[UN]]))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,7])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,7]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,7])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,10])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,10]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,10])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,13])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,13]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,13])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,16])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,16]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,16]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 10
abline(h = -log10(11), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(F[AS]), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

######################## FuniUN ######################## 

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,8])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,8]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,8])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,11])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,11]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,11])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,14])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,14]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,14])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,17])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,17]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,17]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(11), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(F[UN.UN]), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

######################## FuniWE ######################## 

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,9])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,9]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,9])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,12])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,12]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,12])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,15])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,15]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,15])),
              log10(abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,18])+1)*(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,18]/abs(RANDOM_SUB_POLYPED_TwoPopsadmixed_ID10[,18]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(11), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = expression(F[UN.WE]), line = 2, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

#########################################
#### FIGURE SZZZ ID EST WORLD no AFR ####

load("./ID_EST/WORLDnoAFR_1000_genomes_LMM/WORLDnoAFR_ID3_simus_COMPLETE_ID_EST.RData")

#Set colors
colours = c(rep("#3366FF",4))
#Set lables
labs = c("LM", expression(LMM[AS]), expression(LMM[GCTA[WE]]), expression(LMM[GCTA[UN]]))

#Plot the Forth, selfT, selpT, demaT
vioplot(cbind(log10(abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,7])+1)*(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,7]/abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,7])),
              log10(abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,10])+1)*(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,10]/abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,10])),
              log10(abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,13])+1)*(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,13]/abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,13])),
              log10(abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,16])+1)*(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,16]/abs(ID_EST_WORLD_noAFR_LMM_ID3_selhT_selpT_demaT[,16]))),
        col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250),log10(250)))

#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:4, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:4, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 6)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 6, hadj = 1, labels = c(-100,-10,0,10,100))

#Add panel title
mtext(side = 3, outer = F, text = "WORLD no AFR", line = 5, cex = 5, at = 2.5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7, line = -5)
#ADD AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()





###################################################################
#### FIGURE SXX ID EST POLYPED RANDOM SUBSAMPLING ALL scenario ####

#Load polypeds ubsampling AFR
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/AFR_SUB_RANDOM_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
AFR_SUB_RANDOM_POLYPED_ID3 = foreachoutput
#subsample the scenario
AFR_SUB_RANDOM_POLYPED_ID3 = AFR_SUB_RANDOM_POLYPED_ID3[AFR_SUB_RANDOM_POLYPED_ID3$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",]

#Load polypeds ubsampling EAS
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/EAS_SUB_RANDOM_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
EAS_SUB_RANDOM_POLYPED_ID3 = foreachoutput
#subsample the scenario
EAS_SUB_RANDOM_POLYPED_ID3 = EAS_SUB_RANDOM_POLYPED_ID3[EAS_SUB_RANDOM_POLYPED_ID3$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",]

#Load polypeds ubsampling EAS
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/WORLD_SUB_RANDOM_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
WORLD_SUB_RANDOM_POLYPED_ID3 = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",]

#rm
rm(foreachoutput)
rm(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)

pdf("./PLOTvII/FIGSXX_IDest_POLYPED_RANDOM_SUBSAMPLING.pdf", height = 26, width = 57)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9)),
                    ncol = 9, nrow = 7, byrow = T), heights = c(0.4,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

par(mar = c(1,5,5,1))

#### PANNEL A: EAS LM ####

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])+1)*(EAS_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)]/abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  20, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n(SUB: EAS)", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])+1)*(EAS_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)]/abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 20, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])+1)*(EAS_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)]/abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA[WE]]), side = 3, outer = F, line =  18, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])+1)*(EAS_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)]/abs(EAS_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA[WE]]), side = 3, outer = F, line =  18, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: AFR LM ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])+1)*(AFR_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)]/abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n(SUB: AFR)", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])+1)*(AFR_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)]/abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])+1)*(AFR_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)]/abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])+1)*(AFR_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)]/abs(AFR_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(25,5,5,1))

#Set F names
labs = c(expression(F[PED]), expression(F[AS]), expression(F[uniUN]), expression(F[uniWE]),
         expression(F[HBD.100KB]), expression(F[ROH.100KB]),
         expression(F[HBD.1MB]), expression(F[ROH.1MB]))

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n(SUB: WORLD)", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])+1)*(WORLD_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)]/abs(WORLD_SUB_RANDOM_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(6000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()






###################################################################
#### FIGURE SYY ID EST POLYPED RANGED SUBSAMPLING ALL scenario ####

#Load polypeds ubsampling AFR
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/AFR_SUB_RANGED_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
AFR_SUB_RANGED_POLYPED_ID3 = foreachoutput
#subsample the scenario
AFR_SUB_RANGED_POLYPED_ID3 = AFR_SUB_RANGED_POLYPED_ID3[AFR_SUB_RANGED_POLYPED_ID3$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",]

#Load polypeds ubsampling EAS
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/EAS_SUB_RANGED_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
EAS_SUB_RANGED_POLYPED_ID3 = foreachoutput
#subsample the scenario
EAS_SUB_RANGED_POLYPED_ID3 = EAS_SUB_RANGED_POLYPED_ID3[EAS_SUB_RANGED_POLYPED_ID3$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",]

#Load polypeds ubsampling EAS
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/WORLD_SUB_RANGED_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
WORLD_SUB_RANGED_POLYPED_ID3 = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",]

#rm
rm(foreachoutput)
rm(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)

pdf("./PLOTvII/FIGSYY_IDest_POLYPED_RANGED_SUBSAMPLING.pdf", height = 26, width = 57)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9)),
                    ncol = 9, nrow = 7, byrow = T), heights = c(0.4,1,0,1,0,1.5,0), widths = c(0.5,1,.3,1,.3,1,.3,1,0))

par(mar = c(1,5,5,1))

#### PANNEL A: EAS LM ####

#Set colors
colours = c("firebrick",rep("#3366FF",3),rep("darkorange2",4))

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])+1)*(EAS_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)]/abs(EAS_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR MODEL
mtext(text = "LM", side = 3, outer = F, line =  20, cex = 5)
## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n(SUB: EAS)", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: EAS LMM AS ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])+1)*(EAS_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)]/abs(EAS_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line = 20, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: EAS LMM kc0 ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])+1)*(EAS_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)]/abs(EAS_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA[WE]]), side = 3, outer = F, line =  18, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: EAS LMM gcta ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(EAS_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])+1)*(EAS_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)]/abs(EAS_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA[UN]]), side = 3, outer = F, line =  18, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: AFR LM ####

#Plot the first, selfF, selpF, demaF
boxplot(log10(abs(AFR_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])+1)*(AFR_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)]/abs(AFR_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n(SUB: AFR)", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: AFR LMM AS ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(AFR_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])+1)*(AFR_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)]/abs(AFR_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: AFR LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(AFR_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])+1)*(AFR_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)]/abs(AFR_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: AFR LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(AFR_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])+1)*(AFR_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)]/abs(AFR_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(25,5,5,1))

#Set F names
labs = c(expression(F[PED]), expression(F[AS]), expression(F[uniUN]), expression(F[uniWE]),
         expression(F[HBD.100KB]), expression(F[ROH.100KB]),
         expression(F[HBD.1MB]), expression(F[ROH.1MB]))

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(7:8,10:15)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n(SUB: WORLD)", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(16:17,19:24)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(25:26,28:33)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
boxplot(log10(abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])+1)*(WORLD_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)]/abs(WORLD_SUB_RANGED_POLYPED_ID3[,c(34:35,37:42)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = -1, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(6000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression(LOG[10] ~ "(b)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()
