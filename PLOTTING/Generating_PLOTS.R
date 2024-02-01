library("vioplot")

#### LOAD

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
pdf("./PLOTS/FIG1_LM_LMMAS_IDest_POLYPED.pdf", height = 50, width = 32)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

#####################################################
########### FIGURE 2 1KG LM VS LMM GRM AS ###########

#PLOT
pdf("./PLOTS/FIG2_LM_LMMAS_IDest_1KG.pdf", height = 50, width = 32)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

####################################################
#### FIGURE 3 PLOTTING F.UNI.WE for comparisons ####

pdf("./PLOTS/FIG3_FuniWE_ALL.pdf", height = 30, width = 25)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
pdf("./PLOTS/FIGS2_F_COMP_POLYPED.pdf", height = 40, width = 41)

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
pdf("./PLOTS/FIGS3_F_COMP_EAS.pdf", height = 35, width = 35)

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
pdf("./PLOTS/FIGS4_F_COMP_AFR.pdf", height = 35, width = 35)

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
pdf("./PLOTS/FIGS5_F_COMP_WORLD.pdf", height = 40, width = 41)

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

######################################################################
########### FIGURE S6 POLYPED LMM: GCTA WE VS LMM: GCTA UN ###########

#PLOT
pdf("./PLOTS/FIGS6_LMMGCTAWE_LMMGCTAUN_IDest_POLYPED.pdf", height = 50, width = 32)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()

###############################################
#### FIGURE S7 ID EST SUBSAMPLING PEDIGREE ####

pdf("./PLOTS/FIGS7_IDest_PEDIGREEsub.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1,0,1,0,1,0,1.5,0), widths = c(0.6,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: 50 INDVs LM ####

#SET MARGINGS FOR X AXIS
par(mar = c(1,2,5,1))

#Set F names
labs = c(expression(italic(F[PED])), expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c("firebrick", rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(5:12)])+1)*(ID_EST_PGPED_LMM_ID3_SUB50[,c(5:12)]/abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(5:12)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n50 individuals", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: 50 INDVs LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(13:20)])+1)*(ID_EST_PGPED_LMM_ID3_SUB50[,c(13:20)]/abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(13:20)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(600), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: 50 INDVs LMM GCTA we ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(21:28)])+1)*(ID_EST_PGPED_LMM_ID3_SUB50[,c(21:28)]/abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(21:28)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(8000), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: 50 INDVs LMM GCTA un ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(29:36)])+1)*(ID_EST_PGPED_LMM_ID3_SUB50[,c(29:36)]/abs(ID_EST_PGPED_LMM_ID3_SUB50[,c(29:36)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(8000), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: 100 INDVs LM ####

#set margins
par(mar = c(1,2,5,1))

#Set colors

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(5:12)])+1)*(ID_EST_PGPED_LMM_ID3_SUB100[,c(5:12)]/abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(5:12)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n100 individuals", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: 100 INDVs LMM as ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(13:20)])+1)*(ID_EST_PGPED_LMM_ID3_SUB100[,c(13:20)]/abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(13:20)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: 100 INDVs LMM CGTA we ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(21:28)])+1)*(ID_EST_PGPED_LMM_ID3_SUB100[,c(21:28)]/abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(21:28)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: 100 INDVs LMM GCTA un ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(29:36)])+1)*(ID_EST_PGPED_LMM_ID3_SUB100[,c(29:36)]/abs(ID_EST_PGPED_LMM_ID3_SUB100[,c(29:36)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: 250 INDVs LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(5:12)])+1)*(ID_EST_PGPED_LMM_ID3_SUB250[,c(5:12)]/abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(5:12)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n250 individuals", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: 250 INDVs LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(13:20)])+1)*(ID_EST_PGPED_LMM_ID3_SUB250[,c(13:20)]/abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(13:20)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: 250 INDVs LMM GCTA we ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(21:28)])+1)*(ID_EST_PGPED_LMM_ID3_SUB250[,c(21:28)]/abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(21:28)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: 250 INDVs LMM GCTA un ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(29:36),])+1)*(ID_EST_PGPED_LMM_ID3_SUB250[,c(29:36)]/abs(ID_EST_PGPED_LMM_ID3_SUB250[,c(29:36)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: 500 INDVs LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(5:12)])+1)*(ID_EST_PGPED_LMM_ID3_SUB500[,c(5:12)]/abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(5:12)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(2000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "PEDIGREE\n500 individuals", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: 500 INDVs LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(13:20)])+1)*(ID_EST_PGPED_LMM_ID3_SUB500[,c(13:20)]/abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(13:20)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(2000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: 500 INDVs LMM GCTA we ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(21:28)])+1)*(ID_EST_PGPED_LMM_ID3_SUB500[,c(21:28)]/abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(21:28)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: 500 INDVs LMM GCTA un ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(29:36)])+1)*(ID_EST_PGPED_LMM_ID3_SUB500[,c(29:36)]/abs(ID_EST_PGPED_LMM_ID3_SUB500[,c(29:36)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:8, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:8, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()


##################################################################
########### FIGURE S8 1KG LMM: GCTA WE VS LMM: GCTA UN ###########

#PLOT
pdf("./PLOTS/FIGS8_LMMGCTAWE_LMMGCTAUN_IDest_1KG.pdf", height = 50, width = 32)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)


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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 15)

dev.off()


###################################################
#### FIGURE S9 ID EST SUBSAMPLING WORLD Loots ####

pdf("./PLOTS/FIGS9_IDest_WORLDsub.pdf", height = 36, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9),0,5,0,6,0,7,0,8,0,rep(0,9),0,9,0,10,0,11,0,12,0,rep(0,9),0,13,0,14,0,15,0,16,0,rep(0,9)),
                    ncol = 9, nrow = 9, byrow = T), heights = c(0.4,1,0,1,0,1,0,1.5,0), widths = c(0.6,1,.3,1,.3,1,.3,1,0))

#### PANNEL A: 50 INDVs LM ####

#SET MARGINGS FOR X AXIS
par(mar = c(1,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})))

#Set colors
colours = c(rep("#3366FF",3),rep("darkorange2",4))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(7:13)])+1)*(ID_EST_WORLD_LMM_ID3_SUB50[,c(7:13)]/abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(7:13)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "1KG: WORLD\n50 individuals", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL B: 50 INDVs LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(16:22)])+1)*(ID_EST_WORLD_LMM_ID3_SUB50[,c(16:22)]/abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(16:22)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(600), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL C: 50 INDVs LMM GCTA we ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(23:29)])+1)*(ID_EST_WORLD_LMM_ID3_SUB50[,c(23:29)]/abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(23:29)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(8000), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL D: 50 INDVs LMM GCTA un ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(30:36)])+1)*(ID_EST_WORLD_LMM_ID3_SUB50[,c(30:36)]/abs(ID_EST_WORLD_LMM_ID3_SUB50[,c(30:36)])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(8000), labels = NA, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  15, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -1.6, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL E: 100 INDVs LM ####

#set margins
par(mar = c(1,2,5,1))

#Set colors

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB100[,7:13])+1)*(ID_EST_WORLD_LMM_ID3_SUB100[,7:13]/abs(ID_EST_WORLD_LMM_ID3_SUB100[,7:13])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "1KG: WORLD\n100 individuals", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "E", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL F: 100 INDVs LMM as ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB100[,16:22])+1)*(ID_EST_WORLD_LMM_ID3_SUB100[,16:22]/abs(ID_EST_WORLD_LMM_ID3_SUB100[,16:22])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "F", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL G: 100 INDVs LMM CGTA we ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB100[,23:29])+1)*(ID_EST_WORLD_LMM_ID3_SUB100[,23:29]/abs(ID_EST_WORLD_LMM_ID3_SUB100[,23:29])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "G", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL H: 100 INDVs LMM GCTA un ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB100[,30:36])+1)*(ID_EST_WORLD_LMM_ID3_SUB100[,30:36]/abs(ID_EST_WORLD_LMM_ID3_SUB100[,30:36])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "H", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL I: 250 INDVs LM ####

#Plot the first, selfF, selpF, demaF
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB250[,7:13])+1)*(ID_EST_WORLD_LMM_ID3_SUB250[,7:13]/abs(ID_EST_WORLD_LMM_ID3_SUB250[,7:13])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "1KG: WORLD\n250 individuals", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "I", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: 250 INDVs LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB250[,16:22])+1)*(ID_EST_WORLD_LMM_ID3_SUB250[,16:22]/abs(ID_EST_WORLD_LMM_ID3_SUB250[,16:22])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "J", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: 250 INDVs LMM GCTA we ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB250[,23:29])+1)*(ID_EST_WORLD_LMM_ID3_SUB250[,23:29]/abs(ID_EST_WORLD_LMM_ID3_SUB250[,23:29])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "K", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: 250 INDVs LMM GCTA un ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB250[,30:36])+1)*(ID_EST_WORLD_LMM_ID3_SUB250[,30:36]/abs(ID_EST_WORLD_LMM_ID3_SUB250[,30:36])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "L", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL M: 500 INDVs LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB500[,7:13])+1)*(ID_EST_WORLD_LMM_ID3_SUB500[,7:13]/abs(ID_EST_WORLD_LMM_ID3_SUB500[,7:13])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(2000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

## ADD MTEXT FOR POP
mtext(text = "1KG: WORLD\n500 individuals", side = 2, outer = F, line =  30, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "M", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL N: 500 INDVs LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB500[,16:22])+1)*(ID_EST_WORLD_LMM_ID3_SUB500[,16:22]/abs(ID_EST_WORLD_LMM_ID3_SUB500[,16:22])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(250) ,log10(250)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(2000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(100,10)),log10(c(1,10,100))), las = 1, cex.axis = 4, hadj = 1, labels = c(-100,-10,0,10,100))

#ADD MTEXT PANNEL
mtext(text = "N", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL O: 500 INDVs LMM GCTA we ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB500[,23:29])+1)*(ID_EST_WORLD_LMM_ID3_SUB500[,23:29]/abs(ID_EST_WORLD_LMM_ID3_SUB500[,23:29])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "O", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL P: 500 INDVs LMM GCTA un ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_LMM_ID3_SUB500[,30:36])+1)*(ID_EST_WORLD_LMM_ID3_SUB500[,30:36]/abs(ID_EST_WORLD_LMM_ID3_SUB500[,30:36])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:7, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:7, y = -log10(20000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

#ADD MTEXT PANNEL
mtext(text = "P", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

######################################################
#### FIGURE S10 ID EST ALL pops standard scenario ####

pdf("./PLOTS/FIGS10_IDest_ALLpops_standard_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



##################################################
#### FIGURE S11 ID EST ALL pops h scenario #### 

pdf("./PLOTS/FIGS11_IDest_ALLpops_h_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

##################################################
#### FIGURE S12 ID EST ALL pops p scenario #### 

pdf("./PLOTS/FIGS12_IDest_ALLpops_p_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



##################################################
#### FIGURE S13 ID EST ALL pops DEMA scenario ####

pdf("./PLOTS/FIGS13_IDest_ALLpops_dema_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



##################################################
#### FIGURE S14 ID EST ALL pops h & p scenario ####

pdf("./PLOTS/FIGS14_IDest_ALLpops_hp_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



######################################################
#### FIGURE S15 ID EST ALL pops h & dema scenario ####

pdf("./PLOTS/FIGS15_IDest_ALLpops_hdema_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

######################################################
#### FIGURE S16 ID EST ALL pops p & dema scenario ####

pdf("./PLOTS/FIGS16_IDest_ALLpops_pdema_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()



#################################################
#### FIGURE S17 ID EST ALL pops ALL scenario ####

pdf("./PLOTS/FIGS17_IDest_ALLpops_ALL_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

#######################################################################
#### FIGURE S18 ID EST ALL pops ALL scenario filtering on MAF 0.05 ####

#Load data
load("./ID_EST/WORLD_1000_genomes_LMM/WORLDMAF0.05_ID3_simus_COMPLETE_ID_EST.RData")

#PLOT
pdf("./PLOTS/FIGS18_IDest_WORLD_MAF0.05_ALL_scenario.pdf", height = 14, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9)),
                    ncol = 9, nrow = 3, byrow = T), heights = c(0.2,1,0), widths = c(0.1,1,.3,1,.3,1,.3,1,0))

par(mar = c(1,2,5,1))


#### PANNEL A: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,10,5,1))

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()


###############################################
#### FIGURE S19 ID EST WORLD BCFTools Qual ####

#Load data
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_FROHBCFToolsQual_ID3_simus_COMPLETE_ID_EST.RData")
#because of me and names problems
ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT = ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT
#rm(list = ls(pattern = "ID_EST_WORLD_LMM"))
#Load World
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_ID3_simus_COMPLETE_ID_EST.RData")

#Merge WORLD all Fs + BCFTools.Qual
ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT = cbind(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,1:7], FROHs.BCFTools.Qual = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,1],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,8:17], FROHs.BCFTools.Qual.lme.AS = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,2],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,18:27], FROHs.BCFTools.Qual.lme.AS = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,3],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,28:37], FROHs.BCFTools.Qual.lme.AS = ID_EST_WORLD_BCFToolsQualONLY_LMM_ID3_selhT_selpT_demaT[,4],
                                                                ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT[,38:44])

#PLOT
pdf("./PLOTS/FIGS19_IDest_WORLD_BCFTools.Qual_ALL_scenario.pdf", height = 14, width = 54)

#layout WITH POLYPED
layout(mat = matrix(c(rep(0,9),0,1,0,2,0,3,0,4,0,rep(0,9)),
                    ncol = 9, nrow = 3, byrow = T), heights = c(0.2,1,0), widths = c(0.25,1,.3,1,.3,1,.3,1,0))

par(mar = c(1,2,5,1))


#### PANNEL I: WORLD LM ####

#RESET MARGINGS FOR X AXIS
par(mar = c(20,2,5,1))

#Set F names
labs = c(expression(italic(F[AS])), expression(italic(F[UNI]^u)), expression(italic(F[UNI]^w)),
         expression(italic(F[HBD]^{"Qual."})),expression(italic(F[HBD]^{"100KB"})), expression(italic(F[ROH]^{"100KB"})),
         expression(italic(F[HBD]^{"1MB"})), expression(italic(F[ROH]^{"1MB"})),
         expression(italic(F[HBD]^{"5MB"})), expression(italic(F[ROH]^{"5MB"})))

colours = c(rep("#3366FF",3),rep("darkorange2",7))

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,5:14])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,5:14]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,5:14])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(8000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR POP
mtext(text = "WORLD", side = 2, outer = F, line =  30, cex = 5)
## ADD MTEXT FOR LM
mtext(text = "LM", side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "A", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL J: WORLD LMM AS ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,18:27])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,18:27]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,18:27])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(8000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[AS]), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "B", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL K: WORLD LMM kc0 ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,31:40])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,31:40]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,31:40])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(8000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^w), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "C", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

#### PANNEL L: WORLD LMM gcta ####

#Plot the Forth, selfT, selpT, demaT
vioplot(log10(abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,44:53])+1)*(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,44:53]/abs(ID_EST_WORLD_withBCFToolsQual_LMM_ID3_selhT_selpT_demaT[,44:53])), col = colours, yaxt = 'n', xaxt = 'n', bty = 'n', frame.plot = F, ylim = c(-log10(2500) ,log10(2500)))
#Add the line for ID = 3
abline(h = -log10(4), lty = 1, lwd = 9, col = "grey30")
#Add red line for ID = 0
abline(h = log10(1), lty = 2, lwd = 6, col = "red")
#Add x axis
axis(1, at = 1:10, labels = F, line = 0, lwd.ticks = 3.5, tck = -0.02)
#Add x axis labels
text(x = 1:10, y = -log10(8000), labels = labs, srt = 45, xpd = T, adj = 1, cex = 4)
#Add y axis
axis(2, at = c(-log10(c(1000,100,10)),log10(c(1,10,100,1000))), las = 1, cex.axis = 4, hadj = 1, labels = c(-1000,-100,-10,0,10,100,1000))

## ADD MTEXT FOR LM
mtext(text = expression(LMM[GCTA]^u), side = 3, outer = F, line =  14, cex = 5)
#ADD MTEXT PANNEL
mtext(text = "D", side = 3, at = -2, cex = 7)
#ADD LOG AXIS TEXT
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()

#################################################################
#### FIGURE S20 ID EST: Intermediate frequencies causal loci ####

load("./ID_EST/EAS_1000_genomes_LMM/EAS_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")
load("./ID_EST/AFR_1000_genomes_LMM/AFR_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")
load("./ID_EST/POLYPED_20cM_LMM/POLYPED_ID3_simus_intermediateCausalLoci_COMPLETE_ID_EST.RData")

pdf("./PLOTS//FIGS20_IDest_ALLpops_IntermediateCausalLoci_ALL_scenario.pdf", height = 36, width = 54)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

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
mtext(text = expression("b (" ~ log[10] ~ " scale)"), side = 2, at = 0, cex = 4, line = 10)

dev.off()
