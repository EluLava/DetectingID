library(xlsx)

setwd("~/PhD/Humans_Inbreeding_Est/id_est_feb_2022/data/")

#Load the data
load("./ID_EST/EAS_1000_genomes_LMM/EAS_ID3_simus_COMPLETE_ID_EST.RData")
load("./ID_EST/AFR_1000_genomes_LMM/AFR_ID3_simus_COMPLETE_ID_EST.RData")
load("./ID_EST/WORLD_1000_genomes_LMM/WORLD_ID3_simus_COMPLETE_ID_EST.RData")
load("./ID_EST/POLYPED_20cM_LMM/POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/WORLD_SUB_RANDOM_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")
load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/WORLD_SUB_RANGED_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")

#########################################################################################
############ TABLE 1 RMSE both POLYPED (full + both subsampled) four models #############
#########################################################################################

#Create empty dataframe
dataRMSE = as.data.frame(matrix(ncol = 10, nrow = 12))
#Set colnames
colnames(dataRMSE) = c("Model", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB")

#Set model
dataRMSE$Model = rep(c("LM", "LMM: AS", "LMM: GCTA WE", "LMM: GCTA UN"), each = 3)
#Set population
dataRMSE$Population = c("PEDIGREE (complete)","PEDIGREE (random sub-sampling)","PEDIGREE (ranged sub-sampling)")

## FULL PEDIGREE

#LM

#get the dataframe
ID_est = ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT
#Subset only the F we are interested in
ID_est = ID_est[,c(5,6,8:13)]

#Loop through F
for (Finb in 1:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (complete)" & dataRMSE$Model == "LM",(Finb + 2)] = RMSE
  
}

#LMM: AS

#get the dataframe
ID_est = ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT
#Subset only the F we are interested in
ID_est = ID_est[,c(14,15,17:22)]

#Loop through F
for (Finb in 1:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (complete)" & dataRMSE$Model == "LMM: AS",(Finb + 2)] = RMSE
  
}

#LMM: GCTA WE

#get the dataframe
ID_est = ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT
#Subset only the F we are interested in
ID_est = ID_est[,c(23,24,26:31)]

#Loop through F
for (Finb in 1:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (complete)" & dataRMSE$Model == "LMM: GCTA WE",(Finb + 2)] = RMSE
  
}

#LMM: GCTA UN

#get the dataframe
ID_est = ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT
#Subset only the F we are interested in
ID_est = ID_est[,c(32,33,35:40)]

#Loop through F
for (Finb in 1:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (complete)" & dataRMSE$Model == "LMM: GCTA UN",(Finb + 2)] = RMSE
  
}

####################################################################################

# RANDOM SUBSAMPLING PEDIGREE

load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/WORLD_SUB_RANDOM_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")

#LM

ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,6,7,9:14)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (random sub-sampling)" & dataRMSE$Model == "LM",(Finb + 1)] = RMSE
  
}

#LMM: AS

#get the data frame
ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,15,16,18:23)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (random sub-sampling)" & dataRMSE$Model == "LMM: AS",(Finb + 1)] = RMSE
  
}

#LMM: GCTA WE

#get the data frame
ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,24,25,27:32)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (random sub-sampling)" & dataRMSE$Model == "LMM: GCTA WE",(Finb + 1)] = RMSE
  
}

#LMM: GCTA UN

#get the data frame
ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,33,34,36:41)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (random sub-sampling)" & dataRMSE$Model == "LMM: GCTA UN",(Finb + 1)] = RMSE
  
}

####################################################################################

# RANGED SUBSAMPLING PEDIGREE

load("./ID_EST/POLYPED_SUBSMPLING_20cM_LMM/WORLD_SUB_RANGED_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")

#LM

ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,6,7,9:14)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (ranged sub-sampling)" & dataRMSE$Model == "LM",(Finb + 1)] = RMSE
  
}

#LMM: AS

#get the data frame
ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,15,16,18:23)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (ranged sub-sampling)" & dataRMSE$Model == "LMM: AS",(Finb + 1)] = RMSE
  
}

#LMM: GCTA WE

#get the data frame
ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,24,25,27:32)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (ranged sub-sampling)" & dataRMSE$Model == "LMM: GCTA WE",(Finb + 1)] = RMSE
  
}

#LMM: GCTA UN

#get the data frame
ID_est = WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST[WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST$Simuscen == "PGPED_LMM_ID3_selhT_selpT_demaT",2:ncol(WORLD_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST)]
#Subset only the F we are interested in
ID_est = ID_est[,c(1,33,34,36:41)]

#Loop through F
for (Finb in 2:ncol(ID_est)) {
  
  #Extract values
  IDVal = ID_est[,Finb]
  #true value if ID
  trueVal = -3
  
  RMSE = sqrt((mean((IDVal - trueVal)^2)))
  
  #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
  dataRMSE[dataRMSE$Population == "PEDIGREE (ranged sub-sampling)" & dataRMSE$Model == "LMM: GCTA UN",(Finb + 1)] = RMSE
  
}

#Pass into correct notation
for(col in 3:ncol(dataRMSE)){
  #pass into numeric
  dataRMSE[,col] = as.numeric(as.character(dataRMSE[,col]))
  dataRMSE[,col] = formatC(dataRMSE[,col], format = "f", digits = 2)
}

#Save the Table
write.csv(dataRMSE, file = "./Tables/TABLE1.csv", quote = F, row.names = F)

################################################################
############ TABLE 2 RMSE ALL 1KG Pops four models #############
################################################################

#Create empty dataframe
dataRMSE = as.data.frame(matrix(ncol = 9, nrow = 12))
#Set colnames
colnames(dataRMSE) = c("Model", "Population", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB")

#Set model
dataRMSE$Model = rep(c("LM", "LMM: AS", "LMM: GCTA WE", "LMM: GCTA UN"), each = 3)
#Set population
dataRMSE$Population = rep(c("EAS", "AFR", "WORLD"),4)

#LM

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  scenario = "selhT_selpT_demaT"
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  
  ## LM
  
  #Loop through F
  for (Finb in 5:11) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE[dataRMSE$Model == "LM" & dataRMSE$Population == pop,(Finb - 2)] = RMSE
    
  }
  
  ## LMM: AS
  
  #Loop through F
  for (Finb in 15:21) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE[dataRMSE$Model == "LMM: AS" & dataRMSE$Population == pop,(Finb - 12)] = RMSE
    
  }
  
  ## LMM: GCTA WE
  
  #Loop through F
  for (Finb in 25:31) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE[dataRMSE$Model == "LMM: GCTA WE" & dataRMSE$Population == pop,(Finb - 22)] = RMSE
    
  }
  
  ## LMM: GCTA UN
  
  #Loop through F
  for (Finb in 35:41) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2, na.rm = T)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE[dataRMSE$Model == "LMM: GCTA UN" & dataRMSE$Population == pop,(Finb - 32)] = RMSE
    
  }
}

#Pass into correct notation
for(col in 3:ncol(dataRMSE)){
  #pass into numeric
  dataRMSE[,col] = as.numeric(as.character(dataRMSE[,col]))
  dataRMSE[,col] = formatC(dataRMSE[,col], format = "f", digits = 2)
}

#Save the Table
write.csv(dataRMSE, file = "./Tables/TABLE2.csv", quote = F, row.names = F)

###################################################
############ TABLE S1 Samples included ############
###################################################

#Read the F data and keep only samples info
dta = read.table("./F_datasets/WORLD_allchr_allFs.txt", header = T)[,1:3]
colnames(dta) = c("sample","population","super_population")

#Add one column per manuscript population
dta = cbind(dta, EAS_included=vector(mode = "logical", length = nrow(dta)), AFR_included=vector(mode = "logical", length = nrow(dta)), WORLD_included=vector(mode = "logical", length = nrow(dta)))

#add TRUE when samples were incl/uded in the pop
dta$EAS_included[dta$super_population == "EAS"] = TRUE
dta$AFR_included[dta$super_population == "AFR"] = TRUE
dta$WORLD_included = TRUE

write.csv(dta, file = "./Tables/TableS1_1000_Genomes_Samples_description.csv", quote = F, row.names = F)


############################################################################################
############ TABLE S3 RMSE ALL Fs; ALL scenarios; LINEAR MODEL; ALL POPULATIONS ############
############################################################################################

#Create empty dataframe
dataRMSE.LM = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataRMSE.LM) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataRMSE.LM$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataRMSE.LM$Population = rep(c("PGPED","EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
    
    #Loop through scenario
    for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 5:14) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      #true value if ID
      trueVal = -3
      
      
      RMSE = sqrt((mean((IDVal - trueVal)^2)))
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataRMSE.LM[dataRMSE.LM$Scenario == scenario & dataRMSE.LM$Population == pop,(Finb - 1)] = RMSE

    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(5,6,8:13)]
  
  #Loop through F
  #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3

    RMSE = sqrt((mean((IDVal - trueVal)^2)))

    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE.LM[dataRMSE.LM$Scenario == scenario & dataRMSE.LM$Population == pop,(Finb + 2)] = RMSE
    
  }
}

#Set Scenario CLEAN
dataRMSE.LM$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)

#Pass into scientific notation
for(col in 3:ncol(dataRMSE.LM)){
  #pass into numeric
  dataRMSE.LM[,col] = as.numeric(as.character(dataRMSE.LM[,col]))
  dataRMSE.LM[,col] = formatC(dataRMSE.LM[,col], format = "f", digits = 2)
}


#############################################################################################################
############ TABLE S4 RMSE ALL Fs; ALL scenarios; LINEAR MIXED MODEL, GRM = AS ; ALL POPULATIONS ############
#############################################################################################################

#Create empty dataframe
dataRMSE.LMM.AS = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataRMSE.LMM.AS) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataRMSE.LMM.AS$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataRMSE.LMM.AS$Population = rep(c( "PGPED","EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 15:24) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      #true value if ID
      trueVal = -3
      
      
      RMSE = sqrt((mean((IDVal - trueVal)^2)))
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataRMSE.LMM.AS[dataRMSE.LMM.AS$Scenario == scenario & dataRMSE.LMM.AS$Population == pop,(Finb - 11)] = RMSE
      
    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(5,6,8:13)]
  
  #Loop through F
  #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE.LMM.AS[dataRMSE.LMM.AS$Scenario == scenario & dataRMSE.LMM.AS$Population == pop,(Finb + 2)] = RMSE
    
  }
}

#Set Scenario CLEAN
dataRMSE.LMM.AS$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)

#Pass into scientific notation
for(col in 3:ncol(dataRMSE.LMM.AS)){
  #pass into numeric
  dataRMSE.LMM.AS[,col] = as.numeric(as.character(dataRMSE.LMM.AS[,col]))
  dataRMSE.LMM.AS[,col] = formatC(dataRMSE.LMM.AS[,col], format = "f", digits = 2)
}

##############################################################################################################
############ TABLE S5 RMSE ALL Fs; ALL scenarios; LINEAR MIXED MODEL, GRM = kc0 ; ALL POPULATIONS ############
##############################################################################################################

#Create empty dataframe
dataRMSE.LMM.kc0 = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataRMSE.LMM.kc0) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataRMSE.LMM.kc0$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataRMSE.LMM.kc0$Population = rep(c("PGPED", "EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 25:34) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      #true value if ID
      trueVal = -3
      
      
      RMSE = sqrt((mean((IDVal - trueVal)^2)))
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataRMSE.LMM.kc0[dataRMSE.LMM.kc0$Scenario == scenario & dataRMSE.LMM.kc0$Population == pop,(Finb - 21)] = RMSE
      
    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(5,6,8:13)]
  
  #Loop through F
  #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE.LMM.kc0[dataRMSE.LMM.kc0$Scenario == scenario & dataRMSE.LMM.kc0$Population == pop,(Finb + 2)] = RMSE
    
  }
}

#Set Scenario CLEAN
dataRMSE.LMM.kc0$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)

#Pass into scientific notation
for(col in 3:ncol(dataRMSE.LMM.kc0)){
  #pass into numeric
  dataRMSE.LMM.kc0[,col] = as.numeric(as.character(dataRMSE.LMM.kc0[,col]))
  dataRMSE.LMM.kc0[,col] = formatC(dataRMSE.LMM.kc0[,col], format = "f", digits = 2)
}

###############################################################################################################
############ TABLE S6 RMSE ALL Fs; ALL scenarios; LINEAR MIXED MODEL, GRM = GCTA ; ALL POPULATIONS ############
###############################################################################################################

#Create empty dataframe
dataRMSE.LMM.gcta = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataRMSE.LMM.gcta) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataRMSE.LMM.gcta$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataRMSE.LMM.gcta$Population = rep(c("PGPED", "EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 35:44) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      #true value if ID
      trueVal = -3
      
      
      RMSE = sqrt((mean((IDVal - trueVal)^2, na.rm = T)))
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataRMSE.LMM.gcta[dataRMSE.LMM.gcta$Scenario == scenario & dataRMSE.LMM.gcta$Population == pop,(Finb - 31)] = RMSE
      
    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(5,6,8:13)]
  
  #Loop through F
  #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    #true value if ID
    trueVal = -3
    
    RMSE = sqrt((mean((IDVal - trueVal)^2)))
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataRMSE.LMM.gcta[dataRMSE.LMM.gcta$Scenario == scenario & dataRMSE.LMM.gcta$Population == pop,(Finb + 2)] = RMSE
    
  }
}

#Set Scenario CLEAN
dataRMSE.LMM.gcta$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)

#Pass into scientific notation
for(col in 3:ncol(dataRMSE.LMM.gcta)){
  #pass into numeric
  dataRMSE.LMM.gcta[,col] = as.numeric(as.character(dataRMSE.LMM.gcta[,col]))
  dataRMSE.LMM.gcta[,col] = formatC(dataRMSE.LMM.gcta[,col], format = "f", digits = 2)
}

#######################
#write the excel dataset, one Table Sup per onglet
write.xlsx(dataRMSE.LM, file = "./Tables/TableS3-6_ID_EST_RMSE.xlsx", sheetName = "SIMPLE LINEAR MODEL", append = F, col.names = T, row.names = F)
write.xlsx(dataRMSE.LMM.AS, file = "./Tables/TableS3-6_ID_EST_RMSE.xlsx", sheetName = "LINEAR MIXED MODEL GRM AS", append = T, col.names = T, row.names = F)
write.xlsx(dataRMSE.LMM.kc0, file = "./Tables/TableS3-6_ID_EST_RMSE.xlsx", sheetName = "LINEAR MIXED MODEL GRM GCTA WE", append = T, col.names = T, row.names = F)
write.xlsx(dataRMSE.LMM.gcta, file = "./Tables/TableS3-6_ID_EST_RMSE.xlsx", sheetName = "LINEAR MIXED MODEL GRM GCTA UN", append = T, col.names = T, row.names = F)

###########################################################################################
############ TABLE S7 PERCENTAGES OF REPLICATES WHICH DID NOT CONVERGE LMM: AS ############
###########################################################################################

#Create empty dataframe
dataREPCONV.LMM.AS = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataREPCONV.LMM.AS) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataREPCONV.LMM.AS$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataREPCONV.LMM.AS$Population = rep(c("PGPED", "EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 15:24) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      
      #Get # of rep > 1000 OR < -1000 OR NA
      NB.rep.NOTconv = length(IDVal[IDVal > 1000 | IDVal < -1000 | is.na(IDVal)])
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataREPCONV.LMM.AS[dataREPCONV.LMM.AS$Scenario == scenario & dataREPCONV.LMM.AS$Population == pop,(Finb - 11)] = NB.rep.NOTconv
      
    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(14,15,17:22)]
  
  #Loop through F
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    
    #Get # of rep > 1000 OR < -1000 OR NA
    NB.rep.NOTconv = length(IDVal[IDVal > 1000 | IDVal < -1000 | is.na(IDVal)])
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataREPCONV.LMM.AS[dataREPCONV.LMM.AS$Scenario == scenario & dataREPCONV.LMM.AS$Population == pop,(Finb + 2)] = NB.rep.NOTconv
    
  }
}

#Set Scenario CLEAN
dataREPCONV.LMM.AS$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)

############################################################################################
############ TABLE S8 PERCENTAGES OF REPLICATES WHICH DID NOT CONVERGE LMM: kc0 ############
############################################################################################

#Create empty dataframe
dataREPCONV.LMM.kc0 = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataREPCONV.LMM.kc0) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataREPCONV.LMM.kc0$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataREPCONV.LMM.kc0$Population = rep(c("PGPED", "EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 25:34) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      
      #Get # of rep > 1000 OR < -1000 OR NA
      NB.rep.NOTconv = length(IDVal[IDVal > 1000 | IDVal < -1000 | is.na(IDVal)])
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataREPCONV.LMM.kc0[dataREPCONV.LMM.kc0$Scenario == scenario & dataREPCONV.LMM.kc0$Population == pop,(Finb - 21)] = NB.rep.NOTconv
      
    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(23,24,26:31)]
  
  #Loop through F
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    
    #Get # of rep > 1000 OR < -1000 OR NA
    NB.rep.NOTconv = length(IDVal[IDVal > 1000 | IDVal < -1000 | is.na(IDVal)])
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataREPCONV.LMM.kc0[dataREPCONV.LMM.kc0$Scenario == scenario & dataREPCONV.LMM.kc0$Population == pop,(Finb + 2)] = NB.rep.NOTconv
    
  }
}

#Set Scenario CLEAN
dataREPCONV.LMM.kc0$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)


#############################################################################################
############ TABLE S8 PERCENTAGES OF REPLICATES WHICH DID NOT CONVERGE LMM: gcta ############
#############################################################################################

#Create empty dataframe
dataREPCONV.LMM.gcta = as.data.frame(matrix(ncol = 13, nrow = 32))
#Set colnames
colnames(dataREPCONV.LMM.gcta) = c("Scenario", "Population", "Fped", "Fas", "Funi.UN", "Funi.WE", "FROHs.BCFTools.100KB", "FROHs.PLINK.100KB", "FROHs.BCFTools.1MB", "FROHs.PLINK.1MB", "Fas.LDMS", "Funi.UN.LDMS", "Funi.WE.LDMS")

#Set Scenario
dataREPCONV.LMM.gcta$Scenario = rep(c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT"),4)
#Set population
dataREPCONV.LMM.gcta$Population = rep(c("PGPED", "EAS", "AFR", "WORLD"), each = 8)

#Loop through 1 KG populations
for(pop in c("EAS", "AFR", "WORLD")){
  
  #Loop through scenario
  for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
    
    #get the values of ID
    ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
    
    #Loop through F
    #for (Finb in colnames(get(paste0("ID_EST_", pop, "_LMM_ID3_selhF_selpF_demaF")))[5:14]) {
    for (Finb in 35:44) {
      
      #Extract values
      IDVal = ID_est[,Finb]
      
      #Get # of rep > 1000 OR < -1000 OR NA
      NB.rep.NOTconv = length(IDVal[IDVal > 1000 | IDVal < -1000 | is.na(IDVal)])
      
      #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
      dataREPCONV.LMM.gcta[dataREPCONV.LMM.gcta$Scenario == scenario & dataREPCONV.LMM.gcta$Population == pop,(Finb - 31)] = NB.rep.NOTconv
      
    }
  }
}

#ADD PolyPed population
pop = "PGPED"

#Loop through scenario
for (scenario in c("selhF_selpF_demaF", "selhT_selpF_demaF", "selhF_selpT_demaF", "selhF_selpF_demaT", "selhT_selpT_demaF", "selhT_selpF_demaT", "selhF_selpT_demaT", "selhT_selpT_demaT")) {
  
  #get the values of ID
  ID_est = get(paste0("ID_EST_", pop, "_LMM_ID3_", scenario))
  #Subset only the F we are interested in
  ID_est = ID_est[,c(32,33,35:40)]
  
  #Loop through F
  for (Finb in 1:ncol(ID_est)) {
    
    #Extract values
    IDVal = ID_est[,Finb]
    
    #Get # of rep > 1000 OR < -1000 OR NA
    NB.rep.NOTconv = length(IDVal[IDVal > 1000 | IDVal < -1000 | is.na(IDVal)])
    
    #Fill the df FPED EMPTY BECASUE ONLY FOR POLYPED
    dataREPCONV.LMM.gcta[dataREPCONV.LMM.gcta$Scenario == scenario & dataREPCONV.LMM.gcta$Population == pop,(Finb + 2)] = NB.rep.NOTconv
    
  }
}

#Set Scenario CLEAN
dataREPCONV.LMM.gcta$Scenario = rep(c("Standard", "ADD", "DOM", "DEMA", "ADD & DOM", "ADD & DEMA", "DOM & DEMA", "ADD & DOM & DEMA"),4)

#######################
#write the excel dataset, one Table Sup per onglet
write.xlsx(dataREPCONV.LMM.AS, file = "./Tables/TableS7-9_ID_EST_RepConvergence.xlsx", sheetName = "LINEAR MIXED MODEL GRM AS", append = T, col.names = T, row.names = F)
write.xlsx(dataREPCONV.LMM.kc0, file = "./Tables/TableS7-9_ID_EST_RepConvergence.xlsx", sheetName = "LINEAR MIXED MODEL GRM GCTA WE", append = T, col.names = T, row.names = F)
write.xlsx(dataREPCONV.LMM.gcta, file = "./Tables/TableS7-9_ID_EST_RepConvergence.xlsx", sheetName = "LINEAR MIXED MODEL GRM GCTA UN", append = T, col.names = T, row.names = F)

