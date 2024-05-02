#1 Load packages  -----------------------------------------------------
HCBioSIM_SOIL<-c("rpart", "rpart.plot", "Cubist", "caret", "stringi", "plyr", "tidyverse")
lapply(HCBioSIM_SOIL, require, character.only=TRUE)

#2 Getting fingerprints from tsv (Chemotyper)----------------------------------------------------
d1<-read.csv(file.choose(), header = T, sep ="\t")
d1<-as_tibble(d1)

#3 Trimming toxprints ------------------------------------------------------

#Trimming for Sediment 
Fings<-d1 %>% 
  select ("M_SMILES","bond.C.N_nitrile_generic",	"bond.C.N_nitrile",	"bond.CN_amine_aliphatic_generic",	"bond.CS_sulfide_dialkyl",	
          "bond.CS_sulfide",	"chain.alkaneBranch_isopropyl_C3",	"chain.alkaneBranch_isohexyl_pentyl_3.methyl",	
          "chain.alkaneBranch_isooctyl_heptyl_3.methyl",	"chain.alkaneBranch_isononyl_heptyl_2_5.methyl",	"chain.alkaneCyclic_hexyl_C6",
          "chain.alkaneLinear_ethyl_C2_.connect_noZ_CN.4.",	"chain.alkaneLinear_propyl_C3",	"chain.alkaneLinear_butyl_C4",	"chain.alkaneLinear_hexyl_C6",
          "chain.alkaneLinear_octyl_C8",	"chain.alkaneLinear_decyl_C10",	"chain.alkaneLinear_dodedyl_C12",	"chain.alkaneLinear_tetradecyl_C14",	
          "chain.alkaneLinear_hexadecyl_C16",	"chain.alkaneLinear_stearyl_C18",	"chain.alkeneCyclic_diene_cyclopentadiene",	
          "chain.alkeneCyclic_ethene_C_.connect_noZ.",	"chain.alkeneLinear_mono.ene_ehtylene_terminal",	"chain.alkeneLinear_mono.ene_ethylene_generic",	
          "chain.aromaticAlkane_Ar.C_meta",	"chain.aromaticAlkane_Ar.C_ortho",	"chain.aromaticAlkane_Ph.C1_acyclic_generic",	
          "chain.aromaticAlkane_Ph.1_4.C1_acyclic",	"chain.aromaticAlkane_Ph.C2",	"chain.aromaticAlkane_Ph.C1_cyclic",	"chain.aromaticAlkene_Ph.C2_cyclic",	
          "ring.aromatic_benzene",	"ring.aromatic_biphenyl", "ring.aromatic_phenyl",	"ring.fused_.5_6._indane",	ring.fused_.6_6._naphthalene,	ring.fused_PAH_acenaphthylene,	
          "ring.fused_PAH_anthracene",	"ring.fused_PAH_benz.a.anthracene",	"ring.fused_PAH_benzophenanthrene",	
          "ring.fused_PAH_fluorene",	"ring.fused_PAH_phenanthrene",	"ring.fused_PAH_pyrene",	"ring.hetero_.5._N_pyrrole",
          "ring.hetero_.5._N_pyrrole_generic",	"ring.hetero_.5._N_S_thiazole",	"ring.hetero_.5._O_furan",	"ring.hetero_.5._O_oxolane",
          "ring.hetero_.5._S_thiophene",	"ring.hetero_.5._Z_1_3.Z",	"ring.hetero_.5._Z_1.Z",	"ring.hetero_.5_6._N_indole",
          "ring.hetero_.5_6._N_S_benzothiazole_.1_3..",	"ring.hetero_.5_6._O_benzofuran",	"ring.hetero_.5_6._Z_generic",
          "ring.hetero_.6._N_pyrazine",	"ring.hetero_.6._N_pyridine",	"ring.hetero_.6._N_pyridine_generic",	"ring.hetero_.6._Z_1.",
          "ring.hetero_.6._Z_1_4.",	"ring.hetero_.6._Z_generic",	"ring.hetero_.6_5_6._N_carbazole",	"ring.hetero_.6_5_6._O_benzofuran_dibenzo",
          "ring.hetero_.6_6._N_isoquinoline",	"ring.hetero_.6_6._N_quinoline",	"ring.hetero_.6_6._N_quinoxaline",	"ring.hetero_.6_6._Z_generic",
          "ring.hetero_.6_6_6._N_acridine")

#Trimming for Soil
Fings2<-d1 %>% 
  select (M_SMILES,chain.alkaneBranch_isopropyl_C3, chain.alkaneBranch_t.butyl_C4,	chain.alkaneBranch_neopentyl_C5,	
          chain.alkaneBranch_isohexyl_pentyl_3.methyl,	chain.alkaneBranch_isooctyl_heptyl_3.methyl,	
          chain.alkaneBranch_isononyl_heptyl_2_5.methyl,	chain.alkaneCyclic_pentyl_C5,	chain.alkaneCyclic_hexyl_C6,	
          chain.alkaneLinear_ethyl_C2_.connect_noZ_CN.4.,	chain.alkaneLinear_propyl_C3,	chain.alkaneLinear_butyl_C4,	chain.alkaneLinear_hexyl_C6,
          chain.alkaneLinear_octyl_C8,	chain.alkaneLinear_decyl_C10,	chain.alkaneLinear_dodedyl_C12,	chain.alkaneLinear_tetradecyl_C14,	
          chain.alkaneLinear_hexadecyl_C16,	chain.alkaneLinear_stearyl_C18,	chain.alkeneCyclic_diene_cyclohexene,	chain.alkeneCyclic_diene_cyclopentadiene,	
          chain.alkeneCyclic_ethene_C_.connect_noZ.,	chain.aromaticAlkane_Ar.C_meta,	chain.aromaticAlkane_Ar.C_ortho,	
          chain.aromaticAlkane_Ph.C1_acyclic_connect_H_gt_1,	chain.aromaticAlkane_Ph.C1_acyclic_connect_noDblBd,	
          chain.aromaticAlkane_Ph.C1_acyclic_generic,	chain.aromaticAlkane_Ph.1_4.C1_acyclic,	chain.aromaticAlkane_Ph.C2,	
          chain.aromaticAlkane_Ph.C1_cyclic,	chain.aromaticAlkene_Ph.C2_cyclic, ring.aromatic_benzene,	ring.aromatic_biphenyl,	
          ring.aromatic_phenyl,	ring.fused_.6_6._naphthalene,	ring.fused_.6_6._tetralin,	ring.fused_PAH_acenaphthylene,	
          ring.fused_PAH_anthracene,	ring.fused_PAH_benz.a.anthracene,	ring.fused_PAH_benzophenanthrene,	ring.fused_PAH_fluorene,	
          ring.fused_PAH_phenanthrene,	ring.fused_PAH_pyrene,	ring.hetero_.5._O_furan,	ring.hetero_.5._O_oxolane,	ring.hetero_.5._S_thiophene,
          ring.hetero_.5._Z_1.Z,	ring.hetero_.5_6._O_benzofuran,	ring.hetero_.5_6._Z_generic,	ring.hetero_.6_5_6._O_benzofuran_dibenzo)


# Trimming for freshwater/marine 
Fings3<-d1 %>%
  select(M_SMILES, chain.alkaneBranch_isopropyl_C3,	chain.alkaneBranch_t.butyl_C4,	chain.alkaneBranch_neopentyl_C5,	
         chain.alkaneBranch_isohexyl_pentyl_3.methyl,	chain.alkaneBranch_isooctyl_heptyl_3.methyl,	
         chain.alkaneBranch_isooctyl_hexyl_2.ethyl,	chain.alkaneBranch_isooctyl_hexyl_2.methyl,	chain.alkaneBranch_isononyl_heptyl_2_5.methyl,	
         chain.alkaneBranch_isononyl_pentyl_1_1_1_3.metyl,	chain.alkaneBranch_isodecyl_octyl_1_2.methyl,	chain.alkaneCyclic_pentyl_C5,	
         chain.alkaneCyclic_hexyl_C6,	chain.alkaneLinear_ethyl_C2_.connect_noZ_CN.4.,	chain.alkaneLinear_propyl_C3,	
         chain.alkaneLinear_butyl_C4,	chain.alkaneLinear_hexyl_C6,	chain.alkaneLinear_octyl_C8,	chain.alkaneLinear_decyl_C10,	
         chain.alkaneLinear_dodedyl_C12,	chain.alkaneLinear_tetradecyl_C14,	chain.alkaneLinear_hexadecyl_C16,	chain.alkaneLinear_stearyl_C18,	
         chain.aromaticAlkane_Ar.C_meta,	chain.aromaticAlkane_Ar.C_ortho,	chain.aromaticAlkane_Ph.C1_acyclic_generic,	
         chain.aromaticAlkane_Ph.1_4.C1_acyclic,	chain.aromaticAlkane_Ph.C2,	chain.aromaticAlkane_Ph.C4,	chain.aromaticAlkane_Ph.C6,	
         chain.aromaticAlkane_Ph.C8,	chain.aromaticAlkane_Ph.C9_nonylphenyl,	chain.aromaticAlkane_Ph.C10,	chain.aromaticAlkane_Ph.C1_cyclic,	
         ring.aromatic_benzene,	ring.aromatic_biphenyl,	ring.fused_.5_6._indane,	ring.fused_.6_6._naphthalene,	ring.fused_.6_6._tetralin,	
         ring.fused_PAH_acenaphthylene,	ring.fused_PAH_anthracene,	ring.fused_PAH_benz.a.anthracene,	ring.fused_PAH_benzophenanthrene,	
         ring.fused_PAH_fluorene,	ring.fused_PAH_phenanthrene,	ring.fused_PAH_pyrene)

#These default values have been taken from the QMRF FOR Hydrocarbon Biodegradation System-Integrated Model (HC-BioSIM) in Davis et al 2022-2024
#Sediment parameters 
newSed<-Fings %>% 
  mutate(StudyT = 4, logfoc = 0.551450, logWSratio = 0.361728, log_Tdose = 2.30103,log_Dose = 0.69897, 
         Disp = "No", Innoculum = "Marine") %>% 
  rename_at('M_SMILES', ~'SMILES') %>% rowid_to_column(("Data.ID.No"))

#Soil parameters 
newSoil<-Fings2 %>% 
  mutate(StudyT = 20, logfoc = 1.0413926865, logDose = 2, logTdose = 4.301029996, Dispersant = "No") %>% 
  rename_at('M_SMILES', ~'SMILES')%>% rowid_to_column(("Data.ID.No"))

# SW parameters 
Fings3xSW<- Fings3 %>% 
  mutate(Viscosity = 1.023, C_disp = 0, Concentration = 2.5, Temperature = 5, FW = 0, SW = 1, AS = 0) %>%
  rename_at('M_SMILES', ~'Smiles') %>% rowid_to_column(("ChemID"))

# FW parameters 
Fings3xFW<- Fings3 %>% 
  mutate(Viscosity = 1.023, C_disp = 0, Concentration = 0.1, Temperature = 21, FW = 1, SW = 0, AS = 0) %>%
  rename_at('M_SMILES', ~'Smiles') %>% rowid_to_column(("ChemID"))

#4 Sediment module -----------------------------------------------------
sedimentD<-read.csv(file.choose()) #add training dataset 

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(sedimentD$Data.ID.No)){
  if (sedimentD[i,38]==1) {sedimentD[i,30:47]=0}
  if (sedimentD[i,37]==1) {sedimentD[i,30:46]=0}
  if (sedimentD[i,36]==1) {sedimentD[i,30:35]=0}
  if (sedimentD[i,35]==1) {sedimentD[i,30:34]=0}
  if (sedimentD[i,34]==1) {sedimentD[i,30:33]=0}
  if (sedimentD[i,33]==1) {sedimentD[i,30:32]=0}
  if (sedimentD[i,32]==1) {sedimentD[i,30:31]=0}
  if (sedimentD[i,31]==1) {sedimentD[i,30]=0}
}

#Setting same random seed as the pp-LFER models to ensure comparison
set.seed(34353638)
smp_sz = floor(0.80*length(sedimentD$Data.ID.No))
train_indSed = sample(x = seq(1,length(sedimentD$Data.ID.No),1), size = smp_sz)

trainSed = sedimentD[train_indSed,]

#Incorporating the chemical fingerprints and system parameters together
m3_8Sed = cubist(x = trainSed[,c(13, 17,21, seq(24, 82,1))], y=trainSed$logDT50) #Temp + Dose + WS

#arrange data 
newSed<-rbind.fill(sedimentD, newSed) #this merges table with default values for FW 
newSed<-newSed[1234:nrow(newSed),]

#Make new prediction 
#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(newSed$Data.ID.No)){
  if (newSed[i,38]==1) {newSed[i,30:47]=0}
  if (newSed[i,37]==1) {newSed[i,30:46]=0}
  if (newSed[i,36]==1) {newSed[i,30:35]=0}
  if (newSed[i,35]==1) {newSed[i,30:34]=0}
  if (newSed[i,34]==1) {newSed[i,30:33]=0}
  if (newSed[i,33]==1) {newSed[i,30:32]=0}
  if (newSed[i,32]==1) {newSed[i,30:31]=0}
  if (newSed[i,31]==1) {newSed[i,30]=0}
}

SedResults<-data.frame(LogDT50_Sediment = signif(predict(m3_8Sed, newSed), 7))  

#5 Soil module -------------------------------------------------------------
soilD <-read.csv(file.choose()) #add training dataset 

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(soilD$Data.ID.No)){
  if (soilD[i,31]==1) {soilD[i,23:30]=0}
  if (soilD[i,30]==1) {soilD[i,23:29]=0}
  if (soilD[i,29]==1) {soilD[i,23:28]=0}
  if (soilD[i,28]==1) {soilD[i,23:27]=0}
  if (soilD[i,27]==1) {soilD[i,23:26]=0}
  if (soilD[i,26]==1) {soilD[i,23:25]=0}
  if (soilD[i,25]==1) {soilD[i,23:24]=0}
  if (soilD[i,24]==1) {soilD[i,23]=0}
}

#Trimming higher order DAHs and PAHs
for (i in 1:length(soilD$Data.ID.No)){
  if (sum(soilD[i,44:53])>=1) {soilD[i,42]=0}
}

#Setting same random seed as the pp-LFER models to ensure comparison
set.seed(34353638)
smp_sz = floor(0.80*length(soilD$Data.ID.No))
train_indSoil = sample(x = seq(1,length(soilD$Data.ID.No),1), size = smp_sz)

trainSoil = soilD[train_indSoil,]

#Incorporating the chemical fingerprints and system parameters together
m3_7Soil = cubist(x = trainSoil[,c(9, 10, 13, seq(14, 60,1))], y=trainSoil$logDT50)

#To make a new prediction 
#arrange data 
newSoil<-rbind.fill(soilD, newSoil) #this merges table with default values for FW 
newSoil<-newSoil[838:nrow(newSoil),]

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(newSoil$Data.ID.No)){
  if (newSoil[i,31]==1) {newSoil[i,23:30]=0}
  if (newSoil[i,30]==1) {newSoil[i,23:29]=0}
  if (newSoil[i,29]==1) {newSoil[i,23:28]=0}
  if (newSoil[i,28]==1) {newSoil[i,23:27]=0}
  if (newSoil[i,27]==1) {newSoil[i,23:26]=0}
  if (newSoil[i,26]==1) {newSoil[i,23:25]=0}
  if (newSoil[i,25]==1) {newSoil[i,23:24]=0}
  if (newSoil[i,24]==1) {newSoil[i,23]=0}
}

SoilResults<-data.frame(LogDT50_Soil = signif(predict(m3_7Soil, newSoil), 7))  

#6 Fresh Water module ----------------------------------------------
#Load trainset SW'n'FW

Waterdc<-read.csv(file.choose()) #add training dataset 

#Set identical seed to previous models and split dataset for training and validation
set.seed(34353637)
smp_sz2c = floor(0.80*length(Waterdc$Chemical))
train_ind2c = sample(x = seq(1,length(Waterdc$ChemID),1), size = smp_sz2c)

#Trim redundant and conditional ToxPrint values, as described in supplemental information
#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(Waterdc$ChemID)){
  if(Waterdc[i,58]==1){Waterdc[i,50:57]=0}
  if(Waterdc[i,57]==1){Waterdc[i,50:56]=0}
  if(Waterdc[i,56]==1){Waterdc[i,50:55]=0}
  if(Waterdc[i,55]==1){Waterdc[i,50:54]=0}
  if(Waterdc[i,54]==1){Waterdc[i,50:53]=0}
  if(Waterdc[i,53]==1){Waterdc[i,50:52]=0}
  if(Waterdc[i,52]==1){Waterdc[i,50:51]=0}
  if(Waterdc[i,51]==1){Waterdc[i,50]=0}
}

#Trimming redundant n-Ph fragments C2 - C10
for (i in 1:length(Waterdc$ChemID)){
  if(Waterdc[i,68]==1){Waterdc[i,63:67]=0}
  if(Waterdc[i,67]==1){Waterdc[i,63:66]=0}
  if(Waterdc[i,66]==1){Waterdc[i,63:65]=0}
  if(Waterdc[i,65]==1){Waterdc[i,63:64]=0}
  if(Waterdc[i,64]==1){Waterdc[i,63]=0}
}

#Trimming redundant higher aromatics from benzenes (DAH+)
for (i in 1:length(Waterdc$ChemID)){
  if(sum(Waterdc[i,71:81])>=1){Waterdc[i,70]=0}
}

#Trimming redundant higher aromatics from bi/terphenyls (DAH+)
for (i in 1:length(Waterdc$ChemID)){
  if(sum(Waterdc[i,72:81])>=1){Waterdc[i,71]=0}
}

#Define training and test sets
train2c = Waterdc[train_ind2c,]

#Run cubist() model
m3d2 = cubist(x = train2c[,c(13,14,15,16,19,20,21, seq(37, length(names(train2c))))], y=train2c$log_HL_obs)

###To make a new prediction

newFW<-rbind.fill(Waterdc, Fings3xFW) #this merges table with default values for FW 
newFW<-newFW[729:nrow(newFW),]

#Trim redundant and conditional ToxPrint values, as described in supplemental information

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(newFW$ChemID)){
  if(newFW[i,58]==1){newFW[i,50:57]=0}
  if(newFW[i,57]==1){newFW[i,50:56]=0}
  if(newFW[i,56]==1){newFW[i,50:55]=0}
  if(newFW[i,55]==1){newFW[i,50:54]=0}
  if(newFW[i,54]==1){newFW[i,50:53]=0}
  if(newFW[i,53]==1){newFW[i,50:52]=0}
  if(newFW[i,52]==1){newFW[i,50:51]=0}
  if(newFW[i,51]==1){newFW[i,50]=0}
}

#Trimming redundant n-Ph fragments C2 - C10
for (i in 1:length(newFW$ChemID)){
  if(newFW[i,68]==1){newFW[i,63:67]=0}
  if(newFW[i,67]==1){newFW[i,63:66]=0}
  if(newFW[i,66]==1){newFW[i,63:65]=0}
  if(newFW[i,65]==1){newFW[i,63:64]=0}
  if(newFW[i,64]==1){newFW[i,63]=0}
}

#Trimming redundant higher aromatics from benzenes (DAH+)
for (i in 1:length(newFW$ChemID)){
  if(sum(newFW[i,71:81])>=1){newFW[i,70]=0}
}

#Trimming redundant higher aromatics from bi/terphenyls (DAH+)
for (i in 1:length(newFW$ChemID)){
  if(sum(newFW[i,72:81])>=1){newFW[i,71]=0}
}

FWResults<-data.frame(LogDT50_Freshwater = signif(predict(m3d2, newFW), 7)) #Predicted log(DT50). 

#7 Seawater module ---------------------------------------------------------
###To make a new prediction

newSW<-rbind.fill(Waterdc, Fings3xSW)
newSW<-newSW[729:nrow(newSW),]

#Trim redundant and conditional ToxPrint values, as described in supplemental information
#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(newSW$ChemID)){
  if(newSW[i,58]==1){newSW[i,50:57]=0}
  if(newSW[i,57]==1){newSW[i,50:56]=0}
  if(newSW[i,56]==1){newSW[i,50:55]=0}
  if(newSW[i,55]==1){newSW[i,50:54]=0}
  if(newSW[i,54]==1){newSW[i,50:53]=0}
  if(newSW[i,53]==1){newSW[i,50:52]=0}
  if(newSW[i,52]==1){newSW[i,50:51]=0}
  if(newSW[i,51]==1){newSW[i,50]=0}
}

#Trimming redundant n-Ph fragments C2 - C10
for (i in 1:length(newSW$ChemID)){
  if(newSW[i,68]==1){newSW[i,63:67]=0}
  if(newSW[i,67]==1){newSW[i,63:66]=0}
  if(newSW[i,66]==1){newSW[i,63:65]=0}
  if(newSW[i,65]==1){newSW[i,63:64]=0}
  if(newSW[i,64]==1){newSW[i,63]=0}
}

#Trimming redundant higher aromatics from benzenes (DAH+)
for (i in 1:length(newSW$ChemID)){
  if(sum(newSW[i,71:81])>=1){newSW[i,70]=0}
}

#Trimming redundant higher aromatics from bi/terphenyls (DAH+)
for (i in 1:length(newSW$ChemID)){
  if(sum(newSW[i,72:81])>=1){newSW[i,71]=0}
}

SWResults<-data.frame(LogDT50_Seawater = signif(predict(m3d2, newSW), 7))
# add function to turn data into DT50 before getting the .csv file 

SWResults<-SWResults %>%
  mutate(DT50_Seawater = 10^(LogDT50_Seawater))

FWResults<-FWResults%>%
  mutate(DT50_Freshwater = 10^(LogDT50_Freshwater))

SoilResults<-SoilResults %>% 
  mutate(DT50_Soil= 10^(LogDT50_Soil))

SedResults<-SedResults %>% 
  mutate(DT50_Sediment= 10^(LogDT50_Sediment))


write.csv(All_DT50<-data.frame(FWResults, SWResults, SoilResults, SedResults), file = "Predictions_DT50.csv")
