#Model Development, Prediction, and Output

#Loading R packages
library(scales)
library(stringr)
library(rpart)
library(rpart.plot)
library(Cubist)
library(caret)
library(RColorBrewer)

#Reading in Data, parameters, and system information file (tab "R_Read_File_Input" in Supplemental Information File)
d = read.csv(file="R_File_For_Input_v2.csv")


#Setting same random seed for all models to ensure identical training and validation set splits
set.seed(34353637)
smp_sz2a = floor(0.80*length(d$Chemical))
train_ind2a = sample(x = seq(1,length(d$ChemID),1), size = smp_sz2a)

train2 = d[train_ind2a,]
test2 = d[-train_ind2a,]

#Stepwise MLR for BioHCWin and pp-LFER models
bhcw_stepa = step(object = lm((train2$log_HL_obs-train2$log_HL_bhcw)~-1), scope = formula(lm((train2$log_HL_obs-train2$log_HL_bhcw)~train2$TempDiff+train2$Viscosity+train2$FW+train2$SW+train2$AS+train2$Concentration+train2$C_disp,train2)), direction = "forward")
qc_step2c = step(object = lm(train2$log_HL_obs~train2$HLEG-1), scope = formula(lm(train2$log_HL_obs~train2$E_qc+train2$S_qc+train2$B_qc+train2$V_qc+train2$HLEG+train2$Viscosity+train2$TempDiff+(train2$FW+train2$SW)+train2$AS+train2$Concentration+train2$C_disp,train2)), direction = "forward")

#Predicted Values for SA-BioHCWin and pp-LFER Models
bhcw_preda = d$log_HL_bhcw + d$SW*bhcw_stepa$coefficients[[1]] + d$AS*bhcw_stepa$coefficients[[2]] + d$FW*bhcw_stepa$coefficients[[3]] + d$C_disp*bhcw_stepa$coefficients[[4]] + d$TempDiff*bhcw_stepa$coefficients[[5]] + d$Viscosity*bhcw_stepa$coefficients[[6]]
qc_pred2c = d$HLEG*qc_step2c$coefficients[[1]] + d$E_qc*qc_step2c$coefficients[[2]] + d$C_disp*qc_step2c$coefficients[[3]] + d$AS*qc_step2c$coefficients[[4]] + d$TempDiff*qc_step2c$coefficients[[5]] + d$B_qc*qc_step2c$coefficients[[6]] + d$S_qc*qc_step2c$coefficients[[7]] + d$V_qc*qc_step2c$coefficients[[8]] + d$Concentration*qc_step2c$coefficients[[9]] + d$Viscosity*qc_step2c$coefficients[[10]]

#Write out model predictions to CSV files
write.csv(bhcw_preda, file="BHCWSA_Predictions_041921021.csv")
write.csv(qc_pred2c, file="Abraham_QC_Predictions_041921021.csv")

#Re-read in input file for HC-BIOSIM model
dc = read.csv(file = "R_File_For_Input_v2.csv")

#Set identical seed to previous models and split dataset for training and validation
set.seed(34353637)
smp_sz2c = floor(0.80*length(dc$Chemical))
train_ind2c = sample(x = seq(1,length(dc$ChemID),1), size = smp_sz2c)

#Trim redundant and conditional ToxPrint values, as described in supplemental information

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(dc$ChemID)){
  if(dc[i,58]==1){dc[i,50:57]=0}
  if(dc[i,57]==1){dc[i,50:56]=0}
  if(dc[i,56]==1){dc[i,50:55]=0}
  if(dc[i,55]==1){dc[i,50:54]=0}
  if(dc[i,54]==1){dc[i,50:53]=0}
  if(dc[i,53]==1){dc[i,50:52]=0}
  if(dc[i,52]==1){dc[i,50:51]=0}
  if(dc[i,51]==1){dc[i,50]=0}
}

#Trimming redundant n-Ph fragments C2 - C10
for (i in 1:length(dc$ChemID)){
  if(dc[i,68]==1){dc[i,63:67]=0}
  if(dc[i,67]==1){dc[i,63:66]=0}
  if(dc[i,66]==1){dc[i,63:65]=0}
  if(dc[i,65]==1){dc[i,63:64]=0}
  if(dc[i,64]==1){dc[i,63]=0}
}

#Trimming redundant higher aromatics from benzenes (DAH+)
for (i in 1:length(dc$ChemID)){
  if(sum(dc[i,71:81])>=1){dc[i,70]=0}
}

#Trimming redundant higher aromatics from bi/terphenyls (DAH+)
for (i in 1:length(dc$ChemID)){
  if(sum(dc[i,72:81])>=1){dc[i,71]=0}
}

#Define training and test sets
train2c = dc[train_ind2c,]
test2c = dc[-train_ind2c,]

#Run cubist() model
m3d2 = cubist(x = train2c[,c(13,14,15,16,19,20,21, seq(37, length(names(train2c))))], y=train2c$log_HL_obs)

#Predict DT50s for training, test, and complete datasets, write to CSV
P3train = predict(m3d2, train2c)
P3test = predict(m3d2, test2c)
P3all = predict(m3d2, dc)
write.csv(P3all, file="HC-BIOSIM_Predictions_FINAL_04272021.csv")

###To make a new prediction
new = read.csv(file="Predict_HCBioSim.csv")

#Trim redundant and conditional ToxPrint values, as described in supplemental information

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(new$ChemID)){
  if(new[i,58]==1){new[i,50:57]=0}
  if(new[i,57]==1){new[i,50:56]=0}
  if(new[i,56]==1){new[i,50:55]=0}
  if(new[i,55]==1){new[i,50:54]=0}
  if(new[i,54]==1){new[i,50:53]=0}
  if(new[i,53]==1){new[i,50:52]=0}
  if(new[i,52]==1){new[i,50:51]=0}
  if(new[i,51]==1){new[i,50]=0}
}

#Trimming redundant n-Ph fragments C2 - C10
for (i in 1:length(new$ChemID)){
  if(new[i,68]==1){new[i,63:67]=0}
  if(new[i,67]==1){new[i,63:66]=0}
  if(new[i,66]==1){new[i,63:65]=0}
  if(new[i,65]==1){new[i,63:64]=0}
  if(new[i,64]==1){new[i,63]=0}
}

#Trimming redundant higher aromatics from benzenes (DAH+)
for (i in 1:length(new$ChemID)){
  if(sum(new[i,71:81])>=1){new[i,70]=0}
}

#Trimming redundant higher aromatics from bi/terphenyls (DAH+)
for (i in 1:length(new$ChemID)){
  if(sum(new[i,72:81])>=1){new[i,71]=0}
}


newpredict = signif(predict(m3d2, new), 7) #Predicted log(DT50). 
write.csv(newpredict, file = "HCBioSim_predicted.csv")


head(new)

