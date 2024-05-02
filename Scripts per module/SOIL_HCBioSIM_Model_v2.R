library(rpart)
library(rpart.plot)
library(Cubist)
library(caret)
library(stringi)
library(stringr)
library(scales)
library(RColorBrewer)

RMSE <- function(actual, predicted){
  sqrt((sum((actual - predicted)^2))/(length(actual)))
}

d = read.csv(file="Soil_Final_Dataset_06082023_n946_v3.csv")

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(d$Data.ID.No)){
  if (d[i,31]==1) {d[i,23:30]=0}
  if (d[i,30]==1) {d[i,23:29]=0}
  if (d[i,29]==1) {d[i,23:28]=0}
  if (d[i,28]==1) {d[i,23:27]=0}
  if (d[i,27]==1) {d[i,23:26]=0}
  if (d[i,26]==1) {d[i,23:25]=0}
  if (d[i,25]==1) {d[i,23:24]=0}
  if (d[i,24]==1) {d[i,23]=0}
}

#Trimming higher order DAHs and PAHs
for (i in 1:length(d$Data.ID.No)){
  if (sum(d[i,44:53])>=1) {d[i,42]=0}
}

#Setting same random seed as the pp-LFER models to ensure comparison
set.seed(34353638)
smp_sz = floor(0.80*length(d$Data.ID.No))
train_ind = sample(x = seq(1,length(d$Data.ID.No),1), size = smp_sz)

train = d[train_ind,]
test = d[-train_ind,]

#Incorporating the chemical fingerprints and system parameters together
m3 = cubist(x = train[,c(9, 10, 13, seq(14, 60,1))], y=train$logDT50)

#No system params
m3_0 = cubist(x = train[,c(seq(14, 60,1))], y=train$logDT50)

#Single System Parameter
m3_1 = cubist(x = train[,c(10, seq(14, 60,1))], y=train$logDT50) #Temp
m3_2 = cubist(x = train[,c(9, seq(14, 60,1))], y=train$logDT50) #foc
m3_3 = cubist(x = train[,c(13, seq(14, 60,1))], y=train$logDT50) #Dose

#2 system Parameters
m3_4 = cubist(x = train[,c(9, 13, seq(14, 60,1))], y=train$logDT50) #foc + Dose
m3_5 = cubist(x = train[,c(9, 10, seq(14, 60,1))], y=train$logDT50) #foc + Temp
m3_6 = cubist(x = train[,c(10, 13, seq(14, 60,1))], y=train$logDT50) #Dose + Temp

m3_7 = cubist(x = train[,c(9, 10, 13, seq(14, 60,1))], y=train$logDT50)

mat = matrix(nrow = 8, ncol = 6, data = "")

#Predictions on training Dataset
P3_0 = predict(m3_0, train)
P3_0v = predict(m3_0, test)
mat[1, 1] = #number of rules
  mat[1, 2] = length(m3_0$vars$used) # number of variables used
mat[1, 3] = (cor(P3_0, train$logDT50))^2
mat[1, 4] =RMSE(actual = train$logDT50, predicted = P3_0)
mat[1, 5] = (cor(P3_0v, test$logDT50))^2
mat[1, 6] =RMSE(actual = test$logDT50, predicted = P3_0v)

P3_1 = predict(m3_1, train)
P3_1v = predict(m3_1, test)
mat[2, 1] = #number of rules
  mat[2, 2] = length(m3_1$vars$used) # number of variables used
mat[2, 3] = (cor(P3_1, train$logDT50))^2
mat[2, 4] =RMSE(actual = train$logDT50, predicted = P3_1)
mat[2, 5] = (cor(P3_1v, test$logDT50))^2
mat[2, 6] =RMSE(actual = test$logDT50, predicted = P3_1v)

P3_2 = predict(m3_2, train)
P3_2v = predict(m3_2, test)
mat[3, 1] = #number of rules
  mat[3, 2] = length(m3_2$vars$used) # number of variables used
mat[3, 3] = (cor(P3_2, train$logDT50))^2
mat[3, 4] =RMSE(actual = train$logDT50, predicted = P3_2)
mat[3, 5] = (cor(P3_2v, test$logDT50))^2
mat[3, 6] =RMSE(actual = test$logDT50, predicted = P3_2v)

P3_3 = predict(m3_3, train)
P3_3v = predict(m3_3, test)
mat[4, 1] = #number of rules
  mat[4, 2] = length(m3_3$vars$used) # number of variables used
mat[4, 3] = (cor(P3_3, train$logDT50))^2
mat[4, 4] =RMSE(actual = train$logDT50, predicted = P3_3)
mat[4, 5] = (cor(P3_3v, test$logDT50))^2
mat[4, 6] =RMSE(actual = test$logDT50, predicted = P3_3v)

P3_4 = predict(m3_4, train)
P3_4v = predict(m3_4, test)
mat[5, 1] = #number of rules
  mat[5, 2] = length(m3_4$vars$used) # number of variables used
mat[5, 3] = (cor(P3_4, train$logDT50))^2
mat[5, 4] =RMSE(actual = train$logDT50, predicted = P3_4)
mat[5, 5] = (cor(P3_4v, test$logDT50))^2
mat[5, 6] =RMSE(actual = test$logDT50, predicted = P3_4v)

P3_5 = predict(m3_5, train)
P3_5v = predict(m3_5, test)
mat[6, 1] = #number of rules
  mat[6, 2] = length(m3_5$vars$used) # number of variables used
mat[6, 3] = (cor(P3_5, train$logDT50))^2
mat[6, 4] =RMSE(actual = train$logDT50, predicted = P3_5)
mat[6, 5] = (cor(P3_5v, test$logDT50))^2
mat[6, 6] =RMSE(actual = test$logDT50, predicted = P3_5v)

P3_6 = predict(m3_6, train)
P3_6v = predict(m3_6, test)
mat[7, 1] = #number of rules
  mat[7, 2] = length(m3_6$vars$used) # number of variables used
mat[7, 3] = (cor(P3_6, train$logDT50))^2
mat[7, 4] =RMSE(actual = train$logDT50, predicted = P3_6)
mat[7, 5] = (cor(P3_6v, test$logDT50))^2
mat[7, 6] =RMSE(actual = test$logDT50, predicted = P3_6v)

P3_7 = predict(m3_7, train)
P3_7v = predict(m3_7, test)
all_P3 = predict(m3_7, d)
mat[8, 1] = #number of rules
  mat[8, 2] = length(m3_7$vars$used) # number of variables used
mat[8, 3] = (cor(P3_7, train$logDT50))^2
mat[8, 4] =RMSE(actual = train$logDT50, predicted = P3_7)
mat[8, 5] = (cor(P3_7v, test$logDT50))^2
mat[8, 6] =RMSE(actual = test$logDT50, predicted = P3_7v)

write.csv(x = mat, file = "Soil_Model_Performance_Table_02212024.csv")

#Bootstrapping the final selected model
sum_sfdt = matrix(data = "", nrow = 5, ncol = 7)

set.seed(12345678)
f1_train = sample(d$Data.ID.No, size = floor(0.80*length(d$Data.ID.No)), replace = FALSE)
f1_test = sample(d$Data.ID.No, size = ceiling(0.20*length(d$Data.ID.No)), replace = FALSE)
f2_train = sample(d$Data.ID.No, size = floor(0.80*length(d$Data.ID.No)), replace = FALSE)
f2_test = sample(d$Data.ID.No, size = ceiling(0.20*length(d$Data.ID.No)), replace = FALSE)
f3_train = sample(d$Data.ID.No, size = floor(0.80*length(d$Data.ID.No)), replace = FALSE)
f3_test = sample(d$Data.ID.No, size = ceiling(0.20*length(d$Data.ID.No)), replace = FALSE)
f4_train = sample(d$Data.ID.No, size = floor(0.80*length(d$Data.ID.No)), replace = FALSE)
f4_test = sample(d$Data.ID.No, size = ceiling(0.20*length(d$Data.ID.No)), replace = FALSE)
f5_train = sample(d$Data.ID.No, size = floor(0.80*length(d$Data.ID.No)), replace = FALSE)
f5_test = sample(d$Data.ID.No, size = ceiling(0.20*length(d$Data.ID.No)), replace = FALSE)

tr_folds = matrix(ncol = floor(0.80*length(d$Data.ID.No)), nrow = 5, data = c(f1_train, f2_train, f3_train, f4_train, f5_train), byrow = TRUE)

for (i in 1:5){
  tr = data.frame(d[tr_folds[i,],])
  te = data.frame(d[-tr_folds[i,],])
  
  m3 = cubist(x = tr[,c(9, 10, 13, seq(14, 60,1))], y=tr$logDT50)
  summary(m3)
  pr_te = predict(m3, te)
  P3 = predict(m3, d)
  pr_tr = predict(m3, tr)
  
  
  sum_sfdt[i,1] = RMSE(actual = d$logDT50, predicted = P3)
  sum_sfdt[i,2] = cor(d$logDT50, P3)^2
  sum_sfdt[i,3] = RMSE(actual = tr$logDT50, predicted = pr_tr)
  sum_sfdt[i,4] = cor(tr$logDT50, pr_tr)^2
  sum_sfdt[i,5] = RMSE(actual = te$logDT50, predicted = pr_te)
  sum_sfdt[i,6] = cor(te$logDT50, pr_te)^2
  sum_sfdt[i,7] = length(te$Data.ID.No)
  print(i)
}

#Outputting the bootstrap performace table
write.csv(sum_sfdt, file = "Soil_CV_Analysis_02212024.csv")

#Model performance and visualization
class = unique(d$Chemical.Class)
aromatics = class[grep(pattern = "A|a", x = class)]
nonarom = class[grep(pattern = "A|a", x = class, invert = TRUE)]
classv2 = c(as.character(nonarom[c(1)]), as.character(nonarom[c(2,4,3,5)]), as.character(aromatics[c(4,5,2,6,1,3)]))
arom_col = c("darkmagenta", brewer.pal(n = (length(classv2)-1), name = "RdBu"))

#setting up plots
lims = c(-3,5)
axs = c(1:9 %o% 10^(lims[[1]]:lims[[2]]))[c(1:9 %o% 10^(lims[[1]]:lims[[2]])) <= (10^lims[[2]])]
axl = c(rep(c(1, rep(0,8)),8),1)
axl1 = gsub(pattern = "\\<0\\>", replacement = "", x = as.character(c(rep(c(1, rep(0,8)),(lims[[2]]-lims[[1]])),1)*axs))

x11(width = 7, height = 6)
par(mfrow = c(1, 1))
par(cex = 0.8)
par(mar = c(2, 2, 2, 2), oma = c(5, 5, 2, 1))

plot(-1, type="n", xlim=10^(lims), ylim=10^(lims), xaxs="i", yaxs="i", axes=F, log="xy", xlab="", ylab="")
axis(side=1, at = axs, labels = rep("", length(axs)), tick = TRUE, font = 2, cex.axis=1.0, xpd=NA)
axis(side=2, at = axs, labels = rep("", length(axs)), tick = TRUE, font = 2, cex.axis=1.0, xpd=NA, las=2)
text(x = axs, y = rep(10^(lims[[1]]-0.45),length(axs)), labels = axl1, cex=1.4, xpd=NA, font=2)
text(y = axs, x = rep(10^(lims[[1]]-0.5),length(axs)), labels = axl1, cex=1.4, xpd=NA, font=2)
abline(h=axs, col=alpha(colour = "grey20", alpha = .10))
abline(v=axs, col=alpha(colour = "grey20", alpha = .10))
abline(a = 0, b=1, lty=1, lwd=2)
abline(a = 0.3, b=1, lty=6, lwd=2)
abline(a = -0.3, b=1, lty=6, lwd=2)
abline(a = 1, b=1, lty=6, lwd=2)
abline(a = -1, b=1, lty=6, lwd=2)

for (i in 1:length(classv2)){
  points(x=d$DT50[d$Chemical.Class==as.character(classv2[[i]])], y=(10^(all_P3))[d$Chemical.Class==as.character(classv2[[i]])], pch=21, col="black", bg=arom_col[[i]], xpd=NA, cex=1.6)
}

legend("topleft", legend = classv2, pch=rep(21,length(classv2)), pt.bg = arom_col, col=rep("black", 10), bty="n", ncol = 2, pt.cex = 1.6, cex=1.4)

RMSE_train = RMSE(actual = train$logDT50, train_P3)
R2_train= cor(train$logDT50,train_P3)^2
RMSE_test = RMSE(actual = test$logDT50, P3)
R2_test= cor(test$logDT50,P3)^2

mtext(side=1, outer=TRUE, text = "Observed log(DT50) [days]", line = 1, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=2, outer=TRUE, text = "Predicted log(DT50) [days]", line = 2.5, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=3, outer=TRUE, text = "HC-BioSIM - Soil ", line = 0.0, at = 0.5, font = 2, xpd=NA, cex=1.4)

dev.copy2pdf(file="Soil_DT50_Predictions_02212024.pdf")

#Predicted logDT50 values - rows align with rows of input file
all_P3 = predict(m3_7, d)
write.csv(all_P3, file="Soil_DT50_Predictions_02212024.csv")
#End of BioSIM plots and calculations

####
####
#For making predictions from a new input file / dataset:
e = read.csv(file="NAME.csv")

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(e$Data.ID.No)){
  if (e[i,31]==1) {e[i,23:30]=0}
  if (e[i,30]==1) {e[i,23:29]=0}
  if (e[i,29]==1) {e[i,23:28]=0}
  if (e[i,28]==1) {e[i,23:27]=0}
  if (e[i,27]==1) {e[i,23:26]=0}
  if (e[i,26]==1) {e[i,23:25]=0}
  if (e[i,25]==1) {e[i,23:24]=0}
  if (e[i,24]==1) {e[i,23]=0}
}

#Making predictions
predvalues = predict(m3_7, e)
write.csv(predvalues, file="NAME_OUTPUT.csv")
