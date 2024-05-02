library(rpart)
library(rpart.plot)
library(Cubist)
library(caret)
library(stringi)
library(stringr)
library(scales)
library(RColorBrewer)
library(viridis)

RMSE <- function(actual, predicted){
  sqrt((sum((actual - predicted)^2))/(length(actual)))
}

d = read.csv(file="Sediment_Final_Dataset_06072023_n1233.csv")

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(d$Data.ID.No)){
  if (d[i,38]==1) {d[i,30:47]=0}
  if (d[i,37]==1) {d[i,30:46]=0}
  if (d[i,36]==1) {d[i,30:35]=0}
  if (d[i,35]==1) {d[i,30:34]=0}
  if (d[i,34]==1) {d[i,30:33]=0}
  if (d[i,33]==1) {d[i,30:32]=0}
  if (d[i,32]==1) {d[i,30:31]=0}
  if (d[i,31]==1) {d[i,30]=0}
}

#Setting fixed random seed to ensure reproducibility
set.seed(34353638)
smp_sz = floor(0.80*length(d$Data.ID.No))
train_ind = sample(x = seq(1,length(d$Data.ID.No),1), size = smp_sz)

train = d[train_ind,]
test = d[-train_ind,]

#Optimizing Cubist() Model

#Subset 00 - No system parameters
m3_0 = cubist(x = train[,c(seq(24, 82,1))], y=train$logDT50)

#Subsets 01 - 04 - single system parameter
m3_1 = cubist(x = train[,c(13, seq(24, 82,1))], y=train$logDT50) #Temp
m3_2 = cubist(x = train[,c(15, seq(24, 82,1))], y=train$logDT50) #foc
m3_3 = cubist(x = train[,c(17, seq(24, 82,1))], y=train$logDT50) # W:S ratio
m3_4 = cubist(x = train[,c(21, seq(24, 82,1))], y=train$logDT50) #Dose

#Subsets 05 - 7 - T + next best parameter(s)
m3_5 = cubist(x = train[,c(13, 17, seq(24, 82,1))], y=train$logDT50) #Temp + W:S
m3_6 = cubist(x = train[,c(13, 21, seq(24, 82,1))], y=train$logDT50) #Temp + Dose
m3_7 = cubist(x = train[,c(17, 21, seq(24, 82,1))], y=train$logDT50) #Dose + W:S

#Subset 08 - T + Dose + WS
m3_8 = cubist(x = train[,c(13, 17,21, seq(24, 82,1))], y=train$logDT50) #Temp + Dose + WS

#orignal Full Model
m3_9 = cubist(x = train[,c(13, 15, 17, 21, seq(24, 82,1))], y=train$logDT50)

#Evaluating the model(s) performance

mat = matrix(nrow = 10, ncol = 6, data = "")

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
mat[8, 1] = #number of rules
  mat[8, 2] = length(m3_7$vars$used) # number of variables used
mat[8, 3] = (cor(P3_7, train$logDT50))^2
mat[8, 4] =RMSE(actual = train$logDT50, predicted = P3_7)
mat[8, 5] = (cor(P3_7v, test$logDT50))^2
mat[8, 6] =RMSE(actual = test$logDT50, predicted = P3_7v)

P3_8 = predict(m3_8, train)
P3_8v = predict(m3_8, test)
mat[9, 1] = #number of rules
  mat[9, 2] = length(m3_8$vars$used) # number of variables used
mat[9, 3] = (cor(P3_8, train$logDT50))^2
mat[9, 4] =RMSE(actual = train$logDT50, predicted = P3_8)
mat[9, 5] = (cor(P3_8v, test$logDT50))^2
mat[9, 6] =RMSE(actual = test$logDT50, predicted = P3_8v)

P3_9 = predict(m3_9, train)
P3_9v = predict(m3_9, test)
mat[10, 1] = #number of rules
  mat[10, 2] = length(m3_9$vars$used) # number of variables used
mat[10, 3] = (cor(P3_9, train$logDT50))^2
mat[10, 4] =RMSE(actual = train$logDT50, predicted = P3_9)
mat[10, 5] = (cor(P3_9v, test$logDT50))^2
mat[10, 6] =RMSE(actual = test$logDT50, predicted = P3_9v)

#SUmmary performance for all iterations of tested system parameter combination(s)
write.csv(x = mat, file = "Sediment_Model_Performance_Table_02212024.csv")

###
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
  
  m3 = cubist(x = tr[,c(13, 17,21, seq(24, 82,1))], y=tr$logDT50) #Temp + Dose + WS
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

#Bootstrap analysis output table
write.csv(sum_sfdt, file = "Sediment_CV_Analysis_02212024.csv")


#Model performance visualization plot(s)
class = unique(d$Chemical.Class)
aromatics = class[grep(pattern = "A|a", x = class)]
nonarom = class[grep(pattern = "A|a", x = class, invert = TRUE)]
classv2 = c(as.character(nonarom[c(3,1)]), as.character(nonarom[c(2,4,5)]), "DN", as.character(aromatics[c(1,6,2,5,4,3)]))
arom_col = c("darkgreen", "darkmagenta", brewer.pal(n = (length(classv2)-2), name = "RdBu"))

lims = c(-2,6)
axs = c(1:9 %o% 10^(lims[[1]]:lims[[2]]))[c(1:9 %o% 10^(lims[[1]]:lims[[2]])) <= (10^lims[[2]])]
axl = c(rep(c(1, rep(0,8)),8),1)
axl1 = gsub(pattern = "\\<0\\>", replacement = "", x = as.character(c(rep(c(1, rep(0,8)),(lims[[2]]-lims[[1]])),1)*axs))

x11(width = 7, height = 6)
par(mfrow = c(1, 1))
par(cex = 0.8)
par(mar = c(2, 2, 2, 2), oma = c(5, 5, 2, 1))

#Model predictions for complete data set (input file rows - ALL)
all_P3_8 = predict(m3_8, d)

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
  points(x=d$DT50[d$Chemical.Class==as.character(classv2[[i]])], y=(10^(all_P3_8))[d$Chemical.Class==as.character(classv2[[i]])], pch=21, col="black", bg=arom_col[[i]], xpd=NA, cex=1.6)
}

legend("topleft", legend = classv2, pch=rep(21,length(classv2)), col = rep("black", 11), pt.bg=arom_col, bty="n", ncol = 2, pt.cex = 1.6, cex=1.4)

mtext(side=1, outer=TRUE, text = "Observed log(DT50) [days]", line = 1, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=2, outer=TRUE, text = "Predicted log(DT50) [days]", line = 2.5, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=3, outer=TRUE, text = "HC-BioSIM - Sediment ", line = 0.0, at = 0.5, font = 2, xpd=NA, cex=1.4)

#Plot Output
dev.copy2pdf(file="Sed_DT50_Predictions_ALL_02212024.pdf")

#Predicted logDT50 values - rows align with rows of input file
write.csv(all_P3_8, file="Sed_DT50_Predictions_ALL_02212024.csv")

#End of BioSIM plots and calculations

####
####
#For making predictions from a new input file / dataset:
e = read.csv(file="NAME")

#Trimming redundant n-alkyl fragments C3 - C18
for (i in 1:length(e$Data.ID.No)){
  if (e[i,38]==1) {e[i,30:47]=0}
  if (e[i,37]==1) {e[i,30:46]=0}
  if (e[i,36]==1) {e[i,30:35]=0}
  if (e[i,35]==1) {e[i,30:34]=0}
  if (e[i,34]==1) {e[i,30:33]=0}
  if (e[i,33]==1) {e[i,30:32]=0}
  if (e[i,32]==1) {e[i,30:31]=0}
  if (e[i,31]==1) {e[i,30]=0}
}

#Making predictions
predvalues = predict(m3_8, e)
write.csv(predvalues, file="NAME_OUTPUT.csv")
