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

d = read.csv(file="Soil_Final_Dataset_11022021.csv")

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

#Setting same random seed as the pp-LFER models to ensure comparison
set.seed(34353638)
smp_sz = floor(0.80*length(d$Data.ID.No))
train_ind = sample(x = seq(1,length(d$Data.ID.No),1), size = smp_sz)

train = d[train_ind,]
test = d[-train_ind,]

#Incorporating the chemical fingerprints and system parameters together
m3 = cubist(x = train[,c(9, 10, 11, 13, seq(14, 62,1))], y=train$logDT50)
P3 = predict(m3, test)
all_P3 = predict(m3, d)
train_P3 = predict(m3, train)

#Evaluating the model performance
(cor(P3, test$logDT50))^2
RMSE(actual = test$logDT50, predicted = P3)

class = unique(d$Chemical.Class)
aromatics = class[grep(pattern = "A|a", x = class)]
nonarom = class[grep(pattern = "A|a", x = class, invert = TRUE)]
classv2 = c(as.character(nonarom[c(4,1,3)]), as.character(aromatics[c(1,2,3,6,4,5)]), as.character(nonarom[c(2)]))
arom_col = c(brewer.pal(n = (length(classv2)-1), name = "RdBu"), "darkmagenta")

T_miss = as.numeric(is.na(d$StudyT)==TRUE)
D_miss = as.numeric(is.na(d$log_Dose)==TRUE)
OC_miss = as.numeric(is.na(d$logfoc)==TRUE)
Rat_miss = as.numeric(is.na(d$logWSratio)==TRUE)

all = all_P3
noD = all[D_miss==0]
noT = all[T_miss==0 & D_miss==0]
noOC = all[T_miss==0 & D_miss==0 & OC_miss==0]
noWSR = all[T_miss==0 & D_miss==0 & OC_miss==0 & Rat_miss==0]

#Predictive plots by source (similar to other plots for pp-LFER model)
#setting up plots
#lims = c(floor(min(d$log_HL_obs, d$log_HL_bhcw, log10(e$Observed.HL..d.))), ceiling(max(d$log_HL_obs,d$log_HL_bhcw, log10(e$Observed.HL..d.))))
lims = c(-3,5)
axs = c(1:9 %o% 10^(lims[[1]]:lims[[2]]))[c(1:9 %o% 10^(lims[[1]]:lims[[2]])) <= (10^lims[[2]])]
axl = c(rep(c(1, rep(0,8)),8),1)
axl1 = gsub(pattern = "\\<0\\>", replacement = "", x = as.character(c(rep(c(1, rep(0,8)),(lims[[2]]-lims[[1]])),1)*axs))
#refs=unique(d$Ref)
#cols = rainbow(n = length(refs))

x11(width = 7, height = 6)
par(mfrow = c(1, 1))
par(cex = 0.8)
par(mar = c(2, 2, 2, 2), oma = c(5, 5, 2, 1))
 
plot(-1, type="n", xlim=10^(lims), ylim=10^(lims), xaxs="i", yaxs="i", axes=F, log="xy", xlab="", ylab="")
axis(side=1, at = axs, labels = rep("", length(axs)), tick = TRUE, font = 2, cex.axis=0.6, xpd=NA)
axis(side=2, at = axs, labels = rep("", length(axs)), tick = TRUE, font = 2, cex.axis=0.6, xpd=NA, las=2)
text(x = axs, y = rep(10^(lims[[1]]-0.45),length(axs)), labels = axl1, cex=1.4, xpd=NA, font=2)
text(y = axs, x = rep(10^(lims[[1]]-0.45),length(axs)), labels = axl1, cex=1.4, xpd=NA, font=2)
abline(h=axs, col=alpha(colour = "grey20", alpha = .10))
abline(v=axs, col=alpha(colour = "grey20", alpha = .10))
abline(a = 0, b=1, lty=1, lwd=2)
abline(a = 0.3, b=1, lty=6, lwd=2)
abline(a = -0.3, b=1, lty=6, lwd=2)
abline(a = 1, b=1, lty=6, lwd=2)
abline(a = -1, b=1, lty=6, lwd=2)
#legend(x=0.1, y=1000, legend = refs, pch = rep(21,length(refs)), pt.bg = c("grey40", cols), cex = 0.8, bty="n", xpd=NA, ncol=3)
#points(x=d$DT50, y=(10^(all)), pch=21, col="grey70", bg="grey70", xpd=NA, cex=1.2)
for (i in 1:length(classv2)){
  points(x=d$DT50[d$Chemical.Class==as.character(classv2[[i]])], y=(10^(all_P3))[d$Chemical.Class==as.character(classv2[[i]])], pch=21, col="black", bg=arom_col[[i]], xpd=NA, cex=1.6)
}

#points(x=d$DT50[D_miss==0], y=(10^(noD)), pch=21, col="red", bg="red", xpd=NA, cex=1.2)
#points(x=d$DT50[T_miss==0 & D_miss==0], y=(10^(noT)), pch=21, col="goldenrod", bg="goldenrod", xpd=NA, cex=1.1)
#points(x=d$DT50[T_miss==0 & D_miss==0 & OC_miss==0], y=(10^(noOC)), pch=21, col=alpha(colour = "darkgreen", alpha=0.75), bg=alpha(colour = "darkgreen", alpha=0.75), xpd=NA, cex=1.2)
#points(x=d$DT50[T_miss==0 & D_miss==0 & OC_miss==0 & Rat_miss==0], y=(10^(noWSR)), pch=21, col="darkblue", bg="darkblue", xpd=NA, cex=1.2)
legend("topleft", legend = classv2, pch=rep(21,length(classv2)), pt.bg = arom_col, col=arom_col, bty="n", ncol = 2, pt.cex = 1.6, cex=1.4)

RMSE_train = RMSE(actual = train$logDT50, train_P3)
R2_train= cor(train$logDT50,train_P3)^2
RMSE_test = RMSE(actual = test$logDT50, P3)
R2_test= cor(test$logDT50,P3)^2
#mtext(side = 1, at = 10000, line = -5, text = str_c("RMSE (train) = ", signif(RMSE_train,3)), font=2, cex=0.8, xpd=NA)
#mtext(side = 1, at = 10000, line = -4, text = str_c("R-sq (train) = ", signif(R2_train,3)), font=2, cex=0.8, xpd=NA)
#mtext(side = 1, at = 5000, line = -3, text = str_c("RMSE (test) = ", signif(RMSE_test,3)), font=2, cex=1.2, xpd=NA)
#mtext(side = 1, at = 5000, line = -2, text = str_c("R-sq (test) = ", signif(R2_test,3)), font=2, cex=1.2, xpd=NA)
mtext(side=1, outer=TRUE, text = "Observed log(DT50) [days]", line = 1, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=2, outer=TRUE, text = "Predicted log(DT50) [days]", line = 2.5, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=3, outer=TRUE, text = "HC-BioSIM - Soil ", line = 0.0, at = 0.5, font = 2, xpd=NA, cex=1.4)

dev.copy2pdf(file="Soil_DT50_Predictions_ALL_v3.pdf")
write.csv(all_P3, file="Soil_DT50_Predictions_ALL.csv")
write.csv(test$logDT50, file="Soil_pred_test_DT50.csv")
write.csv(P3, file="Soil_exp_test_DT50.csv")

x11(width = 7, height = 6)
par(mfrow = c(1, 1))
par(cex = 0.8)
par(mar = c(2, 2, 2, 2), oma = c(5, 5, 2, 1))

plot(-1, type="n", xlim=10^(lims), ylim=10^(lims), xaxs="i", yaxs="i", axes=F, log="xy", xlab="", ylab="")
axis(side=1, at = axs, labels = rep("", length(axs)), tick = TRUE, font = 2, cex.axis=0.6, xpd=NA)
axis(side=2, at = axs, labels = rep("", length(axs)), tick = TRUE, font = 2, cex.axis=0.6, xpd=NA, las=2)
text(x = axs, y = rep(10^(lims[[1]]-0.45),length(axs)), labels = axl1, cex=1.4, xpd=NA, font=2)
text(y = axs, x = rep(10^(lims[[1]]-0.45),length(axs)), labels = axl1, cex=1.4, xpd=NA, font=2)
abline(h=axs, col=alpha(colour = "grey20", alpha = .10))
abline(v=axs, col=alpha(colour = "grey20", alpha = .10))
abline(a = 0, b=1, lty=1, lwd=2)
abline(a = 0.3, b=1, lty=6, lwd=2)
abline(a = -0.3, b=1, lty=6, lwd=2)
abline(a = 1, b=1, lty=6, lwd=2)
abline(a = -1, b=1, lty=6, lwd=2)
#legend(x=0.1, y=1000, legend = refs, pch = rep(21,length(refs)), pt.bg = c("grey40", cols), cex = 0.8, bty="n", xpd=NA, ncol=3)
#points(x=d$DT50, y=(10^(all)), pch=21, col="grey70", bg="grey70", xpd=NA, cex=1.2)
for (i in 1:length(classv2)){
  points(x=d$DT50[d$Chemical.Class==as.character(classv2[[i]])], y=(d$HL_bhcw[d$Chemical.Class==as.character(classv2[[i]])]), pch=21, col="black", bg=arom_col[[i]], xpd=NA, cex=1.6)
}

#points(x=d$DT50[D_miss==0], y=(10^(noD)), pch=21, col="red", bg="red", xpd=NA, cex=1.2)
#points(x=d$DT50[T_miss==0 & D_miss==0], y=(10^(noT)), pch=21, col="goldenrod", bg="goldenrod", xpd=NA, cex=1.1)
#points(x=d$DT50[T_miss==0 & D_miss==0 & OC_miss==0], y=(10^(noOC)), pch=21, col=alpha(colour = "darkgreen", alpha=0.75), bg=alpha(colour = "darkgreen", alpha=0.75), xpd=NA, cex=1.2)
#points(x=d$DT50[T_miss==0 & D_miss==0 & OC_miss==0 & Rat_miss==0], y=(10^(noWSR)), pch=21, col="darkblue", bg="darkblue", xpd=NA, cex=1.2)
legend("topleft", legend = classv2, pch=rep(21,length(classv2)), pt.bg = arom_col, col=arom_col, bty="n", ncol = 2, pt.cex = 1.6, cex=1.4)

RMSE_train_old = RMSE(actual = train$logDT50[train$logHL_biohcw!=-99], train$logHL_biohcw[train$logHL_biohcw!=-99])
R2_train_old = cor(train$logDT50[train$logHL_biohcw!=-99],train$logHL_biohcw[train$logHL_biohcw!=-99])^2
RMSE_test_old = RMSE(actual = test$logDT50[test$logHL_biohcw!=-99], test$logHL_biohcw[test$logHL_biohcw!=-99])
R2_test_old = cor(test$logDT50[test$logHL_biohcw!=-99],test$logHL_biohcw[test$logHL_biohcw!=-99])^2
#mtext(side = 1, at = 10000, line = -5, text = str_c("RMSE (train) = ", signif(RMSE_train,3)), font=2, cex=0.8, xpd=NA)
#mtext(side = 1, at = 10000, line = -4, text = str_c("R-sq (train) = ", signif(R2_train,3)), font=2, cex=0.8, xpd=NA)
#mtext(side = 1, at = 5000, line = -3, text = str_c("RMSE (test) = ", signif(RMSE_test,3)), font=2, cex=1.2, xpd=NA)
#mtext(side = 1, at = 5000, line = -2, text = str_c("R-sq (test) = ", signif(R2_test,3)), font=2, cex=1.2, xpd=NA)
mtext(side=1, outer=TRUE, text = "Observed log(DT50) [days]", line = 1, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=2, outer=TRUE, text = "Predicted log(DT50) [days]", line = 2.5, at = 0.5, font = 2, xpd=NA, cex=1.4)
mtext(side=3, outer=TRUE, text = "BioHCWin - Soil (1:1 IMEF)", line = 0.0, at = 0.5, font = 2, xpd=NA, cex=1.4)


dev.copy2pdf(file="BIOHCWIN_Soil_DT50_Predictions_ALL_v3.pdf")

#Residual Boxplots as a Function of System Descriptors
x11(width = 16, height = 15)
par(mfrow = c(2, 2))
par(cex = 0.8)
par(mar = c(3, 3, 3, 3), oma = c(4, 3, 2, 1))

boxplot((all_P3 - d$logDT50)~d$StudyT, ylim=c(-2,2), xlab="Test Temperate [oC]", ylab="Logarithmic Residuals [-]", las=1, xpd=NA, varwidth=TRUE)
abline(h=1, lty=3, lwd=1)
abline(h=-1, lty=3, lwd=1)
abline(h=0.3, lty=6, lwd=1)
abline(h=-0.3, lty=6, lwd=1)
abline(h=0, lty=1, lwd=1)

boxplot((all_P3 - d$logDT50)~d$logDose, ylim=c(-2,2), xlab="LOG HC Dose [mg/kg]", ylab="Logarithmic Residuals [-]", las=1, xpd=NA, varwidth=TRUE)
abline(h=1, lty=3, lwd=1)
abline(h=-1, lty=3, lwd=1)
abline(h=0.3, lty=6, lwd=1)
abline(h=-0.3, lty=6, lwd=1)
abline(h=0, lty=1, lwd=1)

boxplot((all_P3 - d$logDT50)~d$logfoc, ylim=c(-2,2), xlab="LOG foc [-]", ylab="Logarithmic Residuals [-]", las=1, xpd=NA, varwidth=TRUE)
abline(h=1, lty=3, lwd=1)
abline(h=-1, lty=3, lwd=1)
abline(h=0.3, lty=6, lwd=1)
abline(h=-0.3, lty=6, lwd=1)
abline(h=0, lty=1, lwd=1)

mtext(text = "Test Temperature [C]", side=1, at=0.25, line=-30, outer=TRUE, xpd=NA, font=2)
mtext(text = "Substance Loading [log(mg/kg)]", side=1, at=0.75, line=-30, outer=TRUE, xpd=NA, font=2)
mtext(text = "log Foc [-]", side=1, at=0.25, line=0, outer=TRUE, xpd=NA, font=2)
mtext(text = "Logarithmic Residual [-]", side=2, at=0.75, line=-34, outer=TRUE, xpd=NA, font=2)
mtext(text = "Logarithmic Residual [-]", side=2, at=0.25, line=0, outer=TRUE, xpd=NA, font=2)
mtext(text = "Logarithmic Residual [-]", side=2, at=0.75, line=0, outer=TRUE, xpd=NA, font=2)
dev.copy2pdf(file="Soil_System_Residuals_v1.pdf")

#Residual Boxplots as a Function of HC Class and C#
x11(width = 12, height = 6)
par(mfrow = c(1, 2))
par(cex = 0.8)
par(mar = c(3, 3, 3, 3), oma = c(4, 3, 2, 1))


boxplot((all_P3 - d$logDT50)~d$Carbon.Number, ylim=c(-2,2), xlab="", ylab="", las=1, xpd=NA, varwidth=TRUE, cex.axis=0.8)
abline(h=1, lty=3, lwd=1)
abline(h=-1, lty=3, lwd=1)
abline(h=0.3, lty=6, lwd=1)
abline(h=-0.3, lty=6, lwd=1)
abline(h=0, lty=1, lwd=1)
#abline(lm((v$log.HL_VEGA. - v$log.HL_obs.)~v$Cn), col="red")
#text(x=seq(1,21,1), y=rep(-2.4, 21), labels = seq(4,24,1), cex=0.8, xpd=NA)


boxplot((all_P3 - d$logDT50)~d$Chemical.Class, ylim=c(-2,2), xlab="", ylab="", las=1, xpd=NA, varwidth=TRUE, cex.axis=0.8)
abline(h=1, lty=3, lwd=1)
abline(h=-1, lty=3, lwd=1)
abline(h=0.3, lty=6, lwd=1)
abline(h=-0.3, lty=6, lwd=1)
abline(h=0, lty=1, lwd=1)

#text(x=seq(1,11,1), y=rep(-2.4, 11), labels = c("nP", "iP", "MN", "DN", "PN", "mAr", "NMAH", "DAH", "NDAH", "PAH", "NPAH"), cex=0.7, xpd=NA)
mtext(text = "Carbon Number [-]", side=1, at=0.25, line=0, outer=TRUE, xpd=NA, font=2)
mtext(text = "Hydrocarbon Class [-]", side=1, at=0.75, line=0, outer=TRUE, xpd=NA, font=2)
mtext(text = "Logarithmic Residual [-]", side=2, at=0.5, line=0, outer=TRUE, xpd=NA, font=2)
mtext(text = "Logarithmic Residual [-]", side=2, at=0.5, line=-36, outer=TRUE, xpd=NA, font=2)


dev.copy2pdf(file="Soil_CnClass_Residuals_v1.pdf")
