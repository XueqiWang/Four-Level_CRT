##################################################
# Collection of results from raw data
# Size and Power
#
# NOTE: This first chunk is not needed if using
# the intermediate results (i.e. the 11 XLSX files)
# currently available in the folder "binResults".
##################################################
library(openxlsx)

binRES<-function(row,case,ind,cv){
  path<-"binResults/binRData/"
  if(case==1){
    name<-"Size"
  } else{
    name<-"Power"
  }
  
  if(cv=="00"){
    cv<-""
    if(ind==1){
      namei<-"_ind"
    } else{
      namei<-""
    }
  } else{
    if(ind==1){
      namei<-"_ind_var_"
    } else{
      namei<-"_var_"
    }
  }
  
  load(paste0(path,row,"row/",name,"_beta",namei,cv,".RData"))
  
  # non-convergence
  err<-(results[,16]==1)
  nerror<-sum(err)
  errSIM<-results[err,1]
  
  # the correct results
  K<-as.numeric(results[1,14])
  p<-as.numeric(results[1,15])
  
  # remove non-convergence
  res1<-results[!err,6:13]
  mcsd1<-sd(res1[,1])
  essd1<-colMeans(res1[,2:8])
  pval1<-colMeans(apply(abs(res1[,1])/res1[,2:8],2,function(x){2*(1-pt(x,df=K-p))})<0.05)
  
  return(data.frame(MB=pval1[[1]], ROB=pval1[[2]], KC=pval1[[3]], MD=pval1[[4]], AVG=pval1[[7]], FG=pval1[[5]], MBN=pval1[[6]], nerror=nerror))
}

# Collect data (Tables 2-3, Web Tables 1-19)
tPower<-read.xlsx("binResults/pred_power_bin.xlsx", colNames=TRUE)

for (ind in 1:2){
  for (cv in c("00", 25, 50, 75, 100)){
    # Size
    sim_size<-NULL
    for (i in 1:30){
      sim_size<-rbind(sim_size,binRES(row=i,case=1,ind=ind,cv=cv))
    }
    
    # Power
    sim_power<-NULL
    for (i in 1:30){
      sim_power<-rbind(sim_power,binRES(row=i,case=2,ind=ind,cv=cv))
    }
    
    if(ind==1){
      namei<-"ind"
    } else{
      namei<-"MAEE"
    }
    
    final_bin<-cbind(tPower, sim_power, sim_size)
    write.xlsx(final_bin, file = paste("binResults/final_bin_",namei,"_",cv,".xlsx", sep=""), row.names = FALSE)
  }
}



##################################################
# Read data
# Size and Power
##################################################

for (ind in 1:2){
  for (cv in c("00", 25, 50, 75, 100)){
    if(ind==1){
      namei<-"_ind_"
    } else{
      namei<-"_MAEE_"
    }
    
    filename <- paste("binResults/final_bin", namei, cv, ".xlsx", sep="")
    data_maee <- read.xlsx(filename, colNames=TRUE)
    colnames(data_maee)[19:25]<-c("MB.1","ROB.1","KC.1","MD.1","AVG.1","FG.1","MBN.1")
    assign(paste("final_bin", namei, cv, sep=""), data_maee)
  }
}



##################################################
# Figures of results for BC1 (KC)
# Size and Power
##################################################

par(mfrow=c(1,2))

# Size (Figure 2 in the manuscript)
# MAEE
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE for BC1",xlim=c(1,30),ylim=c(0.01, 0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0.01, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0.01, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC.1,pch=16,cex = 1.5, col="orange")
points(as.numeric(row.names(final_bin_MAEE_25)),final_bin_MAEE_25$KC.1,pch=16,cex = 1.5,col="limegreen")
points(as.numeric(row.names(final_bin_MAEE_50)),final_bin_MAEE_50$KC.1,pch=16,cex = 1.5,col="mediumpurple2")
points(as.numeric(row.names(final_bin_MAEE_75)),final_bin_MAEE_75$KC.1,pch=16,cex = 1.5,col="lightskyblue")
points(as.numeric(row.names(final_bin_MAEE_100)),final_bin_MAEE_100$KC.1,pch=16,cex = 1.5,col="palevioletred1")

legend(x = "bottomright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25", "CV = 0.50", "CV = 0.75", "CV = 1.00"),
       col=c("orange", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
       cex = 0.8, pch = c(16,16,16,16,16))

# ind
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(b) IND for BC1",xlim=c(1,30),ylim=c(0.01, 0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0.01, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0.01, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC.1,pch=16,cex = 1.5, col="orange")
points(as.numeric(row.names(final_bin_ind_25)),final_bin_ind_25$KC.1,pch=16,cex = 1.5,col="limegreen")
points(as.numeric(row.names(final_bin_ind_50)),final_bin_ind_50$KC.1,pch=16,cex = 1.5,col="mediumpurple2")
points(as.numeric(row.names(final_bin_ind_75)),final_bin_ind_75$KC.1,pch=16,cex = 1.5,col="lightskyblue")
points(as.numeric(row.names(final_bin_ind_100)),final_bin_ind_100$KC.1,pch=16,cex = 1.5,col="palevioletred1")

legend(x = "bottomright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25", "CV = 0.50", "CV = 0.75", "CV = 1.00"),
       col=c("orange", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
       cex = 0.8, pch = c(16,16,16,16,16))


# Power (Figure 3 in the manuscript)
# MAEE
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE for BC1",xlim=c(1,30),ylim=c(-0.15, 0.05),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.15, 0.05, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.15, 0.05, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC-final_bin_MAEE_00$t_test,pch=16,cex = 1.5, col="orange")
points(as.numeric(row.names(final_bin_MAEE_25)),final_bin_MAEE_25$KC-final_bin_MAEE_25$t_test,pch=16,cex = 1.5,col="limegreen")
points(as.numeric(row.names(final_bin_MAEE_50)),final_bin_MAEE_50$KC-final_bin_MAEE_50$t_test,pch=16,cex = 1.5,col="mediumpurple2")
points(as.numeric(row.names(final_bin_MAEE_75)),final_bin_MAEE_75$KC-final_bin_MAEE_75$t_test,pch=16,cex = 1.5,col="lightskyblue")
points(as.numeric(row.names(final_bin_MAEE_100)),final_bin_MAEE_100$KC-final_bin_MAEE_100$t_test,pch=16,cex = 1.5,col="palevioletred1")

legend(x = "bottomright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25", "CV = 0.50", "CV = 0.75", "CV = 1.00"),
       col=c("orange", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
       cex = 0.8, pch = c(16,16,16,16,16))

# ind
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(b) IND for BC1",xlim=c(1,30),ylim=c(-0.15, 0.05),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.15, 0.05, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.15, 0.05, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC-final_bin_ind_00$t_test,pch=16,cex = 1.5, col="orange")
points(as.numeric(row.names(final_bin_ind_25)),final_bin_ind_25$KC-final_bin_ind_25$t_test,pch=16,cex = 1.5,col="limegreen")
points(as.numeric(row.names(final_bin_ind_50)),final_bin_ind_50$KC-final_bin_ind_50$t_test,pch=16,cex = 1.5,col="mediumpurple2")
points(as.numeric(row.names(final_bin_ind_75)),final_bin_ind_75$KC-final_bin_ind_75$t_test,pch=16,cex = 1.5,col="lightskyblue")
points(as.numeric(row.names(final_bin_ind_100)),final_bin_ind_100$KC-final_bin_ind_100$t_test,pch=16,cex = 1.5,col="palevioletred1")

legend(x = "bottomright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25", "CV = 0.50", "CV = 0.75", "CV = 1.00"),
       col=c("orange", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
       cex = 0.8, pch = c(16,16,16,16,16))

par(mfrow=c(1,1))



##################################################
# Figures of results for all variance estimators
# Size and Power
##################################################

par(mfrow=c(1,2))

#################################
# CV = 0.25 compared to CV = 0.00
#################################

# Size (Web Figures 1 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_25
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_25
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5)
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


# Power (Web Figures 2 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_25
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB-final_bin_MAEE_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB-final_bin_MAEE_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC-final_bin_MAEE_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD-final_bin_MAEE_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG-final_bin_MAEE_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG-final_bin_MAEE_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN-final_bin_MAEE_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_25
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB-final_bin_ind_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB-final_bin_ind_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC-final_bin_ind_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD-final_bin_ind_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG-final_bin_ind_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG-final_bin_ind_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN-final_bin_ind_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.25"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


#################################
# CV = 0.50 compared to CV = 0.00
#################################

# Size (Web Figures 3 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_50
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.50"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_50
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5)
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.50"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


# Power (Web Figures 4 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_50
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB-final_bin_MAEE_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB-final_bin_MAEE_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC-final_bin_MAEE_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD-final_bin_MAEE_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG-final_bin_MAEE_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG-final_bin_MAEE_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN-final_bin_MAEE_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.50"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_50
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB-final_bin_ind_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB-final_bin_ind_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC-final_bin_ind_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD-final_bin_ind_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG-final_bin_ind_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG-final_bin_ind_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN-final_bin_ind_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.50"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


#################################
# CV = 0.75 compared to CV = 0.00
#################################

# Size (Web Figures 5 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_75
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.75"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_75
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5)
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.75"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


# Power (Web Figures 6 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_75
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB-final_bin_MAEE_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB-final_bin_MAEE_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC-final_bin_MAEE_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD-final_bin_MAEE_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG-final_bin_MAEE_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG-final_bin_MAEE_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN-final_bin_MAEE_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.75"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_75
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB-final_bin_ind_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB-final_bin_ind_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC-final_bin_ind_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD-final_bin_ind_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG-final_bin_ind_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG-final_bin_ind_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN-final_bin_ind_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 0.75"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


#################################
# CV = 1.00 compared to CV = 0.00
#################################

# Size (Web Figures 7 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_100
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 1.00"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_100
plot(NULL,xlab="Scenario",ylab="Type I Error Rate",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(0,0.09),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, seq(0, 0.09, by = 0.01), las = 1, cex.axis = 0.85) 
axis(2, seq(0, 0.09, by = 0.005), labels=F, tick=T, las = 1, cex.axis = 1.5)
abline(h=0.036,lty=2,col="darkgray",lwd=2)
abline(h=0.064,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB.1,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB.1,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC.1,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD.1,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG.1,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG.1,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN.1,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB.1,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB.1,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC.1,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD.1,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG.1,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG.1,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN.1,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 1.00"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


# Power (Web Figures 8 in the supporting information)
# MAEE
final_bin<-final_bin_MAEE_100
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(a) MAEE",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MB-final_bin_MAEE_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$ROB-final_bin_MAEE_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$KC-final_bin_MAEE_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MD-final_bin_MAEE_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$AVG-final_bin_MAEE_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$FG-final_bin_MAEE_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_MAEE_00)),final_bin_MAEE_00$MBN-final_bin_MAEE_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 1.00"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))

# ind
final_bin<-final_bin_ind_100
plot(NULL,xlab="Scenario",ylab="Power Difference",xaxt = 'n',yaxt = 'n',
     main = "(b) IND",xlim=c(1,30),ylim=c(-0.2,0.1),
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5,cex = 1.5)
axis(1, 1, las = 1, cex.axis = 0.85)
axis(1, seq(5, 30, by = 5), las = 1, cex.axis = 0.85)
axis(1, seq(1, 30, by = 1), labels=F, tick=T, las = 1, cex.axis = 1.5) 
axis(2, round(seq(-0.2, 0.1, by = 0.025),3), las = 1, cex.axis = 0.85) 
axis(2, round(seq(-0.2, 0.1, by = 0.005),3), labels=F, tick=T, las = 1, cex.axis = 1.5) 
abline(h=-0.026,lty=2,col="darkgray",lwd=2)
abline(h=0.026,lty=2,col="darkgray",lwd=2)

points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MB-final_bin_ind_00$t_test,pch=3,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$ROB-final_bin_ind_00$t_test,pch=8,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$KC-final_bin_ind_00$t_test,pch=16,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MD-final_bin_ind_00$t_test,pch=0,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$AVG-final_bin_ind_00$t_test,pch=2,cex = 1.5)
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$FG-final_bin_ind_00$t_test,pch=1,cex = 2) # almost identical to KC
points(as.numeric(row.names(final_bin_ind_00)),final_bin_ind_00$MBN-final_bin_ind_00$t_test,pch=4,cex = 1.5)

points(as.numeric(row.names(final_bin)),final_bin$MB-final_bin$t_test,pch=3,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$ROB-final_bin$t_test,pch=8,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$KC-final_bin$t_test,pch=16,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$MD-final_bin$t_test,pch=0,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$AVG-final_bin$t_test,pch=2,cex = 1.5,col="red")
points(as.numeric(row.names(final_bin)),final_bin$FG-final_bin$t_test,pch=1,cex = 2,col="red") # almost identical to KC
points(as.numeric(row.names(final_bin)),final_bin$MBN-final_bin$t_test,pch=4,cex = 1.5,col="red")

legend(x = "bottomright",inset = 0,
       legend = c("MB", "BC0", "BC1 (KC)", "BC2 (MD)", "AVG (of KC & MD)", "BC3 (FG)", "BC4 (MBN)"), 
       cex = 0.8, pch = c(3,8,16,0,2,1,4))
legend(x = "topright",inset = 0,
       legend = c("CV = 0.00", "CV = 1.00"), col=c("black", "red"),
       cex = 0.8, pch = c(16,16))


par(mfrow=c(1,1))


