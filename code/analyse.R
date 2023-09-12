###############################################
## data
###############################################
data <- read.csv("patients_test.csv")
patientsID <- unique(data$patientsID)
patient1 = data[data$patientsID == patientsID[1],]
patient2 = data[data$patientsID == patientsID[2],]
patient3 = data[data$patientsID == patientsID[3],]
patient4 = data[data$patientsID == patientsID[4],]
## initial data (Ex patient1)
plot(patient1$signals,type = "l",main = patientsID[1],ylab = "SaO2%",xlab = "Time in s")
## Removing error values 
## Delta filter method 
delta_filter <- function(sig,x){
  sig_filtred <- rep(NA,length(sig))
  sig_filtred[1] <- sig[1]
  previous <- 1
  for(i in 2:length(sig)){
    if((abs(sig[previous] - sig[i])/sig[previous])*100<x){
      sig_filtred[i] <- sig[i]
      previous <- i
    }
    else {
      sig_filtred[i]
    }
  }
  return (sig_filtred)
}
## application
filtred_p1 <- delta_filter(patient1$signals,5)
plot(filtred_p1, type = "l", col = "black", ylab = "SaO2%", xlab = "Time in s")
# plot(patient1$signals, col = "red")
# legend("bottomleft", legend = c("Original sig", "Filtered sig"), col = c("black", "red"), lty = 1)
fpatient1 = patient1[!is.na(filtred_p1),]
fpatient2 = patient2[!is.na(delta_filter(patient2$signals,4)),]
fpatient3 = patient3[!is.na(delta_filter(patient3$signals,4)),]
fpatient4 = patient4[!is.na(delta_filter(patient4$signals,4)),]
filtred_data = rbind(fpatient1,fpatient2,fpatient3,fpatient4)
## TODO block method/mean or median replacement
###############################################################
################
## Representing signals
get_apneas_idxs <- function(){
  total_time <- sum(data.apnea$Evènement_Apnée.Centrale) + sum(data.apnea$Evènement_Apnée.Obstructive)
  patients <- unique(data.apnea$patientsID)
  apneas = data.frame(patientID = NA , begin = NA, end = NA,type = NA)
  indexation[1,c(2,3)] = c(1,data.apnea$Durée[1])
  indexation$patientID = data.apnea$patientsID[1]
  s <- 2
  while(sum(indexation$Lasting)<total_time){
    new_idx <- indexation$idxStart[s-1] +  indexation$Lasting[s-1] 
    indexation[s,1] <- data.apnea$patientsID[new_idx]
    indexation[s,2] <- new_idx
    indexation[s,3] <- data.apnea$Durée[new_idx]
    s = s + 1
  }
  curves <- vector("list")
  time <- vector("list")
  for(i in 1:length(patients)){
    patient <- patients[i] 
    apneas = indexation[indexation$patientID == patient,]
    # colors <- rainbow(nrow(apneas))
    curves[[patients[i]]] <- list()
    for (k in 1:nrow(apneas)) {
      time[[patient]][[k]] <- seq(1:apneas$Lasting[k])
      j <- apneas$idxStart[k]
      curves[[patient]][[k]] <- data.apnea[(j:(j+apneas$Lasting[k]-1)),3]
    }
  }
  return(list(curves = curves,time = time))
}
central_sa = get_curves_time("csa")
obstructive_sa = get_curves_time("osa")
display_signals <- function(patient,res){
  y_min = min(sapply(res$curves[[patient]],min))
  t_max = max(sapply(res$time[[patient]],max))
  plot(NULL, xlim = c(1, t_max), ylim = c(y_min, 100), xlab = "Time", ylab = "SaO2%", main = patient)
  colors = rainbow(length(res$time[[patient]]))
  for (i in 1:length(res$time[[patient]])) {
    lines(res$time[[patient]][[i]],res$curves[[patient]][[i]],type='l',lty=1,col = colors[i])
  }
}
###############################################################
######## Display signal for each patient
###############################################################
p1 = "PA1011"
p2 = "PA1089"
p3 = "PA1150"
p4 = "PA1141"
display_signals(p1,central_sa)
display_signals(p2,central_sa)
display_signals(p3,central_sa)
display_signals(p4,central_sa)
patient1 = data[data$patientsID == "PA1011",]
patient2 = data[data$patientsID == "PA1089",]
patient3 = data[data$patientsID == "PA1150",]
patient4 = data[data$patientsID == "PA1141",]
plot(patient1$signals,type = "l")
plot(patient2$signals,type = "l")
plot(patient3$signals,type = "l")
plot(patient4$signals,type = "l")

###############################################################
######## Remove error values <70%
###############################################################
err_idx = data$signals<70
plot(data[err_idx,]$signals,type = "l")

##### delta filter 
delta_filter <- function(sig,x){
  sig_filtred <- rep(NA,length(sig))
  sig_diff <- diff(sig) 
  i = 1
  while(i<(length(sig)-1)){
    if(abs(sig[i]-sig[i+1])<x){
      sig_filtred[i] <- sig[i]
      sig_filtred[i+1] <- sig[i+1]
    }
    if (sig[i] < 70 && sig_diff[i] == 0){
          sig_filtred[i] <- NA
    }
      i <- i+1
  }
  
  return (sig_filtred)
}
pat1_f <- delta_filter(patient3$signals,4)
# pat2_f <- delta_filter(pat2_f,6)
plot(patient3$signals,type='l')
lines(pat1_f,col="red")
length(pat1_f[which(pat1_f<70)])
#abline(h = pat1_f[which(pat1_f<70)],col= "blue")
pat_id = p4
patient = data[data$patientsID == pat_id,]
sig <- patient$signals 
##############################################
######### Derivative of curves
##############################################
plot(diff(patient$signals),type="l")
m = mean(patient$signals[diff(patient$signals)==0])
plot(patient$signals,type = "l",col = "black")
abline(h=m,col = "red")
##############################################
######### Display each apnea + 60s in time 
##############################################
indicatrice <- which(patient$Evènement_Apnée.Centrale ==1)
apnees_centrales <- vector(mode ="list",length=sum(diff(indicatrice)>1))
j<-1
apnee <- patient$signals[indicatrice[1]]
i <- 1
while(i<length(indicatrice)){
  i <- i+1
  if(indicatrice[i]-indicatrice[i-1]==1){
    apnee <- c(apnee,patient$signals[indicatrice[i]])
  }else{
    apnees_centrales[[j]] <-
      c(apnee,patient$signals[indicatrice[i]+(1:60)])
    j <- j+1
    apnee <- patient$signals[i]
  }
}
plot(apnees_centrales[[4]],type='l',lty=1)
apnees_label <- vector(mode = "list",length=sum(diff(indicatrice)>1))
j<-1
apnee <- patient$signals[indicatrice[1]]
i <- 1
while(i<length(indicatrice)){
  i <- i+1
  if(indicatrice[i]-indicatrice[i-1]==1){
    apnee <- c(apnee,patient$signals[indicatrice[i]])
  }else{
    apnees_label[[j]] <- c(apnee,rep(NA,60))
    j <- j+1
    apnee <- patient$signals[i]
  }
}
ix <- 6
plot(apnees_centrales[[ix]],type='l',lty=1)
lines(apnees_label[[ix]],col='red',lwd=2)

##############################################
######### Decalage = valeur constante
##############################################
sig_shift <- rep(NA,length(patient$signals))
shift <- 18
sig <- patient$signals 
for(i in (shift+1):length(sig)){
  if(patient$Evènement[i-shift]=="Apnée Centrale"){
    sig_shift[i] <- patient$signals[i]
  } 
}

par(mfrow=c(1,1))
plot(patient$signals[4000:6000],type='l')
lines(sig_shift[4000:6000],col='red',lwd=2)
##############################################
######### Decalage = duree de l'apnee 
##############################################
sig_shift <- rep(NA,length(patient$signals))
duree <- sapply(central_sa$time[pat_id],lengths)
j <- 0
i <- 1
while(i <length(sig)){
  if(patient$Evènement[i]=="Apnée Centrale"){
    j <- j + 1 
    for (l in 1:duree[[j]]){
      sig_shift[i:(i + duree[[j]] - 1)] <- patient$signals[i: (i + duree[[j]] - 1)]
      } 
    i <- i + duree[[j]]
    } 
  else {
    i <- i + 1
  }
}
print(sum(!is.na(sig_shift)))
print(sum(unlist(duree)))
par(mfrow=c(1,1))
plot(patient$signals[1:1000],type='l')
lines(sig_shift[1:1000],col="red",lwd=2)
## Delta filter method 
delta_filter <- function(sig,x){
  sig_filtred <- rep(NA,length(sig))
  sig_filtred[1] <- sig[1]
  previous <- 1
  for(i in 2:length(sig)){
    if((abs(sig[previous] - sig[i])/sig[previous])*100<x){
      sig_filtred[i] <- sig[i]
      previous <- i
    }
  }
  return (sig_filtred)
}
## application
filtred_p1 <- delta_filter(patient1$signals,4)
plot(patient1$signals, type = "l", col = "black", ylab = "SaO2%", xlab = "Time in s")
lines(filtred_p1, col = "red")
legend("bottomleft", legend = c("Original sig", "Filtered sig"), col = c("black", "red"), lty = 1)
fpatient1 = patient1[!is.na(filtred_p1),]
fpatient2 = patient2[!is.na(delta_filter(patient2$signals,5)),]
fpatient3 = patient3[!is.na(delta_filter(patient3$signals,5)),]
fpatient4 = patient4[!is.na(delta_filter(patient4$signals,4)),]
## TODO block method/mean or median replacement
###############################################################
################

###############################################################
######## Display
###############################################################

plot(fpatient1$signals,type = "l",main = patientsID[1],ylab = "SaO2%",xlab = "Time in s")
plot(fpatient2$signals,type = "l",main = patientsID[2],ylab = "SaO2%",xlab = "Time in s")
plot(fpatient3$signals,type = "l",main = patientsID[3],ylab = "SaO2%",xlab = "Time in s")
plot(fpatient4$signals,type = "l",main = patientsID[4],ylab = "SaO2%",xlab = "Time in s")

##############################################
######### Display each apnea label and 60s after 
##############################################
patient = fpatient1
indicatrice <- which(patient$Evènement_Apnée.Centrale ==1)
apnees_centrales <- vector(mode ="list",length=sum(diff(indicatrice)>1))
apnees_label <- vector(mode = "list",length=sum(diff(indicatrice)>1))
j<-1
apnee <- patient$signals[indicatrice[1]]
i <- 1
while(i<length(indicatrice)){
  i <- i+1
  if(indicatrice[i]-indicatrice[i-1]==1){
    apnee <- c(apnee,patient$signals[indicatrice[i]])
  }else{
    apnees_centrales[[j]] <-c(apnee,patient$signals[indicatrice[i]+(1:60)])
    apnees_label[[j]] <- c(apnee,rep(NA,60))
    j <- j+1
    apnee <- patient$signals[i]
  }
}

par(mfrow = c(4, 4))

for (ix in 1:16) {
  plot(apnees_centrales[[ix]], type = 'l', lty = 1,ylab = "SaO2%",xlab = "time in s")
  lines(apnees_label[[ix]], col = 'blue', lwd = 2)
}


##############################################
######### Shifting test !!
##############################################
patient = fpatient2
patient$app = rep(0,nrow(patient)) 
sig_shift <- rep(NA,length(patient$signals))
## size of shifting
shift1 <- 22
shift2 <- 25
shift3 <- 10
shift4 <- 10
shift <- 20
sig <- patient$signals 
for(i in (shift+1):length(sig)){
  if(patient$Evènement[i-shift]=="Apnée Centrale"){
    sig_shift[i] <- patient$signals[i]
    patient$app[i] = 1
  } 
}
indicatrice <- which(patient$app ==1)
apnees_centrales <- vector(mode ="list",length=sum(diff(indicatrice)>1))
apnees_label <- vector(mode = "list",length=sum(diff(indicatrice)>1))
j<-1
apnee <- patient$signals[indicatrice[1]]
i <- 1
while(i<length(indicatrice)){
  i <- i+1
  if(indicatrice[i]-indicatrice[i-1]==1){
    apnee <- c(apnee,patient$signals[indicatrice[i]])
  }else{
    apnees_centrales[[j]] <-c(apnee,patient$signals[indicatrice[i]+(1:30)])
    apnees_label[[j]] <- c(apnee,rep(NA,30))
    j <- j+1
    apnee <- patient$signals[i]
  }
}
par(mfrow = c(4, 4))
col = "blue" #rainbow(length(apnees_centrales))
for (ix in 1:length(apnees_centrales)) {
  plot(apnees_centrales[[ix]], type = 'l', lty = 1,ylab = "SaO2%",xlab = "time in s",main = length(which(!is.na(apnees_label[[ix]]))))
  lines(apnees_label[[ix]], col = col, lwd = 2)
}
