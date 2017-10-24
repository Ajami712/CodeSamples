#==========================
# This script uses lipidomics data and rate of glucose uptake (Rd) data to make an elasticnet logistic regression model 
# that predicts whether patients are insulin-resistant.
# Rd values are the response variables and lipidomics data comprise the features.
#==========================

library(glmnet)
library(ROCR)

#=================================================================================
#Load Model-Building Data
#=================================================================================

#Assign all important data components to parameters
X <- read.csv("FA_AT.data.csv",header=FALSE)#f.data
fstate <- read.csv("FA_AT.state.csv",header=FALSE,quote="'") #f.state
Patients <- read.csv(file="FA_AT.patients.csv",header=FALSE,quote="'")#f.patients
Lipids <- read.csv("FA_AT.lipids.csv",header=FALSE,quote="'")

YP <- as.vector(unlist(read.csv("FA_AT.Rd.csv",header=FALSE,quote="'")))

X <- as.matrix(X)
Patients <- as.vector(unlist(Patients))
fstate <- as.vector(unlist(fstate))
Lipids <- as.vector(unlist(Lipids))

Y <- 0 #0-vector with same length size as fstate
for(i in 1:(length(fstate) - 1)){Y <- c(Y,0)}

#Correlates each value in Y and fstate to the corresponding patient
rownames(X) <- Patients
colnames(X) <- Lipids
names(fstate) <- Patients
names(Y) <- Patients


#=================================================================================
#Build Model from Non-Screener Model-Building Data
#=================================================================================

#Remove screeners
Xh = X[-grep("S",Patients),]
Yh = Y[grep("h",Patients)]
fstateh = fstate[grep("h",Patients)]

#Set all NRd values to equal 1 within Y
for(i in 1:length(fstateh)){if(fstateh[i] == "NRd"){Yh[i] <- 1}else{Yh[i]<-0}}
		
#Calculate Z-Score
Z <- matrix(0,nrow(Xh),ncol(Xh))


for(i in 1:ncol(Xh)){for(j in 1:nrow(Xh)){Z[j,i] = (Xh[j,i] - mean(Xh[,i]))/sd(Xh[,i])}}
for(i in 1:ncol(Xh)){for(j in 1:nrow(Xh)){if(Z[j,i] == "NaN"){Z[j,i] <- 0}}}
colnames(Z) <- Lipids

LambdaMin <- 1000
for(i in 1:100){
	cverr <- cv.glmnet(Z,Yh,family="binomial",nlambda=1000,alpha=0.5)
	#model = glmnet(Z,Y,family="binomial",nlambda=10,alpha=1,)
	if(LambdaMin > cverr$lambda.min){
		LambdaMin <- cverr$lambda.min
		# Coeffs <- coef(cverr)
    cverrBest <- cverr
	}
}
#The LambdaMin varies a bit every time we run cverr
#To mitigate this, I have cverr run many times, and take the min of LambdaMin

#We can use that recovered LambdaMin to re-make its glmnet model


pred <- predict(cverrBest,newx=Z,s="lambda.min")

pred_object <- prediction(pred,Yh)

performance(pred_object,measure = 'auc')



#=================================================================================
#Use Model On Screener Model-Building Data
#=================================================================================

#Make screener data
Xs = X[grep("S",Patients),]
Ys = Y[-grep("h",Patients)]
fstateS = fstate[-grep("h",Patients)]


#Set all NRd values to equal 1 within Y
for(i in 1:length(fstateS)){if(fstateS[i] == "NRd"){Ys[i] <- 1}else{Ys[i]<-0}}

#Calculate Z-Score
Z_S <- matrix(0,nrow(Xs),ncol(Xs))


for(i in 1:ncol(Xs)){for(j in 1:nrow(Xs)){Z_S[j,i] = (Xs[j,i] - mean(Xh[,i]))/sd(Xh[,i])}}
for(i in 1:ncol(Xs)){for(j in 1:nrow(Xs)){if(Z_S[j,i] == "NaN"){Z_S[j,i] <- 0}}}
colnames(Z_S) <- Lipids


pred <- predict(cverrBest,newx=Z_S,s="lambda.min")
pred_object <- prediction(pred,Ys)
perf <- performance(pred_object,measure = 'auc')
perf_roc <-performance(pred_object,"tpr","fpr")
plot(perf_roc)



#=================================================================================
#Load Post-Treatment Data
#=================================================================================

#Assign all important data components to parameters
XPost <- read.csv("FA_AT_Data_F3.csv",header=FALSE)#f.data
fstatePost <- read.csv("FA_AT_State_F.csv",header=FALSE,quote="'") #f.state
PatientsPost <- read.csv(file="FA_AT_Patients_F.csv",header=FALSE,quote="'")#f.patients
LipidsPost <- read.csv("FA_AT_Lipids_F.csv",header=FALSE,quote="'")

XPost <- as.matrix(XPost)
PatientsPost <- as.vector(unlist(PatientsPost))
fstatePost <- as.vector(unlist(fstatePost))
LipidsPost <- as.vector(unlist(LipidsPost))

YPost <- 0 #0-vector with same length size as fstate
for(i in 1:(length(fstatePost) - 1)){YPost <- c(YPost,0)}

#Correlates each value in Y and fstate to the corresponding patient
rownames(XPost) <- PatientsPost
colnames(XPost) <- LipidsPost
names(fstatePost) <- PatientsPost
names(YPost) <- PatientsPost


#=================================================================================
#Use Model On Screener Model-Building Data
#=================================================================================

#Set all NRd values to equal 1 within Y
for(i in 1:length(fstatePost)){if(fstatePost[i] == "NRd"){YPost[i] <- 1}else{YPost[i]<-0}}

#Calculate Z-Score
Z_Post <- matrix(0,nrow(XPost),ncol(XPost))


for(i in 1:ncol(XPost)){for(j in 1:nrow(XPost)){Z_Post[j,i] = (XPost[j,i] - mean(Xh[,i]))/sd(Xh[,i])}}
for(i in 1:ncol(XPost)){for(j in 1:nrow(XPost)){if(Z_Post[j,i] == "NaN"){Z_Post[j,i] <- 0}}}
colnames(Z_Post) <- LipidsPost


pred <- predict(cverrBest,newx=Z_Post,s="lambda.min")
pred_object <- prediction(pred,YPost)
perf <- performance(pred_object,measure = 'auc')
perf_roc <-performance(pred_object,"tpr","fpr")
plot(perf_roc)
