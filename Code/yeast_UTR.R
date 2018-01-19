
if (!require("pls")) {
  install.packages("pls")
  library(pls)
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("hydroGOF")) {
  install.packages("hydroGOF")
  library(hydroGOF)
}

split_training_test_set<-function(data_ordered,ratio=5,seed=1000){
  nparts = round(nrow(data_ordered)/ratio)
  list_splitted = split_vector_in_equal_parts(data_ordered,nparts,seed=seed)
  set.seed(seed)
  seeds = sample(1:1000,nparts)
  test_set = data_ordered[0,]
  training_set = data_ordered[0,]
  for(i in 1:nparts){
    segment = list_splitted[[i]]
    set.seed(seeds[i])
    select_for_test_set = sample(1:nrow(segment),1)
    test_set = rbind(test_set,segment[select_for_test_set,])
    training_set = rbind(training_set,segment[-select_for_test_set,])
  }
  return(list(training_set,test_set))
}
split_vector_in_equal_parts<-function(dataframe,nparts,seed){
  output_list<-list()
  nrest = nrow(dataframe)%%nparts
  row_dataframe = 1:nrow(dataframe)
  set.seed(seed)
  extra1_added_to_these_parts = sample(1:nparts,nrest)
  n_last_extracted = 1
  for(i in 1:nparts){
    nelements_in_part = round(nrow(dataframe)/nparts)
    if(i %in% extra1_added_to_these_parts){
      nelements_in_part =nelements_in_part+1
    }
    output_list[[i]]<-dataframe[(n_last_extracted:(n_last_extracted+nelements_in_part-1)),]
    n_last_extracted=n_last_extracted+nelements_in_part
  }
  return(output_list)
}


create_interaction_terms<-function(X){
  X_2<-X
  for(col in 1:ncol(X)){
    added_part<-X[,col]*X
    colnames(added_part)<-paste(colnames(X),':','X',col,sep='')
    X_2<-cbind(X_2,added_part)
  }
  return(X_2)
}

# 
# csv='/home/gpeters/Dropbox/Code/YeastUTR/output_analysis.csv'
# data=read.csv(csv)
# Y<-data$protein_abundance/min(data$protein_abundance)
# Xdata_only = data[,-(1:3)]
# 
# data_preprocessed = cbind(Y,Xdata_only)
# list_training_test = split_training_test_set(data_preprocessed,ratio=5)
# training_set = list_training_test[[1]]
# test_set = list_training_test[[2]]
# validation_pls = "CV"
# Y_training=as.matrix(training_set[,1])
# colnames(Y_training)<-c("Y")
# X_training = as.matrix(training_set[,-1])
# #colnames(X_training)<-paste("X",1:ncol(X_training),sep='')
# X_filtered = X_training[,0]
# columns_ok<-c()
# for(col in 1:ncol(X_training)){
#   if(sd(X_training[,col])!=0){ 
#     X_filtered<-cbind(X_filtered,X_training[,col])
#     columns_ok<-c(columns_ok,col)
#   }
# }
# colnames(X_filtered)<-colnames(X_training)[columns_ok]
# 
# predictors = ncol(X_filtered)
# set.seed(1000)
# #LINEAR
# basic_pls = plsr(Y_training~X_filtered,validation=validation_pls,scale=T,ncomp=50)
# summary(basic_pls)
# pdf("/home/gpeters/Dropbox/Code/YeastUTR/pls_plots_without_sticky_w_most_of_all_predictors.pdf")
# biplot(basic_pls,var.axes = TRUE)
# validationplot(basic_pls,ncomp=1:50)
# components=12
# training_predict = predict(basic_pls,X_filtered,ncomp=components)
# training_x = data.frame(training_predict)[,1]
# training_y = training_set[,1]
# training = data.frame(cbind(training_x,training_y))
# rownames(training)<-1:nrow(training)
# colnames(training)<-c("Experimental","Predicted")
# training_ggplot = ggplot(training,aes(x=Predicted,y=Experimental))+geom_point() + geom_abline(intercept = 0, slope = 1)
# nse_training = NSE(training_x,training_y)
# X_test<-as.matrix(test_set[,-(1)])
# X_test<-X_test[,columns_ok]
# test_predict = predict(basic_pls,X_test,ncomp=components)
# test_x = data.frame(test_predict)[,1]
# test_y = test_set[,1]
# test = data.frame(cbind(test_x,test_y))
# rownames(test)<-1:nrow(test)
# colnames(test)<-c("Experimental","Predicted")
# 
# nse_test = NSE(data.frame(test_predict)[,1],test_set[,1])
# r2_test=summary(lm(data.frame(test_predict)[,1]~test_set[,1]))$r.squared
# test_ggplot = ggplot(test,aes(x=Predicted,y=Experimental))+
#   geom_point()  + geom_abline(intercept = 0, slope = 1)+
#   ggtitle(paste("NSE =",round(nse_test,4),"R2 =",round(r2_test,4),"Met. =",validation_pls,"Comps =",components,"Preds =",predictors))
# print(test_ggplot)
# dev.off()
# 
# 
# 
# 
# #basic_pls$loadings
# 
# 
# 





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

#######<<<<<OLD

csv='/home/gpeters/Dropbox/Code/YeastUTR/output_analysis.csv'
data=read.csv(csv)
Y<-data$protein_abundance/min(data$protein_abundance)
Xdata_only = data[,-(1:3)]

data_preprocessed = cbind(Y,Xdata_only[1:13])
set.seed(415324)
data_preprocessed=data_preprocessed[sample(nrow(data_preprocessed)),]

n <- 409
nr <- nrow(data_preprocessed)
r2_vector = c()
data_opslag = data.frame(matrix(ncol=14))
list_data_preprocessed = split(data_preprocessed, rep(1:ceiling(nr/n), each=n, length.out=nr))
for(training_set in list_data_preprocessed){
  print(nrow(training_set))
  validation_pls = "CV"
  Y_training=as.matrix(training_set[,1])
  colnames(Y_training)<-c("Y")
  X_training = as.matrix(training_set[,-1])
  #colnames(X_training)<-paste("X",1:ncol(X_training),sep='')
  predictors = ncol(X_training)
  set.seed(1000)
  #LINEAR SIMPLE
  basic_pls = plsr(Y_training~X_training,validation=validation_pls,scale=T)
#  summary(basic_pls)  
  components=4
  training_predict = predict(basic_pls,X_training,ncomp=components)
  r2_training=summary(lm(data.frame(training_predict)[,1]~training_set[,1]))$r.squared
  print(r2_training)
  coef_data = as.matrix(coef(basic_pls,ncomp=components,intercept = TRUE))
  
  scale_data = cbind(as.matrix(basic_pls$scale),colnames(training_set[,-1]))
  rownames(coef_data) = c("intercept",rownames(scale_data))
  coef_data<-cbind(coef_data,c("intercept",colnames(training_set[,-1])))
  colnames(coef_data) <- c("coefficient","name")
  #print(coef_data)
  colnames(scale_data) <- c("scale","name")
  #print(scale_data)
  data_opslag = rbind(data_opslag,c(r2_training,as.numeric(coef_data[2:14,"coefficient"])/as.numeric(scale_data[,"scale"])))
}
for(column_i in 1:ncol(data_opslag[2:6,])){
  col_i = data_opslag[2:6,column_i]
  print(paste(mean(col_i),' = ',sd(col_i),' : ',sd(col_i)/abs(mean(col_i))))
  
  
}

pdf("/home/gpeters/Dropbox/Code/YeastUTR/pls_plots_without_sticky_simple.pdf")
biplot(basic_pls,var.axes = TRUE)
validationplot(basic_pls)
components=4
training_predict = predict(basic_pls,X_training,ncomp=components)
training_x = data.frame(training_predict)[,1]
training_y = training_set[,1]
training = data.frame(cbind(training_x,training_y))
rownames(training)<-1:nrow(training)
colnames(training)<-c("Experimental","Predicted")
training_ggplot = ggplot(training,aes(x=Predicted,y=Experimental))+geom_point() + geom_abline(intercept = 0, slope = 1)
nse_training = NSE(training_x,training_y)
X_test<-as.matrix(test_set[,-(1)])
#X_test<-X_test[,columns_ok]
test_predict = predict(basic_pls,X_test,ncomp=components)
test_x = data.frame(test_predict)[,1]
test_y = test_set[,1]
test = data.frame(cbind(test_x,test_y))
rownames(test)<-1:nrow(test)
colnames(test)<-c("Experimental","Predicted")
nse_test = NSE(data.frame(test_predict)[,1],test_set[,1])
r2_test=summary(lm(data.frame(test_predict)[,1]~test_set[,1]))$r.squared









########OLD >>>>>
csv='/home/gpeters/Dropbox/Code/YeastUTR/output_analysis.csv'
data=read.csv(csv)
Y<-data$protein_abundance/min(data$protein_abundance)
Xdata_only = data[,-(1:3)]

data_preprocessed = cbind(Y,Xdata_only[1:13])
list_training_test = split_training_test_set(data_preprocessed,ratio=5)
training_set = list_training_test[[1]]
test_set = list_training_test[[2]]
validation_pls = "CV"
Y_training=as.matrix(training_set[,1])
colnames(Y_training)<-c("Y")
X_training = as.matrix(training_set[,-1])
#colnames(X_training)<-paste("X",1:ncol(X_training),sep='')
predictors = ncol(X_training)
set.seed(1000)
#LINEAR SIMPLE
basic_pls = plsr(Y_training~X_training,validation=validation_pls,scale=T)
summary(basic_pls)
pdf("/home/gpeters/Dropbox/Code/YeastUTR/pls_plots_without_sticky_simple.pdf")
biplot(basic_pls,var.axes = TRUE)
validationplot(basic_pls)
components=4
training_predict = predict(basic_pls,X_training,ncomp=components)
training_x = data.frame(training_predict)[,1]
training_y = training_set[,1]
training = data.frame(cbind(training_x,training_y))
rownames(training)<-1:nrow(training)
colnames(training)<-c("Experimental","Predicted")
training_ggplot = ggplot(training,aes(x=Predicted,y=Experimental))+geom_point() + geom_abline(intercept = 0, slope = 1)
nse_training = NSE(training_x,training_y)
X_test<-as.matrix(test_set[,-(1)])
#X_test<-X_test[,columns_ok]
test_predict = predict(basic_pls,X_test,ncomp=components)
test_x = data.frame(test_predict)[,1]
test_y = test_set[,1]
test = data.frame(cbind(test_x,test_y))
rownames(test)<-1:nrow(test)
colnames(test)<-c("Experimental","Predicted")
nse_test = NSE(data.frame(test_predict)[,1],test_set[,1])
r2_test=summary(lm(data.frame(test_predict)[,1]~test_set[,1]))$r.squared
test_ggplot = ggplot(test,aes(x=Predicted,y=Experimental))+
  geom_point()  + geom_abline(intercept = 0, slope = 1)+
  ggtitle(paste("NSE =",round(nse_test,4),"R2 =",round(r2_test,4),"Met. =",validation_pls,"Comps =",components,"Preds =",predictors))
print(test_ggplot)
dev.off()


#LETS TRY MLR



mlr_training = lm(Y~dG_EFE+purineAG_in_min3+U_in_min3+A_in_min1+AA_in_min32+CG_in_min32+AC_in_min21+oof_uAUG+GACA_kmer+GG_kmer+CACC_kmer+CA_in_min76+CC_in_min76,data=training_set)
summary(mlr_training)
test_predict_mlr=predict(mlr_training,test_set)
test_x = test_predict_mlr
test_y = test_set[,1]
test_mlr = data.frame(cbind(test_x,test_y))
rownames(test_mlr)<-1:nrow(test_mlr)
colnames(test_mlr)<-c("Predicted","Experimental")
nse_test_mlr = NSE(test_mlr[,1],test_mlr[,2])
r2_test_mlr=summary(lm(test_mlr[,1]~test_mlr[,2]))$r.squared
predictors = length(coefficients(mlr_training))-1

test_ggplot = ggplot(test_mlr,aes(x=Predicted,y=Experimental))+
  geom_point()  + geom_abline(intercept = 0, slope = 1)+
  ggtitle(paste("NSE =",round(nse_test_mlr,4),"R2 =",round(r2_test_mlr,4),"Preds =",predictors))
print(test_ggplot)


#Y_test<-cbind(data$protein_abundance/min(data$protein_abundance),as.character(data$full_utr))

#full_python_script_test_predict = predict(basic_pls,data_preprocessed[1:100,-1],ncomp=components)

cbind(data_preprocessed[c(4,7),1],as.vector(as.character(data$full_utr))[c(4,7)])
test_x[1:2]
#head(test_set)
# 
# #plot loadings
# row<-1
# for(i in 1:components){
#   loadings = basic_pls$loadings[,paste("Comp",i)]
#   colnames_loadings<-names(loadings)
#   if(i == 1){
#     
#     loadings_cum=loadings
#     loadings_data<-data.frame(nvariable=numeric(length(loadings)*components),variable=character(length(loadings)*components),component=character(length(loadings)*components),loading=numeric(length(loadings)*components),stringsAsFactors = FALSE)
#   }else{
#     loadings_cum=loadings_cum+loadings
#   }
#   for(j in 1:length(loadings)){
#     loadings_data[row,]<-c(j,colnames_loadings[j],i,as.numeric(loadings[j]))
#     row=row+1
#   }
# }
# loadings_data$loading<-as.numeric(loadings_data$loading)
# loadings_data$nvariable<-as.numeric(loadings_data$nvariable)
# pdf("/home/gpeters/Dropbox/Code/YeastUTR/loadings_simple.pdf",paper="a4r")
# ggplot(loadings_data,aes(x=nvariable,y=loading))+geom_bar(position=position_dodge(), width=0.4,stat="identity") +facet_grid(component~.)
# dev.off()
# 
# loading_data<-data.frame(cbind(1:length(loadings_cum),loadings_cum,names(loadings)))
# colnames(loading_data)<-c("n","loadings","variable")
# loading_data$n<-unique(as.numeric(loadings_data$n))
# loading_data$loadings<-as.numeric(as.vector(loading_data$loadings))
# pdf("/home/gpeters/Dropbox/Code/YeastUTR/loadings_cum_simple.pdf",paper="a4r")
# ggplot(loading_data,aes(x=variable,y=loadings))+geom_bar(position=position_dodge(), width=0.4,stat="identity")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="top")
# dev.off()
# 
# 
 manual_matrix = X_test
 manual_scale = as.matrix(basic_pls$scale)
 for(i in 1:ncol(manual_matrix)){
   manual_matrix[,i]=manual_matrix[,i]/manual_scale[i,]
   
 }
 coef7<-coef(basic_pls,ncomp=components,intercept = TRUE)
 coef_matrix = as.matrix(as.matrix(coef7)[-1,])
# 
 test_manual = manual_matrix%*%as.matrix(as.matrix(coef7)[-1,])+as.matrix(coef7)[1,]
# 
# compare_data = cbind(test_manual,test_x)

coef_data = as.matrix(coef(basic_pls,ncomp=components,intercept = TRUE))

scale_data = cbind(as.matrix(basic_pls$scale),colnames(training_set[,-1]))
rownames(coef_data) = c("intercept",rownames(scale_data))
coef_data<-cbind(coef_data,c("intercept",colnames(training_set[,-1])))
colnames(coef_data) <- c("coefficient","name")

colnames(scale_data) <- c("scale","name")

write.csv(coef_data,file="/home/gpeters/Dropbox/Code/YeastUTR/coefficients.csv",quote=FALSE)
write.csv(scale_data,file="/home/gpeters/Dropbox/Code/YeastUTR/scales.csv",quote=FALSE)


# #LOG TRANSFORMATION not great idea...
# basic_pls_log = plsr(log(Y_training)~X_filtered,validation=validation_pls,scale=T,ncomp=50)
# summary(basic_pls_log)
# pdf("/home/gpeters/Dropbox/Code/YeastUTR/pls_log_plots_without_sticky_w_most_of_all_predictors.pdf")
# biplot(basic_pls_log,var.axes = TRUE)
# validationplot(basic_pls_log,ncomp=1:50)
# components=12
# training_predict = predict(basic_pls_log,X_filtered,ncomp=components)
# training_x = exp(data.frame(training_predict)[,1])
# training_y = training_set[,1]
# training = data.frame(cbind(training_x,training_y))
# rownames(training)<-1:nrow(training)
# colnames(training)<-c("Experimental","Predicted")
# training_ggplot = ggplot(training,aes(x=Predicted,y=Experimental))+geom_point() + geom_abline(intercept = 0, slope = 1)
# nse_training = NSE(training_x,training_y)
# X_test<-as.matrix(test_set[,-(1)])
# X_test<-X_test[,columns_ok]
# test_predict = predict(basic_pls_log,X_test,ncomp=components)
# test_x = exp(data.frame(test_predict)[,1])
# test_y = test_set[,1]
# test = data.frame(cbind(test_x,test_y))
# rownames(test)<-1:nrow(test)
# colnames(test)<-c("Experimental","Predicted")
# nse_test = NSE(test_x,test_set[,1])
# test_ggplot = ggplot(test,aes(x=Predicted,y=Experimental))+
#   geom_point()  + geom_abline(intercept = 0, slope = 1)+
#   ggtitle(paste("NSE =",nse_test,"Method =",validation_pls,"Components =",components,"Predictors =",predictors))
# print(test_ggplot)
# dev.off()

# 
# 
# #INCLUDE INTERACTIONS... NO HUGE IMPROVEMENT
# pdf("/home/gpeters/Dropbox/Code/YeastUTR/pls_plots_without_sticky_w_most_of_all_predictors.pdf")
# data_preprocessed = cbind(Y,Xdata_only[1:13])
# set.seed(1000)
# list_training_test = split_training_test_set(data_preprocessed,ratio=5)
# training_set = list_training_test[[1]]
# test_set = list_training_test[[2]]
# validation_pls = "CV"
# Y_training=as.matrix(training_set[,1])
# colnames(Y_training)<-c("Y")
# X_training = as.matrix(training_set[,-1])
# colnames(X_training)<-paste("X",1:ncol(X_training),sep='')
# # 
# # X_training2_filtered = X_training2[,0]
# # columns_ok<-c()
# # for(col in 1:ncol(X_training2)){
# #   if(sd(X_training2[,col])==0){
# #     print("checkenbak")    
# #   }else{
# #     X_training2_filtered<-cbind(X_training2_filtered,X_training2[,col])
# #     columns_ok<-c(columns_ok,col)
# #   }
# # }
# set.seed(1000)
# basic_pls = plsr(Y_training~X_training2,validation=validation_pls,scale=T)
# summary(basic_pls)
# biplot(basic_pls,var.axes = TRUE)
# validationplot(basic_pls)
# components=13
# training_predict = predict(basic_pls,X_training2[,columns_ok],ncomp=components)
# training_x = data.frame(training_predict)[,1]
# training_y = training_set[,1]
# training = data.frame(cbind(training_x,training_y,as.numeric(training_set$sd)/training_set[,1]))
# rownames(training)<-1:nrow(training)
# colnames(training)<-c("Experimental","Predicted")
# training_ggplot = ggplot(training,aes(x=Predicted,y=Experimental))+geom_point() + geom_abline(intercept = 0, slope = 1)
# #ggsave(training_ggplot,file="/home/gpeters/Dropbox/Documents/PhD/Rational_design_translational_riboswitches/Experiments/20151116/pls_ar.pdf")
# #print(training_ggplot)
# nse_training = NSE(training_x,training_y)
# X_test=as.matrix(test_set[,-(1)])
# colnames(X_test)<-paste("X",1:ncol(X_test),sep='')
# 
# X_test2 = cbind(create_interaction_terms(X_test[,thermodynamic_features]),X_test[,-thermodynamic_features])
# test_predict = predict(basic_pls,X_test2,ncomp=components)
# test_x = data.frame(test_predict)[,1]
# test_y = test_set[,1]
# test = data.frame(cbind(test_x,test_y))
# rownames(test)<-1:nrow(test)
# colnames(test)<-c("Experimental","Predicted")
# nse_test = NSE(data.frame(test_predict)[,1],test_set[,1])
# test_ggplot = ggplot(test,aes(x=Predicted,y=Experimental))+
#   geom_point()  + geom_abline(intercept = 0, slope = 1)+
#   # geom_errorbar(aes(x=Predicted,ymin=Experimental-Exp_sd, ymax=Experimental+Exp_sd), width=0.25)
#   ggtitle(paste("NSE =",nse_test,"Method =",validation_pls,"Components =",components))
# print(test_ggplot)
# dev.off()
# 





#USING PLS MODEL TO LINK ORIGINAL VARIABLES TO OUTPUT
#basic_pls$scores
dist_test<-scan("/home/gpeters/Dropbox/Code/YeastUTR/distribution_test.csv")