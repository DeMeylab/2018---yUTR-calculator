
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

#######>>>>> Part of the script used to perform 5-fold cross validation, but not used for the actual model

#path name where you saved 'output_analysis'
#csv='/output_analysis.csv'
#data=read.csv(csv)
#Y<-data$protein_abundance/min(data$protein_abundance)
#Xdata_only = data[,-(1:3)]

#data_preprocessed = cbind(Y,Xdata_only[1:13])
#set.seed(415324)
#data_preprocessed=data_preprocessed[sample(nrow(data_preprocessed)),]

#n <- 409
#nr <- nrow(data_preprocessed)
#r2_vector = c()
#data_opslag = data.frame(matrix(ncol=14))
#list_data_preprocessed = split(data_preprocessed, rep(1:ceiling(nr/n), each=n, length.out=nr))
#for(training_set in list_data_preprocessed){
#  print(nrow(training_set))
#  validation_pls = "CV"
#  Y_training=as.matrix(training_set[,1])
#  colnames(Y_training)<-c("Y")
#  X_training = as.matrix(training_set[,-1])
#  #colnames(X_training)<-paste("X",1:ncol(X_training),sep='')
#  predictors = ncol(X_training)
#  set.seed(1000)
#  #LINEAR SIMPLE
#  basic_pls = plsr(Y_training~X_training,validation=validation_pls,scale=T)
#  summary(basic_pls)  
#  components=4
#  training_predict = predict(basic_pls,X_training,ncomp=components)
#  r2_training=summary(lm(data.frame(training_predict)[,1]~training_set[,1]))$r.squared
#  print(r2_training)
#  coef_data = as.matrix(coef(basic_pls,ncomp=components,intercept = TRUE))
  
#  scale_data = cbind(as.matrix(basic_pls$scale),colnames(training_set[,-1]))
#  rownames(coef_data) = c("intercept",rownames(scale_data))
#  coef_data<-cbind(coef_data,c("intercept",colnames(training_set[,-1])))
#  colnames(coef_data) <- c("coefficient","name")
#  #print(coef_data)
#  colnames(scale_data) <- c("scale","name")
#  #print(scale_data)
#  data_opslag = rbind(data_opslag,c(r2_training,as.numeric(coef_data[2:14,"coefficient"])/as.numeric(scale_data[,"scale"])))
#}
#for(column_i in 1:ncol(data_opslag[2:6,])){
#  col_i = data_opslag[2:6,column_i]
#  print(paste(mean(col_i),' = ',sd(col_i),' : ',sd(col_i)/abs(mean(col_i))))  
#}

#pdf(".../pls_plots.pdf")
#biplot(basic_pls,var.axes = TRUE)
#validationplot(basic_pls)
#components=4
#training_predict = predict(basic_pls,X_training,ncomp=components)
#training_x = data.frame(training_predict)[,1]
#training_y = training_set[,1]
#training = data.frame(cbind(training_x,training_y))
#rownames(training)<-1:nrow(training)
#colnames(training)<-c("Experimental","Predicted")
#training_ggplot = ggplot(training,aes(x=Predicted,y=Experimental))+geom_point() + geom_abline(intercept = 0, slope = 1)
#nse_training = NSE(training_x,training_y)
#X_test<-as.matrix(test_set[,-(1)])
#X_test<-X_test[,columns_ok]
#test_predict = predict(basic_pls,X_test,ncomp=components)
#test_x = data.frame(test_predict)[,1]
#test_y = test_set[,1]
#test = data.frame(cbind(test_x,test_y))
#rownames(test)<-1:nrow(test)
#colnames(test)<-c("Experimental","Predicted")
#nse_test = NSE(data.frame(test_predict)[,1],test_set[,1])
#r2_test=summary(lm(data.frame(test_predict)[,1]~test_set[,1]))$r.squared

########<<<<<

######## Part of the script used to perform the actual PLS-regression

#path name where you saved 'output_analysis'
csv='.../output_analysis.csv'
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
#path name where to save the pls_plots
pdf(".../pls_plots.pdf")
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


coef_data = as.matrix(coef(basic_pls,ncomp=components,intercept = TRUE))
scale_data = cbind(as.matrix(basic_pls$scale),colnames(training_set[,-1]))
rownames(coef_data) = c("intercept",rownames(scale_data))
coef_data<-cbind(coef_data,c("intercept",colnames(training_set[,-1])))
colnames(coef_data) <- c("coefficient","name")
colnames(scale_data) <- c("scale","name")
#path name where you want to save the coefficient and scale data
write.csv(coef_data,file=".../coefficients.csv",quote=FALSE)
write.csv(scale_data,file=".../scales.csv",quote=FALSE)
