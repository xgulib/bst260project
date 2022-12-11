library(tidyverse)
library(readxl)
library(caret)
library(RNOmni)
library(pROC)
library(randomForest) 

#Read in
meta <- read_excel("/Users/xgu/Documents/Harvard/Fall 2022/BST260/bst260project/41591_2022_1688_MOESM3_ESM.xlsx", 
                          sheet = 13, skip = 1, na = "NA", col_types = "guess")
demo <- read_excel("/Users/xgu/Documents/Harvard/Fall 2022/BST260/bst260project/41591_2022_1688_MOESM3_ESM.xlsx", 
                   sheet = 10, skip = 1, na = "NA", col_types = c("text", "text", rep("numeric", 22), "text", "text", "text"))

#Selection and filtering
demo_new <- demo %>%
  filter(Status %in% c("IHD372", "MMC372")) %>%
  mutate(case = case_when(Status == "MMC372" ~ 0, TRUE ~ 1),
         Gender = case_when(Gender == "Male" ~ 1, TRUE ~ 0)) %>%
  rename(age = "Age (years)", tag = "Fasting plasma triglycerides (mmol/L)", 
         adiponectin = "Fasting plasma adiponectin (mg/L)", crp = "Fasting plasma CRP (mg/L)",
         sbp = "Systolic blood pressure (mmHg)", dbp = "Diastolic blood pressure (mmHg)", 
         lvef = "Left ventricular ejection fraction (%)", act = "Physical activity (h/week)") %>%
  select(ID, case, age, tag, adiponectin, crp, sbp, dbp, Gender)
  
meta_new <- meta %>%
  filter(Status %in% c("IHD372", "MMC372")) %>%
  select(-c(Status))

#Merge
main <- demo_new %>%
  left_join(meta_new, by = "ID") %>%
  filter(acetate != "NA", spermidine != "NA")

#Check missing
pctmiss <- function(x){
  pctmiss <- sum(is.na(x))/length(x)
  return(pctmiss)
}
miss <- data.frame(sapply(main, pctmiss))

#replace missing
var <- main %>%
  select(-c("ID", "case")) %>%
  mutate(tag = case_when(is.na(tag) ~ median(tag, na.rm = TRUE), TRUE ~ tag),
         adiponectin = case_when(is.na(adiponectin) ~ median(adiponectin, na.rm = TRUE), TRUE ~ adiponectin),
         crp = case_when(is.na(crp) ~ median(crp, na.rm = TRUE), TRUE ~ crp),
         sbp = case_when(is.na(sbp) ~ median(sbp, na.rm = TRUE), TRUE ~ sbp),
         dbp = case_when(is.na(dbp) ~ median(dbp, na.rm = TRUE), TRUE ~ dbp))

#Preprocessing
nzv <- nearZeroVar(var)
col_index <- setdiff(1:ncol(var), nzv)
length(col_index)

var_proc <- var[,col_index]


#check normality
normality <- data.frame()
for (i in 1:length(colnames(var_proc))){
  normality[i, 1] <- colnames(var_proc)[i]
  normality[i, 2] <- shapiro.test(pull(var_proc[,i]))$p.value
  colnames(normality) <- c("metabolites", "shapiro.p")
}

sum(normality$shapiro.p > 0.05) #only 101 normal
which(normality$shapiro.p > 0.05)
hist(pull(var_proc[,6]))
shapiro.test(pull(var_proc[,6]))

#Log transformation not work
m <- as.matrix(var_proc[,8:1422])
exp_m <- exp(m)
var_proc_exp <- cbind(var_proc[,1:7], as.data.frame(exp_m))

var_proc_int <- as.data.frame(sapply(var_proc_exp, RankNorm))

#check normality again!
normality_int <- data.frame()
for (i in 1:length(colnames(tibble(var_proc_int)))){
  normality_int[i, 1] <- colnames(tibble(var_proc_int))[i]
  normality_int[i, 2] <- shapiro.test(pull(tibble(var_proc_int)[,i]))$p.value
  colnames(normality_int) <- c("metabolites", "shapiro.p")
}

sum(normality_int$shapiro.p > 0.05) #now 840 normal





#Split data
set.seed(34324)
main_new <- cbind(main[,1:2], var_proc_int)
train_index <- createDataPartition(main_new$case, times = 1, p = 0.8, list = FALSE)
train_set <- main_new[train_index,]
test_set <- main_new[-train_index,]

#X and Y matrix
x_train <- as.matrix(train_set[,3:1424])
y_train <- factor(train_set$case)
x_test <- as.matrix(test_set[,3:1424])
y_test <- factor(test_set$case)




#####KNN
set.seed(3245)
b <- 10
control_knn <- trainControl(method = "cv", number = b, p = .9) 
train_knn <- train(x_train, y_train,  
                   method = "knn",  
                   tuneGrid = data.frame(k = seq(5,150,10)), 
                   trControl = control_knn) 
ggplot(train_knn, highlight = TRUE)
train_knn$bestTune 
train_knn$results$Accuracy
fit_knn <- knn3(x_train, y_train, k = train_knn$bestTune$k) 

y_pred_knn <- predict(fit_knn, x_test, type = "class")
y_pred_knn_p <- predict(fit_knn, x_test, type = "prob")
confusionMatrix(y_pred_knn, y_test)
confusionMatrix(y_pred_knn, y_test)$overall["Accuracy"]

#ROC
roc_knn <- roc(as.factor(test_set$case), y_pred_knn_p[, 2])
plot(roc_knn, print.thres="best", type = "line", print.auc = TRUE, grid = TRUE, ylim = c(0,1))

# library(ROCR)
# pred <- prediction(y_pred_knn_p[, 2], test_set$case)
# perf <- performance(pred,"tpr","fpr")
# plot(perf, colorize = T, lwd = 2)




#######PCA
col_means <- colMeans(x_train)

pca <- prcomp(x_train)
summary(pca)
summary(pca)$importance[3,] ##69 pc

pc <- 69
x_train_pc <- pca$x[,1:pc] 


#####PCA + glm
glm_tmp <- as.data.frame(cbind(train_set$case, x_train_pc))
glm_tmp <- glm_tmp %>% rename(case = V1)
fit_glm <- glm(case ~., data = glm_tmp, family = "binomial")
y_prob <- predict(fit_glm, as.data.frame(x_test_pc), type = "response")
y_pred_glm <- factor(ifelse(y_prob > 0.5, 1, 0))
confusionMatrix(y_pred_glm, y_test)$overall[["Accuracy"]] 

#ROC
roc_glm <- roc(as.factor(test_set$case), y_prob)
plot(roc_glm, print.thres="best", type = "line", print.auc = TRUE, grid = TRUE, ylim = c(0,1), col = "Red")


#######PCA + KNN
set.seed(5436)
b <- 10
control_pca <- trainControl(method = "cv", number = b, p = .9) 
train_pcaknn <- train(x_train_pc, y_train,  
                   method = "knn",  
                   tuneGrid = data.frame(k = seq(5,200,20)), 
                   trControl = control_pca) 
ggplot(train_pcaknn, highlight = TRUE)
train_pcaknn$bestTune 
train_pcaknn$results$Accuracy
fit_pcaknn <- knn3(x_train_pc, y_train, k = train_pcaknn$bestTune$k) 

x_test_pc_pre <- sweep(x_test,2,col_means) %*% pca$rotation 
x_test_pc <- x_test_pc_pre[,1:pc]
y_pred_pcaknn <- predict(fit_pcaknn, x_test_pc, type = "class")
y_pred_pcaknn_p <- predict(fit_pcaknn, x_test_pc, type = "prob")
confusionMatrix(y_pred_pcaknn, factor(test_set$case))$overall["Accuracy"]

#ROC
roc_pcaknn <- roc(as.factor(test_set$case), y_pred_pcaknn_p[, 2])
plot(roc_pcaknn, print.thres="best", type = "line", print.auc = TRUE, grid = TRUE, ylim = c(0,1), col = "Red")



##########Random forest
set.seed(3982)
b <- 5
control_rf <- trainControl(method="cv", number = b, p = .9) 
train_rf <-  train(x_train, y_train,  
                   method = "rf",  
                   ntree = 200,
                   trControl = control_rf, 
                   tuneGrid = data.frame(mtry = seq(10, 500, 20))) 
ggplot(train_rf, highlight = TRUE)
train_rf$bestTune
train_rf$results$Accuracy
fit_rf <- randomForest(x_train, y_train,  
                       mtry = train_rf$bestTune$mtry,
                       minNode = 10) 
plot(fit_rf) 

y_pred_rf <- predict(fit_rf, x_test, type = "class")
y_pred_rf_p <- predict(fit_rf, x_test, type = "prob")
confusionMatrix(y_pred_rf, y_test)
confusionMatrix(y_pred_rf, y_test)$overall["Accuracy"]

roc_rf <- roc(as.factor(test_set$case), y_pred_rf_p[, 2])
plot(roc_rf, print.thres="best", type = "line", print.auc = TRUE, grid = TRUE, ylim = c(0,1), col = "blue")



#Ensemble
p_knn <- y_pred_knn_p
p_rf <- y_pred_rf_p / rowSums(y_pred_rf_p) 
p_pcaknn <- y_pred_pcaknn_p
p_glm <- as.matrix(cbind(1-y_prob, y_prob))
colnames(p_glm) <- c(0, 1)
p <- (p_rf + p_knn + p_pcaknn + p_glm)/4
y_pred <- factor(apply(p, 1, which.max)-1) 
confusionMatrix(y_pred, y_test)
confusionMatrix(y_pred, y_test)$overall["Accuracy"] 

roc_es <- roc(as.factor(test_set$case), p[, 2])
plot(roc_es, print.thres="best", type = "line", print.auc = TRUE, grid = TRUE, ylim = c(0,1), col = "green")
