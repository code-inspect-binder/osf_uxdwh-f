# This is R code for the following paper:
# Applying dominance analysis and random forest to improve interpretations of variable importance in multiple regression in L2 research


# Install Required Packages -----------------------------------------------------------------

list.of.packages <- c("remote","rpsychi","MASS","ppcor","lavaan","semPlot","glmnet","relaimpo","yhat","boot","ggplot2","randomForest","Boruta","car","kernlab","xgboost","caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install rpsych from archive
# This is necessary because these packages are not updated for the current version of R.
library(remotes)
install_version("rpsychi", "0.8", repos = "http://cran.us.r-project.org")
install_version("dominanceanalysis", "2.0.0", repos = "http://cran.us.r-project.org")




# Table 1 -----------------------------------------------------------------

correl <- matrix(c(
  1, 0.62, 0.43, 0.56, 0.61,
  0.62, 1, 0.67, 0.55, 0.72,
  0.43, 0.67, 1, 0.67, 0.65,
  0.56, 0.55, 0.67, 1, 0.79,
  0.61, 0.72, 0.65, 0.79, 1),
  nrow=5,
  dimnames=list(c("Speaking","Vocabulary","Grammar","Writing","Reading"),
                c("Speaking","Vocabulary","Grammar","Writing","Reading")))

# Multiple Regression Analysis
library(rpsychi)
multreg.second(Speaking ~ Vocabulary+Grammar+Writing+Reading, corr=correl, n=100)

# Creating Simulated Data 
library(MASS)
set.seed(88)
mu <- rep(c(0), times = length(colnames(correl)))
simdat <- mvrnorm(n=100, mu=mu, Sigma=correl, empirical=TRUE)
colnames(simdat) <- c("Speaking","Vocabulary","Grammar","Writing","Reading")
cor(simdat)

# Calculating p-values
z <- scale(as.data.frame(simdat))
z <- data.frame(z)
summary(lm(Speaking~ ., z)) 

# Partial and Semipartial Correlation (Not in the Paper)
library(ppcor)
pcor(simdat)$estimate # Partial Correlation
spcor(simdat)$estimate # Semipartial Correlation

# Path Tracing (Just for Reference)
# Speaking and Vocabulary
0.47+(-0.19*0.67)+(0.31*0.55)+(0.15*0.72)
# Speaking and Grammer
-0.19+(0.47*0.67)+(0.31*0.67)+(0.15*0.65)
# Speaking and Writing
0.31+(0.47*0.55)+(-0.19*0.67)+(0.15*0.79)
# Speaking and Writing
0.15+(0.47*0.72)+(-0.19*0.65)+(0.31*0.79)




# Figure 1 (Quite Literally)-----------------------------------------------------------------
# https://osf.io/adnj2

correl <- matrix(c(
  1, 0.5, 0.4, 0.4,
  0.5, 1, 0.2, 0.5,
  0.4, 0.2, 1, 0.3,
  0.4, 0.5, 0.3, 1),
  nrow=4,
  dimnames=list(c("X1","X2","X3","Y"),
                c("X1","X2","X3","Y")))

# Multiple Regression Analysis
library(rpsychi)
multreg.second(Y~ X1+X2+X3, corr=correl, n=100)

# SEM
library(lavaan)
regression.model <-'
Y ~ a*X1 + b*X2 + c*X3
# residual variance of Y
Y ~~ z*Y
'
regression.fit <- sem(regression.model, sample.cov=correl, sample.nobs=100)
summary(regression.fit, rsquare=TRUE)

library(semPlot)
semPaths(regression.fit, "std", style="lisrel", 
         mar=c(6,1,3,1), edge.label.cex=.8, fade=F, theme = 'gray') 




# Regularization (Ridge Regression) -----------------------------------------------------------------

correl <- matrix(c(
  1, 0.62, 0.43, 0.56, 0.61,
  0.62, 1, 0.67, 0.55, 0.72,
  0.43, 0.67, 1, 0.67, 0.65,
  0.56, 0.55, 0.67, 1, 0.79,
  0.61, 0.72, 0.65, 0.79, 1),
  nrow=5,
  dimnames=list(c("Speaking","Vocabulary","Grammar","Writing","Reading"),
                c("Speaking","Vocabulary","Grammar","Writing","Reading")))

# Creating Simulated Data 
library(MASS)
set.seed(88)
mu <- rep(c(0), times = length(colnames(correl)))
simdat <- mvrnorm(n=100, mu=mu, Sigma=correl, empirical=TRUE)
colnames(simdat) <- c("Speaking","Vocabulary","Grammar","Writing","Reading")
sim.df <- as.data.frame(simdat)

# Apply the Ridge Regression
library(glmnet)
predictors <- as.matrix(sim.df[-1])
response_variable <- as.matrix(sim.df[1])
lambdas <- 10^seq(2, -2, by = -.1)
fit <- glmnet(predictors, response_variable, alpha = 0, lambda = lambdas)
# alpha: ridge = 0, lasso = 1, elastic net = between 0 and 1
plot(fit, xvar = "lambda", label = TRUE)

# Cross-validation
set.seed(88) # to get the same result
lambda_calc <- cv.glmnet(predictors, response_variable, alpha = 0, lambda = lambdas, grouped = FALSE)
# alpha: ridge = 0, lasso = 1, elastic net = between 0 and 1
plot(lambda_calc)

# 1000 cross-validations and take the average of lambda with minimum average error (lambda.min) 
LambdaCalcs <- NULL
for (i in 1:1000){
  lambda_calc <- cv.glmnet(predictors, response_variable, alpha = 0, lambda = lambdas, grouped = FALSE)
  LambdaCalcs <- cbind(LambdaCalcs, lambda_calc$lambda.min)
}
optlambda <- mean(LambdaCalcs)
# Get the standardized beta coefficients
options(scipen=999) # to turn off scientific notation (e notation)
coef(fit, s = optlambda)




# Table 3 (Dominance Analysis) -----------------------------------------------------------------

# Dominance Analysis with Correlation Matrix (using dominanceanalysis package)
library(dominanceanalysis)
lm.cov <- lmWithCov(Reading~ Vocabulary+Grammar+Writing+Speaking, correl)
da <- dominanceAnalysis(lm.cov)
print(da)

# Dominance Analysis with the Simulated Data
lm.out <- lm(Speaking~., sim.df)

# Relative Importance Analysis (=Dominance Analysis in Larson-Hall, 2016)
library(relaimpo)
calc.relimp(lm.out)

# CI and Statistical Test (Nimon & Oswald, 2013)
library(yhat)
regrOut <- calc.yhat(lm.out)
# Bootstrapping
library(boot)
set.seed(88)
boot.out <- boot(sim.df, boot.yhat, 1000, lmOut=lm.out, regrout0=regrOut)
# Summary Statistics of the Bootstrap Data
result <- booteval.yhat(regrOut, bty= "perc", boot.out)
# See the Results
regrOut
regrOut$PredictorMetrics[,12] #GenDom Weight
regrOut$OrderedPredictorMetrics[,12] #GenDom Order
regrOut$PairedDominanceMetrics #GenDom Comparisons
# (0:Xi dominates Xj / 1:Xj dominates Xi / 0.5:Dominance not established)
result #With 95%CIs
result$combCIpm[,12, drop=FALSE] #DA weights with upper/lower CI
result$combCIpmDiff[,"GenDom", drop=FALSE] #Comparisons of PVs




# Figure 2 -----------------------------------------------------------------

library(ggplot2)
dat <- data.frame(
  dw = regrOut$PredictorMetrics[,12][1:(length(regrOut$PredictorMetrics[,12])-1)],
  low = result$lowerCIpm[,12],
  up = result$upperCIpm[,12])
valnames <- rownames(dat)

ggplot(dat) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  geom_bar(aes(x=reorder(valnames,dw), y=dw), stat="identity", fill="skyblue", alpha=0.6) +
  geom_errorbar(aes(x=rownames(dat), ymin=low, ymax=up), width=0.3, colour="red", alpha=0.9, size=0.8) +
  geom_label(aes(x=valnames, y=dw, label = dw, vjust=0.45), position = position_dodge(width=0.9)) +
  labs(x = "Predictor Variables",
       y = "Dominance Weights",
       caption = paste0("R-squared = ",
                        regrOut$PredictorMetrics[,12][length(regrOut$PredictorMetrics[,12])],
                        " / ", 
                        "N = ", nrow(simdat)
                        )) + 
  ylim(c(-0.001, 0.33)) +
  geom_text(x = 2.5, y = 0.31, label = "*", size = 5) +
  geom_segment(x = 1, xend = 1, y = 0.3, yend = 0.29, size = 0.3) +
  geom_segment(x = 1, xend = 4, y = 0.3, yend = 0.3, size = 0.3) +
  geom_segment(x = 4, xend = 4, y = 0.3, yend = 0.29, size = 0.3) +
  geom_text(x = 2.0, y = 0.27, label = "*", size = 5) +
  geom_segment(x = 1, xend = 1, y = 0.26, yend = 0.25, size = 0.3) +
  geom_segment(x = 1, xend = 3, y = 0.26, yend = 0.26, size = 0.3) +
  geom_segment(x = 3, xend = 3, y = 0.26, yend = 0.25, size = 0.3) +
  coord_flip()




# Figure 4 -----------------------------------------------------------------

# Random Forest
library(randomForest)
forest <- randomForest(Speaking~., data=sim.df)
print(forest)
plot(forest)
forest$importance
varImpPlot(forest)

# Boruta
library(Boruta)
set.seed(88)
boruta <- Boruta(Speaking~., maxRuns = 200, data = sim.df, doTrace = 2)
print(boruta)
plot(boruta)
attStats(boruta)




# Reproduction (Goh et al., 2020) ----------------------------------------------------------

# Correlation Matrix
correl <- matrix(c(1, 0.67, 0.41, 0.02, -0.35, 0.35, 0.45,
                   0.67, 1, 0.18, 0.05, -0.28, 0.24, 0.31,
                   0.41, 0.18, 1, 0.21, -0.32, 0.34, 0.43,
                   0.02, 0.05, 0.21, 1, -0.03, 0.22, 0.47,
                   -0.35, -0.28, -0.32, -0.03, 1, -0.13, -0.22,
                   0.35, 0.24, 0.34, 0.22, -0.13, 1, 0.33,
                   0.45, 0.31, 0.43, 0.47, -0.22, 0.33, 1),
  				nrow=7,
  				dimnames=list(c("Score", "Wordcount", "CLI", "Commas", "Stopwords", "Linking", "WordsSentence"),
                c("Score", "Wordcount", "CLI", "Commas", "Stopwords", "Linking", "WordsSentence")))
                
# Multiple Regression Analysis
library(rpsychi)
multreg.second(Score~ Wordcount+CLI+Commas+Stopwords+Linking+WordsSentence, corr=correl, n=200)

# Dominance Analysis (using dominanceanalysis package)
library(dominanceanalysis)
lm.cov <- lmWithCov(Score ~ Wordcount+CLI+Commas+Stopwords+Linking+WordsSentence, correl)
da <- dominanceAnalysis(lm.cov)
print(da)

# Simulated Dataset Using Means, SDs, and Correlations
mu <- c(25.95, 228.03, 7.94, 0.76, 0.15, 0.007, 15.24)
stddev <- c(12.16, 97.47, 2.01, 0.53, 0.03, 0.007, 5.91)
corMat <- matrix(c(1, 0.67, 0.41, 0.02, -0.35, 0.35, 0.45,
                   0.67, 1, 0.18, 0.05, -0.28, 0.24, 0.31,
                   0.41, 0.18, 1, 0.21, -0.32, 0.34, 0.43,
                   0.02, 0.05, 0.21, 1, -0.03, 0.22, 0.47,
                   -0.35, -0.28, -0.32, -0.03, 1, -0.13, -0.22,
                   0.35, 0.24, 0.34, 0.22, -0.13, 1, 0.33,
                   0.45, 0.31, 0.43, 0.47, -0.22, 0.33, 1),
                 ncol = 7)
covMat <- stddev %*% t(stddev) * corMat
library(MASS)
set.seed(88)
dat1 <- mvrnorm(n = 200, mu = mu, Sigma = covMat, empirical = TRUE)
colnames(dat1) <- c("Score", "Wordcount", "CLI", "Commas", "Stopwords", "Linking", "WordsSentence")
dat1 <- as.data.frame(dat1)
colMeans(dat1) # Means
apply(dat1, 2, sd) # SDs
cor(dat1) # Correlation Matrix

# Multiple Regression Analysis
library(rpsychi)
multreg(Score~ Wordcount+CLI+Commas+Stopwords+Linking+WordsSentence, data=dat1)
lm.out <- lm(Score ~., dat1)

# Calculating p-values
z <- scale(dat1)
z <- data.frame(z)
summary(lm(Score~ ., z)) 

# Check the VIF (variance inflation factor)
library(car)
vif(lm.out) # VIF>10 shows multicollinearity. VIF<2 recommended.

# Relative importance analysis (=Dominance Analysis in Larson-Hall, 2016)
library(relaimpo)
calc.relimp(lm.out)

# Dominance Analysis（Same Result）
library(dominanceanalysis)
da.result <- dominanceAnalysis(lm.out)
print(da.result)

# CI and Statistical Test (Nimon & Oswald, 2013)
library(yhat)
regrOut <- calc.yhat(lm.out)
# Bootstrapping
library(boot)
set.seed(88)
boot.out <- boot(dat1, boot.yhat, 1000, lmOut=lm.out, regrout0=regrOut)
# Summary Statistics of the Bootstrap Data
result <- booteval.yhat(regrOut, bty= "perc", boot.out)
# See the Results
regrOut
regrOut$PredictorMetrics[,14] #GenDom Weight
regrOut$OrderedPredictorMetrics[,14] #GenDom Order
regrOut$PairedDominanceMetrics #GenDom Comparisons
# (0:Xi dominates Xj / 1:Xj dominates Xi / 0.5:Dominance not established)
result #With 95%CIs
result$combCIpm[,14, drop=FALSE] #DA weights with upper/lower CI
result$combCIpmDiff[,"GenDom", drop=FALSE] #Comparisons of PVs




# Figure 5 -----------------------------------------------------------------

library(ggplot2)
dat <- data.frame(
  dw = regrOut$PredictorMetrics[,14][1:(length(regrOut$PredictorMetrics[,14])-1)],
  low = result$lowerCIpm[,14],
  up = result$upperCIpm[,14])
valnames <- rownames(dat)

ggplot(dat) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  geom_bar(aes(x=reorder(valnames,dw), y=dw), stat="identity", fill="skyblue", alpha=0.6) +
  geom_errorbar(aes(x=rownames(dat), ymin=low, ymax=up), width=0.3, colour="red", alpha=0.9, size=0.8) +
  geom_label(aes(x=valnames, y=dw, label = dw, vjust=0.45), position = position_dodge(width=0.9)) +
  labs(x = "Predictor Variables",
       y = "Dominance Weights",
       caption = paste0("R-squared = ",
                        regrOut$PredictorMetrics[,14][length(regrOut$PredictorMetrics[,14])],
                        " / ", 
                        "N = ", nrow(dat1)
                        )) + 
  ylim(c(-0.001, 0.55)) +
  geom_text(x = 5.5, y = 0.43, label = "*", size = 5) +
  geom_segment(x = 6, xend = 6, y = 0.42, yend = 0.40, size = 0.3) +
  geom_segment(x = 6, xend = 5, y = 0.42, yend = 0.42, size = 0.3) +
  geom_segment(x = 5, xend = 5, y = 0.42, yend = 0.40, size = 0.3) +
  geom_text(x = 5, y = 0.46, label = "*", size = 5) +
  geom_segment(x = 6, xend = 6, y = 0.45, yend = 0.43, size = 0.3) +
  geom_segment(x = 6, xend = 4, y = 0.45, yend = 0.45, size = 0.3) +
  geom_segment(x = 4, xend = 4, y = 0.45, yend = 0.43, size = 0.3) + 
  geom_text(x = 4.5, y = 0.49, label = "*", size = 5) +
  geom_segment(x = 6, xend = 6, y = 0.48, yend = 0.46, size = 0.3) +
  geom_segment(x = 6, xend = 3, y = 0.48, yend = 0.48, size = 0.3) +
  geom_segment(x = 3, xend = 3, y = 0.48, yend = 0.46, size = 0.3) + 
  geom_text(x = 4, y = 0.52, label = "*", size = 5) +
  geom_segment(x = 6, xend = 6, y = 0.51, yend = 0.49, size = 0.3) +
  geom_segment(x = 6, xend = 2, y = 0.51, yend = 0.51, size = 0.3) +
  geom_segment(x = 2, xend = 2, y = 0.51, yend = 0.49, size = 0.3) + 
  geom_text(x = 3.5, y = 0.55, label = "*", size = 5) +
  geom_segment(x = 6, xend = 6, y = 0.54, yend = 0.52, size = 0.3) +
  geom_segment(x = 6, xend = 1, y = 0.54, yend = 0.54, size = 0.3) +
  geom_segment(x = 1, xend = 1, y = 0.54, yend = 0.52, size = 0.3) + 
  geom_text(x = 3, y = 0.23, label = "*", size = 5) +
  geom_segment(x = 5, xend = 5, y = 0.22, yend = 0.20, size = 0.3) +
  geom_segment(x = 5, xend = 1, y = 0.22, yend = 0.22, size = 0.3) +
  geom_segment(x = 1, xend = 1, y = 0.22, yend = 0.20, size = 0.3) +
  geom_text(x = 2.5, y = 0.18, label = "*", size = 5) +
  geom_segment(x = 4, xend = 4, y = 0.17, yend = 0.15, size = 0.3) +
  geom_segment(x = 4, xend = 1, y = 0.17, yend = 0.17, size = 0.3) +
  geom_segment(x = 1, xend = 1, y = 0.17, yend = 0.15, size = 0.3) + 
  coord_flip()




# Figure 6 -----------------------------------------------------------------

# Random Forest
library(randomForest)
forest <- randomForest(Score~., data=dat1)
print(forest)
plot(forest)
forest$importance
varImpPlot(forest)

# Boruta
library(Boruta)
set.seed(88)
boruta <- Boruta(Score~., maxRuns = 200, data = dat1, doTrace = 2)
print(boruta)
plot(boruta)
attStats(boruta)




# Cross-validation of Algorithms -----------------------------------------------------------------

library(caret)

# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# control <- trainControl(method="cv", number=10)
# control <- trainControl(method = "LOOCV")

# CART (Classification and Regression Trees)
set.seed(88)
fit.cart <- train(Score~., data=dat1, method="rpart", trControl=control,
                  tuneGrid = data.frame(cp = c(0.01, 0.05, 0.1)), preProcess = c('center', 'scale'))

# SVM (Support Vector Machine with Radial Basis Function)
set.seed(88)
fit.svm <- train(Score~., data=dat1, method="svmRadial", trControl=control, preProcess = c('center', 'scale'))

# kNN (k-Nearest Neighbors)
set.seed(88)
fit.knn <- train(Score~., data=dat1, method="knn", trControl=control, preProcess = c('center', 'scale'))

# Random Forest
set.seed(88)
fit.rf <- train(Score~., data=dat1, method="rf", trControl=control, preProcess = c('center', 'scale'))

# Neural Network
set.seed(88)
fit.nnet <- train(Score~., data = dat1, method="nnet", trControl = control,
                  preProcess = c('center', 'scale'), linout = TRUE, trace=FALSE)

# XgboostLinear（eXtreme Gradient Boosting）
set.seed(88)
fit.xgbLinear <- train(Score~., data=dat1, method="xgbLinear", trControl=control, preProcess = c('center', 'scale'))

# Collect Resamples
results <- resamples(list(CART=fit.cart, SVM=fit.svm, KNN=fit.knn, RF=fit.rf, NNET=fit.nnet, XGB=fit.xgbLinear))

# Summarize the Result
summary(results)

# Box and Whisker Plots to Compare Models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

# Dot plots of Accuracy with 95% Confidence Intervals
scales <- list(x=list(relation="free"), y=list(relation="free"))
dotplot(results, scales=scales)

# Density Plots of Accuracy
# to evaluate the overlap in the estimated behavior of algorithms
scales <- list(x=list(relation="free"), y=list(relation="free"))
densityplot(results, scales=scales, pch = "|", auto.key=T)

# Pair-wise Scatterplots of Predictions
splom(results)

# Difference in Model Predictions
options(scipen=999) # to turn off scientific notation (e notation)
diffs <- diff(results)
# Pair-wise Comparisons
summary(diffs)