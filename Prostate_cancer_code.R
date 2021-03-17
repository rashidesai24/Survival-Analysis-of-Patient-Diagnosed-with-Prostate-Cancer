library(readr)
cancer <- read_csv("/Users/rashidesai/Desktop/training_data.csv")
nrow(cancer) #15,385 rows of 33 variables

str(cancer)

# NAs plot

summary(cancer$age) # Mean age = 76.9; Median = 78
hist(cancer$age, col = 'gold', xlab = "Age", ylab = "Frequency", main = "Plot of Age Distribution" )

#Grouping observations of the variables into different categories
# cancer$age<- ifelse(cancer$age <= 50, 1, ifelse(cancer$age <= 70, 2, ifelse(cancer$age<=90,3,4)))

cancer$age[cancer$age<=60] <- 1
cancer$age[cancer$age>60 & cancer$age<=70] <- 2
cancer$age[cancer$age>70 & cancer$age<=80] <- 3
cancer$age[cancer$age>80 & cancer$age<=90] <- 4
cancer$age[cancer$age>90] <- 5
summary(cancer$age)
cancer$age<-as.factor(cancer$age)
plot(cancer$age, col = "gold", xlab= "Age", ylab = "Frequency", main = "Distribution of Age")
levels(cancer$age)
cancer$age <- within(cancer, age <- relevel(age, ref = 3))

cancer$race<-as.factor(cancer$race)
levels(cancer$race)
cancer$race <- within(cancer, race <- relevel(race, ref = 4))

cancer$tea<- ifelse(cancer$tea <= 3, 1, ifelse(cancer$tea<=7,2,3))
summary(cancer$tea)
cancer$tea<-as.factor(cancer$tea)
plot(cancer$tea)

# Converting data types
cols<-c("stage","race","family_history","first_degree_history","previous_cancer","smoker"
        ,"side","rd_thrpy","h_thrpy","chm_thrpy","cry_thrpy","brch_thrpy","rad_rem","multi_thrpy","survival_7_years")
cancer[,cols] <-  data.frame(apply(cancer[cols],2, as.factor))
cancer$t_score<-as.factor(cancer$t_score)
cancer$n_score<-as.factor(cancer$n_score)
cancer$m_score<-as.factor(cancer$m_score)
cancer$side<-as.factor(cancer$side)
cancer$survival_1_year<-as.factor(cancer$survival_1_year)

str(cancer)

#-------------------------------------------------------------------------------------------------------------------------------------------------
# UNIVARIATE ANALYSIS
#-------------------------------------------------------------------------------------------------------------------------------------------------

# stacked bar chart
library(ggplot2)

# Gleason Score: a measurement of how abnormal the cancer cells look compared to normal cells
ggplot(cancer, aes(x = gleason_score , fill = survival_7_years)) + geom_bar(position = "stack")
# We can see that patients with high gleason score do not survive after 7 years of diagnosis
ggplot(cancer, aes(x = tumor_diagnosis , fill = survival_7_years)) + geom_bar(position = "stack")
# Patients who do not survive after 7 years had a bigger tumor compared to people who survived
ggplot(cancer, aes(x = psa_diagnosis , fill = survival_7_years)) + geom_bar(position = "stack")
ggplot(cancer, aes(x = race , fill = survival_7_years)) + geom_bar(position = "stack")
# Most of the patients are under Race 4. Therefore, this level can be used as reference level.
cancer$gleason_score<-ifelse(cancer$gleason_score<=6,1, 
                             ifelse(cancer$gleason_score==7,2,
                                    ifelse(cancer$gleason_score==8,3,
                                           ifelse(cancer$gleason_score<=10,4,5))))
cancer$gleason_score<-as.factor(cancer$gleason_score)
str(cancer)


#-------------------------------------------------------------------------------------------------------------------------------------
# Subject-matter expertise, create dervied variables and conceptual model
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#Creating BMI Variable
BMI = function(height,weight){
  return(0.45455*weight/(.0254*height)^2)}

cancer$BMI <- BMI(cancer$height,cancer$weight)

cancer$obese <- NA
cancer$obese[cancer$BMI<=25] <- 0
cancer$obese[cancer$BMI>25 & cancer$BMI<30] <- 1
cancer$obese[cancer$BMI>=30] <- 2
hist(cancer$obese, col = "darkgreen", xlab = "Obese", ylab = "Frequency", main = "Distribution of obesity")

cancer$obese <- as.factor(cancer$obese)
cancer<-cancer[!is.na(cancer$obese),]
table(cancer$obese, cancer$survival_7_years)

ggplot(cancer, aes(x = previous_cancer , fill = survival_7_years)) + geom_bar(position = "stack")


# CHI-SQUARED TEST

table_tscore = table(cancer$t_score, cancer$survival_7_years)
chisq.test(table_tscore)

table_nscore = table(cancer$n_score, cancer$survival_7_years)
chisq.test(table_nscore)

table_mscore = table(cancer$m_score, cancer$survival_7_years)
chisq.test(table_mscore)

table_race = table(cancer$race, cancer$survival_7_years)
chisq.test(table_race)

# Not significant: 0.2831
table_previous = table(cancer$previous_cancer, cancer$survival_7_years)
chisq.test(table_previous)

# Not significant: 0.1364
table_smoker = table(cancer$smoker, cancer$survival_7_years)
chisq.test(table_smoker)

table_rd = table(cancer$rd_thrpy, cancer$survival_7_years)
chisq.test(table_rd)

#Not significant = 0.1049
table_h = table(cancer$h_thrpy, cancer$survival_7_years)
chisq.test(table_h)

table_chm = table(cancer$chm_thrpy, cancer$survival_7_years)
chisq.test(table_chm)

table_cry = table(cancer$cry_thrpy, cancer$survival_7_years)
chisq.test(table_cry)

#Not significant = 0.03433
table_brch = table(cancer$brch_thrpy, cancer$survival_7_years)
chisq.test(table_brch)

#Not significant = 0.316
table_rad = table(cancer$rad_rem, cancer$survival_7_years)
chisq.test(table_rad)

table_multi = table(cancer$multi_thrpy, cancer$survival_7_years)
chisq.test(table_multi)

table_1ys = table(cancer$survival_1_year, cancer$survival_7_years)
chisq.test(table_1ys)

# Not significant = 0.748
table_side = table(cancer$side, cancer$survival_7_years)
chisq.test(table_side)

table_stage = table(cancer$stage, cancer$survival_7_years)
chisq.test(table_stage)

# Not significant = 0.625
table_symptoms = table(cancer$symptoms, cancer$survival_7_years)
chisq.test(table_symptoms)

# Not significant = 0.04205 if 0.01 but significant if 0.05
levels(cancer$family_history)
chisq.test(table(cancer$obese, cancer$survival_7_years))


# ----------------------------------------------------------------------------------------------------------
# CORRELATION PLOTS
# ----------------------------------------------------------------------------------------------------------

library(car) # for detailed correlation plot 
library(corrplot) # for correlation plot
library(Hmisc) # for correlation test of multiple variables 
cnum <- cancer[,c( "tumor_diagnosis", "tumor_1_year", "psa_diagnosis", 
                   "psa_1_year")]
cormat <- cor(cnum) # Select only numeric variables
corrplot(cormat, method="circle", addCoef.col="black") # With correlation

# ----------------------------------------------------------------------------------------------------------
# LOGISTIC REGRESSION MODEL 
# ----------------------------------------------------------------------------------------------------------
cancer<-cancer[!is.na(cancer$gleason_score),]
cancer<-cancer[!is.na(cancer$t_score),]
cancer<-cancer[!is.na(cancer$n_score),]
cancer<-cancer[!is.na(cancer$m_score),]
cancer<-cancer[!is.na(cancer$obese),]
cancer<-cancer[!is.na(cancer$smoker),]
cancer<-cancer[!is.na(cancer$age),]
cancer<-cancer[!is.na(cancer$race),]
cancer<-cancer[!is.na(cancer$family_history),]
cancer<-cancer[!is.na(cancer$first_degree_history),]
# cancer<-cancer[!is.na(cancer$previous_cancer),]
# cancer<-cancer[!is.na(cancer$smoker),]
# cancer<-cancer[!is.na(cancer$side),]
# cancer<-cancer[!is.na(cancer$tumor_diagnosis),]
cancer<-cancer[!is.na(cancer$psa_diagnosis),]
# cancer<-cancer[!is.na(cancer$psa_1_year),]
cancer<-cancer[!is.na(cancer$symptoms),]
cancer<-cancer[!is.na(cancer$rd_thrpy),]
# cancer<-cancer[!is.na(cancer$h_thrpy),]
cancer<-cancer[!is.na(cancer$chm_thrpy),]
cancer<-cancer[!is.na(cancer$cry_thrpy),]
# cancer<-cancer[!is.na(cancer$brch_thrpy),]
# cancer<-cancer[!is.na(cancer$rad_rem),]
cancer<-cancer[!is.na(cancer$multi_thrpy),]
# cancer<-cancer[!is.na(cancer$survival_1_year),]
cancer<-cancer[!is.na(cancer$survival_7_years),]
cancer<-cancer[!is.na(cancer$tumor_1_year),]
cancer<-cancer[!is.na(cancer$tea),]


str(cancer)
set.seed(2)
indx = sample(2, nrow(cancer), replace = T, prob = c(0.7, 0.3))
training = cancer[indx == 1,]
testing = cancer[indx == 2,]

# IMPORTANT PREDICTORS
# l1<-glm(survival_7_years~gleason_score+t_score+n_score+m_score+stage+age+tumor_diagnosis+psa_diagnosis+rd_thrpy+chm_thrpy+cry_thrpy+multi_thrpy,data = training,family =binomial(link="logit"))


l1 <- glm(survival_7_years~gleason_score+age+race+stage+
            first_degree_history+family_history+t_score+
            n_score+m_score+smoker+tea+obese+rd_thrpy+
            chm_thrpy+cry_thrpy+brch_thrpy+multi_thrpy+
            tumor_1_year+psa_diagnosis, 
          data=training, family="binomial")
summary(l1)
exp(coef(l1))

predTrain<-predict(l1, newdata=training,type = "response")
results <- ifelse(predTrain > 0.5,1,0)
results <-as.factor(results)
levels(results) <- c("0", "1")

library(caret)
summary(results)
library(e1071)
confusionMatrix(results, training$survival_7_years, positive = "1")


#_---------------------------------------------------------------------------------------------------------------------

library(readr)
test <- read.csv("/Users/rashidesai/Desktop/Rashi_score.csv")
nrow(test) #11531
str(test)
summary(test)

test$age[test$age<=60] <- 1
test$age[test$age>60 & test$age<=70] <- 2
test$age[test$age>70 & test$age<=80] <- 3
test$age[test$age>80 & test$age<=90] <- 4
test$age[test$age>90] <- 5
test$tea<- ifelse(test$tea==0,0,ifelse(test$tea<2,1,ifelse(test$tea<4,2,3)))

test$survival_1_year<-as.factor(test$survival_1_year)
cols<-c("stage","race","family_history","first_degree_history","previous_cancer","smoker"
        ,"side","rd_thrpy","h_thrpy","chm_thrpy","cry_thrpy","brch_thrpy","rad_rem","multi_thrpy","survival_7_years")
test[,cols] <-  data.frame(apply(test[cols],2, as.factor))

test$gleason_score<-ifelse(test$gleason_score<=6,1,
                           ifelse(test$gleason_score==7,2,ifelse(test$gleason_score==8,3,ifelse(test$gleason_score<=10,4,5))))

BMI = function(height,weight){
  return(0.45455*weight/(.0254*height)^2)}

test$BMI <- BMI(test$height,test$weight)

test$obese <- NA
test$obese[test$BMI<=25] <- 0
test$obese[test$BMI>25 & test$BMI<30] <- 1
test$obese[test$BMI>=30] <- 2

test$race<-as.factor(test$race)
levels(test$race)
test$race <- within(test, race <- relevel(race, ref = 4))

test$t_score<-as.factor(test$t_score)
test$n_score<-as.factor(test$n_score)
test$m_score<-as.factor(test$m_score)
test$age<-as.factor(test$age)
test$tea<-as.factor(test$tea)
test$symptoms<-as.factor(test$symptoms)
test$gleason_score<-as.factor(test$gleason_score)
test$obese<-as.factor(test$obese)


write.csv(test, "/Users/rashidesai/Desktop/test.csv")

predTest<-predict(object=l1,newdata = test, type="response")
results_test<-ifelse(predTest > 0.5,1,0)
results_test<-as.factor(results_test)
levels(results_test)<-c("0","1")

library(caret)
summary(results)
library(e1071)
confusionMatrix(results, test$survival_7_years, positive = "1")

summary(results_test)
summary(results) 

test$results<-results_test
levels(results)
levels(cancer$survival_7_years)
typeof(results)
typeof(test$survival_7_years)

table(results$survival_7_years)
View(results)

write.csv(results_test, "/Users/rashidesai/Desktop/round2.csv")










