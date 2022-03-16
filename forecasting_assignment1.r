# Preparing Dataset for the Time Series Analysis

library(readr)
library(tseries)
library(TSA)
library(dLagM)
library(forecast)
library(car)
library(x12)
library(kableExtra)


# Descriptive Analysis function
descriptive_analysis <- function(ts, object)
{
  plot(ts,
       ylab = c(paste0(toString(object))),
       main = c(paste0("Monthly Time Series Plot of ",toString(object))),
       type="o")
  points(y=ts,x=time(ts), pch=as.vector(season(ts)))
  
  par(mfrow=c(2,1))
  acf(ts,
      lag.max = 48,
      main = c(paste0("ACF plot of ",toString(object))))
  
  pacf(ts, 
       lag.max = 48,
       main = c(paste0("PACF plot of ",toString(object))))
  par(mfrow=c(1,1))
  
  print(adf.test(ts))
  print(pp.test(ts))
}

# Function for Decomposition
decom <- function(ts)
{
  decom.x12 = x12(ts)
  plot(decom.x12 , sa=TRUE , trend=TRUE)
  
  
  plotSeasFac(decom.x12)
  
  
  decomposition <- stl(ts, t.window=15, s.window="periodic", robust=TRUE)
  plot(decomposition)
}


ASX_data <- read_csv("ASX_data.csv")
head(ASX_data)


#### ASX ordinaries Price Index

ASX_TS <- ts(ASX_data$`ASX price`, start = c(2004,1), frequency = 12)
class(ASX_TS)
head(ASX_TS)


#### Gold Price (AUD)

gold_TS <- ts(ASX_data$`Gold price`, start = c(2004,1), frequency = 12)
class(gold_TS)
head(gold_TS)


#### Crude Oil (Brent, USD/bbl)

oil_TS <- ts(ASX_data$`Crude Oil (Brent)_USD/bbl`, start = c(2004,1), frequency = 12)
class(oil_TS)
head(oil_TS)


#### Copper (USD/tonne)

copper_TS <- ts(ASX_data$`Copper_USD/tonne`, start = c(2004,1), frequency = 12)
class(copper_TS)
head(copper_TS)


#### Whole Dataset

ASX_data_TS = ts(ASX_data, start = c(2004,1), frequency = 12)
class(ASX_data_TS)
head(ASX_data_TS)




ASX_TS_lambda = BoxCox.lambda(ASX_TS) 
ASX_TS_lambda

ASX_TSBC = ((ASX_TS^(ASX_TS_lambda)) - 1) / ASX_TS_lambda
descriptive_analysis(ASX_TSBC,"ASX ords price Index(Box-cox Transformed)")



gold_TS_lambda = BoxCox.lambda(gold_TS) 
gold_TS_lambda

gold_TSBC = ((gold_TS^(gold_TS_lambda)) - 1) / gold_TS_lambda
descriptive_analysis(gold_TSBC,"gold price(Box-cox Transformed) in AUD")



oil_TS_lambda = BoxCox.lambda(oil_TS) 
oil_TS_lambda

oil_TSBC = ((oil_TS^(oil_TS_lambda)) - 1) / oil_TS_lambda
descriptive_analysis(oil_TSBC, "Crude oil price(Box-cox Transformaed) in USD/bbl")


copper_TS_lambda = BoxCox.lambda(copper_TS) 
copper_TS_lambda

copper_TSBC = ((copper_TS^(copper_TS_lambda)) - 1) / copper_TS_lambda
descriptive_analysis(copper_TSBC,"Copper price(Box-cox Transformed) in USD/Tonne")


## Differencing

ASX_TSBC_diff = diff(ASX_TSBC)
descriptive_analysis(ASX_TSBC_diff,"ASX price(Bc Transformed-1st diff)")

gold_TSBC_diff = diff(gold_TSBC)
descriptive_analysis(gold_TSBC_diff,"gold price(BC Transformed-1st diff)")

oil_TSBC_diff = diff(oil_TSBC)
descriptive_analysis(oil_TSBC_diff, "Crude oil price(BC Transformed 1st diff)")

copper_TSBC_diff = diff(copper_TSBC)
descriptive_analysis(copper_TSBC_diff,"Copper price(BC Transformed-1st diff)")


## Decomposition
decom(ASX_TS)

decom(gold_TS)

decom(oil_TS)

decom(copper_TS)


# The most accurate and suitable distributed lag model

ASX_scaled= scale(ASX_data_TS)
plot(ASX_scaled, plot.type="s",
     col =	c("blue", "red","green","black"), 
     main= "ASX ordinaries and gold, oil, copper") 
legend("topleft",lty=1,col=c("blue","red","green","black"), c("ASX(Y)", "Gold(X1)","Oil(X2)","Copper(X3)"))

cor <- as.data.frame(cor(ASX_data_TS)) %>% round(3)
kbl(cor) %>%
  kable_paper()


## Finite Distributed Lag Models.
  
ASX_m = as.data.frame(ASX_data[,c(1,2,3,4)])
colnames(ASX_m) <- c("asx", "gold", "oil", "copper")
for ( i in 1:12){
  model1.1 = dlm(formula = asx ~ gold + copper, data = data.frame(ASX_m), q = i )
  cat("q = ", i, "AIC = ", AIC(model1.1$model), "BIC = ", BIC(model1.1$model),"\n")
}


model1 = dlm(formula = asx ~ gold + copper, data = data.frame(ASX_m), q = 12)
summary(model1$model)
checkresiduals(model1$model)
vif(model1$model)


## Polynomial Distributed Lag Model

finiteDLMauto(x = as.vector(ASX_m$copper), y = as.vector(ASX_m$asx), q.min = 1, q.max = 12, k.order = 1,
              model.type = "poly", error.type ="AIC", trace = TRUE)



model2 = polyDlm(x = as.vector(ASX_m$copper), y = as.vector(ASX_m$asx), q = 10, k = 1, show.beta = TRUE)
summary(model2$model)
checkresiduals(model2$model)
vif(model2$model)


## Koyck Distributed Lag Model

model3 = koyckDlm(x = as.vector(ASX_m$copper), y = as.vector(ASX_m$asx))
summary(model3,diagnostics = TRUE)
checkresiduals(model3$model)
vif(model3$model)


## Autoregressive Distributed Lag Model

for (i in 1:5)
{
  for(j in 1:5)
  { 
    model4 = ardlDlm(formula = asx ~ gold + copper, data = data.frame(ASX_m),p = i , q = j )
    cat("p = ", i, "q =" , j,"AIC = ", AIC(model4$model), "BIC = ", BIC(model4$model), "\n")
  }
}


# ardlDLM(1,5)
  
model4_1 = ardlDlm(formula = asx ~ gold + copper, data = data.frame(ASX_m), p = 1 , q = 5)
summary(model4_1)
checkresiduals(model4_1$model)
vif(model4_1$model)

# ardlDLM(2,5)
  
model4_2 = ardlDlm(formula = asx ~ gold + copper, data = data.frame(ASX_m), p = 2 , q = 5)
summary(model4_2)
checkresiduals(model4_2$model)
vif(model4_2$model)

# ardlDLM(3,5)
  
model4_3 = ardlDlm(formula = asx ~ gold + copper, data = data.frame(ASX_m), p = 3 , q = 5)
summary(model4_3)
checkresiduals(model4_3$model)
vif(model4_3$model)