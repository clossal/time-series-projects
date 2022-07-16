library(ggplot2)
library(astsa)
library(forecast)
#install.packages("MTS")
library(MTS)


MeanTemp <- na.omit(read.csv("MonthlyMeanTemp.csv", header = TRUE))
Rainfall <- na.omit(read.csv("MonthlyRainfall.csv", header = TRUE))
Sunshine <- na.omit(read.csv("MonthlySunshine.csv", header = TRUE))

######################################################################
#---------Q1---------------
class(MeanTemp)
MeanTemp.vec <- as.vector(t(as.matrix(MeanTemp[,2:13])))
head(MeanTemp.vec)
tail(MeanTemp.vec)
MeanTemp.ts <- ts(MeanTemp.vec, start = c(1884,1),end = c(2020,10),frequency = 12)
tail(MeanTemp.ts)
start(MeanTemp.ts);end(MeanTemp.ts);frequency(MeanTemp.ts)
plot(MeanTemp.ts, xlab = "time(years)", ylab = "MeanTemp")
title = ("MonthlyMeanTemp")
acf2(MeanTemp.ts, na.action = na.pass)
acf2(diff(diff(MeanTemp.ts),12))

class(Rainfall)
Rainfall.vec <- as.vector(t(as.matrix(Rainfall[,2:13])))
head(Rainfall.vec)
tail(Rainfall.vec)
Rainfall.ts <- ts(Rainfall.vec, start = c(1884,1),end = c(2020,10),frequency = 12)
tail(Rainfall.ts)
start(Rainfall.ts);end(Rainfall.ts);frequency(Rainfall.ts)
plot(Rainfall.ts, xlab = "time(years)", ylab = "Rainfall")
title = ("MonthlyRainFall")
acf2(Rainfall.ts, na.action = na.pass)
acf2(diff(diff(Rainfall.ts),12))

class(Sunshine)
Sunshine.vec <- as.vector(t(as.matrix(Sunshine[,2:13])))
head(Sunshine.vec)
tail(Sunshine.vec)
Sunshine.ts <- ts(Sunshine.vec, start = c(1884,1),end = c(2020,10),frequency = 12)
str(Sunshine.ts)
start(Sunshine.ts);end(Sunshine.ts);frequency(Sunshine.ts)
plot(Sunshine.ts, xlab = "time(years)", ylab = "Sunshine")
title = ("MonthlySunshine")
acf2(Sunshine.ts, na.action = na.pass)
acf2(diff(diff(Sunshine.ts),12))

######################################################################
#-------Q2------------

##For the MeanTemp case
autofitMeanTemp <- auto.arima(MeanTemp.ts)
forecast(auto.arima(MeanTemp.ts), h=2)
checkresiduals(autofitMeanTemp)

acf2(diff(MeanTemp.ts,12))

#Testimonial 1 
arima.fitMeanTemp1 <- sarima(MeanTemp.ts,1,1,2,1,1,0,12)
summary(arima.fitMeanTemp1$fit)

fit.ar <- arima(MeanTemp.ts,order=c(1,1,2),
                seasonal=list(order=c(1,1,0), period=12))

checkresiduals(fit.ar)
accuracy(fit.ar)

#par(mar=c(1,1,1,1))
tsdiag(fit.ar)

#predict the values 

plot(forecast(auto.arima(MeanTemp.ts),12),include = 18)
plot(forecast(auto.arima(MeanTemp.ts),12),include = 20)


##For the Rainfall case

autofitRainfall <- auto.arima(Rainfall.ts)
forecast(auto.arima(Rainfall.ts), h=2)

summary(autofitRainfall)
tsdiag(autofitRainfall)

acf2(diff(Rainfall.ts,12))
#testimonial 1
arima.fitRainfall1 <- sarima(Rainfall.ts,1,0,1,1,0,2,12)
arima.fitRainfall1$fit

#testimonial 2 
arima.fitRainfall2 <- sarima(Rainfall.ts,1,0,1,1,0,1,12)
arima.fitRainfall2$fit

##For the sunshine case
autofitSunshine <- auto.arima(Sunshine.ts)
forecast(auto.arima(Sunshine.ts), h=2)

summary(autofitSunshine)
tsdiag(autofitSunshine)
dev.off()
checkresiduals(autofitSunshine)
#Testmonial 1 

temp.diff.EA <- diff(MeanTemp.ts, lag=12, differences=1)

acf(na.omit(temp.diff.EA), lag.max=36,xlab="Lag", ylab="ACF",
    main="Differenced autocorrelation plot")

pacf(na.omit(temp.diff.EA), lag.max=36,xlab="Lag", ylab="ACF",
    main="Differenced autocorrelation plot")



######################################################################
#Q3 ------ We are going to apply auto.arima because with rolling window,
# it is expected that p, q, d value will change quite drastically,
# so auto.arima can do this for us.
# After that we will extract error metrics with our forecast value

prediction.roll = function(data, M = M, A=12){
  n = length(data)
  N = ceiling((n-M)/A) + 1
  d = list()
  prediction_error = list()
  fit = list()
  forecast <- list()
  predictions <- list()
  Observations <- list()
  for (i in 1:N){
    range = c(((i-1)*A+1):(min((i-1)*A+M,n)))
    name <- paste0("rolling_Window", sep="_", (i-1)*A+1,
                   sep="_", min((i-1)*A+M))
    d[[name]] <- ts(data[range], frequency = 12)
    set.seed(1)
    m <- auto.arima(data[range], D=12, stepwise=FALSE,
                    seasonal=TRUE,approximation=FALSE)
    fit[[name]] <- m
    fcst = predict(object=m,n.ahead=1)$pred
    forecast[[name]] <- fcst
    so <- data[i+1]
    Observations[[name]] <- so
    prediction_error[[name]] = mean((as.numeric(fcst)-as.numeric(so))^2, na.rm=T)
    print(prediction_error)
    
  }
  return(list(fit = fit, pred_Squared_Error = prediction_error,
              forecasts = forecast, real_observations=Observations,
              dataframe = d))
}

# split data in 2 dataframe
split_data <- function(df, train_window, test_window){
  n = length(df)
  
  # divide data in train/validation: 75/25 ratio
  train_length <- n*0.75
  train <- ts(df[1:train_length], frequency=12)
  test <- ts(df[n-train_length:(n)], frequency = 12)
  
  # apply rolling window function
  train.fit <- prediction.roll(train, M=train_window)
  test.fit <- prediction.roll(test, M=test_window)
  
  # apply another function so we have nice table for errors
  train_error <- list_to_df(train.fit)
  test_error <- list_to_df(test.fit)
  return(list(train_data = train, train = train.fit,
              validation_data = test, validation = test.fit,
              train_error= train_error, validation_error=test_error))
}


# function to get dataframe for errors,
# observations and forecast values
list_to_df <- function(df){
  
  # extract list from dataframe
  l1 <- df$pred_Squared_Error
  l2 <- df$forecasts
  l3 <- df$real_observations
  
  # bind them in rows
  d1 <- do.call(rbind, l1)
  d2 <- do.call(rbind, l2)
  d3 <- do.call(rbind, l3)
  
  # proper names for columns
  colnames(d1) <- "Error"
  colnames(d2) <- "forecast"
  colnames(d3) <- "Real Observations"
  
  # create dataframe with row.names and column
  d <- data.frame(rolling_window = row.names(d1), d1)
  d_2 <- data.frame(rolling_window = row.names(d2), d2)
  d_3 <- data.frame(rolling_window = row.names(d3), d3)
  
  # we don't need row.names anymore
  rownames(d) <- NULL
  rownames(d_2) <- NULL
  rownames(d_3) <- NULL
  
  # merge all three dataframe
  dd<- merge(d, d_2, on="rolling_window")
  ddd <- merge(dd, d_3, on="rolling_window")
  final <- ddd[order(ddd[, 2], decreasing = F), ]
  
  # only return first row, which has least square error
  return(head(final, 1))
 
}

# apply on datasets
data_list = list(Sunshine.ts, MeanTemp.ts, Rainfall.ts)

# give names accordingly
names(data_list) <- c("Sunshine", "Temperature", "Rainfall")

# apply function on the list
datasets <- lapply(data_list, split_data, 180, 120)

# lets apply different window lengths
datasets1 <- lapply(data_list, split_data, 120, 72)

# apply another window
datasets2 <- lapply(data_list, split_data, 144, 96)

# from validation data, we select the data which has 
# least squared error and then forecast Nov and Dec Values

# look at Sunshine data validation error
datasets$Sunshine$validation_error
datasets1$Sunshine$validation_error
datasets2$Sunshine$validation_error

datasets$Sunshine$train_error
datasets1$Sunshine$train_error
datasets2$Sunshine$train_error
# check residuals
checkresiduals(datasets$Sunshine$validation$fit$rolling_Window_301_420)

# convert to ts object
sunshine.forecast <- ts(datasets$Sunshine$validation$dataframe$rolling_Window_301_420,
                        start = c(2008,1),end = c(2020,10), frequency = 12)

# comparison of forecast with question-2
forecast(auto.arima(sunshine.forecast), h=2)
forecast(auto.arima(Sunshine.ts), h=2) # forecast from q-2

# Validation Error on Temperature data
datasets$Temperature$validation_error
datasets1$Temperature$validation_error
datasets2$Temperature$validation_error

datasets$Temperature$train_error
datasets1$Temperature$train_error
datasets2$Temperature$train_error
# check residuals
checkresiduals(datasets1$Temperature$validation$fit$rolling_Window_145_216)

# create time series object
temp.forecast <- ts(datasets1$Temperature$validation$dataframe$rolling_Window_145_216,
        start = c(2008,1),end = c(2020,10), frequency = 12)

forecast(auto.arima(temp.forecast), h=2)
forecast(auto.arima(MeanTemp.ts), h=2) # forecast from q-2

# Validation Error on Rainfall data
datasets$Rainfall$validation_error
datasets1$Rainfall$validation_error
datasets2$Rainfall$validation_error

datasets$Rainfall$train_error
datasets1$Rainfall$train_error
datasets2$Rainfall$train_error
# the least error is in 96 window, lets check 
# residuals and make forecast
checkresiduals(datasets$Rainfall$validation$fit$rolling_Window_229_348)

Rainfall.forecast <- ts(datasets2$Rainfall$validation$dataframe$rolling_Window_229_324,
                        start = c(2008,1),end = c(2020,10), frequency = 12)

forecast(auto.arima(Rainfall.forecast), h=2)
forecast(auto.arima(Rainfall.ts), h=2) # forecast from q-2

# It seems that in most of cases, 120 window length is good for validation data
# except, temperatue data, whose window length is 72
# and 180 window length is working excellently on train data
# so we have made our forecasts with best validation window

######################################################################
#Q-4 ------ VARIMA model


all_data <- cbind(Sunshine.ts, MeanTemp.ts, Rainfall.ts)

# lets analyze correlation of all time series
corrplot::corrplot(cor(all_data))


# It is clear from the correlation plot that,
# sunshine and temperature has strong correlation with each
# other which is understanable, because increase in sunshine also
# increases temperature

# Now build VARMA
mod <- VARMA(all_data)
varma_predictions <- VARMApred(mod, h = 2, orig = 0)

p_varma <- ts(varma_predictions$pred, start = c(2020, 11),
              end = c(2020, 12), frequency = 12)

# lets also apply VARMA on our rolling window datasets
window_data <- cbind(sunshine.forecast, temp.forecast, Rainfall.forecast)
var_model2 <- VARMA(window_data)

var_pred2 <- VARMApred(var_model2, h = 2, orig = 0)
var_pred2$pred
