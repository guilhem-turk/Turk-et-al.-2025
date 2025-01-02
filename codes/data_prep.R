library(readxl)
library(ggpubr)
library(scales)
library(dplyr)
library(tidyr)
library(stats)
library(gridExtra)
library(purrr)
library(viridis)
library(imputeTS)
library(lubridate)
library(data.table)
library(zoo)


isotopes <- read_excel("~/MUSES/Isotope_model/data/Precip&isotopes_Belvaux_2016-2023_update_02022023.xlsx", 
                       col_types = c("numeric", "date", "date", 
                                     "text", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric"))

date_corrections <- read_excel("~/MUSES/Isotope_model/data/Belval_rain isotopes.xlsx", 
                               sheet = "isotope samples", col_types = c("date", 
                                                                        "date", "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "text", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "text"))

df_corr <- data.frame(
  date = date_corrections$`Original date`,
  time = date_corrections$`original time`,
  duration = date_corrections$`original duration (h)`
)

df_corr <- df_corr %>% drop_na(date)
df <- isotopes %>% drop_na(`Sampling start`)

df$d_excess <- df$dD - 8*df$d18O
df$date <- date(df$`Sampling start`)
df$time <- hour(df$`Sampling start`)


df$duration <- ceiling(as.numeric(df$`Sampling end` - df$`Sampling start`)/3600)


df$date[23:1465] <- df_corr$date
df$time[23:1465] <- df_corr$time
df$duration[23:1465] <- df_corr$duration

df[df$date >= "2022-03-27" & df$date <= "2022-10-30", 11] <- df[df$date >= "2022-03-27" & df$date <= "2022-10-30", 11] -1

start <- min(df$`Sampling start`)
end <- max(df$`Sampling end`)

df_helper <- data.frame(timestamp = seq.POSIXt(start, end, by="hour"))


gwt <- read_excel("~/MUSES/Isotope_model/data/GWL_1954-2022.xlsx", 
                  col_types = c("date", "text", "numeric", 
                                "text", "text", "skip", "skip", "skip", 
                                "skip", "skip", "skip", "skip", "skip", 
                                "skip"))


rain <- read_excel("~/MUSES/Isotope_model/data/CR10X_rainfall_Belvaux_LIST_2016-2021_Guilhem.xlsx", 
                   col_types = c("date", "numeric"))

rain.oberkorn <- read_excel("~/MUSES/Isotope_model/data/local_parameters.xlsx", 
                            sheet = "rain_oberkorn", col_types = c("date", "numeric"))

rain$timestamp <- round(rain$timestamp, units = "hours")

df.rain <- merge(rain, df_helper, by="timestamp", all.y=T)

df.rain <- merge(df.rain, rain.oberkorn, by="timestamp", all.x=T)

df.rain$date <- as.Date(df.rain$timestamp)


dff.rain <- df.rain %>%
  mutate(date=as.Date(timestamp)) %>%
  group_by(date) %>%
  reframe(rain=sum(rain, na.rm=T),
          rain_alt=sum(value2, na.rm=T))


ggplot(dff.rain) + geom_point(aes(x=rain, y=rain_alt)) + theme_bw()


ddf.rain <- dff.rain %>% 
  reframe(rain=rollapply(rain, width = 3, FUN = sum, fill = NA, align="center"),
          rain_alt=rollapply(rain_alt, width = 3, FUN = sum, fill = NA, align="center"))

ddf.rain$date <- dff.rain$date

ggplot(ddf.rain) + geom_point(aes(x=rain, y=rain_alt)) + theme_bw()


ddf.rain <- ddf.rain %>% mutate(diff=rain_alt-rain)

ggplot(ddf.rain) + geom_point(aes(x=date, y=diff)) + 
  geom_point(data=subset(ddf.rain, rain<=0.5), aes(x=date, y=diff), color="red") +
  geom_point(data=subset(ddf.rain, rain!=0 & rain_alt/rain >=2), aes(x=date, y=diff), color="blue") + theme_bw()


ddf.rain$index <- 0
ddf.rain$index[ddf.rain$rain<=0.5 & ddf.rain$diff >= 1] <- 1

ddf.rain <- dff.rain %>% 
  reframe(index=rollapply(index, width = 3, FUN = sum, fill = NA, align="center"))


df.rain <- merge(df.rain, ddf.rain[,c(3,5)], by="date")

df.rain$rain[df.rain$index==1] <- df.rain$value2[df.rain$index==1]


dff.rain <- df.rain %>%
  mutate(date=as.Date(timestamp)) %>%
  group_by(date) %>%
  reframe(rain=sum(rain, na.rm=T),
          rain_alt=sum(value2, na.rm=T))


ggplot(dff.rain) + geom_point(aes(x=rain, y=rain_alt)) + theme_bw()





temp <- read_excel("~/MUSES/Isotope_model/data/local_parameters.xlsx", 
                   sheet = "temp")

rh <- read_excel("~/MUSES/Isotope_model/data/local_parameters.xlsx", 
                 sheet = "rh", col_types = c("date", "numeric"))

gph.roodt <- read_excel("~/MUSES/Isotope_model/data/local_parameters.xlsx", 
                        sheet = "gph", col_types = c("date", "skip", "skip", "numeric"))

gph.findel <- read_excel("~/MUSES/Isotope_model/data/local_parameters.xlsx", 
                         sheet = "gph2")

gph.findel <- gph.findel[!duplicated(gph.findel$Timestamp),]

gir <- read_excel("~/MUSES/Isotope_model/data/local_parameters.xlsx", 
                  sheet = "gir")




gph2 <- merge(gph.findel , gph.roodt, by="Timestamp", all=T)


gph2$gph <- gph2$Value
gph2$gph[is.na(gph2$gph)==T] <- gph2$sea_level[is.na(gph2$gph)==T]

gph <- gph2[,c(1,4)]

local_param <- merge(temp, rh, by="Timestamp", all=T)
local_param <- merge(local_param, gph, by="Timestamp", all=T)
local_param <- merge(local_param, gir, by="Timestamp", all=T)

colnames(local_param) <- c("timestamp","temp","rh","gph","gir")


local_param <- merge(local_param, df.rain[,c(1:3)], by="timestamp", all.y=T)

local_param$time <- hour(local_param$timestamp)



df <- merge(df, gwt[,c(1,4)], by="date")
local_param <- merge(local_param, gwt[,c(1,4)], by="date")




write.table(df, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\data\\clean_isotopes.csv', sep = ",", quote = FALSE)
write.table(local_param, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\data\\clean_meteo.csv', sep = ",", quote = FALSE)






