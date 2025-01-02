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
library(gghighlight)
library(hrbrthemes)
library(MASS)
library(moments)
library(caret)
library(broom)
library(lme4)
library(data.table)
library(zoo)
library(sf)     # For spatial data handling
library(elevatr)



gwt <- read_excel("~/MUSES/Isotope_model/data/GWL_1954-2022.xlsx", 
                  col_types = c("date", "text", "numeric", 
                                "text", "text", "skip", "skip", "skip", 
                                "skip", "skip", "skip", "skip", "skip", 
                                "skip"))

df_meteo <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/output_data/monthly_data.csv", header = TRUE, sep = ",", dec = ".")

dff_param <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/output_data/model_param.csv", header = TRUE, sep = ",", dec = ".")




gwt$gwt <- factor(gwt$gwt, levels=c("West","North West","South West","High CE","Low CE","North","South","East"))

levels(gwt$gwt) <- list(W = "West", NW = "North West", SW = "South West", 
                        HCE = "High CE", LCE = "Low CE",
                        N= "North", S = "South", E = "East")



## extrapolation to GNIP stations



df_gnip <- read_excel("~/MUSES/Isotope_model/data/gnip_stations.xlsx")



dff_gnip <- df_gnip[,c(9,12:15,20,35,39)]
dff_gnip$date <- as.Date(substr(dff_gnip$`Sample Date`, 1, 10), format = "%Y-%m-%d")

dff_gnip <- pivot_wider(dff_gnip, names_from = "Measurand Symbol", values_from = "Measurand Amount")
names(dff_gnip) <- c("region","station","lat","long","elevation","timestamp","date","P","dO")

dff_gnip$region <- factor(dff_gnip$region)
levels(dff_gnip$region)

dff_gnip <- dff_gnip %>% mutate(year=year(date), month=month(date)) 



dff_meteo <- drop_na(df_meteo) %>% mutate(year=year(date), month=month(date))

dff_meteo <- rename(dff_meteo, station = site)



df_gwt <- merge(gwt[,c(1,4)], dff_param)

dff_gwt <- df_gwt %>%
  mutate(year=year(date), month=month(date)) %>%
  group_by(month, year) %>%
  reframe(date = as.Date(min(date)),
          a = mean(a),
          b = mean(b))


dff <- merge(dff_meteo, dff_gwt[,-c(3)], by=c("year","month"))

dff$temp <- dff$t2m - 273.15

dff <- dff %>% mutate(dO.rec = a + b*temp)


dff <- dff %>%
  group_by(station) %>%  
  distinct(date, .keep_all = TRUE)





# Download DEM data for Europe with specified resolution


europe_dem <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/output_data/meteo&elevation_raster.csv", header = TRUE, sep = ",", dec = ".")
europe_dem <- europe_dem %>% distinct(lat, long, .keep_all = TRUE)


df_europe_dem <- rename(europe_dem[,-c(3:6)], avg_alt = elevation)

df_europe_dem <- df_europe_dem %>%
  mutate(long = round(long * 2) / 2,
         lat = round(lat * 2) / 2)

dff_gnip <- dff_gnip %>%
  mutate(long = round(long * 2) / 2,
         lat = round(lat * 2) / 2)


ddf_gnip <- merge(dff_gnip, df_europe_dem, by=c("long", "lat"), all=T)




data <- merge(ddf_gnip[,-c(3,6)], dff[,-c(4:7)], by=c("station", "year", "month"))

data <- data[is.na(data$dO)==FALSE,]

data <- data %>% 
  group_by(station) %>%
  mutate(N = n())



data <- arrange(data, station, year, month)


data$avg_alt[data$avg_alt< (-200)] <- NA
data$avg_alt[is.na(data$avg_alt)==T] <- data$elevation[is.na(data$avg_alt)==T]




data.mm <- data %>%
  group_by(station, month, lat, long, avg_alt) %>%
  reframe(dO = mean(dO, na.rm=T),
          temp = mean(t2m, na.rm=T),
          dO.rec = mean(dO.rec, na.rm=T),
          N = n())


data.mm <- data.mm %>%
  group_by(station) %>%
  mutate(n = n()) %>%
  ungroup()


data.mm <- subset(data.mm, n == 12)

n <- length(unique(data.mm$station))




data.mm <- data.mm %>% 
  mutate(index = row_number(),
         type = ifelse(index<max(index)*3/4, "training","testing"))



train_set <- subset(data.mm, type == "training")
test_set <- subset(data.mm, type == "testing")

length(unique(test_set$station))



mlrm <- lm(dO~lat+long+avg_alt+dO.rec, train_set)
stats <- summary(mlrm)
stats







coefs <- stats$coefficients[,1]

test_set <- test_set %>% mutate(fits = coefs[1] + coefs[2]*lat + coefs[3]*long + 
                                  coefs[4] * avg_alt + coefs[5]*dO.rec,
                                residuals = dO - fits)




model_eval.a <- test_set %>% 
  group_by(station, lat, long) %>%
  summarise(RMSE = sqrt(sum(residuals^2)/n()),
            mean = round(mean(dO, na.rm=T), 2),
            N = n())



data.mm_sd <- data %>%
  group_by(station) %>%
  reframe(sd = sd(dO, na.rm=T))


model_eval.a <- merge(model_eval.a, data.mm_sd)




model_eval.b <- test_set %>% 
  summarise(tss = sum((dO - mean(dO))^2),
            rss = sum(residuals^2),
            rsq = 1 - (rss / tss),
            RMSE = sqrt(sum(residuals^2)/n()),
            N = n())




sf_eval <- st_as_sf(model_eval.a, coords = c("long", "lat"), crs = 4326)



shape_europe <- st_read("~/MUSES/Isotope_model/data/NUTS_RG_20M_2021_4326.shp")
shape_europe <- subset(shape_europe, LEVL_CODE==0)





spatial_model <- ggplot() + 
  geom_sf(data = shape_europe, colour = "grey30", fill = "white") +
  geom_sf(data = subset(sf_eval, RMSE >= sd & N == 12), aes(), color="red", size=1, shape=4) + 
  geom_sf(data = subset(sf_eval, RMSE < sd & N == 12), aes(color=RMSE/sd), size=1) + 
  theme_bw() + labs(color="RMSE/SD [-]") +
  scale_colour_gradientn(colours = c("yellow", "green", "darkred")) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))

spatial_model



ggplot() + 
  geom_sf(data = shape_europe, colour = "grey30", fill = "white") +
  geom_sf(data = sf_eval, aes(), color="darkgrey", size=1) + 
  geom_sf(data = subset(sf_eval, RMSE < range), aes(color=RMSE), size=1) + 
  theme_bw() + labs(color="RMSE [\u2030]") +
  scale_colour_gradientn(colours = c("yellow", "green", "darkred")) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))







## now apply the spatial model on the temp & precip data from the ERA network for all the GNIP stations


data <- data %>% mutate(fits = coefs[1] + coefs[2]*lat + coefs[3]*long + 
                          coefs[4] * avg_alt + coefs[5]*dO.rec,
                        residuals = dO - fits)



model_eval <- data %>% 
  group_by(station, lat, long) %>%
  summarise(tss = sum((dO - mean(dO))^2),
            rss = sum(residuals^2),
            rsq = 1 - (rss / tss),
            RMSE = sqrt(sum(residuals^2) / n()),
            median = round(median(dO, na.rm=T), 2),
            IQR = round(IQR(dO, na.rm=T), 2),
            sd = round(sd(dO, na.rm=T), 2),
            ratio = RMSE/sd,
            N = n())




model_eval.all <- data.frame(
  rsq = weighted.mean(model_eval$rsq, model_eval$N, na.rm=T),
  RMSE = weighted.mean(model_eval$RMSE, model_eval$N, na.rm=T))



sf_eval <- st_as_sf(model_eval, coords = c("long", "lat"), crs = 4326)



shape_europe <- st_read("~/MUSES/Isotope_model/data/NUTS_RG_20M_2021_4326.shp")
shape_europe <- subset(shape_europe, LEVL_CODE==0)






general_model <- ggplot() + 
  geom_sf(data = shape_europe, colour = "grey30", fill = "white") +
  geom_sf(data = subset(sf_eval, RMSE >= sd), aes(), color="red", size=1, shape=4) + 
  geom_sf(data = subset(sf_eval, RMSE < sd), aes(color=RMSE/sd), size=1) + 
  theme_bw() + labs(color="RMSE/SD [-]") +
  scale_colour_gradientn(colours = c("yellow", "green", "darkred")) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))


general_model





ggplot() + 
  geom_sf(data = shape_europe, colour = "grey30", fill = "white") +
  geom_sf(data = subset(sf_eval, RMSE/sd > 1.5), aes(), color="red", size=1, shape=4) + 
  geom_sf(data = subset(sf_eval, RMSE/sd < 1.5), aes(color=RMSE/sd), size=1) + 
  theme_bw() + labs(color="RMSE/SD [-]") +
  scale_colour_gradientn(colours = c("yellow", "green", "darkred")) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))




ggplot() + 
  geom_sf(data = shape_europe, colour = "grey30", fill = "white") +
  geom_sf(data = sf_eval, aes(), color="darkgrey", size=1) + 
  geom_sf(data = subset(sf_eval, RMSE < IQR), aes(color=RMSE), size=1) + 
  theme_bw() + labs(color="RMSE [\u2030]") +
  scale_colour_gradientn(colours = c("yellow", "green", "darkred")) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))









sf_eval2 <- sf_eval[sf_eval$N >= 96,] %>%  arrange(desc(RMSE)) ## select the stations with over 8 years of data and arrange them in descending RMSE order
#sf_eval2 <- sf_eval2[seq(7,95,8),]  ## select 12 stations at equal RMSE steps
sf_eval2 <- sf_eval2[c(8, 15, 23, 31, 39, 47, 56, 63, 71, 79, 87, 95),]

data2 <- merge(data, sf_eval2[,1])


data2 <- arrange(data2, desc(date)) %>% 
  group_by(station) %>%
  mutate(index = row_number(),
         type = ifelse(index>min(sf_eval2$N), "no","yes"))

data2 <- subset(data2, type == "yes")


# Create a new data frame for facet labels
facet_labels <- data2 %>%
  group_by(station, avg_alt) %>%
  summarise(date=min(date), dO = max(data2$dO) * 0.95)


facet_labels <- merge(facet_labels, sf_eval2, by="station")

facet_labels <- facet_labels %>% mutate(rsq = round(rsq, digits=2),
                                        RMSE = round(RMSE, digits=1))

facet_labels$country <- c("(MD)","(ES)","(DE)","(AT)","(DE)","(CH)","(CH)","(CH)","(PT)","(FI)","(DE)","(IT)")

facet_labels <- facet_labels %>% mutate(label = paste(station, country, ", ", avg_alt, "masl"))




plot_eval <- ggplot(data2) +
  geom_line(aes(x=date, y=dO, color="GNIP data")) +
  geom_line(aes(x=date, y=fits, color="Model data")) +
  geom_text(data = facet_labels, aes(x=date, y=dO, label=label), 
            hjust = 0, vjust = 1, size = 2) + 
  geom_text(data = facet_labels, aes(x=date, y=-17, label=paste("R² =", rsq)), 
            hjust = 0, vjust = 1, size = 2) + 
  geom_text(data = facet_labels, aes(x=date, y=-20, label=paste("RMSE =", RMSE)), 
            hjust = 0, vjust = 1, size = 2) + 
  facet_wrap(~station, scales="free_x", ncol = 3) + theme_bw() + 
  labs(x="", y=expression(paste(delta^{18}, "O [\u2030]"))) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(10),
        legend.title=element_blank(),
        legend.text = element_text(size = 8),
        legend.position = "top")  


plot_eval









setwd("C:/Users/turk/Documents/MUSES/Isotope_model")


ggsave(spatial_model, path="figures", filename = "rmse_iqr(sp)_new.png",
       device = ragg::agg_png, dpi=350,
       width = 12, height = 12, units = "cm",
       bg="white")


ggsave(general_model, path="figures", filename = "rmse_iqr(gg).png",
       device = ragg::agg_png, dpi=350,
       width = 12, height = 12, units = "cm",
       bg="white")


ggsave(plot_eval, path="figures", filename = "stations_eval(new).png",
       device = ragg::agg_png, dpi=350,
       width = 15, height = 18, units = "cm",
       bg="white")



df_coef <- t(data.frame(coefs))


write.table(df_coef, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\model_param2.csv', sep = ",", quote = FALSE)










data2 <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/output_data/meteo&elevation_raster.csv", header = TRUE, sep = ",", dec = ".")


data2 <- subset(data2, elevation > - 50)
data2 <- data2 %>% mutate(date = as.Date(date), year = year(date), month=month(date))


data2 <- merge(data2, dff_gwt[,-c(3)], by=c("year","month"))

data2$temp <- data2$t2m - 273.15

data2 <- data2 %>% mutate(dO.rec = a + b*temp)


data2 <- data2 %>% mutate(dO.pred = coefs[1] + coefs[2]*lat + coefs[3]*long + 
                            coefs[4] * elevation + coefs[5]*dO.rec)

data2$date <- as.Date(data2$date)


date_interest <- "1973-01-01"
data2_sf <- st_as_sf(subset(data2, date == date_interest), coords = c("long", "lat"), crs = 4326)




plot_isoscape <- ggplot() + ggtitle(paste(months(as.Date(date_interest)), year(date_interest))) +
  geom_sf(data = data2_sf, aes(color=dO.pred), shape = 15, size=2) + 
  geom_sf(data = shape_europe, colour = "grey30", fill = "transparent") +
  theme_bw() + labs(color=expression(paste(delta^{18}, "O [\u2030]"))) +
  scale_colour_gradientn(colours = c("yellow", "orange", "darkred")) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))





setwd("C:/Users/turk/Documents/MUSES/Isotope_model")


ggsave(plot_isoscape, path="figures", filename = "model_isoscape.png",
       device = ragg::agg_png, dpi=350,
       width = 12, height = 12, units = "cm",
       bg="white")






