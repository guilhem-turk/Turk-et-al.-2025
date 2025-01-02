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
library(Hmisc)
library(ggrepel)


ddf <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/output_data/ddf.csv", header = TRUE, sep = ",", dec = ".")

ddf$P <- as.numeric(ddf$P)

ddf$date <- as.Date(ddf$date)
ddf$year <- year(ddf$date)



ddf.mm <- ddf %>%
  mutate(month = month(date)) %>%
  group_by(year, month, gwt, gwt2) %>%
  reframe(date = min(date),
          dD = weighted.mean(dD, P, na.rm=T),
          dO = weighted.mean(dO, P, na.rm=T),
          d_excess = weighted.mean(d_excess, P, na.rm=T),
          P = sum(P, na.rm=T),
          temp = mean(temp, na.rm=T))





ddf2 <- subset(ddf.mm, year < 2021) %>% 
  group_by(gwt) %>% 
  reframe(a =  round(summary(lm(dO ~ temp, weights = P))$coefficients[1,1], 2),
          b = round(summary(lm(dO ~ temp, weights = P))$coefficients[2,1], 2),
          rmse = round(sqrt(sum(summary(lm(dO ~ temp, weights = P))$residuals^2) / sum(P)), digits = 2))


ddf2$gwt <- factor(ddf2$gwt, levels=c("W","NW","SW","HCE","LCE","N","S","E"))


dd <- merge(ddf, ddf2, by="gwt")
dd <- dd %>% mutate(dO.rec = a + b*temp)


dd.weekly <- dd %>%
  mutate(week = week(date)) %>%
  group_by(year, week) %>%
  reframe(date = min(date),
          dD = weighted.mean(dD, P, na.rm=T),
          dO = weighted.mean(dO, P, na.rm=T),
          d_excess = weighted.mean(d_excess, P, na.rm=T),
          dO.rec = weighted.mean(dO.rec, P, na.rm=T),
          P = sum(P, na.rm=T),
          temp = mean(temp, na.rm=T))

# dd.weekly <- drop_na(dd.weekly)


test_set <- subset(dd.weekly, year >= 2021)
test_set <- test_set %>% mutate(residuals = dO - dO.rec)


model_eval <- test_set %>%
  summarise(
    rmse = round(sqrt(sum(summary(lm(dO ~ temp, weights = P))$residuals^2) / sum(P)), digits = 2),
    rss=sum(summary(lm(dO ~ temp, weights = P))$residuals^2),
    tss=sum(P*(dO-weighted.mean(dO, P))^2),
    rsq = round(1-rss/tss, 2),
    N = n()
  )




test_set$date <- as.Date(test_set$date)
test_set <- arrange(test_set, date)

model_eval$date <- test_set$date[1] + months(2)
model_eval$dO <- max(dd$dO, na.rm=T) *0.9 




plot_test <- ggplot(dd.weekly) + 
  geom_rect(aes(xmin = as.Date("2021-01-01"), xmax = as.Date(max(date)), ymin = -Inf, ymax = Inf),
            fill = "skyblue", alpha=0.3) +
  geom_smooth(aes(x=date, y=dO, color="Observation", weight = P), 
              method = "loess", linetype="twodash", span=0.075, linewidth=0.5, se=F) +
  geom_smooth(aes(x=date, y=dO, weight = P), color="transparent", 
              method = "loess", linetype="twodash", span=0.075) +
  geom_line(aes(x=date, y=dO, color="Observation"), linewidth=0.5) +
  geom_line(aes(x=date, y=dO.rec, color="Model"), linewidth=0.5) +
  geom_vline(aes(xintercept=as.Date("2021-01-01")), linetype="dashed") +
  geom_text(data=model_eval, aes(x=date, y=dO, label=paste("R² =", rsq, ", RMSE =", rmse)), 
            hjust = 0, vjust = 1, size = 3, color="darkred") + 
  scale_color_manual(values=c("darkred", "black")) +
  scale_x_date(date_breaks="1 year", date_labels="%Y") + 
  theme_bw() + labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), colour="") +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.position = c(0.2,0.9),
        legend.direction = "horizontal")  


plot_test







dd.monthly <- dd %>%
  mutate(month = month(date)) %>%
  group_by(year, month) %>%
  reframe(date = min(date),
          dD = weighted.mean(dD, P, na.rm=T),
          dO = weighted.mean(dO, P, na.rm=T),
          d_excess = weighted.mean(d_excess, P, na.rm=T),
          dO.rec = weighted.mean(dO.rec, P, na.rm=T),
          P = sum(P, na.rm=T),
          temp = mean(temp, na.rm=T))

# dd.weekly <- drop_na(dd.weekly)


test_set <- subset(dd.monthly, year >= 2021)
test_set <- test_set %>% mutate(residuals = dO - dO.rec)


model_eval.monthly <- test_set %>%
  summarise(
    rmse = round(sqrt(sum(summary(lm(dO ~ temp, weights = P))$residuals^2) / sum(P)), digits = 2),
    rss=sum(summary(lm(dO ~ temp, weights = P))$residuals^2),
    tss=sum(P*(dO-weighted.mean(dO, P))^2),
    rsq = round(1-rss/tss, 2),
    N = n()
  )









setwd("C:/Users/turk/Documents/MUSES/Isotope_model")


ggsave(plot_test, path="figures", filename = "model_test(Aug 24).png",
       device = ragg::agg_png, dpi=350,
       width = 16, height = 8, units = "cm",
       bg="white")




ddf2 <- arrange(ddf2, gwt)

write.table(ddf2, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\model_param.csv', sep = ",", quote = FALSE)



