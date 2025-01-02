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
  reframe(a =  round(summary(lm(dO ~ temp))$coefficients[1,1], 2),
          b = round(summary(lm(dO ~ temp))$coefficients[2,1], 2),
          rmse = round(sqrt(sum(summary(lm(dO ~ temp))$residuals^2) / sum(P)), digits = 1))


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
    rmse = round(sqrt(sum(residuals^2) / n()), digits = 1),
    sd = sd(dO, na.rm=T),
    rss=sum(residuals^2),
    tss=sum((dO-mean(dO))^2),
    rsq = round(1-rss/tss, 2),
    N = n()
  )





test_set$date <- as.Date(test_set$date)
test_set <- arrange(test_set, date)

model_eval$date <- test_set$date[1] + months(2)
model_eval$dO <- max(dd$dO, na.rm=T) *0.9 




plot_test.a <- ggplot(dd.weekly) + 
  geom_rect(aes(xmin = as.Date("2021-01-01"), xmax = as.Date(max(date)), ymin = -Inf, ymax = Inf),
            fill = "skyblue", alpha=0.3) +
  geom_line(aes(x=date, y=dO, color="Observation"), linewidth=0.5) +
  geom_line(aes(x=date, y=dO.rec, color="Model"), linewidth=0.5) +
  geom_vline(aes(xintercept=as.Date("2021-01-01")), linetype="dashed") +
  scale_color_manual(values=c("darkred", "black")) +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + 
  theme_bw() + labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), colour="") +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="darkgrey"),
        panel.grid.minor = element_blank(),
        text = element_text(size=10),
        legend.text = element_text(size = 8),
        legend.position = c(0.3,0.9),
        legend.direction = "horizontal")  


plot_test.a






plot_test.b <- ggplot(test_set) +
  geom_abline(aes(intercept=0, slope=1)) +
  geom_smooth(aes(x=dO, y=dO.rec), method="lm", se=F, linetype="dashed", color="black") +
  geom_point(aes(x=dO, y=dO.rec, size=P), shape=21) +
  geom_text(data=model_eval, aes(x=-17, y=-1.5, label=paste("R² =", rsq, ", RMSE =", rmse)), 
            hjust = 0, vjust = 1, size = 2.5) + 
  scale_size_continuous(breaks=c(40,80,120), range=c(1,5)) +
  theme_bw() + labs(x=expression(paste(delta^{18}, "O [\u2030] (observation)")),
                    y=expression(paste(delta^{18}, "O [\u2030] (model)")), size="P [mm]") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=10),
        legend.text = element_text(size = 8),
        legend.position = c(0.8,0.25))  



plot_test.b





plot_test <- ggarrange(plot_test.a, plot_test.b, 
                      align="hv", ncol = 2, nrow = 1,
                      labels = c("(a)", "(b)"), widths = c(1.8,1) ,
                      font.label = list(size = 10, face = "plain"))

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
    rmse = round(sqrt(sum(residuals^2) / n()), digits = 1),
    rss=sum(residuals^2),
    tss=sum((dO-mean(dO))^2),
    rsq = round(1-rss/tss, 2),
    N = n()
  )




test_set$date <- as.Date(test_set$date)
test_set <- arrange(test_set, date)

model_eval.monthly$date <- test_set$date[1] + months(2)
model_eval.monthly$dO <- max(dd$dO, na.rm=T) *0.9 




ggplot(dd.monthly) + 
  geom_rect(aes(xmin = as.Date("2021-01-01"), xmax = as.Date(max(date)), ymin = -Inf, ymax = Inf),
            fill = "skyblue", alpha=0.3) +
  geom_line(aes(x=date, y=dO, color="Observation"), linewidth=0.5) +
  geom_line(aes(x=date, y=dO.rec, color="Model"), linewidth=0.5) +
  geom_vline(aes(xintercept=as.Date("2021-01-01")), linetype="dashed") +
  geom_text(data=model_eval.monthly, aes(x=date, y=dO, label=paste("R² =", rsq, ", RMSE =", rmse)), 
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












ddf2_alt <- subset(ddf.mm, year < 2021) %>% 
  reframe(a =  round(summary(lm(dO ~ temp))$coefficients[1,1], 2),
          b = round(summary(lm(dO ~ temp))$coefficients[2,1], 2),
          rmse = round(sqrt(sum(summary(lm(dO ~ temp))$residuals^2) / sum(P)), digits = 1))


dd <- dd %>% mutate(dO.rec = ddf2_alt$a + ddf2_alt$b*temp)


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
    rmse = round(sqrt(sum(residuals^2) / n()), digits = 1),
    sd = sd(dO, na.rm=T),
    rss=sum(residuals^2),
    tss=sum((dO-mean(dO))^2),
    rsq = round(1-rss/tss, 2),
    N = n()
  )







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
    rmse = round(sqrt(sum(residuals^2) / n()), digits = 1),
    rss=sum(residuals^2),
    tss=sum((dO-mean(dO))^2),
    rsq = round(1-rss/tss, 2),
    N = n()
  )











setwd("C:/Users/turk/Documents/MUSES/Isotope_model")


ggsave(plot_test, path="figures", filename = "model_test(Sep 24).png",
       device = ragg::agg_png, dpi=350,
       width = 20, height = 8, units = "cm",
       bg="white")




ddf2 <- arrange(ddf2, gwt)

write.table(ddf2, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\model_param.csv', sep = ",", quote = FALSE)
