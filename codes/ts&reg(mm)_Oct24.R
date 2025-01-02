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


df <- read_excel("~/MUSES/Isotope_model/data/Precip&isotopes_Belvaux_2016-2023_update_02022023.xlsx", 
                 col_types = c("numeric", "date", "date", 
                               "text", "numeric", "numeric", "numeric", 
                               "numeric", "numeric"))

gwt <- read_excel("~/MUSES/Isotope_model/data/GWL_1954-2022.xlsx", 
                  col_types = c("date", "text", "numeric", 
                                "text", "text", "skip", "skip", "skip", 
                                "skip", "skip", "skip", "skip", "skip", 
                                "skip"))

temp <- read.csv("~/MUSES/1st publication/data/Schifflange Air Temperature Mean (4727010).csv", sep=",", dec=".", skip=17)


rain <- read_excel("~/MUSES/Isotope_model/data/CR10X_rainfall_Belvaux_LIST_2016-2021_Guilhem.xlsx", 
                   col_types = c("date", "numeric"))


colnames(df) <- c("index","start","end","type","P","dD","dD_dev","dO","dO_dev")
colnames(temp) <- c("timestamp","temp","code")



df$d_excess <- df$dD - 8*df$dO
df$date <- date(df$end)
df$duration <- ceiling(as.numeric(df$end - df$start)/3600)

df <- drop_na(df)

df <- df %>% mutate(index=row_number())


df[duplicated(df$start)==T,]




temp <- temp %>% mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M"))

rain <- rain %>% mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M"),
                        timestamp = round_date(timestamp, unit = "minute"))

df_meteo <- merge(temp[,-3], rain, by="timestamp", all.y=T)
df_meteo <- drop_na(df_meteo)
df_meteo$temp <- as.numeric(df_meteo$temp)





## index all the obseration times to calculate event-based meteorological variables



dff <- pivot_longer(df[,c(1:3)],!index, values_to = "timestamp", names_to = "start/end")

dff <- dff %>% mutate(timestamp = round_date(timestamp, unit = "hour"),
                      control = as.numeric(lead(timestamp)-timestamp)/3600)

dff$control[is.na(dff$control)==T] <- 0

dff$timestamp[dff$control<0] <- lead(dff$timestamp)[dff$control<0]
dff <- dff %>% mutate(control = as.numeric(lead(timestamp)-timestamp)/3600)


df_helper <- dff %>% 
  group_by(index) %>%
  reframe(timestamp = seq.POSIXt(min(timestamp), max(timestamp), by="hour"))


ddf <- merge(dff, df_helper, by=c("index","timestamp"), all.y=T)



dff_meteo <- merge(ddf[,c(1,2)], df_meteo)

dff_meteo.dd <- dff_meteo %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(index, date) %>%
  reframe(precip = sum(rain, na.rm=T)) %>%
  group_by(index) %>%
  mutate(max = max(precip, na.rm=T))


dates <- data.frame(index = dff_meteo.dd$index[dff_meteo.dd$precip==dff_meteo.dd$max],
                    date = dff_meteo.dd$date[dff_meteo.dd$precip==dff_meteo.dd$max])


dates <- dates %>% distinct(index, .keep_all = T)


ddf_meteo <- dff_meteo %>%
  group_by(index) %>%
  reframe(min_temp =  min(temp, na.rm=T),
          max_temp =  max(temp, na.rm=T),
          temp = median(temp, na.rm=T),
          precip = sum(rain, na.rm=T),
          duration = n())

ddf_meteo <- merge(dates, ddf_meteo)





ddf <- merge(df[,-c(11,12)], ddf_meteo, by="index")
ddf <- merge(ddf, gwt[,c(1,4)])


ddf$gwt <- factor(ddf$gwt, levels=c("West","North West","South West","High CE","Low CE","North","South","East"))

levels(ddf$gwt) <- list(W = "West", NW = "North West", SW = "South West", 
                        HCE = "High CE", LCE = "Low CE",
                        N= "North", S = "South", E = "East")

ddf <- ddf %>%  mutate(gwt2=recode(gwt, "W"="Zonal", "NW"="Mixed", "SW"="Mixed", "HCE"="Mixed",
                                   "LCE"="Mixed", "N"="Meridional", "S"="Meridional", "E"="Meridional"))


ggplot(ddf) +
  geom_point(aes(x=P, y=precip)) + theme_bw()






ddf %>%
  reframe(span_dD = round(max(dD, na.rm=T) - min(dD, na.rm=T),1),
          span_dO = round(max(dO, na.rm=T) - min(dO, na.rm=T),1),
          span_d = max(d_excess, na.rm=T) - min(d_excess, na.rm=T),
          dD = round(weighted.mean(dD, P, na.rm=T), 1),
          dO = round(weighted.mean(dO, P, na.rm=T), 1),
          d_excess = round(weighted.mean(d_excess, P, na.rm=T), 1),
          N = n())


ddf %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(year, month) %>%
  reframe(dD = sd(dD, na.rm=T),
          dO = sd(dO, na.rm=T),
          d_excess = sd(d_excess, na.rm=T)) %>%
  reframe(dD = round(mean(dD, na.rm=T),1),
          dO = round(mean(dO, na.rm=T),1),
          d_excess = round(mean(d_excess, na.rm=T),1)) 



ddf %>% 
  mutate(month_nr = month(date),
         season = ifelse(month_nr %in% c(12,1,2), "winter",
                         ifelse(month_nr %in% c(6,7,8), "summer","other"))) %>% 
  group_by(season) %>%
  reframe(dD = round(weighted.mean(dD, P, na.rm=T), 1),
          dO = round(weighted.mean(dO, P, na.rm=T), 1),
          d_excess = round(weighted.mean(d_excess, P, na.rm=T), 1))







## calculate event-based variables and plot results


ddf_event <- ddf %>%
  mutate(Group = cumsum(gwt != lag(gwt, default = first(gwt))))

ddf_event <- ddf_event %>%
  group_by(Group, gwt, gwt2) %>%
  reframe(date = min(date),
          dD = weighted.mean(dD, P, na.rm=T),
          dO = weighted.mean(dO, P, na.rm=T),
          d_excess = weighted.mean(d_excess, P, na.rm=T),
          P = sum(P, na.rm=T),
          temp = median(temp, na.rm=T),
          duration = n())

ddf_event <- drop_na(ddf_event)


ggplot(ddf_event) +
  geom_point(aes(x=temp, y=dO)) + theme_bw()





ggplot(ddf) + 
  geom_smooth(aes(x=temp, y=dO), method = lm, formula = y ~ splines::ns(x, 2), color="black",  linetype="dashed", fullrange = T) +
  geom_point(aes(x=temp, y=dO, size=max_temp-min_temp, color=gwt), shape=21, alpha=0.7, stroke=1.1) +
  stat_regline_equation(aes(x=temp, y=dO, label = ..eq.label..), formula = y ~ splines::ns(x, 2), label.y = 1.5, size=3) + 
  stat_regline_equation(aes(x=temp, y=dO, label = ..rr.label..), formula = y ~ splines::ns(x, 2), label.y = 0, size=3) +
  scale_color_viridis(discrete=T) + theme_bw() +  
  labs(x="Temperature [°C]", y=expression(paste(delta^{18}, "O [\u2030]")), 
       color="Atmospheric\ncirculation type", size="Temperature range [°C]") +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size=8),
        legend.text = element_text(size = 8))   




ddf.mm <- ddf %>%
  mutate(month = month(date)) %>%
  group_by(gwt, month) %>%
  reframe(date = min(date),
          dD = weighted.mean(dD, precip, na.rm=T),
          dO = weighted.mean(dO, precip, na.rm=T),
          d_excess = weighted.mean(d_excess, precip, na.rm=T),
          precip = sum(precip, na.rm=T),
          temp = mean(temp, na.rm=T))





plot_dO <- ggplot(ddf_event, aes(x=date, y=dO)) + 
  geom_line(data=df, aes(x=date, y=dO), color="grey") +
  geom_point(aes(color=gwt, size=P), shape=21, alpha=0.8, stroke=0.9) +
  geom_smooth(method="loess", span=0.1, se=F, color="black", linetype="twodash") +
  scale_color_viridis(discrete=T) +
  theme_bw() + labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), 
                    color="", size="P [mm]") +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="darkgrey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1),
        text = element_text(size=10))


boxplot_dO <- ggplot(ddf_event) + 
  geom_boxplot(aes(x=gwt, y=dO, fill=gwt), alpha=0.7) +
  theme_bw() + labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), fill="") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        plot.margin = margin(1,1,1,1),
        text=element_text(size=8)) + 
  scale_fill_viridis(discrete=T) 




plot_dexcess <- ggplot(ddf_event, aes(x=date, y=d_excess)) + 
  geom_line(data=df, aes(x=date, y=d_excess), color="grey") +
  geom_point(aes(color=gwt, size=P), shape=21, alpha=0.8, stroke=0.9) +
  geom_smooth(method="loess", span=0.1, se=F, color="black", linetype="twodash") +
  scale_color_viridis(discrete=T) +
  theme_bw() + labs(x="", y=expression(paste("d-excess [\u2030]")), 
                    color="", size="P [mm]") +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="darkgrey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1),
        text = element_text(size=10))


boxplot_dexcess <- ggplot(ddf_event) + 
  geom_boxplot(aes(x=gwt, y=d_excess, fill=gwt), alpha=0.7) +
  theme_bw() + labs(x="", y=expression(paste("d-excess [\u2030]")), fill="") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        plot.margin = margin(1,1,1,1),
        text=element_text(size=8)) + 
  scale_fill_viridis(discrete=T) 



ddf.weekly <- ddf %>%
  mutate(year=year(date), week = week(date)) %>%
  group_by(year, week) %>%
  reframe(date = min(date),
          precip = sum(precip, na.rm=T),
          temp = mean(temp, na.rm=T))


plot_meteo <- ggplot(ddf.weekly) + 
  geom_col(aes(x=date, y=-(precip)/4), color="lightblue") +
  geom_point(data=ddf_event, aes(x=date, y=temp-40), shape=4) +
  scale_y_continuous(labels = function(y) y + 40, 
                     sec.axis = sec_axis(~.*(-4), name="Precipitation [mm]")) + 
  theme_bw() + labs(x="", y="Temperature [°C]") +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="darkgrey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1),
        text = element_text(size=10))






iso_plot <- ggarrange(plot_dO, plot_dexcess,  plot_meteo,
                      align="hv", common.legend = T,  ncol = 1, nrow = 3,
                      labels = c("(a)", "(b)",  "(c)"), font.label = list(size = 10, face = "plain"),
                      legend = "top")

iso_plot










ggplot(ddf) +  
  geom_smooth(aes(x=dO, y=dD, group=gwt), method = lm, color="transparent",  formula = y ~ x, fullrange = T) +
  geom_smooth(aes(x=dO, y=dD), method = lm, formula = y ~ x, color="black",  
              linetype="dashed", linewidth=0.7, fullrange = T, se=F) +
  geom_abline(aes(intercept=10, slope=8)) +
  geom_point(aes(x=dO, y=dD), size=0.7) + 
  theme_bw() + labs(x=expression(paste(delta^{18}, "O [\u2030]")), 
                    y=expression(paste(delta^{2}, "H [\u2030]"))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))





dual_eq <- ddf %>% 
  reframe(a =  round(summary(lm(dD ~ dO, weights = P))$coefficients[1,1], 2),
          b = round(summary(lm(dD ~ dO, weights = P))$coefficients[2,1], 2),
          rss=sum(P*summary(lm(dD ~ dO, weights = P))$residuals^2),
          tss=sum(P*(dD-weighted.mean(dD, P))^2),
          rsq = round(1-rss/tss, 2), 
          equation = paste("y = ", a, " + ", b, " * x"))


ddf %>% 
  group_by(gwt) %>% 
  reframe(a =  round(summary(lm(dD ~ dO, weights = P))$coefficients[1,1], 2),
          b = round(summary(lm(dD ~ dO, weights = P))$coefficients[2,1], 2),
          rss=sum(P*summary(lm(dD ~ dO, weights = P))$residuals^2),
          tss=sum(P*(dD-weighted.mean(dD, P))^2),
          rsq = round(1-rss/tss, 2), 
          equation = paste("y = ", a, " + ", b, " * x"))











## temperature-d18O regression part


ddf.mm <- ddf %>%
  mutate(month = month(date)) %>%
  group_by(gwt, gwt2, month) %>%
  reframe(date = min(date),
          dD = weighted.mean(dD, P, na.rm=T),
          dO = weighted.mean(dO, P, na.rm=T),
          d_excess = weighted.mean(d_excess, P, na.rm=T),
          P = sum(P, na.rm=T),
          temp = mean(temp, na.rm=T))


ddf.mm2 <- ddf %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  reframe(P_tot = sum(P, na.rm=T))


ddf.mm <- merge(ddf.mm, ddf.mm2, by="month")
ddf.mm <- ddf.mm %>% mutate(contr = P/P_tot)



ddf2 <- ddf %>% 
  group_by(gwt) %>% 
  reframe(a =  round(summary(lm(dO ~ temp, weights = P))$coefficients[1,1], 2),
          b = round(summary(lm(dO ~ temp, weights = P))$coefficients[2,1], 2),
          rmse = round(sqrt(sum(summary(lm(dO ~ temp, weights = P))$residuals^2) / sum(P)), digits = 2),
          rss=sum(summary(lm(dO ~ temp, weights = P))$residuals^2),
          tss=sum(P*(dO-weighted.mean(dO, P))^2),
          rsq = round(1-rss/tss, 2), 
          median = round(median(dO, na.rm=T), 2),
          IQR= round(IQR(dO, na.rm=T), 2),
          equation = paste("y = ", a, " + ", b, " * x"))


ddf2 <- ddf.mm %>% 
  group_by(gwt) %>% 
  reframe(a =  round(summary(lm(dO ~ temp))$coefficients[1,1], 1),
          b = round(summary(lm(dO ~ temp))$coefficients[2,1], 2),
          rmse = round(summary(lm(dO ~ temp))$sigma, 1),
          rss=sum(summary(lm(dO ~ temp))$residuals^2),
          tss=sum((dO-mean(dO))^2),
          rsq = round(1-rss/tss, 2), 
          median = round(median(dO, na.rm=T), 1),
          IQR= round(IQR(dO, na.rm=T), 1),
          equation = paste("y = ", a, " + ", b, " * x"))



ddf_stats <- ddf2[,c(4,7:10)] 




ddf2 %>%
  reframe(a = mean(a),
          b= mean(b),
          rmse = mean(rmse),
          rsq = mean(rsq))





# Create a new data frame for facet labels
facet_labels <- ddf.mm %>%
  group_by(gwt2) %>%
  summarise(temp=min(ddf_event$temp), dO = max(ddf_event$dO) * 0.95)



plot_dO <- ggplot(ddf.mm) + 
  geom_smooth(data=ddf.mm[,-1], aes(x=temp, y=dO), method = lm, formula = y ~ x, color="black",  linetype="dashed", fullrange = T) +
  geom_point(aes(x=temp, y=dO, size=contr, color=gwt), shape=21, alpha=0.7, stroke=1.1) +
  geom_smooth(aes(x=temp, y=dO, color=gwt), method=lm, formula = y ~ x, se=F, fullrange = T) +
  geom_text(data = facet_labels, aes(x=temp, y=dO, label=gwt2), 
            hjust = 0, vjust = 1, size = 3) + 
  scale_color_viridis(discrete=T) + theme_bw() + facet_wrap(~gwt2) +
  labs(x="Temperature [°C]", y=expression(paste(delta^{18}, "O [\u2030]")), 
       color="CP", size="Contribution [%]") +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size=10),
        legend.text = element_text(size = 8), 
        legend.position = "top")   






# Create a new data frame for facet labels
facet_labels <- ddf.mm %>%
  group_by(gwt2) %>%
  summarise(temp=min(ddf_event$temp), d_excess = max(ddf_event$d_excess) * 0.95)


plot_dexcess <- ggplot(ddf.mm) + 
  geom_smooth(data=ddf.mm[,-1], aes(x=temp, y=d_excess), method = lm, formula = y ~ x, color="black",  linetype="dashed", fullrange = T) +
  geom_point(aes(x=temp, y=d_excess, size=contr, color=gwt), shape=21, alpha=0.7, stroke=1.1) +
  geom_smooth(aes(x=temp, y=d_excess, color=gwt), method=lm, formula = y ~ x, se=F, fullrange = T) +
  geom_text(data = facet_labels, aes(x=temp, y=d_excess, label=gwt2), 
            hjust = 0, vjust = 1, size = 3) + 
  scale_color_viridis(discrete=T) + theme_bw() + facet_wrap(~gwt2) +
  labs(x="Temperature [°C]", y=expression(paste("d-excess [\u2030]")), 
       color="Atmospheric\ncirculation type", size="Total precipitation [mm]") +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size=8),
        legend.title=element_blank(),
        legend.text = element_text(size = 8))  





reg_plot <- ggarrange(plot_dO, plot_dexcess, 
                      align="hv", common.legend = T,  ncol = 1, nrow = 2,
                      labels = c("a", "b"), font.label = list(size = 10, face = "plain"),
                      legend = "top")

reg_plot











setwd("C:/Users/turk/Documents/MUSES/Isotope_model")


ggsave(iso_plot, path="figures", filename = "iso_ts.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 18, units = "cm",
       bg="white")


ggsave(reg_plot, path="figures", filename = "iso_temp_reg(new).png",
       device = ragg::agg_png, dpi=350,
       width = 17, height = 17, units = "cm",
       bg="white")


ggsave(plot_dO, path="figures", filename = "dO-T_plot.png",
       device = ragg::agg_png, dpi=350,
       width = 17, height = 9, units = "cm",
       bg="white")





write.table(ddf, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\ddf.csv', sep = ",", quote = FALSE)
write.table(ddf2[,c(1:4,7)], file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\dO-T.csv', sep = ",", quote = FALSE)

