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


ddf <- ddf %>% 
  mutate(P = as.numeric(P),
         date = as.Date(date),
         year = year(date),
         month = months(date, abbr = T))

ddf$gwt <- factor(ddf$gwt, levels=c("W","NW","SW","HCE","LCE","N","S","E"))


df <- subset(ddf, year > 2016 & year < 2023)


df1 <- drop_na(df) %>%                                      
  group_by(gwt, year) %>% 
  summarise(P=sum(P, na.rm=T))

df1 <- df1 %>% group_by(year) %>%complete(gwt,fill=list(P=0))

df2 <- drop_na(df) %>%                                      
  group_by(year) %>% 
  summarise(tot_rain=sum(P, na.rm=T))

dff <- merge(df1, df2)
dff <- arrange(dff, year, gwt)

dff$contr <- dff$P/dff$tot_rain



plot1 <- ggplot(dff, aes(x=year, y=contr, fill=gwt)) + 
  geom_area(alpha=0.7, size=1, color="white") + 
  scale_fill_viridis(discrete = T) + theme_ipsum() +
  labs(x="Year", y="Fraction of total rain amount [-]", fill="") +
  theme(panel.grid.minor.x = element_blank(),
        text = element_text(size=10),
        legend.text = element_text(size = 8))




df1 <- drop_na(df) %>%                                      
  group_by(gwt, month) %>% 
  summarise(P=sum(P, na.rm=T))

df2 <- drop_na(df) %>%                                      
  group_by(month) %>% 
  summarise(tot_rain=sum(P, na.rm=T))

df2$month_ <- 1:nrow(df2)

dff <- merge(df1, df2)
dff <- arrange(dff, month, gwt)

dff$contr <- dff$P/dff$tot_rain



dte_formatter <- function(x) { 
  #formatter for axis labels: J12, F12, M12, etc... 
  mth <- substr(format(x, "%b"),1,1) 
  #yr <- format(x, "%y") 
  #paste0(mth, yr) 
  mth
} 


dff$month <- factor(dff$month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))



plot2 <- ggplot(dff) + 
  geom_col(aes(x=month, y=contr, fill=gwt), position = "stack", color = "white", width = 0.7, alpha=0.7) +
  scale_fill_viridis(discrete = T) + theme_ipsum() + labs(fill="") +
  labs(x="Month", y="Fraction of total rain amount [-]") + 
  scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=10),
        legend.text = element_text(size = 8))


plot_stacked <- ggarrange(plot1, plot2 , 
                          labels=c("(a)","(b)"), font.label = list(size = 10, face = "plain"), 
                          common.legend = T, legend="right")

plot_stacked








setwd("C:/Users/turk/Documents/MUSES/Isotope_model")
getwd()


ggsave(plot_stacked, path="figures", filename = "gwt_contr.png",
       device = ragg::agg_png, dpi=350,
       width = 20, height = 10, units = "cm",
       bg="white")














