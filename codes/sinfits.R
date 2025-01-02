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


df <-  read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/output_data/ddf.csv", header = TRUE, sep = ",", dec = ".")

df$P <- as.numeric(df$P)
df$date <- as.Date(df$date)


dff <- drop_na(df) %>%
  mutate(week = week(date), month = month(date), year= year(date)) %>%
  group_by(gwt, month) %>%
  reframe(de.sd = sd(d_excess),
          dO.sd = sd(dO),
          de = weighted.mean(d_excess, P),
          dO = weighted.mean(dO, P),
          precip = sum(P),
          count = n())


dff$month <- factor(dff$month, levels=c("10","11","12","1","2","3","4","5","6","7","8","9"))
dff$gwt <- factor(dff$gwt, levels=c("W","NW","SW","HCE","LCE","N","S","E"))



plot_dO <- ggplot(dff) + geom_point(aes(x=month, y=dO, fill=gwt, size=count), shape=21, color="black", alpha=0.7) + 
  geom_smooth(aes(x=as.numeric(month), y=dO), se=F, linetype="dashed", color="black") +
  scale_size(range=c(2,10)) + scale_fill_viridis(discrete=T) + theme_minimal() +  
  scale_x_discrete(labels=c("O","N","D", "J","F","M","A","M","J","J","A","S")) +
  labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), size="Number of\nsamples", fill="CP") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(text = element_text(size=10),
        plot.margin = margin(0.5,0.5,0.5,0.5, unit="cm"))



plot_dO.sd <- ggplot(dff) + geom_point(aes(x=month, y=dO.sd, fill=gwt, group=gwt, size=count), shape=21, color="black", alpha=0.7) + 
  scale_size(range=c(1,12)) + scale_fill_viridis(discrete=T) + theme_minimal() +  
  scale_x_discrete(labels=c("O","N","D", "J","F","M","A","M","J","J","A","S")) +
  labs(x="", y=expression(paste(delta^{18}, "O [\u2030] standard deviation")), 
       size="Number of\nsamples", fill="CP") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(text = element_text(size=10),
        plot.margin = margin(0.5,0.5,0.5,0.5, unit="cm"))


plot_de <- ggplot(dff) + geom_point(aes(x=month, y=de, fill=gwt, size=count), shape=21, color="black", alpha=0.7) + 
  geom_smooth(aes(x=as.numeric(month), y=de), se=F, linetype="dashed", color="black") +
  scale_size(range=c(2,10)) + scale_fill_viridis(discrete=T) + theme_minimal() +  
  scale_x_discrete(labels=c("O","N","D", "J","F","M","A","M","J","J","A","S")) +
  labs(x="", y="d-excess [\u2030]", size="Number of\nsamples", fill="CP") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(text = element_text(size=10),
        plot.margin = margin(0.5,0.5,0.5,0.5, unit="cm"))



plot_de.sd <- ggplot(dff) + geom_point(aes(x=month, y=de.sd, fill=gwt, group=gwt, size=count), shape=21, color="black", alpha=0.7) + 
  scale_size(range=c(1,12)) + scale_fill_viridis(discrete=T) + theme_minimal() + 
  scale_x_discrete(labels=c("O","N","D", "J","F","M","A","M","J","J","A","S")) +
  labs(x="", y="d-excess [\u2030] standard deviation", size="Number of\nsamples", fill="CP") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(text = element_text(size=10),
        plot.margin = margin(0.5,0.5,0.5,0.5, unit="cm"))


plot_gwt <- ggarrange(plot_dO , plot_dO.sd, plot_de, plot_de.sd,
                      labels = c("(a)", "(b)", "(c)", "(d)"), font.label = list(size = 10, face = "plain"),
                      common.legend = T, legend="right",
                      ncol = 2, nrow = 2)

plot_gwt










ddf <- drop_na(df) %>%
  mutate(week = week(date), month = month(date), year= year(date)) %>%
  group_by(gwt, year, month) %>%
  reframe(de.sd = sd(d_excess),
          dO.sd = sd(dO),
          de = weighted.mean(d_excess, P),
          dO = weighted.mean(dO, P),
          P = sum(P),
          count = n())


ddf$gwt <- factor(ddf$gwt, levels=c("W","NW","SW","HCE","LCE","N","S","E"))


df_helper.a <- data.frame(gwt=rep(c("W","NW","SW","HCE","LCE","N","S","E"), 12),
                          month=rep(c(1:12), each=8))

df_helper.b <- data.frame(month=rep(c(1:12), 6),
                          year=rep(c(2017:2022),each=12))

df_helper <- merge(df_helper.a, df_helper.b, by="month")



ddf <- merge(ddf, df_helper, by=c("gwt","year","month"), all=T)







sin.fit <- function (df) {
  
  T = 12
  x <- seq(1, length(df$dO))
  
  df$dO[is.na(df$dO)==TRUE] <- 999
  df$P[is.na(df$P)==TRUE] <- 0
  
  fit.lm <- nls(dO~A*sin(2*pi*x*(1/T)-P)+k, data=df, start=list(A=2,P=1, k=0),
                weights = P)
  
  fits <- data.frame(fits=fitted(fit.lm))
  stats <- summary(fit.lm)
  A <- round(stats$parameters[1,1], digits=1)
  phi <- round(stats$parameters[2,1], digits=0)
  k <- round(stats$parameters[3,1], digits=1)
  
  dd <- fits
  dd$dO <- df$dO
  dd$P <- df$P
  dd$dO[dd$dO==0] <- NA
  dd <- dd %>% drop_na(dO)
  
  tss =  sum(dd$P*(dd$dO - weighted.mean(dd$dO, dd$P))^2 )
  rss =  sum(dd$P*(dd$dO-dd$fits)^2)
  rsq  = round(1 - (rss/tss), digits=2)
  
  return(list(fits, stats, A, phi, k, rsq))
  
} # fits a sine wave curve by applying a non-linear fitting algorithm




df_rsq <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[6]]) %>%
  setNames(unique(sort(ddf$gwt)))

df_A <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[3]]) %>%
  setNames(unique(sort(ddf$gwt)))

df_phi <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[4]]) %>%
  setNames(unique(sort(ddf$gwt)))

df_k <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[5]]) %>%
  setNames(unique(sort(ddf$gwt)))


df_results.dO <- data.frame(
  A = unlist(df_A), 
  phi = unlist(df_phi),
  k = unlist(df_k),
  rsq = unlist(df_rsq)
)








sin.fit <- function (df) {
  
  T = 12
  x <- seq(1, length(df$de))
  
  df$de[is.na(df$de)==TRUE] <- 999
  df$P[is.na(df$P)==TRUE] <- 0
  
  fit.lm <- nls(de~A*sin(2*pi*x*(1/T)-P)+k, data=df, start=list(A=2,P=1, k=0),
                weights = P)
  
  fits <- data.frame(fits=fitted(fit.lm))
  stats <- summary(fit.lm)
  A <- round(stats$parameters[1,1], digits=1)
  phi <- round(stats$parameters[2,1], digits=0)
  k <- round(stats$parameters[3,1], digits=1)
  
  dd <- fits
  dd$de <- df$de
  dd$P <- df$P
  dd$de[dd$de==0] <- NA
  dd <- dd %>% drop_na(de)
  
  tss =  sum(dd$P*(dd$de - weighted.mean(dd$de, dd$P))^2 )
  rss =  sum(dd$P*(dd$de-dd$fits)^2)
  rsq  = round(1 - (rss/tss), digits=2)
  
  return(list(fits, stats, A, phi, k, rsq))
  
} # fits a sine wave curve by applying a non-linear fitting algorithm




df_rsq <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[6]]) %>%
  setNames(unique(sort(ddf$gwt)))

df_A <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[3]]) %>%
  setNames(unique(sort(ddf$gwt)))

df_phi <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[4]]) %>%
  setNames(unique(sort(ddf$gwt)))

df_k <- ddf %>%                                      
  group_by(gwt) %>% 
  group_map(~sin.fit(.x)[[5]]) %>%
  setNames(unique(sort(ddf$gwt)))


df_results.de <- data.frame(
  A = unlist(df_A), 
  phi = unlist(df_phi),
  k = unlist(df_k),
  rsq = unlist(df_rsq)
)





drop_na(df) %>%
  mutate(week = week(date), month = month(date), year= year(date)) %>%
  group_by(month) %>%
  reframe(de.sd = sd(d_excess),
          dO.sd = sd(dO)) %>%
  reframe(de.sd = mean(de.sd),
          dO.sd = mean(dO.sd))


drop_na(ddf) %>%
  group_by(gwt) %>%
  reframe(de.sd = mean(de.sd, na.rm=T),
          dO.sd = mean(dO.sd, na.rm=T))



ddf2 <- ddf %>%
  group_by(year, month) %>%
  reframe(de = weighted.mean(de, P, na.rm=T),
          dO = weighted.mean(dO, P, na.rm=T),
          P = sum(P, na.rm=T),
          count = n())


T = 12
x <- seq(1, length(ddf2$dO))


nls(dO~A*sin(2*pi*x*(1/T)-P)+k, data=ddf2, start=list(A=2,P=1, k=0), weights = P)         
nls(de~A*sin(2*pi*x*(1/T)-P)+k, data=ddf2, start=list(A=2,P=1, k=0), weights = P)





ggplot(ddf) +
  geom_point(aes(x=dO, y=de)) +
  theme_bw() + facet_wrap(~gwt)


cor.test(ddf$dO, ddf$de)


df_cor <- drop_na(ddf) %>% 
  group_by(gwt) %>%
  reframe(R = round(cor.test(dO, de)$estimate, 2),
          p = round(cor.test(dO, de)$p.value, 3))



df_results.de <- cbind(df_results.de, df_cor[,c(2,3)])







setwd("C:/Users/turk/Documents/MUSES/Isotope_model")


ggsave(plot_gwt, path="figures", filename = "sinfits.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")




write.table(df_results.dO, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\sinfits_dO.csv', sep = ",", quote = FALSE)
write.table(df_results.de, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\output_data\\sinfits_dex.csv', sep = ",", quote = FALSE)








