
setwd("D:/Users/delan/files/GitHub/OPumilio_Aim2_CodeSharing")


#setwd("C:/Users/delco/Documents/GitHub/OPumilio_Aim2_CodeSharing")

library(ggplot2)
library(dplyr)
library(plotly)


## MD data 
# results1 = [repmat([sigF, sigM, optM, tMax, n, nLoci, nLociP, nBins, mutRt, rngSeed, condDens],[nPoints,1]), initCoords, results0];

# [sigF, sigM, optM, tMax, n, nLoci, nLociP, nBins, mutRt, rngSeed,...
#  condDens+1, initCoordT, initCoordP, corr_inner, path_length, end_T, end_P ]

# my_files <- list.files(pattern = "\\.csv$")
# my_data <- lapply(my_files, read.csv)

#df <- read.csv("results2_4x4_08-17-2025 07-13.csv", header=FALSE)

df <- read.csv("results2_4x4_08-17-2025 07-13.csv", header=FALSE) %>% 
  mutate_at(c('V1','V2','V3','V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11','V12','V13'), as.numeric) %>% 
  rename(
    sigF = V1,
    sigM = V2,
    optM = V3,
    tMax = V4,
    n = V5,
    nLociT = V6,
    nLociP = V7,
    nBins = V8,
    mutRt = V9, 
    condDens = V10,
    rngSeed = V11, # swapped w conddens in next iter
    initT = V12,
    initP = V13,
    corrIn = V14,
    pathL = V15,
    endT = V16,
    endP = V17
  )  %>% 
  mutate(
    x0 = 0,
    x1 = 1,
    y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
    y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
  as.data.frame()

vnBins = as.numeric(levels(as.factor(df$nBins))) #c(2, 3, 4, 5, 10, 100)

dflist <- list()
for (i in 1:length(vnBins)){
  dflist[[i]] <- df %>% filter(nBins %in% vnBins[i]) 
}
# df_bins2 <- df %>% filter(nBins %in% 2)
# df_bins3 <- df %>% filter(nBins %in% 3)
# df_bins4 <- df %>% filter(nBins %in% 4)
# df_bins5 <- df %>% filter(nBins %in% 5)
# df_bins10 <- df %>% filter(nBins %in% 10)
# df_binsL <- df %>% filter(nBins %in% 100)

#   dflist[[1]] <- df_bins2
#   dflist[[2]] <- df_bins3
#   dflist[[3]] <- df_bins4
#   dflist[[4]] <- df_bins5
#   dflist[[5]] <- df_bins10
#   dflist[[6]] <- df_binsL

# making sure all data exists:  
df %>% #mutate(
      #x0 = 0,
      #x1 = 1,
      #y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
      #y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
    ggplot(aes(x=initT, y=initP, shape = factor(mutRt))) + 
    geom_jitter(width=0.05, height=0.05)+
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1))


# ----------------
# vnBins = c(2, 3, 4, 5, 10, 100)
index <- 1 #of above
dfp <- dflist[[index]] %>% mutate(
  x0 = 0,
  x1 = 1,
  y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
  y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
  as.data.frame()
  
ggplot(dfp, aes(x=initT, y=initP, color=pathL)) + 
  ggtitle(paste( "nBins = ", vnBins[index]))+
  scale_color_continuous(limits=c(3,15))+
  coord_cartesian(xlim= c(0, 1),ylim = c(0, 1), expand=FALSE) +
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1), color='gray')+
  geom_point()+
  geom_point(aes(x=endT,y=endP), shape=4)+
  geom_segment(aes(xend=endT,yend=endP), linetype="dotted")+
  facet_grid(rows=vars(dfp$mutRt), cols=vars(dfp$sigF))+ # , as.table = FALSE
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  # trick ggplot into using same num-scale on axes with major, while using minor for the bin lines:
  # need to move minor lines slightly so not deleted when overlap with major (intended behavior)
  # this is just less work now vs. making new data set to plot grid
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color = "gray", linetype = 2))+ # 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))




# dfp %>% ggplot(aes(x=initT, y=initP, shape = factor(mutRt))) + 
#   geom_jitter(width=0.05, height=0.05)+
#   geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1))#df_bins[1,2,3]
# -------------

# check patl L more granular at mutrt = 0
# full data is in df
df %>% 
  filter(mutRt == 0.001) %>% 
  mutate(
  x0 = 0,
  x1 = 1,
  y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
  y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
  as.data.frame() %>%
  ggplot(aes(x=initT, y=initP, color=pathL)) + 
  ggtitle("mutRt = 0; sigF vs nBins")+
  scale_color_continuous(limits=c(3.5,7))+
  coord_cartesian(xlim= c(0, 1),ylim = c(0, 1), expand=FALSE) +
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1), color='gray')+
  geom_point()+
  geom_point(aes(x=endT,y=endP), shape=4)+
  geom_segment(aes(xend=endT,yend=endP), linetype="dotted")+
  facet_grid(rows=vars(mutRt), cols=vars(sigF))+ # , as.table = FALSE
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  # trick ggplot into using same num-scale on axes with major, while using minor for the bin lines:
  # need to move minor lines slightly so not deleted when overlap with major (intended behavior)
  # this is just less work now vs. making new data set to plot grid
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color = "gray", linetype = 2))+ # 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))

# all path lengths per mutation rate
ggplot(df, aes(x=corrIn, color=as.factor(mutRt), fill=as.factor(mutRt))) + 
  geom_density(alpha=0.5)+ 
  ylim(c(0,25))+
  xlim(c(-0.5,0.5))+
  facet_grid(rows=vars(nBins), cols=vars(sigF))


df %>% 
  mutate(
    x0 = 0,
    x1 = 1,
    y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
    y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
  as.data.frame() %>%
  ggplot(aes(x=initT, y=initP, color=as.factor(nBins) )) + 
  ggtitle("mutRt = 0; sigF vs nBins")+
  coord_cartesian(xlim= c(0, 1),ylim = c(0, 1), expand=FALSE) +
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1), color='gray')+
  geom_point()+
  geom_point(aes(x=endT,y=endP), shape=4)+
  geom_segment(aes(xend=endT,yend=endP), linetype="dotted")+
  facet_grid(rows=vars(mutRt), cols=vars(sigF))+ # , as.table = FALSE
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  # trick ggplot into using same num-scale on axes with major, while using minor for the bin lines:
  # need to move minor lines slightly so not deleted when overlap with major (intended behavior)
  # this is just less work now vs. making new data set to plot grid
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color = "gray", linetype = 2))+ # 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))

dfsig01 %>% 
  filter(mutRt == 0.001) %>% 
  as.data.frame() %>%
  ggplot(aes(x=initT, y=initP, color=corrIn)) + 
  scale_color_continuous(limits=c(-0.2,0.5))+
  coord_cartesian(xlim= c(0, 1),ylim = c(0, 1), expand=FALSE) +
  geom_point()+
  geom_point(aes(x=endT,y=endP), shape=4)+
  geom_segment(aes(xend=endT,yend=endP), linetype="dotted")+
  facet_grid(~nBins)+ # , as.table = FALSE
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  # trick ggplot into using same num-scale on axes with major, while using minor for the bin lines:
  # need to move minor lines slightly so not deleted when overlap with major (intended behavior)
  # this is just less work now vs. making new data set to plot grid
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color = "gray", linetype = 2))+ # 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.001,1.001, 1/vnBins[index]))

# -------------------
#correlation
ggplot(dfsig01, aes(x=as.factor(initT), y=as.factor(initP), fill=corrIn)) + 
  geom_tile()+
  scale_fill_continuous(limits=c(-0.2, 0.5))+ #df instead of dfp is INTENTIONAL
  facet_grid(rows=vars(dfsig01$mutRt), cols=vars(dfsig01$nBins))+ # , as.table = FALSE
  theme_bw()

# path L - going to have a hard time comparing these? may not be good measure
ggplot(dfp, aes(x=as.factor(initT), y=as.factor(initP), fill=pathL)) + 
  geom_tile()+
  scale_fill_continuous(limits=range(df$pathL))+ #df instead of dfp is INTENTIONAL
  facet_grid(rows=vars(dfp$mutRt), cols=vars(dfp$sigF))+ # , as.table = FALSE
  theme_bw()

ggplot(dfsig01, aes(x=as.factor(initT), y=as.factor(initP), fill=corrIn))+ 
  geom_point()+
  geom_point(aes(x=endT,y=endP), shape=4)+
  geom_segment(aes(xend=endT,yend=endP), linetype="dotted")+
  facet_grid(rows=vars(dfsig01$mutRt), cols=vars(dfsig01$nBins))


