


setwd("C:/Users/delco/Documents/GitHub/OPumilio_Aim2_CodeSharing")

library(ggplot2)
library(dplyr)
library(plotly)


## MD data 
# results1 = [repmat([sigF, sigM, optM, tMax, n, nLoci, nLociP, nBins, mutRt, rngSeed, condDens],[nPoints,1]), initCoords, results0];

# [sigF, sigM, optM, tMax, n, nLoci, nLociP, nBins, mutRt, rngSeed,...
#  condDens+1, initCoordT, initCoordP, corr_inner, path_length, end_T, end_P ]
df <- read.csv("results2_4x4_08-17-2025 07-13.csv", header=FALSE)

df <- df %>% mutate_at(c('V1','V2','V3','V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11','V12','V13'), as.numeric) %>% #,'V14','V15','V16','V17'
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
    rngSeed = V10, # swapped w conddens in next iter
    condDens = V11,
    initT = V12,
    initP = V13,
    corrIn = V14,
    pathL = V15,
    endT = V16,
    endP = V17
  )  
df$condDens <- df$condDens - 1 #changed how this was input, but not output correcting for that 
df_bins2 <- df %>% filter(nBins %in% 2)
df_bins3 <- df %>% filter(nBins %in% 3)
df_bins4 <- df %>% filter(nBins %in% 4)
df_bins5 <- df %>% filter(nBins %in% 5)
df_bins10 <- df %>% filter(nBins %in% 10)
df_binsL <- df %>% filter(nBins %in% 100)
  

df %>% mutate(
      x0 = 0,
      x1 = 1,
      y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
      y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
    ggplot(aes(x=initT, y=initP, shape = factor(mutRt))) + 
    geom_jitter(width=0.05, height=0.05)+
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1))

dfp <- df_bins10 %>% mutate(
  x0 = 0,
  x1 = 1,
  y0 = (((sigF^2)/100) + 1)*(0 - (0.5*((sigF^2)/100))),
  y1 = (((sigF^2)/100) + 1)*(1 - (0.5*((sigF^2)/100)))) %>%
  as.data.frame()
  
dfp %>% ggplot(aes(x=initT, y=initP, shape = factor(mutRt))) + 
  geom_jitter(width=0.05, height=0.05)+
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1))#df_bins[1,2,3]
  

ggplot(dfp, aes(x=initT, y=initP)) + 
  # ylim(0,1)+
  # xlim(0,1)+
  coord_cartesian(xlim= c(0, 1),ylim = c(0, 1), expand=FALSE) +
  geom_segment(aes(x=x0, xend=x1, y = y0, yend = y1), color='gray')+
  geom_point()+
  geom_point(aes(x=endT,y=endP), shape=4)+
  geom_segment(aes(xend=endT,yend=endP), linetype="dotted")+
  facet_grid(rows=vars(dfp$mutRt), cols=vars(dfp$sigF))+ # , as.table = FALSE
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))
# maybe color as path length? cont scale

#correlation
ggplot(dfp, aes(x=as.factor(initT), y=as.factor(initP), fill=corrIn)) + 
  geom_tile()+
  scale_fill_continuous(limits=range(df$corrIn))+ #df instead of dfp is INTENTIONAL
  facet_grid(rows=vars(dfp$mutRt), cols=vars(dfp$sigF))+ # , as.table = FALSE
  theme_bw()

# path L - going to have a hard time comparing these? may not be good measure
ggplot(dfp, aes(x=as.factor(initT), y=as.factor(initP), fill=pathL)) + 
  geom_tile()+
  scale_fill_continuous(limits=range(df$pathL))+ #df instead of dfp is INTENTIONAL
  facet_grid(rows=vars(dfp$mutRt), cols=vars(dfp$sigF))+ # , as.table = FALSE
  theme_bw()

