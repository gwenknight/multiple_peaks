#### Analysis of 8 strains that are not in MACOTRA
# Use to look at standard optical density (OD) data vs calorimeter data (CS)

### Explored these to look at linear relationship 

#################**************** (1) Libraries, data and code needed *******************###############
library(tidyverse) 
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(ggplot2)
library(patchwork) # for combining plots
library(zoo)
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

#####*************************** READ IN DATA *******************###############
ddm_orig <- read.csv("output/cut_all_time_series_fit_params.csv")[,-1]

ddm <- ddm_orig %>% filter(source == "Other")

param <- read.csv("output/cut_all_model_fit_params.csv")[,-1]

strains <- unique(ddm$strain)

ggplot(param %>% filter(strain_name %in% strains), aes(x= inocl, y = timepeak, group = strain_name)) + 
  geom_point(aes(col = strain_name)) + 
  facet_wrap(~strain_name, nrow = 2) + theme(legend.position = "none") + 
  scale_x_continuous(breaks = seq(2:7), labels = function(x) parse(text=paste("10^",x)),"Inoculum", limits = c(1.5,6.5)) + 
  scale_y_continuous(expression(paste("Time to first peak (", italic(tmax(h)),")")))
ggsave("proof_of_principle/time_peak_vs_inoculum_strain_col.png", width = 12, height = 12)

ggplot(param %>% filter(strain_name %in% strains), aes(x= inocl, y = v_m_h_flow, group = strain_name)) + 
  geom_point(aes(col = strain_name)) + 
  facet_wrap(drytime~strain_name) + theme(legend.position = "none") + 
  geom_smooth(method = lm) + 
  scale_x_continuous(breaks = seq(2:7), labels = function(x) parse(text=paste("10^",x)),"Inoculum", limits = c(1.5,6.5)) + 
  scale_y_continuous("Height of highest peak")
ggsave("plots/height_peak_vs_inoculum_strain_col.png", width = 12, height = 12)


ga <- ggplot(param %>% filter(strain_name %in% strains), aes(x = inocl, y = timepeak, group = inocl)) + 
  geom_boxplot() + 
  scale_x_continuous(breaks = seq(2:7), labels = function(x) parse(text=paste("10^",x)),"Inoculum", limits = c(1.5,6.5)) + 
  scale_y_continuous(expression(paste("Time to first peak (", italic(tmax(h)),")")))
ggsave("proof_of_principle/time_peak_vs_inoculum.pdf")


## Data
data_od <- read_csv("proof_of_principle/data/ddm_OD.csv")[,-1]
data_cs <- read_csv("proof_of_principle/data/ddm_CS.csv")[,-1]


### Look at data
g1 <- ggplot(data_cs %>% filter(strain == "SA2704"), aes(x=Time, y = value_J, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Heat flow (mW)") +
  scale_x_continuous("Time (h)", lim = c(0,25)) + 
  #labs(subtitle = "A") +
  scale_color_discrete("Inoculum", labels = c(expression(paste(10^2)), expression(paste(10^3)), 
                                              expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6))))

### Smooth OD values
data_od <- data_od %>% group_by(rep, exp, variable, strain, inoc) %>%
  mutate(ma_value = rollapply(value, 10, mean,fill = NA),
         differ = c(0,diff(ma_value))) %>%
  ungroup()

g2 <- ggplot(data_od %>% filter(strain == "SA2704"), aes(x=Time, y = ma_value, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("OD600")+ 
  scale_x_continuous("Time (h)")+
  #labs(subtitle = "B") +
  scale_color_discrete("Inoculum", labels = c(expression(paste(10^2)), expression(paste(10^3)), 
                                              expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6))))

g3 <- ggplot(data_od %>% filter(strain == "SA2704"), aes(x=Time, y = differ, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("OD600")+ 
  scale_x_continuous("Time (h)")+
  #labs(subtitle = "D") +
  scale_color_discrete("Inoculum", labels = c(expression(paste(10^2)), expression(paste(10^3)), 
                                              expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6))))
ggsave("proof_of_principle/OD_data_difference.png")

### Cumulative heat flow values
g4 <- ggplot(data_cs %>% filter(strain == "SA2704"), aes(x=Time, y = csum, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Heat (mJ)") + 
  scale_x_continuous("Time (h)", lim = c(0,25)) + 
  #labs(subtitle = "C") +
  scale_color_discrete("Inoculum", labels = c(expression(paste(10^2)), expression(paste(10^3)), 
                                              expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6))))

gb <- ((g1 + g4) / (g3 + g2)) + plot_layout(guides = 'collect') #+ plot_annotation(tag_levels = 'A') & theme(legend.position='bottom')
(gb |  ga) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(legend.position='bottom')
ggsave("proof_of_principle/figure2.pdf",width = 16, height = 10)
ggsave("proof_of_principle/figure2.png",width = 16, height = 10)
