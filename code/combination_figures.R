#### Combination figures
var_name <- c(
  auc = "AUC",
  valpeak = "Max value 1st peak",
  exp_gr = "Max exp growth rate",
  second_peak_h = "Max value 2nd peak",
  timepeak = "Time max peak", 
  v_max_h_flow = "Max value 1st peak",
  exp_gr = "Max exp growth rate",
  v_min_h_flow = "Max value 2nd peak",
  mp_h2 = "Max value 2nd peak",
  t_max_h_flow = "Time max peak"
)



####### Figure 2
para <- read_csv("output/clustered_parameters.csv")[,c(-1,-2)]

w<-which(is.na(para$cluster))
para[w,"cluster"] <- "unclustered"

para <- para %>% ungroup() %>% mutate(second_peak_h = ifelse((t_m_h_flow > (timepeak + 2)), v_m_h_flow, ifelse(mp_t2 > (timepeak + 4), mp_h2, 0))) 

cluster_data <- para %>% dplyr::select(c("cluster", "inocl","drytime","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                         "shoulder_point_past_v","shoulder_point_past_t", "second_peak_h")) %>% 
  pivot_longer(cols = c(valpeak:second_peak_h)) 
cluster_data$cluster <- factor(cluster_data$cluster, levels = c("normal","double","spike","post_shoulder","wide","unclustered"))


g1a <- ggplot(cluster_data %>% filter(inocl %in% c(3,5), name %in% c("timepeak")), aes(x = cluster, y = value)) +
  geom_boxplot(aes(fill = cluster), show.legend = FALSE) +
  scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous("Value") +
  facet_grid(drytime + inocl ~ name, scales = "free_x", labeller = labeller(name = var_name)) +
  labs(fill = "Cluster") 

g2a <- ggplot(cluster_data %>% filter(inocl %in% c(3,5), name %in% c("auc")), aes(x = cluster, y = value)) +
  geom_boxplot(aes(fill = cluster), show.legend = FALSE) +
  scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous("Value") +
  facet_grid(drytime + inocl ~ name, scales = "free_x", labeller = labeller(name = var_name)) +
  labs(fill = "Cluster") 

g3a <- ggplot(cluster_data %>% filter(inocl %in% c(3,5), name %in% c("valpeak")), aes(x = cluster, y = value)) +
  geom_boxplot(aes(fill = cluster), show.legend = FALSE) +
  scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous("Value") +
  facet_grid(drytime + inocl ~ name, scales = "free_x", labeller = labeller(name = var_name)) +
  labs(fill = "Cluster") 

g4a <- ggplot(cluster_data %>% filter(inocl %in% c(3,5), name %in% c("second_peak_h")), aes(x = cluster, y = value)) +
  geom_boxplot(aes(fill = cluster), show.legend = FALSE) +
  scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous("Value") +
  facet_grid(drytime + inocl ~ name, scales = "free_x", labeller = labeller(name = var_name)) +
  labs(fill = "Cluster") 

g1a +  g2a  + g3a + g4a + plot_layout(ncol = 4)
ggsave("plots/final/figure2_col.png", width = 13, height = 13)


############### Glucose + dehydration 
#### No correction 
param <- read_csv("output/simple_extract_glucose_allinoc_hf.csv")[,-1]
glc_palette = c("#d7301f","#abd9e9","#2c7bb6","#238b45")


i = "no"
j = "all"
#### Plot over glucose concentrations
param_long <- param %>% pivot_longer(cols = t_max_h_flow:auc) %>% filter(baseline == i) %>% filter(name %in% c("t_max_h_flow" ,"auc", "v_max_h_flow","v_min_h_flow"))
param_long$name <- factor(param_long$name, levels = c("t_max_h_flow" ,"auc", "v_max_h_flow","v_min_h_flow"))

g1g <- ggplot(param_long %>% filter(name == "t_max_h_flow"), aes(x=inoc, y = value, group = interaction(inoc, glucose))) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 

g2g <- ggplot(param_long %>% filter(name == "auc"), aes(x=inoc, y = value, group = interaction(inoc, glucose))) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 

g3g <- ggplot(param_long %>% filter(name == "v_max_h_flow"), aes(x=inoc, y = value, group = interaction(inoc,glucose))) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 

g4g <- ggplot(param_long %>% filter(name == "v_min_h_flow"), aes(x=inoc, y = value, group = interaction(inoc,glucose))) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 


g1gs <- g1g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction(glucose), fill = factor(glucose)),method='loess', formula= y~x) 
g2gs <- g2g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction(glucose), fill = factor(glucose)),method='loess', formula= y~x) 
g3gs <- g3g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction( glucose), fill = factor(glucose)),method='loess', formula= y~x) 
g4gs <- g4g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction(glucose), fill = factor(glucose)),method='loess', formula= y~x)

g1gs + g2gs + g3gs + g4gs + plot_layout(ncol = 4, guides = "collect")
ggsave("plots/glucose_fig_smooth.png")

g1gb <- g1g + geom_boxplot(aes(col = factor(glucose)))
g2gb <- g2g + geom_boxplot(aes(col = factor(glucose)))
g3gb <- g3g + geom_boxplot(aes(col = factor(glucose)))
g4gb <- g4g + geom_boxplot(aes(col = factor(glucose)))

g1gb + g2gb + g3gb + g4gb + plot_layout(ncol = 4, guides = "collect")
ggsave("plots/glucose_fig_boxplot.png")

###### Figure 4
param <- read.csv("output/param_multiple_peaks.csv")

param_long <- param %>% filter(strain %in% c(11016, 11051, 11210, 11257)) %>% pivot_longer(cols = t_m_h_flow:mp_h2) %>% filter(name %in% c("t_m_h_flow" ,"auc", "v_m_h_flow","mp_h2"))
param_long$name <- factor(param_long$name, levels = c("t_m_h_flow" ,"auc", "v_m_h_flow","mp_h2"))

ggplot(param_long, aes(x=inocl, y = value, group = strain)) +
  facet_wrap(~name,  labeller = labeller(name = var_name), scales = "free", ncol = 4) +
  geom_point(aes(col = factor(strain))) + 
  geom_smooth(aes(col = factor(strain), group = strain, fill = factor(strain)),method='loess', formula= y~x) 

g1g <- ggplot(param_long %>% filter(name == "t_m_h_flow"), aes(x=inocl, y = value, group = strain)) +
  facet_grid( ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 

g2g <- ggplot(param_long %>% filter(name == "auc"), aes(x=inocl, y = value)) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 

g3g <- ggplot(param_long %>% filter(name == "v_m_h_flow"), aes(x=inocl, y = value, group = interaction(inocl,glucose))) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 

g4g <- ggplot(param_long %>% filter(name == "mp_h2"), aes(x=inocl, y = value, group = interaction(inocl,glucose))) +
  facet_grid(drytime ~ name, scales = "free_x",  labeller = labeller(name = var_name)) +
  scale_color_manual("Glucose\nconcentration", values = glc_palette) + scale_fill_manual("Glucose\nconcentration", values = glc_palette) + 
  scale_x_continuous("Inoculum") 


g1gs <- g1g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction(glucose), fill = factor(glucose)),method='loess', formula= y~x) 
g2gs <- g2g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction(glucose), fill = factor(glucose)),method='loess', formula= y~x) 
g3gs <- g3g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction( glucose), fill = factor(glucose)),method='loess', formula= y~x) 
g4gs <- g4g + geom_point(aes(col = factor(glucose))) + geom_smooth(aes(col = factor(glucose), group = interaction(glucose), fill = factor(glucose)),method='loess', formula= y~x)

g1gs + g2gs + g3gs + g4gs + plot_layout(ncol = 4, guides = "collect")
ggsave("plots/inoc_4_fig_smooth.png")

g1gb <- g1g + geom_boxplot(aes(col = factor(glucose)))
g2gb <- g2g + geom_boxplot(aes(col = factor(glucose)))
g3gb <- g3g + geom_boxplot(aes(col = factor(glucose)))
g4gb <- g4g + geom_boxplot(aes(col = factor(glucose)))

g1gb + g2gb + g3gb + g4gb + plot_layout(ncol = 4, guides = "collect")
ggsave("plots/inoc_4_fig_boxplot.png")
