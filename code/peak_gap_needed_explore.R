### need a gap between double peaks? 

param_far <- read_csv("output/param_multiple_peaks_farapart.csv")[,-1]
param <- read_csv("output/param_multiple_peaks.csv")[,-1]


param_far <- read_csv("output/clustered_parameters_farapart.csv")
ddm_far <- read_csv("output/clustered_time_series_farapart.csv")

param <- read_csv("output/clustered_parameters.csv")
ddm <- read_csv("output/clustered_time_series.csv")

# How does forcing double peaks to be far apart change the distribution of clusters? 
table(param$cluster)
table(param_far$cluster)
# Time series
ggplot(ddm %>% filter(inoc == 5, drytime == 0), aes(x=Time, y = value, group = interaction(rep, strain))) + geom_line(aes(col = cluster)) + facet_wrap(~strain)
ggplot(ddm_far %>% filter(inoc == 5, drytime == 0), aes(x=Time, y = value, group = interaction(rep, strain))) + geom_line(aes(col = cluster)) + facet_wrap(~strain)


### second peakers
## in wide cluster
secondp <- unlist(param %>% filter(cluster == "wide", inoc == 5, drytime == 0, mp_h2 > 0) %>% dplyr::select(strain) %>% unique())

ggplot(ddm %>% filter(inoc == 5, drytime == 0, strain %in% secondp), aes(x=Time, y = value, group = interaction(rep, strain))) + geom_line(aes(col = cluster)) + 
  facet_wrap(~strain) + 
  geom_point(data = param %>% filter(inoc == 5, drytime == 0, mp_h2 > 0, strain %in% secondp), aes(x=mp_t2, y = mp_h2))

## in double
secondp <- unlist(param %>% filter(cluster == "spike", inoc == 5, drytime == 0, mp_h2 > 0) %>% dplyr::select(strain) %>% unique())

ggplot(ddm %>% filter(inoc == 5, drytime == 0, strain %in% secondp), aes(x=Time, y = value, group = interaction(rep, strain))) + geom_line(aes(col = cluster)) + 
  facet_wrap(~strain) + 
  geom_point(data = param %>% filter(inoc == 5, drytime == 0, mp_h2 > 0, strain %in% secondp), aes(x=mp_t2, y = mp_h2))

ddm %>% filter(inoc == 5, drytime == 0, strain %in% secondp)


