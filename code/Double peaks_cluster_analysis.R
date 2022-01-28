##Can't do PCA on binary data
## Try this clustering from: https://towardsdatascience.com/hierarchical-clustering-on-categorical-data-in-r-a27e578f2995
library(tidyverse)

#----- Dummy Data -----#
# the data will be sterile clean in order to not get distracted with other issues that might arise, but I will also write about some difficulties I had, outside the codelibrary(dplyr)# ensuring reproducibility for sampling
set.seed(40)# generating random variable set
# specifying ordered factors, strings will be converted to factors when using data.frame()# customer ids come first, we will generate 200 customer ids from 1 to 200
id.s <- c(1:200) %>%
  factor()
budget.s <- sample(c("small", "med", "large"), 200, replace = T) %>%
  factor(levels=c("small", "med", "large"), ordered = TRUE)

origins.s <- sample(c("x", "y", "z"), 200, replace = T, prob = c(0.7, 0.15, 0.15))

area.s <- sample(c("area1", "area2", "area3", "area4"), 200, replace = T, prob = c(0.3, 0.1, 0.5, 0.2))

source.s <- sample(c("facebook", "email", "link", "app"), 200, replace = T, prob = c(0.1,0.2, 0.3, 0.4))## day of week - probabilities are mocking the demand curve

dow.s <- sample(c("mon", "tue", "wed", "thu", "fri", "sat", "sun"), 200, replace = T,
                prob = c(0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2)) %>%
  factor(levels=c("mon", "tue", "wed", "thu", "fri", "sat", "sun"), 
         ordered = TRUE)# dish 

dish.s <- sample(c("delicious", "the one you don't like", "pizza"), 200, replace = T)

# by default, data.frame() will convert all the strings to factors - DOESN"T seem to anymore (perhaps within tidyverse there is a clash?)
synthetic.customers <- data.frame(id.s, budget.s, origins.s, area.s, source.s, dow.s, dish.s)

synthetic.customers[sapply(synthetic.customers, is.character)] <- lapply(synthetic.customers[sapply(synthetic.customers, is.character)], as.factor) # need this to ensure factors

library(cluster) 
# to perform different types of hierarchical clustering
# package functions used: daisy(), diana(), clusplot()

gower.dist <- daisy(synthetic.customers[ ,2:7], metric = c("gower")) # class(gower.dist) 

## dissimilarity , dist

# The main input for the code below is dissimilarity (distance matrix)
# After dissimilarity matrix was calculated, the further steps will be the same for all data types
# I prefer to look at the dendrogram and fine the most appealing one first - in this case, I was looking for a more balanced one - to further continue with assessment
#------------ DIVISIVE CLUSTERING ------------#
divisive.clust <- diana(as.matrix(gower.dist), 
                        diss = TRUE, keep.diss = TRUE)
plot(divisive.clust, main = "Divisive")

##### how get the acutal clusters??

#Loading data
param$odd_peaks <- as.numeric(param$odd_peaks)
param$odd_shoulder <- as.numeric(param$odd_shoulder)
param$odd_shoulder_past <- as.numeric(param$odd_shoulder_past)
param$odd_width <- as.numeric(param$odd_width)
param$odd_minor_peak <- as.numeric(param$odd_minor_peak)
data<-param[,c("odd_peaks","odd_width","odd_shoulder","odd_shoulder_past","odd_minor_peak")]#iris[,-c(5)]

# To standarize the variables
data = scale(data)

# Assessing cluster tendency
if(!require(clustertend)) install.packages("clustertend")
library(clustertend)
# Compute Hopkins statistic for the dataset
set.seed(123)
hopkins(data, n = nrow(data)-1)
#Since the H value = 0.2809 which is far below the threshold 0.5, it is highly clusterable

###########################################################################
####################### K Means clustering ################################
###########################################################################

# K-mean - Determining optimal number of clusters
# NbClust Package : 30 indices to determine the number of clusters in a dataset
# If index = 'all' - run 30 indices to determine the optimal no. of clusters
# If index = "silhouette" - It is a measure to estimate the dissimilarity between clusters.
# A higher silhouette width is preferred to determine the optimal number of clusters

if(!require(NbClust)) install.packages("NbClust")
library(NbClust)
nb <- NbClust(data,  distance = "euclidean", min.nc=2, max.nc=15, method = "kmeans",
              index = "silhouette")
nb$All.index
nb$Best.nc # 15 clusters...

#Method II : Same Silhouette Width analysis with fpc package
library(fpc)
pamkClus <- pamk(data, krange = 2:15, criterion="multiasw", ns=2, critout=TRUE)
pamkClus$nc
cat("number of clusters estimated by optimum average silhouette width:", pamkClus$nc, "\n")

#Method III : Scree plot to determine the number of clusters
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(data,centers=i)$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(data,pamkClus$nc)

# get cluster means
aggregate(data,by=list(fit$cluster),FUN=mean)

# append cluster assignment
data <- data.frame(data, clusterid=fit$cluster)

### 15 = too many clusters

###########################################################################
####################### Hierarchical clustering############################
###########################################################################

# Hierarchical clustering - Determining optimal number of clusters
library(NbClust)
res<-NbClust(data, diss=NULL, distance = "euclidean", min.nc=2, max.nc=20,
             method = "ward.D2", index = "kl")
res$All.index
res$Best.nc # 2 clusters? 

# Ward Hierarchical Clustering
d <- dist(data, method = "euclidean")
fit <- hclust(d, method="ward.D2")
plot(fit) # display dendogram

# cluster assignment (members)
groups <- cutree(fit, k=2)
data = cbind(data,groups)

# draw dendogram with red borders around the 6 clusters
rect.hclust(fit, k=2, border="red")
