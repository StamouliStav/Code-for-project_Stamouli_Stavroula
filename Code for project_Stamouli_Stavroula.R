rm(list=ls(all=TRUE)) #empty the environment
setwd("C:/Users/stavroula stamoulis/Downloads") #set the working directory
#Load all the packages needed
require(vegan)
require(ggpubr)
require(ggeffects)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(Hmisc)
library(cluster)
library(ggrepel)

####Perform Hierarchical Clustering
data1<- read.csv("Terrestrial_Community.csv") #load the data
sum(is.na(data1)) #check for missing values
rownames(data1) <- data1$Sample
bray  <- vegdist(data1[,-c(1)], method="bray", binary=FALSE) #Compute dissimilarity matrix

hc <- hclust(bray, method = "average") #hierarchical clustering with average linkage
num_clusters <- 3 # Three clusters are obvious

# Cut the dendrogram
clusters <- cutree(hc, k = num_clusters)

# Plot the dendrogram
plot(hc, main = "Dendrogram of Hierarchical Clustering using Bray-Curtis", xlab = "", sub = "")

# Add rectangles to the plot to indicate the clusters
rect.hclust(hc, k = num_clusters, border = 2:5)

#Calculate the Within Cluster Sum of Squares (WCSS), to see the heterogeneity of each cluster
#First create this function
wcss_per_cluster <- function(data, clusters) {
  unique_clusters <- unique(clusters)
  wcss_values <- numeric(length(unique_clusters))
  
  for (i in seq_along(unique_clusters)) {
    cluster_points <- data[clusters == unique_clusters[i], ]
    cluster_center <- colMeans(cluster_points)
    wcss_values[i] <- sum(rowSums((cluster_points - cluster_center)^2))
  }
  
  return(wcss_values)
}

# Then calculate WCSS for each cluster
wcss_values <- wcss_per_cluster(data1[,-c(1)], clusters)
print(wcss_values)

#Save the histogram in a png file using RMarkdown and the following code
#png("Hist.png", width=1800, height=1200, res=200)
#plot(hc, main = "Dendrogram of Hierarchical Clustering using Bray-Curtis", xlab = "", sub = "")
#rect.hclust(hc, k = num_clusters, border = 2:5)
#dev.off()

####Perform K-means analysis
#Analysis for the terrestrial community
data1<- read.csv("Terrestrial_Community.csv") 
sum(is.na(data1)) 
rownames(data1) <- data1$Sample
#Turn raw abundances to relative abundances
rownames(data1) <- data1$Sample
data1$Total_Abundance <- rowSums(data1[, -1])  

# Calculate relative abundances
relative_abundance_data <- data1[, -1] / data1$Total_Abundance
rownames(relative_abundance_data) <- data1$Sample
# Set the number of clusters
set.seed(123) # Setting seed for reproducibility
k <- 3 #What will be the composition of three clusters?

# Perform K-means clustering
kmeans_result <- kmeans(relative_abundance_data, centers = k)

# Print the results
print(kmeans_result)


####Perform PCoA
data1<- read.csv("Terrestrial_Community.csv") 
sum(is.na(data1)) 
rownames(data1) <- data1$Sample

bray  <- vegdist(data1[,-c(1)], method="bray", binary=FALSE) #dissimilarity matrix
pcoa <- cmdscale(bray, k=8, eig=TRUE) #compute the first 8 axes (less will be retained)
# Extract the scores 
pcoa_axes <- pcoa$points
colnames(pcoa_axes) <- paste0('bray_raw_pcoa_', 1:8)
# Convert the pcoa axis values to a data frame and label by sample
pcoa_axes_df <- data.frame(pcoa_axes)
pcoa_axes_df$Sample <- rownames(pcoa_axes)
# Add an index column to preserve original order
data1$Index <- seq_along(data1$Sample)
pcoa_axes_df$Index <- seq_along(pcoa_axes_df$Sample)

# Merge the data
merged_data <- merge(data1, pcoa_axes_df, by = "Sample")

# Reorder the merged dataset to match the original order of data1
merged_data <- merged_data[order(merged_data$Index.x), ]

# Remove the index columns -no longer needed
merged_data <- merged_data[, !grepl("Index", colnames(merged_data))]
Ful <-merged_data


#Add centroids in the plot
#Separate both dataframes into clusters
Ful$Clusters <- c("II", "II", "II", "Ib", "Ib", "Ib", "Ia", "Ib", "Ib", 
                  "Ia", "Ia", "Ia", "Ia", "Ia", "Ia", "Ia", "Ia",
                  "II", "II", "III", "II", "III", "II", "II", "II",
                  "III", "III", "III", "III", "III","III", "III")
pcoa_axes_df$Clusters <- c("II", "II", "II", "Ib", "Ib", "Ib", "Ia", "Ib", "Ib", 
                           "Ia", "Ia", "Ia", "Ia", "Ia", "Ia", "Ia", "Ia",
                           "II", "II", "III", "II", "III", "II", "II", "II",
                           "III", "III", "III", "III", "III","III", "III")
# Calculate centroids for each cluster
centroids <- pcoa_axes_df%>%
group_by(Clusters) %>%
  summarize(bray_raw_pcoa_1 = mean(bray_raw_pcoa_1),
            bray_raw_pcoa_2 = mean(bray_raw_pcoa_2))

# Print centroids
print(centroids)

plot <-ggplot(Ful, aes(bray_raw_pcoa_1, bray_raw_pcoa_2)) +
  geom_point(colour = "black") + 
  geom_point(data = centroids, aes(x = bray_raw_pcoa_1, y = bray_raw_pcoa_2, colour = Clusters), size = 2) + 
  stat_ellipse(data = Ful, aes(bray_raw_pcoa_1, y = bray_raw_pcoa_2, colour = Clusters, fill= Clusters), alpha = 0.2, geom = "polygon") +  # Add ellipses with 95% confidence level
  theme_classic() + 
  geom_text_repel(aes(label = Sample), 
                  size = 4,
                  box.padding = 0.1,  # Space between text and point
                  point.padding = 0.1,  # Space around points
                  nudge_y = 0.01,  # Adjust vertical positioning (if needed)
                  nudge_x = 0.001) +  # Adjust horizontal positioning (if needed)
  scale_colour_manual(values = c("Ia" = "forestgreen", "Ib" = "green", "II" = "red", "III" = "royalblue1" )) +  
  scale_fill_manual(values = c("Ia" = "forestgreen", "Ib" = "green", "II" = "red", "III" = "royalblue1" )) + 
  labs(title = "PCoA Plot with Centroids and Ellipses",
       x = "PCoA1",
       y = "PCoA2") +
  theme_classic()
plot  
par(mfrow=c(1,2))
eigenvalues <- pcoa$eig[pcoa$eig >0]
print(eigenvalues)
barplot(eigenvalues / sum(eigenvalues), main='Axis variation')
barplot(cumsum(eigenvalues)/ sum(eigenvalues), main='Cumulative variation')
# Print the percentage variation of the first 8 
head(sprintf('%0.2f%%', (eigenvalues / sum(eigenvalues)) * 100), n=8) # five axes should be retained
# Compute correlations between the original variables and the principal coordinates
loadings <- cor(data1[,-c(1)], pcoa_axes)
print(loadings)
loadings <- as.data.frame(loadings)
write.csv(loadings, file = "Loadings_pcoa.csv", row.names = TRUE) #store the outcome
dev.off()

####Perform NMDS
data1<- read.csv("Terrestrial_Community.csv") 
sum(is.na(data1)) #check for missing values
rownames(data1) <- data1$Sample
bray  <- vegdist(data1[,-c(1)], method="bray", binary=FALSE)
# Perform NMDS
nmds_result <- metaMDS(bray, k = 5)
# Convert NMDS result to data frame
nmds_data <- data.frame(nmds_result$points)

# Extract site scores
site_scores <- scores(nmds_result, display = "sites")

# Extract species scores
species_scores <- scores(nmds_result, display = "species")
print(species_scores)
ggplot(nmds_data, aes(x = MDS1, y = MDS2)) +
  geom_point() +
  geom_text(label = rownames(nmds_data), vjust = -0.5, hjust = -0.5, size = 3) +  
  ggtitle("NMDS Plot") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() #Same result as in PCoA


####Create taxa's plots
#The terrestrail community
data1<- read.csv("Ter_dist.csv") #dataset with raw abundances and distances (relative heights in meters) for each sample
sum(is.na(data1)) #check for missing values
data1$Total_Abundance <- rowSums(data1[, -1]) 

# Calculate relative abundances
relative_abundance_data <- data1[, -1] / data1$Total_Abundance

# Multiply each column by 100 to convert to percentages
relative_abundance_data <- relative_abundance_data * 100

# Add back the first column 
relative_abundance_data$Distances <- data1$Distances 

# Display the first few rows of the relative abundance data frame
head(relative_abundance_data)
#Groupings
relative_abundance_data$Cordaites <- relative_abundance_data$Cladaitina.veteadensis +relative_abundance_data$Inaperturopollenites + relative_abundance_data$Densipollenites + relative_abundance_data$Plicatipollenites
relative_abundance_data$Gnetopsida <- relative_abundance_data$Gnetaceaepollenites + relative_abundance_data$Ephedripites
relative_abundance_data$Cycads <- relative_abundance_data$Cycadopites + relative_abundance_data$Eucomiidites




####Create the plots
#Greens: the survivors
plot1 <- ggplot(relative_abundance_data, aes(x = Distances, y = Protohaploxypinus)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Protohaploxypinus", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) + 
  scale_y_continuous(
    expand = c(0, 0)  
  )
plot1
plot2 <- ggplot(relative_abundance_data, aes(x = Distances, y = Vittatina)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Vittatina", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0)  
  )
plot2
plot5 <- plot1 <- ggplot(relative_abundance_data, aes(x = Distances, y = Lunatisporites)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Lunatisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0)  
  )
plot5
plot21 <- ggplot(relative_abundance_data, aes(x = Distances, y = Striatoabietes)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Striatoabietes", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Striatoabietes), by = 2),  
    expand = c(0, 0)  # Ensure y-axis starts at 0
  )
plot21
plot3<- ggplot(relative_abundance_data, aes(x = Distances, y = Alisporites)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Alisporites", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Alisporites), by = 1.5),  
    expand = c(0, 0)  
  )
plot3
plot18 <-ggplot(relative_abundance_data, aes(x = Distances, y = Kraeuselisporites)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Kraeuselisporites", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Kraeuselisporites), by = 3),  
    expand = c(0, 0)  
  )
plot18
plot19 <- ggplot(relative_abundance_data, aes(x = Distances, y = Endosporites)) +
  geom_col(color = "#24ff24", fill = "#24ff24") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Endosporites", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Endosporites), by = 1.5),  
    expand = c(0, 0)  
  )
plot19
grid.arrange(plot19, plot18, plot3, plot21, plot5, plot2, plot1, nrow=7)
#Store this plot in a png file using Rmarkdown and the following code
#png("greens.png", width=2600, height=2200, res=200)
#grid.arrange(plot19, plot18, plot3, plot21, plot5, plot2, plot1, nrow=7)
#dev.off()

#Movs- Purple: the opportunistic 
plot14 <- ggplot(relative_abundance_data, aes(x = Distances, y = Cycads)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Cycads", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Cycads), by = 5),  
    expand = c(0, 0)  
  )
plot14
plot15<- ggplot(relative_abundance_data, aes(x = Distances, y = Gnetopsida)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Gnetopsida", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) + 
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Gnetopsida), by = 5), 
    expand = c(0, 0)  
  )
plot15
plot13 <-ggplot(relative_abundance_data, aes(x = Distances, y = Dictyotriletes)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Dictyotriletes", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Dictyotriletes), by = 1.5),  
    expand = c(0, 0)  
  )
plot13
plot6 <-ggplot(relative_abundance_data, aes(x = Distances, y = Indotriradites)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Indotriradites", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) + 
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Indotriradites), by = 2),  
    expand = c(0, 0)  
  )
plot6
plot12 <- ggplot(relative_abundance_data, aes(x = Distances, y = Lundbladispora)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Lundbladispora", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0)  
  )
plot12
plot4 <- ggplot(relative_abundance_data, aes(x = Distances, y = Densosporites..Densoisporites)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Densosporites/Desnoisp.", x = "Relative Height (m)"  
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Densosporites..Densoisporites), by = 5), 
    expand = c(0, 0)  
  )
plot4
plot10 <- ggplot(relative_abundance_data, aes(x = Distances, y = Secarisporites)) +
  geom_col(color = "#490092", fill = "#490092") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),  
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Secarisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0)  
  )
plot10
grid.arrange(plot10, plot4, plot12, plot6, plot13, plot15, plot14, nrow=7)
#Store this plot in a png file using Rmarkdown and the following code
#png("movs.png", width=2600, height=2200, res=200)
#grid.arrange(plot10, plot4, plot12, plot6, plot13, plot15, plot14, nrow=7)
#dev.off()

#Red: the victims
plot8 <- ggplot(relative_abundance_data, aes(x = Distances, y = Limitisporites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Limitisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Limitisporites), by = 1.5), 
    expand = c(0, 0) 
  )
plot8
plot9 <-ggplot(relative_abundance_data, aes(x = Distances, y = Weylandites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Weylandites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Weylandites), by = 1.5), 
    expand = c(0, 0) 
  )
plot9
plot17 <-ggplot(relative_abundance_data, aes(x = Distances, y = Jayantisporites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Jayantisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Jayantisporites), by = 2), 
    expand = c(0, 0) 
  )
plot17
plot11 <- ggplot(relative_abundance_data, aes(x = Distances, y = Cordaites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Cordaites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Cordaites), by = 5), 
    expand = c(0, 0) 
  )
plot11
plot7 <- ggplot(relative_abundance_data, aes(x = Distances, y = Sulcatisporites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Sulcatisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Sulcatisporites), by = 1.5), 
    expand = c(0, 0) 
  )
plot7
plot16 <- ggplot(relative_abundance_data, aes(x = Distances, y = Gondisporites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Gondisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Gondisporites), by = 2.5), 
    expand = c(0, 0) 
  )
plot16
plot23 <- ggplot(relative_abundance_data, aes(x = Distances, y = Lueckisporites)) +
  geom_col(color = "#920000", fill ="#920000") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Lueckisporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data$Lueckisporites), by = 1.5), 
    expand = c(0, 0) 
  )
plot23
grid.arrange(plot23, plot16, plot7, plot11, plot17, plot9, plot8,  nrow=7)
#Store this plot in a png file using Rmarkdown and the following code
#png("reds.png", width=2600, height=2200, res=200)
#grid.arrange(plot23, plot16, plot7, plot11, plot17, plot9, plot8,  nrow=7)
#dev.off()

####Create  taxa's plots
#The aquatic community
data_a<- read.csv("Aquatics_dist.csv") #dataset with raw abundances and distances (relative heights in meters) for each sample
sum(is.na(data_a)) #check for missing values
data_a$Total_Abundance <- rowSums(data_a[, -1])  

# Calculate relative abundances
relative_abundance_data_a <- data_a[, -1] / data_a$Total_Abundance

# Multiply each column by 100 to convert to percentages
relative_abundance_data_a <- relative_abundance_data_a * 100

# Add back the first column 
relative_abundance_data_a$Distances <- data_a$Distances 
#Groupings
relative_abundance_data_a$Other_Micr <- relative_abundance_data_a$Barathrispaeridium+ relative_abundance_data_a$Buedingiisphaeridium+relative_abundance_data_a$Micrhystridium.subgenus.Brachiprojectidium+ relative_abundance_data_a$Filisphaeridium

plot1 <- ggplot(relative_abundance_data_a, aes(x = Distances, y = Micrhystridium.breve)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Micrhystridium breve", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0) 
  )
plot1
plot2 <- ggplot(relative_abundance_data_a, aes(x = Distances, y = Micrhystridium.pentagonale)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Micrhystridium pentagonale", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data_a$Micrhystridium.pentagonale), by = 2), 
    expand = c(0, 0) 
  )
plot2
plot3 <- ggplot(relative_abundance_data_a, aes(x = Distances, y = Other_Micr)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Other Micr. Complex", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0) 
  )
plot3
plot4 <-ggplot(relative_abundance_data_a, aes(x = Distances, y = Veryhachium.laidii)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Veryhachium laidii", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data_a$Veryhachium.laidii), by = 1.5), 
    expand = c(0, 0) 
  )
plot4
plot5 <-ggplot(relative_abundance_data_a, aes(x = Distances, y = Verychachium.trispinosum)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Veryhachium trispinosum", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data_a$Verychachium.trispinosum), by = 1.5), 
    expand = c(0, 0) 
  )
plot5
grid.arrange(plot5, plot4, plot3, plot2, plot1,  nrow=5)
#Store this plot in a png file using Rmarkdown and the following code
#png("Ac1.png", width=2600, height=2200, res=200)
#grid.arrange(plot5, plot4, plot3, plot2, plot1,  nrow=5)
#dev.off()

plot6 <- ggplot(relative_abundance_data_a, aes(x = Distances, y = Leiosphaeridia)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Leiosphaeridia", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0) 
  )
plot6
plot7 <- ggplot(relative_abundance_data_a, aes(x = Distances, y = Tasmanites)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Tasmanites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data_a$Tasmanites), by = 1.5), 
    expand = c(0, 0) 
  )
plot7
plot8 <- ggplot(relative_abundance_data_a, aes(x = Distances, y = Maculatasporites)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Maculatasporites", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    expand = c(0, 0) 
  )
plot8
plot9 <-ggplot(relative_abundance_data_a, aes(x = Distances, y = Reduviasporonites)) +
  geom_col(color = "#006ddb", fill ="#006ddb") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1), 
    axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),  
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    y = "Reduviasporonites chalastus", x = "Relative Height (m)" 
  ) +
  scale_x_reverse(breaks = unique(relative_abundance_data_a$Distances)) +  
  scale_y_continuous(
    breaks = seq(0, max(relative_abundance_data_a$Reduviasporonites), by = 5), 
    expand = c(0, 0) 
  )
plot9
grid.arrange(plot9, plot8, plot7, plot6,  nrow=4)
#Store this plot in a png file using Rmarkdown and the following code
#png("Ac2.png", width=2600, height=2000, res=200)
#grid.arrange(plot9, plot8, plot6, plot5, nrow=4)
#dev.off()



####Shannon Indices
#Calculate Shannon Index for terrestrial community
#Analysis for the terrestrial community
data<- read.csv("Terrestrial_Community.csv")
sum(is.na(data))
data1 <- data
rownames(data1) <- data1$Sample

#Calculate Shannon Index 
data1 <- select(data1, -Sample)
data1$H <- diversity(data1, index='shannon')
data1$Sample <- rownames(data1)

#Store the results in order to merge with the indices calculated straight from the counts
write.csv(data1, file = "Shannon_Index_terrestrial_r.csv", row.names = FALSE)


###Calculate Shannon Index for pollen producing taxa
data_p<-read.csv("Revised_data_af.csv") #Dataset that has inside the type of palynomorph
data_p <- subset(data_p, Type == "Pollen")
sum(is.na(data_p))
data_p <- select(data_p, -Type)
#Pivot the data
long_data_p <- pivot_longer(data_p, cols = -Genus, names_to = "Sample", values_to="Value")
wide_data_p <- pivot_wider(long_data_p, names_from = Genus, values_from = Value)
data_p1 <- as.data.frame(wide_data_p)
rownames(data_p1) <- data_p1$Sample

##Calculate Shannon only for Pollen
data_p1 <- select(data_p1, -Sample)
data_p1$Hp <- diversity(data_p1, index='shannon')
data_p1$Sample <- rownames(data_p1)
write.csv(data_p1, file = "Shannon_Index_pollen_r.csv", row.names = FALSE) #store the outcome

###Calculate Shannon Index for spores producing taxa
data_s<- read.csv("Revised_data_af.csv") #same dataset with type of palynomoprhs
data_s <- subset(data_s, Type == "Spore")
sum(is.na(data_s))
data_s <- select(data_s, -Type)
#Pivot the data
long_data_s <- pivot_longer(data_s, cols = -Genus, names_to = "Sample", values_to="Value")
wide_data_s <- pivot_wider(long_data_s, names_from = Genus, values_from = Value)
data_s1 <- as.data.frame(wide_data_s)
rownames(data_s1) <- data_s1$Sample

##Calculate Shannon only for Spores
data_s1 <- select(data_s1, -Sample)
data_s1$Hs <- diversity(data_s1, index='shannon')
data_s1$Sample <- rownames(data_s1)
write.csv(data_s1, file = "Shannon_Index_spores_r.csv", row.names = FALSE) #store the outcome

###Calculate Shannon Index for marine plant community- Proven to be unreliable. High values were the artefact of low abundances
Mar<- read.csv("Marine.csv", heade=TRUE)
Mar <- select(Mar, -Type) #type is not needed here since he have only the marine palynomorphs
sum(is.na(Mar))
#Pivot the data
long_Mar <- pivot_longer(Mar, cols = -Genus, names_to = "Sample", values_to="Value")
wide_Mar <- pivot_wider(long_Mar, names_from = Genus, values_from = Value)
Mar1 <- as.data.frame(wide_Mar)
rownames(Mar1) <- Mar1$Sample

#Calculate Shannon Index 
Mar1 <- select(Mar1, -Sample)
Mar1$Hm <- diversity(Mar1, index='shannon')
Mar1$Sample <- rownames(Mar1)
write.csv(Mar1, file = "Shannon_Index_marine_r.csv", row.names = FALSE) #store the outcome



###Plot Indices all together
data2 <- read.csv("Indices_1.csv")
data2 <- dplyr::select(data2, -Sample)
plot1 <- ggplot(data2, aes(x = Distances, y = H))+
  geom_line(linetype = "solid", color = "#24ff24", linewidth = 1.5)+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + labs(x= "Relative Height (m)",
           y = "H"     # Set the y-axis title
  ) +  scale_x_continuous(
    limits = c(-3.72, 3.41),    # Set the x-axis limits
    breaks = seq(-4, 3.5, by = 0.5)  # Set the x-axis breaks
  ) +
  coord_flip()  # Flip coordinates

plot2 <- ggplot(data2, aes(x = Distances, y = Hp)) +
  geom_line(linetype = "solid", color ="#009292", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Hp"     # Set the y-axis title
  ) +
  coord_flip()

plot3 <- ggplot(data2, aes(x = Distances, y = Hs)) +
  geom_line(linetype = "solid", color ="#004949", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Hs"     # Set the y-axis title
  ) + scale_y_continuous(
    breaks = seq(1.50, max(data2$Hs), by = 0.25)) +
  coord_flip()
plot3

plot4 <- ggplot(data2, aes(x = Distances, y = Hm)) +
  geom_line(linetype = "solid", color ="#b6dbff", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + labs(x= "Relative Height (m)",
           y = "Hm"     # Set the y-axis title
  ) +  scale_x_continuous(
    limits = c(-3.72, 3.41),    # Set the x-axis limits
    breaks = seq(-4, 3.5, by = 0.5)  # Set the x-axis breaks
  ) +
  coord_flip()  # Flip coordinates
plot5 <- ggplot(data2, aes(x = Distances, y = MvsT)) +
  geom_line(linetype = "solid", color ="#490092",linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Marine/Terrestrial"     # Set the y-axis title
  ) +
  coord_flip()

plot6 <- ggplot(data2, aes(x = Distances, y = FvsT)) +
  geom_line(linetype = "solid", color ="#006ddb", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Freshwater/Terrestrial"     # Set the y-axis title
  ) +
  coord_flip()
plot7 <- ggplot(data2, aes(x = Distances, y = UnseparatedvsTotal)) +
  geom_line(linetype = "solid", color ="#ff6db6", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Unseparated/Total"     # Set the y-axis title
  ) +
  coord_flip()
plot8 <- ggplot(data2, aes(x = Distances, y = Bisaccate.Indet)) +
  geom_line(linetype = "solid", color ="#924900", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Bisaccate Indet"     # Set the y-axis title
  ) +
  coord_flip()

plot9 <- ggplot(data2, aes(x = Distances, y = Pollen.Indet)) +
  geom_line(linetype = "solid", color ="#db6d00", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Pollen Indet"     # Set the y-axis title
  ) +
  coord_flip()
plot10 <- ggplot(data2, aes(x = Distances, y = Spore.Indet)) +
  geom_line(linetype = "solid", color ="#b66dff", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Spores Indet"     # Set the y-axis title
  ) +
  coord_flip()
plot11 <- ggplot(data2, aes(x = Distances, y = Acritarchs.Indet)) +
  geom_line(linetype = "solid", color ="#b66dff", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ 
  coord_flip()

plot12 <- ggplot(data2, aes(x = Distances, y = SporesvsPollen)) +
  geom_line(linetype = "solid", color ="#920000", linewidth = 1.5) + theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+ labs(
    y = "Spores/Pollen"     # Set the y-axis title
  ) +
  coord_flip()

# Arrange plots in a grid
grid.arrange(plot1, plot2, plot3,plot12, plot7, ncol = 5)
grid.arrange(plot4, plot5, plot6, plot8, plot9,  plot10, ncol=6)

#Save those plots in png images using RMarkdown and the following code
#png("Ind1.png", width=1800, height=1200, res=200)
#grid.arrange(plot1, plot2, plot3,plot12, plot7,  ncol = 5)
#dev.off()

#png("Ind2.png", width=1800, height=1200, res=200)
#grid.arrange(plot5, plot6, plot8, plot9, plot10, ncol=5) #plot together with Pollen Indet and Spores Indet to get the same size and the crop
#dev.off()



####Test the correlations of the indices
hist(data2$Hs) #First quickly check normal distribution for the indices. Some of them deviate from normallity but Pearson is robust enough and results were not unexpected after plotting the indices.
# Calculate correlation matrix with p-values
correlation_matrix <- rcorr(as.matrix(data2))

# Extracting the correlation coefficients
correlations <- correlation_matrix$r
print(correlations)

# Extracting the p-values
p_values <- correlation_matrix$P
print(p_values)



