#### ORGANIZATION OF THE TABLE ####

# Load packages
library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)


# Data Reading
microData <- read.csv2("microData.csv")
head(microData)

# Column Selection
microDataHab <- microData[, c(5, 6, 8, 9, 10, 11)]

# Conversion to Binary Values
microDataHabPA <- ifelse(microDataHab == "A" | microDataHab == "" | microDataHab == "a", 0, 1)
microDataHabPA[(microDataHabPA[, 4] > 0) & (microDataHabPA[, 1] > 0), 1] <- 0

# Creation of the "Habitat" Column
microData$Habitat <- unlist(apply(microDataHabPA, 1, function(x) {
  resu <- c(1:6)[x > 0]
  if (length(resu) == 0) {
    resu <- NA
  }
  resu
}))

# Factorization of the "HabitatFac" Column
microData$HabitatFac <- factor(colnames(microDataHab)[microData$Habitat])
microData$N <- 1

# Removal of Rows with NA Values in "HabitatFac"
microData <- microData[!is.na(microData$HabitatFac), ]

# Calculation of Abundance and Richness per Habitat
microSpAb <- tapply(microData$N, list(paste0(microData$Parcela, "_", microData$HabitatFac), microData$Taxon), sum, default = 0)

microEnv <- tapply(microData$HabitatFac, paste0(microData$Parcela, "_", microData$HabitatFac), function(x) as.character(x[1]))

#### MICROSCALE ANALYSIS ####

# Abundance and Richness
N <- rowSums(microSpAb)
S <- rowSums(microSpAb > 0)

# PCOA
# Calculating Bray-Curtis dissimilarities
bray <- vegdist(microSpAb, "bray")

# Performing PCoA
pcoa <- cmdscale(bray,eig = TRUE,add=TRUE)

# Percentage of variance in each ordination axis
round(pcoa$eig/sum(pcoa$eig),2)

pcoaScores<-scores(pcoa)


## Data analysis

# N
summary(aovN <- aov(N ~ microEnv)) 


# S
summary(aovS <- aov(S ~ microEnv)) 

# PCOA
summary(aovPCOA <- aov(pcoaScores[,1] ~ microEnv)) 


# Graphics

# N

# Plotting Richness
new_names <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying wood")
custom_colors <- c("#1f7a1f", "#2e8b57", "#66cdaa", "#8b4513", "#daa520", "#694004") 


ggplot(data.frame(N = N, microEnv = factor(microEnv)), aes(x = microEnv, y = N, fill = microEnv)) +
  geom_violin(alpha = 0.6, trim = FALSE, scale = "area", width = 1.9) +  
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = custom_colors, labels = new_names) +  
  scale_x_discrete(labels = new_names) +
  labs(
    x = "Habitat Type",
    y = "Abundance",
    title = "Abundance Distribution by Habitat Type",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    legend.position = "none", 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor = element_blank()
  )



# S

# Plotting Richness
new_names <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying wood")
custom_colors <- c("#1f7a1f", "#2e8b57", "#66cdaa", "#8b4513", "#daa520", "#694004") 


ggplot(data.frame(S = S, microEnv = factor(microEnv)), aes(x = microEnv, y = S, fill = microEnv)) +
  geom_violin(alpha = 0.6, trim = FALSE, scale = "area", width = 1.9) +  
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = custom_colors, labels = new_names) +  
  scale_x_discrete(labels = new_names) +
  labs(
    x = "Habitat Type",
    y = "Species Richness",
    title = "Richness Distribution by Habitat Type",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    legend.position = "none", 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor = element_blank()
  )



# PCOA

# Plotting Richness
new_names <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying wood")
custom_colors <- c("#1f7a1f", "#2e8b57", "#66cdaa", "#8b4513", "#daa520", "#694004") 


ggplot(data.frame(pcoa = pcoaScores[,1], microEnv = factor(microEnv)), aes(x = microEnv, y = pcoa, fill = microEnv)) +
  geom_violin(alpha = 0.6, trim = FALSE, scale = "area", width = 1.9) +  
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = custom_colors, labels = new_names) +  
  scale_x_discrete(labels = new_names) +
  labs(
    x = "Habitat Type",
    y = "PCOA",
    title = "PCOA Distribution by Habitat Type",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    legend.position = "none", 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor = element_blank()
  )



#### ANALYSIS MESOSCALE ####

# Reading environmental data
envHar <- read.csv("mesoEscala.csv")
envPPbio <- read.csv("originalData/PPBio2014_Environment.RFAD.csv")

# Manipulating parcel data
envPPbio$Parcela <- paste(envPPbio$Trail, envPPbio$Plot)
envHar$Trail <- gsub("DU_(.)O([0-9])_(.*)", "\\1\\2", envHar$Parcela)
envHar$Plot <- as.numeric(gsub("DU_(.)O([0-9])_(.*)", "\\3", envHar$Parcela))
envHar$Parcela2 <- paste(envHar$Trail, envHar$Plot)


# Matching parcel data
envPPbioReordered <- envPPbio[match(envHar$Parcela2, envPPbio$Parcela),]



# Subsetting data
mesoSpAb0 <- envHar[, grepl("_", colnames(envHar))]
mesoSpAb <- mesoSpAb0[!is.na(envPPbioReordered$Location),] # Getting only data from Ducke forest
mesoEnv <- mesoEnv0[!is.na(mesoEnv0$Location),]

# N and S calculation
N <- rowSums(mesoSpAb)
S <- rowSums(mesoSpAb > 0)


# PCOA
# Calculating Bray-Curtis dissimilarities
bray <- vegdist(mesoSpAb, "bray")

# Performing PCoA
pcoa <- cmdscale(bray,eig = TRUE,add=TRUE)

# Percentage of variance in each ordination axis
round(pcoa$eig/sum(pcoa$eig),2)

pcoaScores<-scores(pcoa)

## Data analysis

# N
# Multiple linear regression for abundance
summary(lmN <- lm(N ~  Sand + P + Vegetation.structure + Vegetaion.structure2, data = mesoEnv)) 

# S
summary(lmS <- lm(S ~  Sand + P + Vegetation.structure + Vegetaion.structure2, data = mesoEnv)) 


# PCOA
summary(lmPCOA <- lm(pcoaScores[,1] ~  Sand + P + Vegetation.structure + Vegetaion.structure2, data = mesoEnv)) 


# Graphics
#N
ggNSand <- ggplot(mesoEnv, aes(x = Sand, y = N, color = "Sand")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#f4a261") +
  labs(title = "A)",
       x = "Amount of Sand",
       y = "N") +
  guides(color = FALSE)
  

ggNP <- ggplot(mesoEnv, aes(x = P, y = N, color = "Phosphorus")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#669bbc") +
  labs(title = "D)",
       x = "Phosphorus Quantity",
       y = "N") +
  guides(color = FALSE)

ggNTree <- ggplot(mesoEnv, aes(x = Vegetation.structure, y = N, color = "Trees")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#6a994e") +
  labs(title = "G)",
       x = "Number of Trees",
       y = "N") +
  guides(color = FALSE)


ggNPalm <- ggplot(mesoEnv, aes(x = Vegetaion.structure2, y = N, color = "Palms")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#344e41") +
  labs(title = "J)",
       x = "Quantity of Palm Trees",
       y = "N") +
  guides(color = FALSE)


#S
ggSSand <- ggplot(mesoEnv, aes(x = Sand, y = S, color = "Sand")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#f4a261") +
  labs(title = "B)",
       x = "Amount of Sand",
       y = "S") +
  guides(color = FALSE)


ggSP <- ggplot(mesoEnv, aes(x = P, y = S, color = "Phosphorus")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#669bbc") +
  labs(title = "E)",
       x = "Phosphorus Quantity",
       y = "S") +
  guides(color = FALSE)

ggSTree <- ggplot(mesoEnv, aes(x = Vegetation.structure, y = S, color = "Trees")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#6a994e") +
  labs(title = "H)",
       x = "Number of Trees",
       y = "S") +
  guides(color = FALSE)


ggSPalm <- ggplot(mesoEnv, aes(x = Vegetaion.structure2, y = S, color = "Palms")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#344e41") +
  labs(title = "K)",
       x = "Quantity of Palm Trees",
       y = "S") +
  guides(color = FALSE)


#pcoa
ggPcoaSand <- ggplot(mesoEnv, aes(x = Sand, y = pcoaScores[,1], color = "Sand")) +
  geom_point() +
  geom_smooth(method="lm")+
  scale_color_manual(values = "#f4a261") +
  labs(title = "C)",
       x = "Amount of Sand",
       y = "PCoA 1") +
  guides(color = FALSE)


ggPcoaP <- ggplot(mesoEnv, aes(x = P, y = pcoaScores[,1], color = "Phosphorus")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#669bbc") +
  labs(title = "F)",
       x = "Phosphorus Quantity",
       y = "PCoA 1") +
  guides(color = FALSE)

ggPcoaTree <- ggplot(mesoEnv, aes(x = Vegetation.structure, y = pcoaScores[,1], color = "Trees")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#6a994e") +
  labs(title = "I)",
       x = "Number of Trees",
       y = "PCoA 1") +
  guides(color = FALSE)


ggPcoaPalm <- ggplot(mesoEnv, aes(x = Vegetaion.structure2, y = pcoaScores[,1], color = "Palms")) +
  geom_point() +
  #geom_smooth(method="lm")+
  scale_color_manual(values = "#344e41") +
  labs(title = "L)",
       x = "Quantity of Palm Trees",
       y = "PCoA 1") +
  guides(color = FALSE)


title1 <- textGrob("Abundance (N)", gp = gpar(fontsize = 12, fontface = "bold"))
title2 <- textGrob("Richness (S)", gp = gpar(fontsize = 12, fontface = "bold"))
title3 <- textGrob("First Axis of the PCoA (PCoA 1)", gp = gpar(fontsize = 12, fontface = "bold"))

combined_plots <- grid.arrange(ggNSand, ggSSand, ggPcoaSand,ggNP, ggSP, ggPcoaP,ggNTree, ggSTree, ggPcoaTree, ggNPalm, ggSPalm, ggPcoaPalm, ncol = 3)

combined_plots <- grid.arrange(title1, title2, title3,  
  ggNSand, ggSSand, ggPcoaSand,
  ggNP, ggSP, ggPcoaP,
  ggNTree, ggSTree, ggPcoaTree,
  ggNPalm, ggSPalm, ggPcoaPalm,
  ncol = 3,
  layout_matrix = rbind( c(1, 2, 3),  
    c(4, 5, 6),  
    c(7, 8, 9),  
    c(10, 11, 12) 
)
