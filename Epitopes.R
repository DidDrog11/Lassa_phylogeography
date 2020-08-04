library('ggplot2')
library('tidyverse')
library('magrittr')
library('gridExtra')
library('TTR')
library('wesanderson')
library('RColorBrewer')

path <- getwd()
Bepipred <-read_excel("C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx")
Chou_Fasman <- readxl::read_xlsx("C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx", sheet = 2)
Emini <- readxl::read_xlsx(path = "C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx", sheet = 3)
Karplus_Schulz <- readxl::read_xlsx(path = "C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx", sheet = 4)
Kolaska_Tongaonkar <- readxl::read_xlsx(path = "C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx", sheet = 5)
Parker <- readxl::read_xlsx(path = "C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx", sheet = 6)
Sequence_identity <- readxl::read_xlsx(path = "C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/AA_identity.xlsx")
pal <- wes_palette("Zissou1", 100, type = "continuous")
pal2 <- brewer.pal(9,"YlOrRd")

Bepipred$Identity <- (1:length(Bepipred$Position))
Bepipred$Identity[1:length(Bepipred$Position)] <- Sequence_identity$Identity[1:length(Bepipred$Position)]
Bepipred1 <- print(ggplot(data = Bepipred, mapping = aes(x = Position,
                                      y = `Bepipred Score`,
                                      fill = Identity))+
  geom_col()+
  xlab('Amino Acid Residue')+
  ylab('Bepipred Score')+
  scale_fill_gradientn(colours = pal2))+
  theme_bw()+
  theme(legend.position="none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

ggsave("Bepipred1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Bepipred2 <- print(ggplot(data = Bepipred, mapping = aes(x = Position,
                                                   y = `Bepipred-2 Score`,
                                                   fill = Identity))+
  geom_col()+
  xlab('Amino Acid Residue')+
  ylab('Bepipred-2 Score')+
  coord_cartesian(ylim = c(0.35,0.65))+
  scale_fill_gradientn(colours = pal2))+
  theme(legend.position="none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank()
       )


Chou_Fasman$Chou_Fasman_Assignment <- Chou_Fasman$`Score (highlighted top 20%)`
Chou_Fasman$Chou_Fasman_Assignment <- ifelse(Chou_Fasman$Chou_Fasman_Assignment >= 1, 'E', '.')

Chou_Fasman1 <- print(ggplot(data = Chou_Fasman, mapping = aes(x = Position,
                                                               y = `Score (highlighted top 20%)`,
                                                               fill = Chou_Fasman_Assignment))+
                        geom_col()+
                        xlab('Amino Acid Residue')+
                        ylab('Chou Fasman Score')+
                        theme_bw()
                      )

Emini$Emini_Assignment <- Emini$`Score (highlighted top 20%)`
Emini$Emini_Assignment <- ifelse(Emini$Emini_Assignment >= 1.5, 'E','.')

Emini1 <- print(ggplot(data = Emini, mapping = aes(x = Position,
                                                   y = `Score (highlighted top 20%)`,
                                                   fill = Emini_Assignment))+
                  geom_col()+
                  xlab('Amino Acid Residue')+
                  ylab('Emini Score')+
                  theme_bw()
                )

Karplus_Schulz1 <- print(ggplot(data = Karplus_Schulz, mapping = aes(x = Position,
                                                            y = `Score (highlighted top 20%)`))+
                           geom_col()+
                           xlab('Amino Acid Residue')+
                           ylab('Karplus Schulz Score')+
                           theme_bw()
                          )

Kolaska_Tongaonkar$Identity <- (1:length(Kolaska_Tongaonkar$Position))
Kolaska_Tongaonkar$Identity[1:length(Kolaska_Tongaonkar$Position)] <- Sequence_identity$Identity[(1:length(Kolaska_Tongaonkar$Position))]
Kolaska_Tongaonkar1 <- print(ggplot(data = Kolaska_Tongaonkar, mapping = aes(x = Position,
                                                                             y = `Score (highlighted top 20%)`,
                                                                             fill = Identity))+
                               geom_col()+
                               xlab('Amino Acid Residue')+
                               ylab('K-T Score')+
                               coord_cartesian(ylim = c(0.75,1.2))+
                               scale_fill_gradientn(colours = pal2)+
                               theme(legend.position="none"))

#, guide = guide_colourbar(barwidth = 0.5, barheight = 10))
Sequence_identity1 <- print(ggplot(data = Sequence_identity, mapping = aes(x = Position,
                                                                           y = Identity
                                                                           ))+
                              geom_col()+
                              xlab('Amino Acid Residue')+
                              ylab('Sequence Identity')+
                              coord_cartesian(ylim = c(0,1))+
                              theme(legend.position="none"))


Sequence_identity2 <- print(ggplot(data = Sequence_identity, mapping = aes(x = Position,
                                                                           y = Identity))+
                              coord_cartesian(ylim = c(0,0.1))+
                              geom_ribbon(aes(ymin = 0, ymax = Identity))+
                              geom_col(aes(fill = Identity))+
                              scale_fill_gradient2(position = "bottom", low = "red", mid = muted("red"), high = "blue",
                                                   midpoint = 0.4)+
                              theme(legend.position = "none", axis.text.y = element_blank())
                              )


Parker1 <- print(ggplot(data = Parker, mapping = aes(x = Position,
                                                     y = `Score (highlighted top 20%)`))+
                               geom_col()+
                               xlab('Amino Acid Residue')+
                               ylab('Parker Score')+
                               theme_bw()
)

combined_epi <- readxl::read_xlsx(path = "C:/Users/ucbtds4/Google Drive/PhD/LIPS/Sequences/Epitopes/Potential_Epitopes.xlsx", sheet = 7)

data <- combined_epi %>%
  filter(Algorithm %in% c("Bepipred", "Bepipred_2", "Kolaska_Tongaonkar"))
data %>%
  ggplot(aes(x = Position, y = Score, group = Algorithm, fill = Algorithm))+
  geom_area()+
  theme(legend.position = "none")+
  ggtitle("Potential Epitope Region")+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1,"lines"),
    strip.text.x = element_text(size = 8),
    plot.title = element_text(size = 14)
  )+
  facet_wrap(~Algorithm, ncol = 1, scale = "free")

cut_y_axis <- grid.arrange(Bepipred1, Bepipred2, Kolaska_Tongaonkar1)

epitope_1_KT <- subset(Kolaska_Tongaonkar, `Score (highlighted top 20%)` >= 1.10)
epitope_2_B <- Bepipred %>%
  filter(`Bepipred Score` >= 0.33 | `Bepipred-2 Score` >= 0.5)
start <- 544
end <- 556
t_bl <- end-start
epitope_test <- print(epitope_2_B %>%
                        filter(Position >= start & Position <= end),
                      n = t_bl)
mean(epitope_test$Identity)

ggplot(combined_epi, aes(x = Position, y = Algorithm)) +
  geom_raster(aes(fill = Score)) +
  scale_fill_gradientn(colours = terrain.colors(10))
  
combined_epi$Scale_score <- scale(combined_epi$Score)
ggplot(combined_epi, aes(Position, Algorithm['Bepipred_2'], colour = Algorithm))+
  geom_raster(aes(fill = Score))+
  scale_fill_gradientn(colours = terrain.colors(10))


# Exploring as a moving average
aa_length <- 105
Bepipred$Bepipred_SMA <- ifelse(Bepipred$Position < aa_length, 'NA', Bepipred$`Bepipred Score` %>% TTR::SMA(., n = aa_length) %>% na.omit())
Bepipred$Bepipred_SMA <-Bepipred$`Bepipred Score` %>% TTR::SMA(., n = aa_length)

## best indexes
(which.max(aa_region)):(which.max(aa_region) + aa_length - 1)

plot(aa_region)
#sapply(1:569, function(x){sum(diff([x:(x+14)])*rollmean(scores[x:(x+14)],2))})


