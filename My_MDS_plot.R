library(ggplot2)
library(dplyr)
library(ggrepel)
#sessionInfo()
data <- read.table("MDS_merge2.mds", header = T)
race<- read.table(file="racefile.txt", header=T)
datafile<- merge(data,race,by=c("IID","FID"))

png("MDS_plot2.png", height = 8, width = 10, units = "in", res = 300)
datafile %>%
  arrange(race) %>%
  ggplot(aes(x = C1, y = C2, color = race, shape = race)) +
  geom_point(size = 2) +
  scale_shape_manual(values=c(1, 1, 1, 1, 3)) +
  geom_text_repel(data=subset(datafile, race == "OWN" & !grepl("HG", datafile$IID)),
            aes(C1, C2, label=IID)) + 
  theme(text = element_text(size = 20))
dev.off()





png("MDS_plot3.png", height = 8, width = 10, units = "in", res = 300)
datafile %>%
  arrange(race) %>%
  ggplot(aes(x = C1, y = C2, color = race, shape = race)) +
  geom_point(size = 2) +
  scale_shape_manual(values=c(1, 1, 1, 1, 3)) +
  geom_text_repel(data = subset(datafile, race == "OWN" & !grepl("HG", datafile$IID) & !grepl("NA", datafile$IID)),
                  aes(C1, C2, label = IID),
                  color = "black",            # Set label color here
                  segment.color = "gray",   # Line color to connect label to point
                  segment.size = 0.5,
                  nudge_x = .01,
                  nudge_y = .05) +     # Thickness of the line
  theme(text = element_text(size = 20))
dev.off()


png("MDS_plot4.png", height = 8, width = 10, units = "in", res = 300)
datafile %>%
  arrange(race) %>%
  ggplot(aes(x = C3, y = C2, color = race, shape = race)) +
  geom_point(size = 2) +
  scale_shape_manual(values=c(1, 1, 1, 1, 3)) +
  geom_text_repel(data = subset(datafile, race == "OWN" & !grepl("HG", datafile$IID) & !grepl("NA", datafile$IID)),
                  aes(C3, C2, label = IID),
                  color = "black",            # Set label color here
                  segment.color = "gray",   # Line color to connect label to point
                  segment.size = 0.5,
                  nudge_x = .01,
                  nudge_y = .05) +     # Thickness of the line
  theme(text = element_text(size = 20))
dev.off()
