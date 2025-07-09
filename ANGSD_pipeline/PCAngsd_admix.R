# Create admixture plots for each admix_K value


library(ggplot2)
library(tidyr)
library(dplyr)
library(pals)
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/pcangsd/admix")

# Load the admixture data for K=2
#admix_data <- read.table("Nleu_autos_rm_inversions_filt_excl_1st_thin_pcangsd_K2_out.admix.2.Q", header = FALSE)
admix_data <- read.table("Nleu_autos_rm_inversions_filt_excl_1st_thin_39scaf_pcangsd_K2_out.admix.2.Q", header = FALSE)


pop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
                   "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee",
                   "Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom",
                   "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
                   "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")


admix_data_1 <- admix_data %>% mutate(Individual = paste0("Ind", row_number()))

admix.id = as.data.frame(cbind(pop, admix_data_1))
names(admix.id) = c("pop","q1","q2","Individual")

admix.id.sort <- admix.id[order(-admix.id$q1), ] 


# Sort by population, then by q1
admix.id.sort <- admix.id[order(-admix.id$q1), ]  


plot = barplot(t(as.matrix(subset(admix.id.sort, select=q1:q2))), col=c("firebrick","royalblue"), border=NA)

#admix_sorted <- admix.id %>%
#  arrange(desc(q1))

admix_data <- admix.id %>%
  arrange(desc(q2)) %>%
  mutate(
    group_id = as.numeric(factor(pop, levels = unique(pop))))%>%
  arrange(group_id)%>%
  mutate(x_pos = row_number() + (group_id - 1) * 1.1)  # Add space between groups
  #mutate(x_pos = row_number())


# Pivot to long format
admix_long <- admix_data %>%
  pivot_longer(
    cols = starts_with("q"),  # Columns to pivot (q1, q2, etc.)
    names_to = "population",   # Name for the variable storing column names
    values_to = "proportion"  # Name for the variable storing values
  )

# View the long-format data
print(admix_long)
admix_long$Individual <- factor(admix_long$Individual, levels = unique(admix_data$Individual))


pop_labels <- admix_long %>%
  group_by(pop) %>%
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    midpoint = (min(x_pos) + max(x_pos)) / 2,  # Midpoint for labeling
    .groups = "drop"
  )%>% as.data.frame()

pop_labels$pop_f <- factor(pop_labels$pop)


colors <- unname(glasbey()[1:14])

custom_colors<- c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300", "#FFD300",
                  "#009FFF","#9A4D42", "#00FFBE", "#783FC1","#1F9698" ,"#FFACFD","#B1CC71","#555555","lightgrey")

#custom_colors<- c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300", "#FFD300",
#                  "#009FFF","#9A4D42", "#00FFBE", "#783FC1","#1F9698" ,"#FFACFD","#B1CC71")

#custom_colors <- c("blue","lightblue","pink","orange","yellow","purple","maroon",
#                   "red","green","darkgray","black","magenta","brown","darkorange","darkblue","darkgreen")


admix_autos <- ggplot(admix_long, aes(x = x_pos, y = proportion)) +
  geom_bar(aes(fill=population),stat = "identity", width = 0.9,show.legend=FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text for individuals
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Location", y = "Admixture Proportion", fill = "Location") +
  ggtitle("Admixture Analysis (K=2)")+
  ggtitle("A.")+
  geom_rect(data = pop_labels, aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -0.05, ymax = 0,fill=pop_f), 
            inherit.aes = FALSE,alpha = 1,show.legend=FALSE)+
  scale_fill_manual(values = custom_colors)+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=16),
        )+
  geom_text(
    data = pop_labels,
    aes(
      x = (start + end) / 2,
      y = -.1, 
      label = pop
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface="bold"
   # fill = "white"   # a label's background color
  )

admix_autos

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Figure5_admixture_autos.png", 
       plot = admix_autos, dpi = 300, width = 10, height = 8, units = "in")


admix_long$population







############################################################################################################

###### Load the admixture data for K=3
admix_data <- read.table("Nleu_autos_rm_inversions_filt_excl_1st_thin_pcangsd_K3_out.admix.3.Q", header = FALSE)
admix_data <- read.table("Nleu_autos_rm_inversions_filt_excl_1st_thin_39scaf_pcangsd_K3_out.admix.3.Q", header = FALSE)


pop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
         "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mall","Mall","Mall","Mall","Mall","Mall","Mall",
         "Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom",
         "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
         "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")


admix_data_1 <- admix_data %>% mutate(Individual = paste0("Ind", row_number()))

admix.id = as.data.frame(cbind(pop, admix_data_1))
names(admix.id) = c("pop","q1","q2","q3","Individual")

admix.id.sort <- admix.id[order(-admix.id$q1), ] 


# Sort by population, then by q1
admix.id.sort <- admix.id[order(-admix.id$q1), ]  

plot = barplot(t(as.matrix(subset(admix.id.sort, select=q1:q3))), col=c("firebrick","royalblue","purple"), border=NA)


admix_data <- admix.id %>%
  arrange(desc(q1)) %>%
  mutate(
    group_id = as.numeric(factor(pop, levels = unique(pop))))%>%
  arrange(group_id)%>%
  mutate(x_pos = row_number() + (group_id - 1) * 1.5)  # Add space between groups
#mutate(x_pos = row_number())


# Pivot to long format
admix_long <- admix_data %>%
  pivot_longer(
    cols = starts_with("q"),  # Columns to pivot (q1, q2, etc.)
    names_to = "population",   # Name for the variable storing column names
    values_to = "proportion"  # Name for the variable storing values
  )

# View the long-format data
print(admix_long)
admix_long$Individual <- factor(admix_long$Individual, levels = unique(admix_data$Individual))


pop_labels <- admix_long %>%
  group_by(pop) %>%
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    midpoint = (min(x_pos) + max(x_pos)) / 2,  # Midpoint for labeling
    .groups = "drop"
  )


custom_colors <- c("blue","lightblue","pink","orange","yellow","purple","maroon",
                   "red","green","darkgray","black","tan","magenta","brown","darkorange","darkblue","darkgreen")

ggplot(admix_long, aes(x = x_pos, y = proportion, fill = population)) +
  geom_bar(stat = "identity", width = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text for individuals
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Individuals", y = "Admixture Proportion", fill = "Population") +
  ggtitle("Admixture Analysis (K=3)")+
  geom_rect(data = pop_labels, aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -0.05, ymax = 0,fill=pop), 
            inherit.aes = FALSE,alpha = 1)+
  scale_fill_manual(values = custom_colors)



###############################################################################################################
######  Load the admixture data for K=2 neo-Z
###### neo-Z 

admix_data <- read.table("Nleu_neoZ_males_rm_filt_excl_1st_order_thin_pcangsd_K2_out.admix.2.Q", header = FALSE)


pop <- c("Nom","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Bog","Weddin","Gund","Gund","Gund","Mallee","Mallee",
                   "Mallee","Mallee","Mallee","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom",
                   "Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Walch","Walch","Walch","Walch","Zost","Zost")


admix_data_1 <- admix_data %>% mutate(Individual = paste0("Ind", row_number()))

admix.id = as.data.frame(cbind(pop, admix_data_1))
names(admix.id) = c("pop","q1","q2","Individual")

admix.id.sort <- admix.id[order(-admix.id$q1), ] 


# Sort by population, then by q1
admix.id.sort <- admix.id[order(-admix.id$q1), ]  


plot = barplot(t(as.matrix(subset(admix.id.sort, select=q1:q2))), col=c("firebrick","royalblue"), border=NA)

#admix_sorted <- admix.id %>%
#  arrange(desc(q1))

admix_data <- admix.id %>%
  arrange(q1) %>%
  mutate(
    group_id = as.numeric(factor(pop,levels=c("Mallee","Binya","Nom","Weddin","Zost","Gund","Pill","Mull","Reedy","Moon","Bog","Walch"))))%>%
  arrange(group_id)%>%
  mutate(x_pos = row_number() + (group_id - 1) * 1.5)  # Add space between groups
#mutate(x_pos = row_number())


# Pivot to long format
admix_long <- admix_data %>%
  pivot_longer(
    cols = starts_with("q"),  # Columns to pivot (q1, q2, etc.)
    names_to = "population",   # Name for the variable storing column names
    values_to = "proportion"  # Name for the variable storing values
  )

# View the long-format data
print(admix_long)
admix_long$Individual <- factor(admix_long$Individual, levels = unique(admix_data$Individual))


pop_labels <- admix_long %>%
  group_by(pop) %>%
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    midpoint = (min(x_pos) + max(x_pos)) / 2,  # Midpoint for labeling
    .groups = "drop"
  )%>% as.data.frame()


pop_labels$pop_f <- factor(pop_labels$pop)

#custom_colors <- c("blue","lightblue","pink","yellow","purple","maroon",
#                   "red","green","darkgray","black","magenta","darkorange","darkblue","darkgreen")

Mcustom_colors <- c("#0000FF","#FF0000","#00FF00","#FF00B6","#005300","#FFD300","#009FFF","#9A4D42",
                    "#00FFBE","#1F9698","#FFACFD","#B1CC71","#555555","lightgrey")

admix_neoZ <- ggplot(admix_long, aes(x = x_pos, y = proportion)) +
  geom_bar(aes(fill=population),stat = "identity", width = 0.9,show.legend = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text for individuals
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Location", y = "Admixture Proportion", fill = "Population") +
  ggtitle("Admixture Analysis (K=2)")+
  ggtitle("B.")+
  geom_rect(data = pop_labels, aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -0.05, ymax = 0,fill=pop_f), 
            inherit.aes = FALSE,alpha = 1,show.legend=FALSE)+
  scale_fill_manual(values = Mcustom_colors)+
  geom_text(
    data = pop_labels,
    aes(
      x = (start + end) / 2,
      y = -.1, 
      label = pop
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface="bold"
   # fill = "white"   # a label's background color
  )+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=16),
  )

p <- admix_autos / admix_neoZ

p


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Figure5_admixture_autos_neoZ.png", 
       plot = p, dpi = 300, width = 8.5, height = 10, units = "in")


######################################################################################################
#### Load the admixture for k=3 neo-Z


admix_data <- read.table("Nleu_neoZ_males_rm_filt_excl_1st_order_thin_pcangsd_K3_out.admix.3.Q", header = FALSE)


pop <- c("Nom","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Bog","Weddin","Gund","Gund","Gund","Mallee","Mallee",
         "Mallee","Mallee","Mallee","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom",
         "Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Walch","Walch","Walch","Walch","Zost","Zost")


admix_data_1 <- admix_data %>% mutate(Individual = paste0("Ind", row_number()))

admix.id = as.data.frame(cbind(pop, admix_data_1))
names(admix.id) = c("pop","q1","q2","q3","Individual")

admix.id.sort <- admix.id[order(-admix.id$q1), ] 


# Sort by population, then by q1
admix.id.sort <- admix.id[order(-admix.id$q1), ]  


plot = barplot(t(as.matrix(subset(admix.id.sort, select=q1:q3))), col=c("firebrick","royalblue","tan"), border=NA)

#admix_sorted <- admix.id %>%
#  arrange(desc(q1))

admix_data <- admix.id %>%
  arrange(desc(q1)) %>%
  mutate(
    group_id = as.numeric(factor(pop,levels=c("Mallee","Binya","Nom","Weddin","Zost","Gund","Pill","Mull","Reedy","Moon","Bog","Walch"))))%>%
  arrange(group_id)%>%
  mutate(x_pos = row_number() + (group_id - 1) * 1.5)  # Add space between groups
#mutate(x_pos = row_number())


# Pivot to long format
admix_long <- admix_data %>%
  pivot_longer(
    cols = starts_with("q"),  # Columns to pivot (q1, q2, etc.)
    names_to = "population",   # Name for the variable storing column names
    values_to = "proportion"  # Name for the variable storing values
  )

# View the long-format data
print(admix_long)
admix_long$Individual <- factor(admix_long$Individual, levels = unique(admix_data$Individual))


pop_labels <- admix_long %>%
  group_by(pop) %>%
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    midpoint = (min(x_pos) + max(x_pos)) / 2,  # Midpoint for labeling
    .groups = "drop"
  )


custom_colors <- c("blue","lightblue","pink","yellow","purple","maroon",
                   "red","green","darkgray","black","tan","magenta","darkorange","darkblue","darkgreen")

ggplot(admix_long, aes(x = x_pos, y = proportion, fill = population)) +
  geom_bar(stat = "identity", width = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text for individuals
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Individuals", y = "Admixture Proportion", fill = "Population") +
  ggtitle("Admixture Analysis (K=3)")+
  geom_rect(data = pop_labels, aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -0.05, ymax = 0,fill=pop), 
            inherit.aes = FALSE,alpha = 1)+
  scale_fill_manual(values = custom_colors)
