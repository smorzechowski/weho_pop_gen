# Add lat long to inversion cluster files for ArcGIS Pro

library(dplyr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/local_pcangsd/tables")
pops <- read.csv('~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/all_samples_pops_lat_long.csv',header=T)


#####################################################################################
table <- read.table('scaffold_4_pop_cluster1_labels.tsv',header=T)
Chr4.1_combined <- left_join(pops,table)
Chr4.1_combined$Chr <- "Chr_4.1"

# Group by Pop and summarize cluster counts
result <- Chr4.1_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_4_pop_cluster1_updated_sum_labels_latlong.csv')



#####################################################################################
table <- read.table('scaffold_4_pop_cluster2_labels.tsv',header=T)
Chr4.2_combined <- left_join(pops,table)
Chr4.2_combined$Chr <- "Chr_4.2"

# Group by Pop and summarize cluster counts
result <- Chr4.2_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_4_pop_cluster2_updated_sum_labels_latlong.csv')



#####################################################################################
table <- read.table('scaffold_6_pop_cluster_labels.tsv',header=T)
Chr6_combined <- left_join(pops,table)
Chr6_combined$Chr <- "Chr_6"

# Group by Pop and summarize cluster counts
result <- Chr6_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_6_pop_cluster_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_8_pop_cluster1_labels.tsv',header=T)
Chr8.1_combined <- left_join(pops,table)
Chr8.1_combined$Chr <- "Chr_8.1"

# Group by Pop and summarize cluster counts
result <- Chr8.1_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_8_pop_cluster1_updated_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_8_pop_cluster2_labels.tsv',header=T) # NEED to recheck the clustering
Chr8.2_combined <- left_join(pops,table)
Chr8.2_combined$Chr <- "Chr_8.2"

# Group by Pop and summarize cluster counts
result <- Chr8.2_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_8_pop_cluster2_updated_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_9_pop_cluster_labels.tsv',header=T)
Chr9_combined <- left_join(pops,table)
Chr9_combined$Chr <- "Chr_9"

# actually the very end of Chr 9 pseudochromosome, just note!

# Group by Pop and summarize cluster counts
result <- Chr9_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_9_pop_cluster_sum_labels_latlong.csv')

#####################################################################################
table <- read.table('scaffold_10_pop_cluster_labels.tsv',header=T)
Chr10_combined <- left_join(pops,table)
Chr10_combined$Chr <- "Chr_10"

# Group by Pop and summarize cluster counts
result <- Chr10_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_10_pop_cluster_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_13_pop_cluster_labels.tsv',header=T)
Chr13_combined <- left_join(pops,table)
Chr13_combined$Chr <- "Chr_13"

# Group by Pop and summarize cluster counts
result <- Chr13_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_13_pop_cluster_sum_labels_latlong.csv')



#####################################################################################
table <- read.table('scaffold_14_pop_cluster_labels.tsv',header=T)
Chr14_combined <- left_join(pops,table)
Chr14_combined$Chr <- "Chr_14"


# Group by Pop and summarize cluster counts
result <- Chr14_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_14_pop_cluster_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_15pop_cluster_labels.tsv',header=T)
Chr15_combined <- left_join(pops,table)
Chr15_combined$Chr <- "Chr_15"

# Group by Pop and summarize cluster counts
result <- Chr15_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_15_pop_cluster_sum_labels_latlong.csv')



#####################################################################################
table <- read.table('scaffold_17_pop_cluster_labels.tsv',header=T)
Chr17_combined <- left_join(pops,table)
Chr17_combined$Chr <- "Chr_17"

# Group by Pop and summarize cluster counts
result <- Chr17_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_17_pop_cluster_sum_labels_latlong.csv')




#####################################################################################
table <- read.table('scaffold_18_pop_cluster_labels.tsv',header=T)
Chr18_combined <- left_join(pops,table)
Chr18_combined$Chr <- "Chr_18"

# Group by Pop and summarize cluster counts
result <- Chr18_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'scaffold_18_pop_cluster_sum_labels_latlong.csv')





#####################################################################################
table <- read.table('pri_ptg000045l_pop_cluster_labels.tsv',header=T)
pri45_combined <- left_join(pops,table)
pri45_combined$Chr <- "pri#ptg000045l"

# Group by Pop and summarize cluster counts
result <- pri45_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'pri_ptg000045l_pop_cluster_sum_labels_latlong.csv')




#####################################################################################
table <- read.table('pri_ptg000015l_pop_cluster_labels.tsv',header=T)
pri15_combined <- left_join(pops,table)
pri15_combined$Chr <- "pri#ptg000015l"

# Group by Pop and summarize cluster counts
result <- pri15_combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
#print(result)


#write.csv(result,'pri_ptg000015l_pop_cluster_sum_labels_latlong.csv')


