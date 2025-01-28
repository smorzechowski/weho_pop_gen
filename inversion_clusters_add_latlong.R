# Add lat long to inversion cluster files for ArcGIS Pro

library(dplyr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/local_pcangsd/tables")
pops <- read.csv('~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/all_samples_pops_lat_long.csv',header=T)


#####################################################################################
table <- read.table('scaffold_4_pop_cluster1_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_4_pop_cluster1_updated_sum_labels_latlong.csv')




#####################################################################################
table <- read.table('scaffold_4_pop_cluster2_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_4_pop_cluster2_updated_sum_labels_latlong.csv')









#####################################################################################
table <- read.table('scaffold_6_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_6_pop_cluster_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_8_pop_cluster1_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_8_pop_cluster1_updated_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_8_pop_cluster2_labels.tsv',header=T) # NEED to recheck the clustering
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_8_pop_cluster2_updated_sum_labels_latlong.csv')






#####################################################################################
table <- read.table('scaffold_9_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_9_pop_cluster_sum_labels_latlong.csv')

#####################################################################################
table <- read.table('scaffold_10_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_10_pop_cluster_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_13_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_13_pop_cluster_sum_labels_latlong.csv')



#####################################################################################
table <- read.table('scaffold_14_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_14_pop_cluster_sum_labels_latlong.csv')


#####################################################################################
table <- read.table('scaffold_15pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_15_pop_cluster_sum_labels_latlong.csv')



#####################################################################################
table <- read.table('scaffold_18_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)


# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'scaffold_18_pop_cluster_sum_labels_latlong.csv')





#####################################################################################
table <- read.table('pri_ptg000045l_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)


# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'pri_ptg000045l_pop_cluster_sum_labels_latlong.csv')




#####################################################################################
table <- read.table('pri_ptg000015l_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)


# Group by Pop and summarize cluster counts
result <- combined %>%
  group_by(Pop, Latitude, Longitude) %>%
  summarise(
    Cluster_0 = sum(Cluster == 0),
    Cluster_1 = sum(Cluster == 1),
    Cluster_2 = sum(Cluster == 2),
    .groups = "drop"
  )

# View the result
print(result)


write.csv(result,'pri_ptg000015l_pop_cluster_sum_labels_latlong.csv')


