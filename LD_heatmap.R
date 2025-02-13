# Make LD heatmaps from ngsLD

library(ggplot2)
library(reshape2)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/ngsLD")


data <- read.table('windows_sorted_reformatted_priptg000045l_lea_thin1kb_100kb.ld',header=T,comment.char = '')
head(data)
ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)


ggplot(ld_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="lightgreen", high="red") +
  labs(x="Position", y="Position", fill="LD (r²)") +
  theme_minimal()

###########################################################################
#data <- read.table('windows_sorted_reformatted_priptg000015l_lea_thin1kb_100kb.ld',header=T,comment.char = '')
#ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
#ld_melt <- melt(ld_matrix, na.rm=TRUE)

#ggplot(ld_melt[ld_melt$value<0.2,], aes(x=Var1, y=value)) +
#  geom_point()

###########################################################################
data <- read.table('priptg000015l_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="#3F0071") +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  theme_minimal()


###########################################################################
data <- read.table('priptg000045l_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
#  scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.3),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  theme_minimal()

###########################################################################


data <- read.table('scaffold_15_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.3),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  theme_minimal()

###########################################################################



data <- read.table('scaffold_18_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
   #scale_fill_gradient(low="white", high="#3F0071") +
   scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.3),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  theme_minimal()


