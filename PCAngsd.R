#Load the RcppCNPy package
library(RcppCNPy)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/pcangsd")
samples<-read.table('~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/ngsRelate/allbamspath.txt',sep="/")
samples$samples <- str_remove_all(samples$V12,'_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam')

#### Autosomal structure, including inversions
C <- as.matrix(read.table("Nleu_autos_ANGSD_PCApcangsd.cov")) # Reads estimated covariance matrix


Apop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
          "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mall","Mall","Mall","Mall","Mall","Mall","Mall","Mall",
          "Moon","Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom","Nom",
          "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
          "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")


mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Apop, data.frame(eigenvectors))) #combine with our population assignments


pca <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Apop)) +
  geom_point(size = 3) +
  geom_text(aes(label = samples$samples), vjust = -1, size = 3) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

print(pca)


#### New PAR structure
C <- as.matrix(read.table("Nleu_newPAR_pcangsd_unpruned.cov")) # Reads estimated covariance matrix


Apop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
          "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mall","Mall","Mall","Mall","Mall","Mall","Mall","Mall",
          "Moon","Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom","Nom",
          "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
          "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")



mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Apop, data.frame(eigenvectors))) #combine with our population assignments

pca <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Apop)) +
  geom_point(size = 3) +
  geom_text(aes(label = samples$samples), vjust = -1, size = 3) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

print(pca)

#### Old code

#C <- as.matrix(read.table("Nleu_neoW_pcangsd_unpruned.cov")) # Reads estimated covariance matrix

#C <- as.matrix(read.table("Nleu_neoZ_pcangsd_unpruned.cov")) # Reads estimated covariance matrix
C <- as.matrix(read.table("Nleu_neoZ_males_pcangsd_unpruned.cov")) # Reads estimated covariance matrix
#D <- as.matrix(read.table("pcangsd.selection")) # Reads PC based selection statistics

# Plot PCA plot
e <- eigen(C)
plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd")

# Obtain p-values from PC-based selection scan
#p <- pchisq(D, 1, lower.tail=FALSE)


#Load the covariance matrix

#We will also add a column with population assignments
Wpop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Bog","Ingalba","Ingalba","Ingalba","Mallee","Mallee","Moonbi","Moonbi"
         ,"Mull","Nom","Nom","Pill","Reedy","Talla","Walch","Walch","Walch","Walch","Zost")


Mpop <- c("Nom","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Bog","Weddin","Gund","Gund","Gund","Mallee","Mallee",
          "Mallee","Mallee","Mallee","Mallee","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom",
          "Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Walch","Walch","Walch","Walch","Zost","Zost")





Msam <- c("Nom-93","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Bog","Weddin","Gund","Gund","Gund","Mallee","Mallee",
          "Mallee","Mallee","Mallee","Mallee","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Nom-94","Nom-96","Nom-100","Nom-103","Nom-105",
          "Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Walch","Walch","Walch","Walch","Zost","Zost")

mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Apop, data.frame(eigenvectors))) #combine with our population assignments
pca.vectors = as_tibble(cbind(Wpop, data.frame(eigenvectors))) #combine with our population assignments
pca.vectors = as_tibble(cbind(Mpop, data.frame(eigenvectors))) #combine with our population assignments


pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = Mpop)) + geom_point()
pca



pca <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Apop)) +
  geom_point(size = 3) +
  geom_text(aes(label = Asam), vjust = -1, size = 3) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

print(pca)

pca <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Mpop)) +
  geom_point(size = 3) +
  geom_text(aes(label = Msam), vjust = -1, size = 3) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

print(pca)


#ggsave(filename = paste0(basedir, "/results/pca_LDpruned_pcangsd_plot.pdf"), plot = pca) #change file path if data on your own computer

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

varPC1
varPC2
varPC3
varPC4