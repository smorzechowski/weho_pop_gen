# Chromosome limits


# Install packages if not already installed
# install.packages("vcfR")
# install.packages("dplyr")
# install.packages("tidyr")

library(vcfR)
library(dplyr)
library(tidyr)

# Read your gzipped VCF file
vcf <- read.vcfR("/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/17-beagle/Nleu_autos_lea_depth_filt_imputed_thin.vcf.gz")

# Convert the fixed fields (which include CHROM and POS) to a data frame
fix_df <- as.data.frame(vcf@fix)

# Convert the POS column to numeric
fix_df$POS <- as.numeric(fix_df$POS)

# For each chromosome, compute the minimum and maximum positions
dummy_limits <- fix_df %>%
  group_by(CHROM) %>%
  summarize(min_pos = min(POS), max_pos = max(POS)) %>%
  pivot_longer(cols = c(min_pos, max_pos),
               names_to = "limit",
               values_to = "x_dummy") %>%
  mutate(y_dummy = 0)  # y_dummy is arbitrary and not used in plotting


dummy_limits$x_dummy[dummy_limits$limit=="min_pos"] <- 0

# Print the dummy limits data frame
#print(dummy_limits)