# Contingency test of inversion genes involved in climate adaptation
# Sophie MacRae Orzechowski
# February 2025

# Bird candidate genes

# Fill in your total gene count here:
#G <- 23414   # Total number of genes in the Nleucotis gff file from liftoff
G <- 23414 - (3057-2098) # excludes ambiguous MSTRG and reg_ in the inversion region

# Known numbers for overlap between inversion genes and climate genes
overlap           <- 96
inversion_total   <- 2098 # excludes ambiguous MSTRG and reg_
#inversion_total   <- 3057 
climate_total     <- 628 # excludes ambiguous LOC

# Derive counts:
inversion_not_climate <- inversion_total - overlap
climate_not_inversion <- climate_total - overlap
not_climate_not_inversion <- G - (overlap + inversion_not_climate + climate_not_inversion)

# Build the table:
contingency_table <- matrix(
  c(overlap, inversion_not_climate,
    climate_not_inversion, not_climate_not_inversion),
  nrow = 2, byrow = TRUE
)
# OPTIONAL: Name the rows/columns for clarity
dimnames(contingency_table) <- list(
  Inversion  = c("Inversion", "Not_Inversion"),
  Climate    = c("Climate", "Not_Climate")
)

# Perform Fisher's exact test:
fisher.test(contingency_table)



#################################################################
# Vertebrate candidate genes

# Fill in your total gene count here:
#G <- 23414   # Total number of genes in the Nleucotis gff file from liftoff
G <- 23414 - (3057-2098) # excludes ambiguous MSTRG and reg_ in the inversion region

# Known numbers for overlap between inversion genes and climate genes
overlap           <- 165
inversion_total   <- 2098 # excludes ambiguous MSTRG and reg_
#inversion_total   <- 3057 
climate_total     <- 1101 # excludes ambiguous LOC

# Derive counts:
inversion_not_climate <- inversion_total - overlap
climate_not_inversion <- climate_total - overlap
not_climate_not_inversion <- G - (overlap + inversion_not_climate + climate_not_inversion)

# Build the table:
contingency_table <- matrix(
  c(overlap, inversion_not_climate,
    climate_not_inversion, not_climate_not_inversion),
  nrow = 2, byrow = TRUE
)
# OPTIONAL: Name the rows/columns for clarity
dimnames(contingency_table) <- list(
  Inversion  = c("Inversion", "Not_Inversion"),
  Climate    = c("Climate", "Not_Climate")
)

# Perform Fisher's exact test:
fisher.test(contingency_table)


#################################################################
# Vertebrate candidate genes + Avian candidate genes

# Fill in your total gene count here:
#G <- 23414   # Total number of genes in the Nleucotis gff file from liftoff
G <- 23414 - (3057-2098) # excludes ambiguous MSTRG and reg_ in the inversion region

# Known numbers for overlap between inversion genes and climate genes
overlap           <- 250
inversion_total   <- 2098 # excludes ambiguous MSTRG and reg_
#inversion_total   <- 3057 
climate_total     <- 1648 # excludes ambiguous LOC

# Derive counts:
inversion_not_climate <- inversion_total - overlap
climate_not_inversion <- climate_total - overlap
not_climate_not_inversion <- G - (overlap + inversion_not_climate + climate_not_inversion)

# Build the table:
contingency_table <- matrix(
  c(overlap, inversion_not_climate,
    climate_not_inversion, not_climate_not_inversion),
  nrow = 2, byrow = TRUE
)
# OPTIONAL: Name the rows/columns for clarity
dimnames(contingency_table) <- list(
  Inversion  = c("Inversion", "Not_Inversion"),
  Climate    = c("Climate", "Not_Climate")
)

# Perform Fisher's exact test:
fisher.test(contingency_table)











# Hypergeometric test
# 146 IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#phyper(146 - 1, 632, G - 632, 3057, lower.tail = FALSE)


#phyper(69 - 1, 632, G - 632, 1040, lower.tail = FALSE)


######################################################################################
# I could run a contingency test to see if hits are significantly more likely to be in inversions or not
# Of course, this is flawed because of the issue of linkage within inversions. But still, it does show that 
# significant GEA hits are more likely to be within inversions rather than outside of it.
# Something is always flawed, it seems like. 

G <- 23414 - (3057-2098) - (1040 - 757) # excludes ambiguous MSTRG and reg_ 

# Known numbers for overlap between GEA results and climate genes
overlap           <- 695
GEA_total   <- 757
inversion_total     <- 2098


# Derive counts:
GEA_not_inversion <- GEA_total - overlap
inversion_not_GEA <- inversion_total - overlap
not_GEA_not_inversion <- G - (overlap + inversion_not_GEA + GEA_not_inversion)


# Build the table:
contingency_table <- matrix(
  c(overlap, inversion_not_GEA,
    GEA_not_inversion, not_GEA_not_inversion),
  nrow = 2, byrow = TRUE
)
# OPTIONAL: Name the rows/columns for clarity
dimnames(contingency_table) <- list(
  Inversion  = c("Inversion", "Not_Inversion"),
  GEA    = c("GEA", "Not_GEA")
)

# Perform Fisher's exact test:
fisher.test(contingency_table)


######################################################################################
# I could run a contingency test to see if climate-related GEA hits are more likely to be within inversions
# Perhaps this is more comparable?? I dunno. I guess it's still complicated by linkage

G <- 632 # excludes LOC

# Known numbers for overlap between GEA results and climate genes
overlap           <- 70
GEA_total   <- 78
inversion_total     <- 96


# Derive counts:
GEA_not_inversion <- GEA_total - overlap
inversion_not_GEA <- inversion_total - overlap
not_GEA_not_inversion <- G - (overlap + inversion_not_GEA + GEA_not_inversion)


# Build the table:
contingency_table <- matrix(
  c(overlap, inversion_not_GEA,
    GEA_not_inversion, not_GEA_not_inversion),
  nrow = 2, byrow = TRUE
)
# OPTIONAL: Name the rows/columns for clarity
dimnames(contingency_table) <- list(
  Inversion  = c("Inversion", "Not_Inversion"),
  GEA    = c("GEA", "Not_GEA")
)

# Perform Fisher's exact test:
fisher.test(contingency_table)



