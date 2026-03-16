### Préparer les databases pour les utiliser pour l'annotation 


###########" HelixMTdb

data <- read.table("HelixMTdb_20200327.tsv", 
                   sep = "\t",       # specify tab as the separator
                   header = TRUE,    # assuming the file has a header row
                   quote = "",       # disable quoting if needed
                   stringsAsFactors = FALSE)  # keep strings as strings


# Remove the "chrM:" prefix from the locus column
data$locus <- sub("^chrM:", "", data$locus)


# Remove the brackets and quotes from the alleles column
data$alleles_clean <- gsub("\\[|\\]|\"", "", data$alleles)

# Split the cleaned string at the comma
allele_split <- strsplit(data$alleles_clean, ",")

# Create new columns for the reference (first allele) and alternate (second allele)
data$Ref <- sapply(allele_split, `[`, 1)
data$Alt <- sapply(allele_split, `[`, 2)

# Optionally, remove the temporary cleaned column
data$alleles_clean <- NULL

# Remove the 'alleles' column
data$alleles <- NULL

# Reorder columns:
data <- data[, c("locus", "Ref", "Alt", "feature", "gene", "counts_hom",
                 "AF_hom", "counts_het", "AF_het", "mean_ARF", "max_ARF",
                 "haplogroups_for_homoplasmic_variants", "haplogroups_for_heteroplasmic_variants")]

# Check the new structure
str(data)


## delete mean_ARF and max_ARF
data$mean_ARF <- NULL
data$max_ARF <- NULL

# change the name of the first column from locus to Pos
names(data)[1] <- "Pos"



write.csv(data, file = "modified_HelixMTdb.csv", row.names = FALSE)

################################## 
## MitImpact database

mitimpact_data <- read.table("MitImpact_db_3.1.3.txt", header = TRUE, sep = "\t")

mitimpact_data$Mithril_id <- NULL
mitimpact_data$MitImpact_id <- NULL
mitimpact_data$Chr <- NULL
mitimpact_data$Molecule_type <- NULL
mitimpact_data$Gene_start <- NULL
mitimpact_data$Gene_position <- NULL
mitimpact_data$Gene_end <- NULL
mitimpact_data$Gene_strand <- NULL
names (mitimpact_data) [1] <- "Pos"
mitimpact_data$Extended_annotation <- NULL


write.csv(mitimpact_data, file = "MitImpact-db-3.1.3_modified.csv", row.names = FALSE)



####################################
### tApoGEE

t_APOGEE_data <- read.table("t-APOGEE_2024.0.1.txt", header = TRUE, sep = "\t", fill = TRUE)


write.csv(t_APOGEE_data, file = "t_APOGEE_data_modified.csv", row.names = FALSE)


###############################""
### MitoMAP RNS/tRNAs

Mitomap <- read.csv ("MutationsRNA MITOMAP Foswiki.csv")


names(Mitomap)[1] <- "Pos"

# Extract Ref (first character) and Alt (last character) from Allele
Mitomap$Ref <- substring(Mitomap$Allele, 1, 1)
Mitomap$Alt <- substring(Mitomap$Allele, nchar(Mitomap$Allele), nchar(Mitomap$Allele))

# Reorder columns so that Pos, Ref, and Alt come first
Mitomap <- Mitomap[, c("Pos", "Ref", "Alt", "Locus", "Disease", "Allele", "RNA",
                       "Homoplasmy", "Heteroplasmy", "Status", "MitoTIP.",
                       "GB.Freq..FL..CR...", "GB.Seqs.FL..CR..", "References")]

# Check the new structure
str(Mitomap)

write.csv(Mitomap, file = "MitoMAP_RNA_tRNA.csv", row.names = FALSE)

##############################################################################################"
## Comparaiosn de gnomAD et mes data
# savoir si on a utilisé la meme seq de ref 

# lire gnomAD 
gnomAd_data <- read.table ("gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv", 
                           sep = "\t",       # specify tab as the separator
                           header = TRUE,    # assuming the file has a header row
                           quote = "",       # disable quoting if needed
                           stringsAsFactors = FALSE)  # keep strings as strings)

gnomAd_data$chromosome <-NULL
gnomAd_data$filters <-NULL
gnomAd_data$AN <-NULL
gnomAd_data$AC_hom <-NULL
gnomAd_data$AC_het <-NULL


# changer le nom des colonnes
colnames (gnomAd_data) <- c ("Pos", "Ref", "Alt", "AF_hom_gnomAD_V3", "AF_het_gnomAD_V3", "max_observed_heteroplasmy")

## save la data filtred 
write.csv(gnomAd_data, file = "gnomAD_filtred.csv", row.names = FALSE)



### comparer gnomAD et HelixMTdb ----------------------------------

# lire les data
Helix_data <- read.csv ("HelixMTdb_modified.csv")
gnomAd_data <- read.csv ("gnomAD_filtred.csv")



# Harmoniser les noms de colonnes pour faciliter la jointure
names(gnomAd_data)[names(gnomAd_data) == "position"] <- "Pos"
names(gnomAd_data)[names(gnomAd_data) == "ref"] <- "Ref"

# Extraire uniquement les colonnes nécessaires pour la comparaison
helix_subset <- Helix_data[, c("Pos", "Ref")]
gnomad_subset <- gnomAd_data[, c("Pos", "Ref")]

# Faire une jointure sur Pos
merged <- merge(helix_subset, gnomad_subset, by = "Pos", suffixes = c("_Helix", "_gnomAD"))

# Comparer les Ref
diff_refs <- merged[merged$Ref_Helix != merged$Ref_gnomAD, ]

# Afficher les résultats
if (nrow(diff_refs) > 0) {
  print("Positions avec des Ref différents :")
  print(diff_refs)
} else {
  print("Tous les Ref sont identiques pour les positions communes.")
}





nrow(diff_refs)


















