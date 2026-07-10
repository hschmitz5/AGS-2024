agglom_genus <- function(ps) {
 
  # If Genus is NA, then replace with lowest classified taxon 
  taxonomy <- data.frame(ps@tax_table) %>%
    rownames_to_column("OTU") %>%
    mutate(
      Genus = case_when(
        !is.na(Genus)  ~ Genus,
        !is.na(Family) ~ ifelse(startsWith(Family, "midas"), paste0("Unk_", Family), paste0("Unk_f_", Family)),
        !is.na(Order)  ~ ifelse(startsWith(Order, "midas"), paste0("Unk_", Order), paste0("Unk_o_", Order)),
        !is.na(Class)  ~ ifelse(startsWith(Class, "midas"), paste0("Unk_", Class), paste0("Unk_c_", Class)),
        !is.na(Phylum) ~ ifelse(startsWith(Phylum, "midas"), paste0("Unk_", Phylum), paste0("Unk_p_", Phylum)),
        is.na(Phylum)  ~ paste0("Unknown_Phylum"),
        .default = Genus
      )
    ) %>%
    column_to_rownames("OTU") %>%
    as.matrix()
  
  # Add new taxonomy to phyloseq object
  tax_table(ps) <- tax_table(taxonomy)
  
  # uses updated taxonomy to keep NA values
  ps_genus = tax_glom(ps, "Genus")
}

agglom_species <- function(ps) {
  
  # If Species is NA, then replace with lowest classified taxon 
  taxonomy <- data.frame(ps@tax_table) %>%
    rownames_to_column("OTU") %>%
    mutate(
      Species = case_when(
        !is.na(Species) ~ Species,
        !is.na(Genus)   ~ ifelse(startsWith(Genus, "midas"), paste0("Unk_", Genus), paste0("Unk_g_", Genus)),
        !is.na(Family)  ~ ifelse(startsWith(Family, "midas"), paste0("Unk_", Family), paste0("Unk_f_", Family)),
        !is.na(Order)   ~ ifelse(startsWith(Order, "midas"), paste0("Unk_", Order), paste0("Unk_o_", Order)),
        !is.na(Class)   ~ ifelse(startsWith(Class, "midas"), paste0("Unk_", Class), paste0("Unk_c_", Class)),
        !is.na(Phylum)  ~ ifelse(startsWith(Phylum, "midas"), paste0("Unk_", Phylum), paste0("Unk_p_", Phylum)),
        is.na(Phylum)   ~ paste0("Unknown_Phylum"),
        .default = Species
      )
    ) %>%
    column_to_rownames("OTU") %>%
    as.matrix()
  
  # Add new taxonomy to phyloseq object
  tax_table(ps) <- tax_table(taxonomy)
  
  # uses updated taxonomy to keep NA values
  ps_species = tax_glom(ps, "Species")
}