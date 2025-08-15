#### Libraries ####
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(tidyr)
library(decontam)
library(microbiome)
library(vegan)
library(tidyverse)
library(dplyr)
library(microViz)
library(ggpubr)
library(cowplot)
library(gridGraphics)
library(ggh4x)
library(stringr)
library(pairwiseAdonis)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(Biostrings)
library(reshape2)
library(ggnewscale)

#### Importing Data ####
asvs = read_qza("deepsea-18s-merged-asv-table-05022025.qza")
metadata = read.csv("combined_metadata.csv")
taxonomy = read_qza("18S-rep-sequences-taxonomy.qza")
tree <- read_qza("rooted-18S-tree.qza")$data

# Extracting data
asv_df = asvs$data
taxonomy_df = taxonomy$data

#### Formatting Data ####
taxonomy_fixed_df <- taxonomy_df %>% separate_wider_delim(Taxon, delim = ";", names_sep = "", too_few = "align_start")
taxonomy_fixed_df[is.na(taxonomy_fixed_df)] <- "Unassigned" # rename NAs into unassigned

#### Make Feature ID rownames ####
taxonomy_fixed_df = as.data.frame(taxonomy_fixed_df)
row.names(taxonomy_fixed_df) = taxonomy_fixed_df$Feature.ID
taxonomy_fixed_df$Feature.ID = NULL

#### Remove prefix from taxonomy assignments and replace NAs with "unassigned" ####
taxonomy_fixed_df[] = lapply(taxonomy_fixed_df, function(x) gsub("D_\\d+__", "", x))
taxonomy_fixed_df[taxonomy_fixed_df==""] = NA
taxonomy_fixed_df = taxonomy_fixed_df %>% replace(is.na(.), "Unassigned")

taxonomy_matrix = as.matrix(taxonomy_fixed_df)

#### Remove typo entry in asv table ####
asv_df <- asv_df[, !(colnames(asv_df) %in% "GordaRdige_Vent090_T24_2019")]

#### Merge Into Phyloseq ####
physeq_asv = otu_table(asv_df, taxa_are_rows = T)
physeq_tax = tax_table(taxonomy_matrix)
physeq_meta = sample_data(metadata)

#### Fix Row Names in Metadata ####
row.names(physeq_meta) = physeq_meta$Sample_Name
physeq_meta$Sample_Name = NULL

phylo_object <- phyloseq(physeq_asv, physeq_tax, physeq_meta)

#### Decontamination ####
sample_data(phylo_object)$is.neg <- sample_data(phylo_object)$Sample_or_Control == "Control" # create a sample-variable for contaminants
phylo_object_contaminants = isContaminant(phylo_object, method = "prevalence", neg = "is.neg", threshold = 0.5, detailed = TRUE, normalize = TRUE) # Detect contaminants based on control samples and ASV prevalence
table(phylo_object_contaminants$contaminant)

#### Make phyloseq object of presence-absence in negative controls and true samples ####
phylo_object_contaminants.pa = transform_sample_counts(phylo_object, function(abund) 1 * (abund > 0)) # Convert phyloseq to presence-absence
ps.pa.neg = subset_samples(phylo_object_contaminants.pa, Sample_or_Control == "Control")
ps.pa.pos = subset_samples(phylo_object_contaminants.pa, Sample_or_Control == "Sample")
df.pa = data.frame(pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg), contaminant = phylo_object_contaminants$contaminant)

#### Make phyloseq object of presence-absence in negative controls and true samples ####
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#### Remove contaminants from dataset ####
phylo_obj_sans_contam = prune_taxa(!phylo_object_contaminants$contaminant, phylo_object)
phylo_obj_sans_contam_sans_controls = subset_samples(phylo_obj_sans_contam, Sample_or_Control != "Control")

#### Remove troublesome samples ####
to_remove = "Axial_Dependable_FS900_2013"
phylo_obj_sans_contam_sans_controls = prune_samples(!(sample_names(phylo_obj_sans_contam_sans_controls) %in% to_remove), phylo_obj_sans_contam_sans_controls)

to_remove = "GordaRidge_BSW020_sterivex_2019_REPa"
phylo_obj_sans_contam_sans_controls = prune_samples(!(sample_names(phylo_obj_sans_contam_sans_controls) %in% to_remove), phylo_obj_sans_contam_sans_controls)

#### Subset samples to remove control samples and Unassgined taxa at domain level ####
phylo_obj_sans_contam_sans_controls = subset_samples(phylo_obj_sans_contam_sans_controls, Sample_or_Control == "Sample")

phylo_obj_sans_contam_sans_controls = phylo_obj_sans_contam_sans_controls %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)

phylo_obj_sans_contam_sans_controls_sans_empty = subset_taxa(phylo_obj_sans_contam_sans_controls, Taxon1 != "Unassigned")

#### Normalize phyloseq object ####
#### Use tax_fix to fill in unassigned taxa for names at higher ranks ####
phylo_notnorm = phylo_obj_sans_contam_sans_controls_sans_empty

phylo_notnorm_taxfix = tax_fix(phylo_notnorm, unknowns = "Unassigned")

phylo_normalized = microbiome::transform(phylo_obj_sans_contam_sans_controls_sans_empty, "compositional")

phylo_norm_taxfix = tax_fix(phylo_normalized, unknowns = "Unassigned")

phylo_not_norm_df = ps_melt(phylo_notnorm)

#### Subset by Region ####
phylo_norm_oldaxial = subset_samples(phylo_normalized, Region == "Old Axial")
phylo_norm_newaxial = subset_samples(phylo_normalized, Region == "New Axial")
phylo_norm_gorda = subset_samples(phylo_normalized, Region == "Gorda Ridge")
phylo_norm_mcr = subset_samples(phylo_normalized, Region == "MCR")
phylo_norm_siders = subset_samples(phylo_normalized, Region == "Siders")

phylo_notnorm_oldaxial = subset_samples(phylo_notnorm, Region == "Old Axial")
phylo_notnorm_newaxial = subset_samples(phylo_notnorm, Region == "New Axial")
phylo_notnorm_gorda = subset_samples(phylo_notnorm, Region == "Gorda Ridge")
phylo_notnorm_mcr = subset_samples(phylo_notnorm, Region == "MCR")

#### Tax bar plot at rank  7 (Metazoa rank) ####
ggplot(phylo_zoa, aes(x=Sample, y=Abundance, fill = Taxon7)) +
  ggtitle("Zoa") +
  geom_bar(stat = "identity") +
  facet_wrap(~ Region, scales = "free_x") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

plot_bar(phylo_zoa, fill = "Taxon7") +
  facet_wrap(~ Region, scales = "free_x") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#### Remove Possible Outlier ####
to_remove = "GordaRidge_BSW020_sterivex_2019_REPa"
phylo_no_outlier = prune_samples(!(sample_names(phylo_notnorm) %in% to_remove), phylo_notnorm)

#### Ordination by region and sample type ####
phylo_norm_ord = ordinate(phylo_normalized, "PCoA", "bray")

plot_ordination(phylo_normalized, phylo_norm_ord, type = "samples", color = "Region", shape = "sample_type") +
  geom_point(size = 3) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  annotate("text",
           x = Inf, y = Inf,
           label = "p = 0.001",
           hjust = 1, vjust = 1,
           size = 4)

#### PERMANOVA (Adonis test) ####
phylo_norm_dist = distance(phylo_normalized, method = "bray")
adonis_data = as(sample_data(phylo_normalized), "data.frame")

#Adonis and disper test by region
phylo_norm_adonis_region = adonis2(phylo_norm_dist ~ Region, data = adonis_data)
write.csv(phylo_norm_adonis_region, "phylo_norm_adonis_region.csv")
res <- betadisper(phylo_norm_dist, adonis_data$Region)
region_disper = permutest(betadisper(phylo_norm_dist, adonis_data$Region))
region_disper_df = as.data.frame(region_disper$tab)
write.csv(region_disper_df, "region_disper_df.csv")

#Adonis and disper test by sample type
phylo_norm_adonis_sampletype = adonis2(phylo_norm_dist ~ sample_type, data = adonis_data)
write.csv(phylo_norm_adonis_sampletype, "phylo_norm_adonis_sampletype.csv")
res2 <- betadisper(phylo_norm_dist, adonis_data$sample_type)
sampletype_disper = permutest(betadisper(phylo_norm_dist, adonis_data$sample_type))
sampletype_disper_df = as.data.frame(sampletype_disper$tab)
write.csv(sampletype_disper_df, "sampletype_disper_df.csv")

#Adonis and disper tests by sample type within each region
phylo_oldaxial_dist = distance(phylo_norm_oldaxial, method = "bray")
adonis_data_oldaxial = as(sample_data(phylo_norm_oldaxial), "data.frame")
phylo_oldaxial_adonis = adonis2(phylo_oldaxial_dist ~ sample_type, data = adonis_data_oldaxial)
write.csv(phylo_oldaxial_adonis, "phylo_oldaxial_adonis.csv")
oldaxial_disper = permutest(betadisper(phylo_oldaxial_dist, adonis_data_oldaxial$sample_type))
oldaxial_disper_df = as.data.frame(oldaxial_disper$tab)
write.csv(oldaxial_disper_df, "oldaxial_disper_df.csv")

phylo_newaxial_dist = distance(phylo_norm_newaxial, method = "bray")
adonis_data_newaxial = as(sample_data(phylo_norm_newaxial), "data.frame")
phylo_newaxial_adonis = adonis2(phylo_newaxial_dist ~ sample_type, data = adonis_data_newaxial)
write.csv(phylo_newaxial_adonis, "phylo_newaxial_adonis.csv")
newaxial_disper = permutest(betadisper(phylo_newaxial_dist, adonis_data_newaxial$sample_type))
newaxial_disper_df = as.data.frame(newaxial_disper$tab)
write.csv(newaxial_disper_df, "newaxial_disper_df.csv")

phylo_gorda_dist = distance(phylo_norm_gorda, method = "bray")
adonis_data_gorda = as(sample_data(phylo_norm_gorda), "data.frame")
phylo_gorda_adonis = adonis2(phylo_gorda_dist ~ sample_type, data = adonis_data_gorda)
write.csv(phylo_gorda_adonis, "phylo_gorda_adonis.csv")
gorda_disper = permutest(betadisper(phylo_gorda_dist, adonis_data_gorda$sample_type))
gorda_disper_df = as.data.frame(gorda_disper$tab)
write.csv(gorda_disper_df, "gorda_disper_df.csv")

phylo_mcr_dist = distance(phylo_norm_mcr, method = "bray")
adonis_data_mcr = as(sample_data(phylo_norm_mcr), "data.frame")
phylo_mcr_adonis = adonis2(phylo_mcr_dist ~ sample_type, data = adonis_data_mcr)
write.csv(phylo_mcr_adonis, "phylo_mcr_adonis.csv")
mcr_disper = permutest(betadisper(phylo_mcr_dist, adonis_data_mcr$sample_type))
mcr_disper_df = as.data.frame(mcr_disper$tab)
write.csv(mcr_disper_df, "mcr_disper_df.csv")

phylo_siders_dist = distance(phylo_norm_siders, method = "bray")
adonis_data_siders = as(sample_data(phylo_norm_siders), "data.frame")
phylo_siders_adonis = adonis2(phylo_siders_dist ~ sample_type, data = adonis_data_siders)
write.csv(phylo_siders_adonis, "phylo_siders_adonis.csv")
siders_disper = permutest(betadisper(phylo_siders_dist, adonis_data_siders$sample_type))
siders_disper_df = as.data.frame(siders_disper$tab)
write.csv(sideres_disper_df, "siders_disper_df.csv")

#### Pairwise Adonis ####
grouping_var <- adonis_data$Region

sample_names <- labels(phylo_norm_dist)
meta <- data.frame(SampleID = sample_names,
                   Group = grouping_var,
                   stringsAsFactors = FALSE)

group_levels <- unique(meta$Group)
pairwise_combos <- combn(group_levels, 2, simplify = FALSE)

results <- list()

# Loop through each pairwise combination
for (pair in pairwise_combos) {
  g1 <- pair[1]
  g2 <- pair[2]
  
  # Subset metadata and distance matrix for just those two groups
  keep_samples <- meta$SampleID[meta$Group %in% c(g1, g2)]
  sub_meta <- meta %>% filter(SampleID %in% keep_samples)
  sub_dist <- as.dist(as.matrix(phylo_norm_dist)[keep_samples, keep_samples])
  
  # Run adonis2
  ad <- adonis2(sub_dist ~ Group, data = sub_meta)
  
  # Save result
  results[[paste(g1, "vs", g2)]] <- data.frame(
    Comparison = paste(g1, "vs", g2),
    F.Model = ad$F[1],
    R2 = ad$R2[1],
    p.value = ad$`Pr(>F)`[1]
  )
}

pairwise_results <- bind_rows(results)

# Adjust p-values for multiple testing
pairwise_results$p.adjusted <- p.adjust(pairwise_results$p.value, method = "fdr")

# View results
print(pairwise_results)

# Write to CSV
write.csv(pairwise_results, "pairwise_adonis_results.csv", row.names = FALSE)

#### Subsetting by taxa ####
phylo_norm_animal = subset_taxa(phylo_normalized, Taxon8 == "D_7__Animalia")

phylo_norm_nematode = subset_taxa(phylo_normalized, Taxon14 == "D_13__Nematoda")

#### Alpha Diversity ####
alpha_div_obs = phyloseq::estimate_richness(phylo_notnorm, measures = "Observed")

alpha_div_observed_metadata = data.frame(sample_data(phylo_notnorm), 
                                         "Reads" = sample_sums(phylo_notnorm), 
                                         "Observed" = estimate_richness(phylo_notnorm, measures = "Observed"), 
                                         "Shannon" = estimate_richness(phylo_notnorm, measures = "Shannon"), 
                                         "InvSimpson" = estimate_richness(phylo_notnorm, measures = "InvSimpson"))

#### Evenness & Alpha Diversity Plot ####
alpha_div_observed_metadata$Evenness = alpha_div_observed_metadata$Shannon / log(alpha_div_observed_metadata$Observed)

ggplot(alpha_div_observed_metadata, aes(x=Region, y=Reads)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test") +
  geom_boxplot() +
  theme_minimal() +
  labs(y = "Reads", x = "Region", title = "Reads by Region") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

#### Multi-index graph ####
alpha_long = alpha_div_observed_metadata %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson, Evenness), names_to = "Metric", values_to = "Value")

metric_labels = c(Observed = "Observed ASVs", Shannon = "Shannon Diversity", InvSimpson = "Inverse Simpson", Evenness = "Evenness")

comparisons <- list(
  c("Old Axial", "New Axial"),
  c("Old Axial", "Gorda Ridge"),
  c("Old Axial", "MCR"),
  c("New Axial", "Gorda Ridge"),
  c("New Axial", "MCR"),
  c("Gorda Ridge", "MCR"),
  c("Siders", "Old Axial"),
  c("Siders", "New Axial"),
  c("Siders", "Gorda Ridge"),
  c("Siders", "MCR")
)

alpha_long$Region <- factor(
  alpha_long$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

ggplot(alpha_long, aes(x = Region, y = Value, fill = Region)) +  # Add 'fill = Region'
  geom_boxplot() +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test") +
  facet_wrap(~ Metric, scales = "free_y", labeller = as_labeller(metric_labels)) +
  scale_fill_manual(values = c(
    "Old Axial" = "red",
    "New Axial" = "orange",
    "Gorda Ridge" = "green",
    "MCR" = "blue",
    "Siders" = "violet"
  )) +
  theme_minimal() +
  labs(y = "Diversity Metric", x = "Sample Type", title = "Alpha Diversity w/ Stat Tests") +
  theme(panel.border = element_rect(color = "black", fill = NA), 
        strip.text = element_text(size = 14, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + guides(fill = "none")

ggplot(data = alpha_long, aes(x = sample_type, y = Value)) + 
  geom_boxplot() + 
  facet_nested(Metric ~ Region, scales = "free_y", drop = TRUE) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.5))
  
 
#### Nematode Genera Summary ####
nematode_raw = subset_taxa(phylo_notnorm, Taxon14 == "Nematoda")
nematode_raw_df = psmelt(nematode_raw)

genus_summary_raw = nematode_raw_df %>%
  group_by(Taxon22) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

genus_summary_raw$Taxon22[14] = "Unclassified"

genus_summary_raw$Taxon22 = gsub("^D_21__", "", genus_summary_raw$Taxon22)

write.csv(genus_summary_raw, "NematodeSummary.csv", row.names = FALSE)

genus_summary_raw_condensed$Taxon22 <- factor(genus_summary_raw_condensed$Taxon22, levels = genus_summary_raw_condensed$Taxon22[order(-genus_summary_raw_condensed$total_reads)])

ggplot(genus_summary_raw_condensed, aes(x = Taxon22, y = total_reads)) +
  geom_bar(stat = "identity") +
  xlab("Nematode Genus") +
  ylab("Total Reads") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

genus_summary_raw_condensed = genus_summary_raw %>%
  mutate(Taxon22 = if_else(row_number() > 9, "Others", Taxon22)) %>%
  group_by(Taxon22) %>%
  summarise(total_reads = sum(total_reads)) %>%
  ungroup()

#### Top 10 Taxa Function ####
merge_top10_phylo <- function(phylo_region_metazoa, top=9){
  transformed <- transform_sample_counts(phylo_region_metazoa, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

# Run function on phyloseq object and make into data frame
phy_18_ludox_top10_V14 <- merge_top20_18s_ludox(glom_18s_ludox, top=9)
phy_18_ludox_top10_V14_df <- psmelt(phy_18_ludox_top10_V14)
# Add common factors to use for plotting
phy_18_ludox_top10_agr = aggregate(Abundance~Sample+Site+Habitat+V14, data=phy_18_ludox_top10_V14_df, FUN=mean) 
unique(phy_18_ludox_top10_agr$V14)
# Put "Others" to the final of the Phylum list - top 10
phy_18_ludox_top10_agr$V14 <- factor(phy_18_ludox_top10_agr$V14,
                                     levels = c("Crustacea", "Gonyaulax spinifera V9","Lecudina phyllochaetopteri V9",
                                                "Maullinia V6","Polychaeta", "Protoperidinium conicum V9", "Protoperidinium leonis V9", "Nematoda",
                                                "Rhabditophora","Others"))
# Reorder Site levels
phy_18_ludox_top10_agr$Site = factor(phy_18_ludox_top10_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))
# Plot by site - Class level top 10
taxonomy_bar_18s_ludox_top10_phylum <- ggplot(phy_18_ludox_top10_agr, aes(x = Sample, y = Abundance, fill = V14)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top10_2, name = "Ludox Phylum") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12)) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) 
taxonomy_bar_18s_ludox_top10_phylum

#### Merge Attempt####
phylo_phylum = phylo_notnorm %>%
  tax_fix(unknowns = "Unassigned") %>%
  tax_glom(taxrank = "Taxon14")

phylo_merged = merge_top10_phylo(phylo_phylum, top = 29)
phylo_merged_df = psmelt(phylo_merged)

phylo_merged_df_agr = aggregate(Abundance~Sample+Region+Taxon14, data=phylo_merged_df, FUN=mean) 

unique(phylo_merged_df$Taxon14)

top_phyla <- sort(taxa_sums(phy_phylum), decreasing = TRUE)
head(top_phyla, 20)

#### Metazoa ####
phylo_zoa = phylo_notnorm %>%
  tax_fix(unknowns = "Unassigned") %>%
  tax_glom(taxrank = "Taxon7") %>%
  merge_top10_phylo(top = 9) %>%
  psmelt()

phylo_zoa_taxfix = tax_fix(phylo_notnorm)
phylo_zoa_taxglom = tax_glom(phylo_zoa_taxfix, taxrank = "Taxon7")

#### All metazoa Taxbar ####
merge_top10_phylo_everything_zoa <- function(phylo_notnorm_taxfix, top=9){
  transformed <- transform_sample_counts(phylo_notnorm_taxfix, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_notnorm_taxfix_glom = tax_glom(phylo_notnorm_taxfix, taxrank = "Taxon7")
phylo_notnorm_taxfix_glom_merged = merge_top10_phylo_everything_zoa(phylo_notnorm_taxfix_glom)
phylo_notnorm_taxfix_glom_merged_df = psmelt(phylo_notnorm_taxfix_glom_merged)

#Plot
ggplot(phylo_notnorm_taxfix_glom_merged_df, aes(x=Sample, y=Abundance, fill = Taxon7)) +
  ggtitle("Zoa Taxbar by Region") +
  geom_bar(stat = "identity") +
  facet_wrap(~ Region, scales = "free_x") + theme(axis.text.x = element_blank(), 
                                                  axis.ticks.x = element_blank(), 
                                                  plot.title = element_text(size = 10, face = "bold", hjust = 0.5))


#### Old Axial Taxbar ####
colors_top10_merged <- c("#9A6324", "#46F0F0", "#1F78B4", "#AAFFC3", "#E6194B", 
                  "#33A02C", "#4363D8", "#008080", "#FB9A99", "grey")

merge_top10_phylo_oldaxial <- function(phylo_norm_oldaxial_metazoa, top=9){
  transformed <- transform_sample_counts(phylo_norm_oldaxial_metazoa, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_norm_oldaxial = subset_samples(phylo_notnorm_taxfix, Region == "Old Axial")
phylo_norm_oldaxial_metazoa = subset_taxa(phylo_norm_oldaxial, Taxon7 == "Metazoa")

phylo_norm_oldaxial_metazoa_glom = tax_glom(phylo_norm_oldaxial_metazoa, taxrank = "Taxon14")
phylo_norm_oldaxial_metazoa_glom_merged = merge_top10_phylo_oldaxial(phylo_norm_oldaxial_metazoa_glom)
phylo_norm_oldaxial_metazoa_glom_merged_df = psmelt(phylo_norm_oldaxial_metazoa_glom_merged)

phylo_norm_oldaxial_metazoa_glom_merged_df$Taxon14 <- phylo_norm_oldaxial_metazoa_glom_merged_df$Taxon14 %>%
  replace_na('Others')

phylo_norm_oldaxial_metazoa_glom_merged_df_agr = aggregate(Abundance~sample_type+location_name+Taxon14, data=phylo_norm_oldaxial_metazoa_glom_merged_df, FUN=mean) 

phylo_norm_oldaxial_metazoa_glom_merged_df_agr$Taxon14 <- factor(
  phylo_norm_oldaxial_metazoa_glom_merged_df_agr$Taxon14,
  levels = c("Beroe forskalii Taxon13", "Crustacea", "Cunina frugifera", "Cydippida Taxon12", "Gastropoda",
             "Ophiuroidea", "Phacellophora camtschatica (eggyolk jelly) Taxon13", "Polychaeta", "Siphonophorae Taxon13", "Others")
)

#Plot
ggplot(phylo_norm_oldaxial_metazoa_glom_merged_df_agr, aes(x =location_name, y = Abundance, fill = Taxon14)) +
  ggtitle("Old Axial Taxbar by Location") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top10_merged) +
  facet_nested(. ~ sample_type + location_name, scales = "free") +
  labs(x = "Location", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### New Axial Taxbar ####
merge_top10_phylo_newaxial <- function(phylo_norm_newaxial_metazoa, top=9){
  transformed <- transform_sample_counts(phylo_norm_newaxial_metazoa, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_norm_newaxial = subset_samples(phylo_notnorm_taxfix, Region == "New Axial")
phylo_norm_newaxial_metazoa = subset_taxa(phylo_norm_newaxial, Taxon7 == "Metazoa")

phylo_norm_newaxial_metazoa_glom = tax_glom(phylo_norm_newaxial_metazoa, taxrank = "Taxon14")
phylo_norm_newaxial_metazoa_glom_merged = merge_top10_phylo_newaxial(phylo_norm_newaxial_metazoa_glom)
phylo_norm_newaxial_metazoa_glom_merged_df = psmelt(phylo_norm_newaxial_metazoa_glom_merged)

phylo_norm_newaxial_metazoa_glom_merged_df$Taxon14 <- phylo_norm_newaxial_metazoa_glom_merged_df$Taxon14 %>%
  replace_na('Others')

phylo_norm_newaxial_metazoa_glom_merged_df_agr = aggregate(Abundance~sample_type+location_name+Taxon14, data=phylo_norm_newaxial_metazoa_glom_merged_df, FUN=mean) 

phylo_norm_newaxial_metazoa_glom_merged_df_agr$Taxon14 <- factor(
  phylo_norm_newaxial_metazoa_glom_merged_df_agr$Taxon14,
  levels = c("Annelida Taxon13", "Crustacea", "Cunina frugifera", "Cydippida Taxon12", "Gastropoda",
             "Phacellophora camtschatica (eggyolk jelly) Taxon13", "Polychaeta", "Siphonophorae Taxon13", "Thaliacea", "Others")
)

#Plot
ggplot(phylo_norm_newaxial_metazoa_glom_merged_df_agr, aes(x =location_name, y = Abundance, fill = Taxon14)) +
  ggtitle("New Axial Taxbar by Location") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top10_merged) +
  facet_nested(. ~ sample_type + location_name, scales = "free") +
  labs(x = "Location", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Gorda Taxbar ####
merge_top10_phylo_gorda <- function(phylo_norm_gorda_metazoa, top=9){
  transformed <- transform_sample_counts(phylo_norm_gorda_metazoa, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}


#Prepare dataframe
phylo_norm_gorda = subset_samples(phylo_notnorm_taxfix, Region == "Gorda Ridge")
phylo_norm_gorda_metazoa = subset_taxa(phylo_norm_gorda, Taxon7 == "Metazoa")

phylo_norm_gorda_metazoa_glom = tax_glom(phylo_norm_gorda_metazoa, taxrank = "Taxon14")
phylo_norm_gorda_metazoa_glom_merged = merge_top10_phylo_gorda(phylo_norm_gorda_metazoa_glom, top = 20)
phylo_norm_gorda_metazoa_glom_merged_df = psmelt(phylo_norm_gorda_metazoa_glom_merged)

phylo_norm_gorda_metazoa_glom_merged_df$Taxon14 <- phylo_norm_gorda_metazoa_glom_merged_df$Taxon14 %>%
  replace_na('Others')

phylo_norm_gorda_metazoa_glom_merged_df_agr = aggregate(Abundance~sample_type+location_name+Taxon14, data=phylo_norm_gorda_metazoa_glom_merged_df, FUN=mean) 

phylo_norm_gorda_metazoa_glom_merged_df_agr$Taxon14 <- factor(
  phylo_norm_gorda_metazoa_glom_merged_df_agr$Taxon14,
  levels = c("Annelida Taxon13", "Bargmannia elongata", "Bdelloidea", "Chelicerata", "Crustacea",
             "Cunina frugifera", "Cydippida Taxon12", "Earleria quadrata", "Eleutheria claparedei", "Enopla", "Gastropoda",
             "Gnathostomata", "Nematoda", "Pandea sp. AGC-2005", "Pantachogon haeckeli", "Pleurobrachia pileus Taxon13", 
             "Polychaeta", "Siphonophorae Taxon13", "Terrazoanthus onoi", "Others")
)

#Plot
ggplot(phylo_norm_gorda_metazoa_glom_merged_df_agr, aes(x =location_name, y = Abundance, fill = Taxon14)) +
  ggtitle("Gorda Ridge Taxbar by Location") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top20) +
  facet_nested(. ~ sample_type + location_name, scales = "free") +
  labs(x = "Location", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### MCR Taxbar ####
merge_top10_phylo_mcr <- function(phylo_norm_mcr_metazoa, top=9){
  transformed <- transform_sample_counts(phylo_norm_mcr_metazoa, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_norm_mcr = subset_samples(phylo_notnorm_taxfix, Region == "MCR")
phylo_norm_mcr_metazoa = subset_taxa(phylo_norm_mcr, Taxon7 == "Metazoa")

phylo_norm_mcr_metazoa_glom = tax_glom(phylo_norm_mcr_metazoa, taxrank = "Taxon14")
phylo_norm_mcr_metazoa_glom_merged = merge_top10_phylo_mcr(phylo_norm_mcr_metazoa_glom, top = 19)
phylo_norm_mcr_metazoa_glom_merged_df = psmelt(phylo_norm_mcr_metazoa_glom_merged)

phylo_norm_mcr_metazoa_glom_merged_df$Taxon14 <- phylo_norm_mcr_metazoa_glom_merged_df$Taxon14 %>%
  replace_na('Others')

phylo_norm_mcr_metazoa_glom_merged_df_agr = aggregate(Abundance~sample_type+location_name+Taxon14, data=phylo_norm_mcr_metazoa_glom_merged_df, FUN=mean) 

phylo_norm_mcr_metazoa_glom_merged_df_agr$Taxon14 <- factor(
  phylo_norm_mcr_metazoa_glom_merged_df_agr$Taxon14,
  levels = c("Aglantha digitale", "Annelida Taxon13", "Aphragmophora", "Beroe abyssicola Taxon13", "Crustacea", 
             "Ctenophora sp. Bahamas Taxon11", "Cydippida Taxon12", "Gastropoda", "Gnathostomata", "Jasonactis erythraios", 
             "Millepora sp. AGC-2001", "Monobrachium parasiticum", "Obeliida sp. USNM 1420678", "Polychaeta", "Porifera Taxon8", 
             "Siphonophorae Taxon13", "Tiaropsis multicirrata", "Trachymedusae Taxon13", "Zanclea sp. KA037", "Others")
)

#Plot
ggplot(phylo_norm_mcr_metazoa_glom_merged_df_agr, aes(x =location_name, y = Abundance, fill = Taxon14)) +
  ggtitle("MCR Taxbar by Location") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top20) +
  facet_nested(. ~ sample_type + location_name, scales = "free") +
  labs(x = "Location", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )


#### Cowplot Ordination ####
#Individual ordinations
phylo_norm_oldaxial_ord = ordinate(phylo_norm_oldaxial, "PCoA", "bray")
phylo_norm_newaxial_ord = ordinate(phylo_norm_newaxial, "PCoA", "bray")
phylo_norm_gorda_ord = ordinate(phylo_norm_gorda, "PCoA", "bray")
phylo_norm_mcr_ord = ordinate(phylo_norm_mcr, "PCoA", "bray")

sample_colors = c("Background" = "green", "Plume" = "orange", "Vent" = "red", "Incubation" = "blue", "Microcolonizer" = "violet")

p1 = plot_ordination(phylo_norm_oldaxial, phylo_norm_oldaxial_ord, type = "samples", color = "sample_type", title = "Old Axial") +
  geom_point(size = 3) +
  scale_color_manual(values = sample_colors) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), legend.position = "none") +
  annotate("text",
             x = Inf, y = Inf,
             label = "p = 0.015",
             hjust = 1, vjust = 1,
             size = 4)

p2 = plot_ordination(phylo_norm_newaxial, phylo_norm_newaxial_ord, type = "samples", color = "sample_type", title = "New Axial") +
  geom_point(size = 3) +
  scale_color_manual(values = sample_colors) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), legend.position = "none") +
  annotate("text",
             x = Inf, y = Inf,
             label = "p = 0.001",
             hjust = 1, vjust = 1,
             size = 4)

p3 = plot_ordination(phylo_norm_gorda, phylo_norm_gorda_ord, type = "samples", color = "sample_type", title = "Gorda Ridge") +
  geom_point(size = 3) +
  scale_color_manual(values = sample_colors) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) + 
  annotate("text",
             x = Inf, y = Inf,
             label = "p = 0.001",
             hjust = 1, vjust = 1,
             size = 4)

#Optionally remove or add legend to edit in photoshop on the final graph
 p3 = p3 + 
  labs(color = "Sample Type") + 
  theme(
    legend.title = element_text(size = 14, face = "bold"),
   legend.text = element_text(size = 12),              
    legend.key.size = unit(1.5, "lines"),                 
    legend.spacing = unit(0.5, "lines")         
  )

# gorda_legend = get_legend(p3)
# p3 = p3 + theme(legend.position = "none")

p4 = plot_ordination(phylo_norm_mcr, phylo_norm_mcr_ord, type = "samples", color = "sample_type", title = "MCR") +
  geom_point(size = 3) +
  scale_color_manual(values = sample_colors) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), legend.position = "none") +
  annotate("text",
             x = Inf, y = Inf,
             label = "p = 0.125",
             hjust = 1, vjust = 1,
             size = 4)

plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol = 4) 

#### Metazoa Top 20 Taxbar ####
colors_top20 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#000075",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "gray")

merge_top10_phylo_metazoa <- function(phylo_notnorm_taxfix_metazoa_glom, top=9){
  transformed <- transform_sample_counts(phylo_notnorm_taxfix_metazoa_glom, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_notnorm_taxfix_metazoa = subset_taxa(phylo_notnorm_taxfix, Taxon7 == "Metazoa")
phylo_notnorm_taxfix_metazoa_glom = tax_glom(phylo_notnorm_taxfix_metazoa, taxrank = "Taxon14")
phylo_notnorm_taxfix_metazoa_glom_merged = merge_top10_phylo_metazoa(phylo_notnorm_taxfix_metazoa_glom)
phylo_notnorm_taxfix_metazoa_glom_merged_df = psmelt(phylo_notnorm_taxfix_metazoa_glom_merged)


phylo_notnorm_taxfix_metazoa_glom_merged_df$Taxon14 <- phylo_notnorm_taxfix_metazoa_glom_merged_df$Taxon14 %>%
  replace_na('Others')

phylo_notnorm_taxfix_metazoa_glom_merged_df = aggregate(Abundance~Region+Sample+Taxon14, data=phylo_notnorm_taxfix_metazoa_glom_merged_df, FUN=mean) 


phylo_notnorm_taxfix_metazoa_glom_merged_df$Sample <- str_trunc(phylo_notnorm_taxfix_metazoa_glom_merged_df$Sample, 38)

phylo_notnorm_taxfix_metazoa_glom_merged_df$Region <- factor(
  phylo_notnorm_taxfix_metazoa_glom_merged_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

phylo_notnorm_taxfix_metazoa_glom_merged_df$Taxon14 <- factor(
  phylo_notnorm_taxfix_metazoa_glom_merged_df$Taxon14,
  levels = c("Annelida Taxon13", "Crustacea", "Cunina frugifera", "Cydippida Taxon12", "Gastropoda", 
             "Monogononta", "Pleurobrachia pileus Taxon13", "Polychaeta", "Siphonophorae Taxon13", "Others")
)

phylo_notnorm_taxfix_metazoa_glom_merged_df$Taxon14 <- factor(
  phylo_notnorm_taxfix_metazoa_glom_merged_df$Taxon14,
  levels = c("Annelida Taxon13", "Beroe forskalii Taxon13", "Crustacea", "Ctenophora sp. Bahamas Taxon11", "Cunina frugifera", 
             "Cydippida Taxon12", "Eleutheria claparedei", "Enopla", "Gastropoda", "Gnathostomata", 
             "Monogononta", "Pantachogon haeckeli", "Phacellophora camtschatica (eggyolk jelly) Taxon13", "Pleurobrachia pileus Taxon13", "Polychaeta", 
             "Porifera Taxon8", "Siphonophorae Taxon13", "Terrazoanthus onoi", "Thaliacea", "Others")
)

#Plot by sample type
ggplot(phylo_notnorm_taxfix_metazoa_glom_merged_df, aes(x = sample_type, y = Abundance, fill = Taxon14)) +
  ggtitle("Metazoa Taxbar by Sample Type") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top20) +
  facet_nested(. ~ Region + sample_type, scales = "free") +
  labs(x = "Sample Type", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#Plot by sample
ggplot(phylo_notnorm_taxfix_metazoa_glom_merged_df, aes(x = Sample, y = Abundance, fill = Taxon14)) +
  ggtitle("Metazoa Taxbar by Region") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top10_merged) +
  facet_wrap(. ~ Region, scales = "free", nrow = 1) +
  labs(x = "Sample", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )
               
               #_x", nrow = 1) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                 # plot.title = element_text(size = 10, face = "bold", hjust = 0.5))


#### Nematoda Taxbar ####
#Prepare dataframe
phylo_notnorm_taxfix_nosiders = subset_samples(phylo_notnorm_taxfix, Region != "Siders")
phylo_notnorm_taxfix_nematoda = subset_taxa(phylo_notnorm_taxfix_nosiders, Taxon14 == "Nematoda")
phylo_notnorm_taxfix_nematoda_glom = tax_glom(phylo_notnorm_taxfix_nematoda, taxrank = "Taxon22")
phylo_norm_taxfix_nematoda_glom = transform_sample_counts(phylo_notnorm_taxfix_nematoda_glom, function(x) x/sum(x))

nematoda_indices <- taxa_names(phylo_norm_taxfix_nematoda_glom)[tax_table(phylo_norm_taxfix_nematoda_glom)[,"Taxon14"] == "Nematoda"]
phylo_norm_taxfix_nematoda_glom_pruned <- prune_samples(sample_sums(prune_taxa(nematoda_indices, phylo_norm_taxfix_nematoda_glom)) > 0, phylo_norm_taxfix_nematoda_glom)

phylo_norm_taxfix_nematoda_glom_pruned_df = psmelt(phylo_norm_taxfix_nematoda_glom_pruned)

phylo_norm_taxfix_nematoda_glom_pruned_df$Taxon22 <- phylo_norm_taxfix_nematoda_glom_pruned_df$Taxon22 %>%
  replace_na('Others')

phylo_norm_taxfix_nematoda_glom_pruned_df$Sample <- str_trunc(phylo_norm_taxfix_nematoda_glom_pruned_df$Sample, 21)

phylo_norm_taxfix_nematoda_glom_pruned_df = aggregate(Abundance~Region+Sample+sample_type+Taxon22, data=phylo_norm_taxfix_nematoda_glom_pruned_df, FUN=mean) 

phylo_norm_taxfix_nematoda_glom_pruned_df$Region <- factor(
  phylo_norm_taxfix_nematoda_glom_pruned_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR")
)

#Plot
ggplot(phylo_norm_taxfix_nematoda_glom_pruned_df, aes(x = Sample, y = Abundance, fill = Taxon22)) +
  ggtitle("Nematoda Taxbar by Region") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top20) +
  facet_nested(. ~ Region+sample_type, scales = "free") +
  labs(x = "Sample", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # rotate x labels 90° and center on ticks
    axis.ticks.x = element_line(color = "black"),  # restore x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Ultimate Alpha Diversity Graph ####
#Prepare Old Axial dataframe
phylo_notnorm_oldaxial = subset_samples(phylo_notnorm, Region == "Old Axial")

alpha_div_observed_metadata_oldaxial = data.frame(sample_data(phylo_notnorm_oldaxial), 
                                         "Reads" = sample_sums(phylo_notnorm_oldaxial), 
                                         "Observed" = estimate_richness(phylo_notnorm_oldaxial, measures = "Observed"), 
                                         "Shannon" = estimate_richness(phylo_notnorm_oldaxial, measures = "Shannon"), 
                                         "InvSimpson" = estimate_richness(phylo_notnorm_oldaxial, measures = "InvSimpson"))

alpha_div_observed_metadata_oldaxial$Evenness = alpha_div_observed_metadata_oldaxial$Shannon / log(alpha_div_observed_metadata_oldaxial$Observed)

alpha_long_oldaxial = alpha_div_observed_metadata_oldaxial %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson, Evenness), names_to = "Metric", values_to = "Value")

metric_labels = c(Observed = "Observed ASVs", Shannon = "Shannon Diversity", InvSimpson = "Inverse Simpson", Evenness = "Evenness")

comparisons1 <- list(
  c("Background", "Plume"),
  c("Background", "Vent"),
  c("Plume", "Vent")
)


#Plot old axial
q1 = ggplot(data = alpha_long_oldaxial, aes(x = sample_type, y = Value)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comparisons1, method = "wilcox.test") +
  facet_wrap(~Metric, scales = "free_y", ncol = 1,  labeller = labeller(Metric = metric_labels)) + 
  ggtitle("Old Axial") + 
  labs(x = "Sample Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", color = "red", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "red", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#Prepare new axial dataframe
phylo_notnorm_newaxial = subset_samples(phylo_notnorm, Region == "New Axial")

alpha_div_observed_metadata_newaxial = data.frame(sample_data(phylo_notnorm_newaxial), 
                                                  "Reads" = sample_sums(phylo_notnorm_newaxial), 
                                                  "Observed" = estimate_richness(phylo_notnorm_newaxial, measures = "Observed"), 
                                                  "Shannon" = estimate_richness(phylo_notnorm_newaxial, measures = "Shannon"), 
                                                  "InvSimpson" = estimate_richness(phylo_notnorm_newaxial, measures = "InvSimpson"))

alpha_div_observed_metadata_newaxial$Evenness = alpha_div_observed_metadata_newaxial$Shannon / log(alpha_div_observed_metadata_newaxial$Observed)

alpha_long_newaxial = alpha_div_observed_metadata_newaxial %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson, Evenness), names_to = "Metric", values_to = "Value")

#Plot new axial
q2 = ggplot(data = alpha_long_newaxial, aes(x = sample_type, y = Value)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comparisons1, method = "wilcox.test") +
  facet_wrap(~Metric, scales = "free_y", ncol = 1,  labeller = labeller(Metric = metric_labels)) + 
  ggtitle("New Axial") + 
  labs(x = "Sample Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", color = "orange", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "orange", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#Prepare gorda dataframe
phylo_notnorm_gorda = subset_samples(phylo_notnorm, Region == "Gorda Ridge")

alpha_div_observed_metadata_gorda = data.frame(sample_data(phylo_notnorm_gorda), 
                                               "Reads" = sample_sums(phylo_notnorm_gorda), 
                                               "Observed" = estimate_richness(phylo_notnorm_gorda, measures = "Observed"), 
                                               "Shannon" = estimate_richness(phylo_notnorm_gorda, measures = "Shannon"), 
                                               "InvSimpson" = estimate_richness(phylo_notnorm_gorda, measures = "InvSimpson"))

alpha_div_observed_metadata_gorda$Evenness = alpha_div_observed_metadata_gorda$Shannon / log(alpha_div_observed_metadata_gorda$Observed)

alpha_long_gorda = alpha_div_observed_metadata_gorda %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson, Evenness), names_to = "Metric", values_to = "Value")

comparisons2 <- list(
  c("Background", "Plume"),
  c("Background", "Vent"),
  c("Background", "Incubation"),
  c("Background", "Microcolonizer"),
  c("Plume", "Vent"),
  c("Plume", "Incubation"),
  c("Plume", "Microcolonizer"),
  c("Vent", "Incubation"),
  c("Vent", "Microcolonizer"),
  c("Incubation", "Microcolonizer")
)

#Plot gorda
q3 = ggplot(data = alpha_long_gorda, aes(x = sample_type, y = Value)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comparisons2, method = "wilcox.test") +
  facet_wrap(~Metric, scales = "free_y", ncol = 1,  labeller = labeller(Metric = metric_labels)) + 
  ggtitle("Gorda Ridge") + 
  labs(x = "Sample Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", color = "green", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "green", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#Prepare MCR dataframe
phylo_notnorm_mcr = subset_samples(phylo_notnorm, Region == "MCR")

alpha_div_observed_metadata_mcr = data.frame(sample_data(phylo_notnorm_mcr), 
                                             "Reads" = sample_sums(phylo_notnorm_mcr), 
                                             "Observed" = estimate_richness(phylo_notnorm_mcr, measures = "Observed"), 
                                             "Shannon" = estimate_richness(phylo_notnorm_mcr, measures = "Shannon"), 
                                             "InvSimpson" = estimate_richness(phylo_notnorm_mcr, measures = "InvSimpson"))

alpha_div_observed_metadata_mcr$Evenness = alpha_div_observed_metadata_mcr$Shannon / log(alpha_div_observed_metadata_mcr$Observed)

alpha_long_mcr = alpha_div_observed_metadata_mcr %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson, Evenness), names_to = "Metric", values_to = "Value")

comparisons3 <- list(
  c("Background", "Plume"),
  c("Background", "Vent"),
  c("Background", "Incubation"),
  c("Plume", "Vent"),
  c("Plume", "Incubation"),
  c("Vent", "Incubation")
)

#Plot MCR
q4 = ggplot(data = alpha_long_mcr, aes(x = sample_type, y = Value)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comparisons3, method = "wilcox.test") +
  facet_wrap(~Metric, scales = "free_y", ncol = 1,  labeller = labeller(Metric = metric_labels)) + 
  ggtitle("MCR") + 
  labs(x = "Sample Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", color = "blue", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "blue", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

plot_grid(q1, q2, q3, q4, labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol = 4) 

#### Locations Ordination ####
colors_top10 <- c("#9A6324", "#46F0F0", "#1F78B4", "#AAFFC3", "#E6BEFF", 
                  "#33A02C", "#4363D8", "#008080", "#FB9A99", "#E6194B")

colors_top11 <- c("#9A6324", "#46F0F0", "#1F78B4", "#AAFFC3", "#E6BEFF", 
                  "#33A02C", "#4363D8", "#008080", "#FB9A99", "#E6194B", "#F032E6")

colors_top12 <- c("#9A6324", "#46F0F0", "#1F78B4", "#AAFFC3", "#E6BEFF", 
                  "#33A02C", "#4363D8", "#008080", "#FB9A99", "#E6194B", "#F032E6", "grey")

#Old axial ordination
phylo_norm_oldaxial_ord = ordinate(phylo_norm_oldaxial, "PCoA", "bray")

plot_ordination(phylo_norm_oldaxial, phylo_norm_oldaxial_ord, type = "samples", color = "location_name", title = "Old Axial by Location") +
  geom_point(size = 3) +
  scale_color_manual(values = colors_top11) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

#New axial ordination
phylo_norm_newaxial_ord = ordinate(phylo_norm_newaxial, "PCoA", "bray")

plot_ordination(phylo_norm_newaxial, phylo_norm_newaxial_ord, type = "samples", color = "location_name", title = "New Axial by Location") +
  geom_point(size = 3) +
  scale_color_manual(values = colors_top11) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

#Gorda ordination
phylo_norm_gorda_ord = ordinate(phylo_norm_gorda, "PCoA", "bray")

plot_ordination(phylo_norm_gorda, phylo_norm_gorda_ord, type = "samples", color = "location_name", title = "Gorda Ridge by Location") +
  geom_point(size = 3) +
  scale_color_manual(values = colors_top12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

#MCR ordination
phylo_norm_mcr_ord = ordinate(phylo_norm_mcr, "PCoA", "bray")

plot_ordination(phylo_norm_mcr, phylo_norm_mcr_ord, type = "samples", color = "location_name", title = "MCR by Location") +
  geom_point(size = 3) +
  scale_color_manual(values = colors_top11) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

#### Adonis by Location ####
#Old axial tests
phylo_oldaxial_dist = distance(phylo_norm_oldaxial, method = "bray")
adonis_data_oldaxial = as(sample_data(phylo_norm_oldaxial), "data.frame")
phylo_oldaxial_adonis = adonis2(phylo_oldaxial_dist ~ location_name, data = adonis_data_oldaxial)
write.csv(phylo_oldaxial_adonis, "phylo_oldaxial_adonis_location.csv")
oldaxial_disper = permutest(betadisper(phylo_oldaxial_dist, adonis_data_oldaxial$sample_type))
oldaxial_disper_df = as.data.frame(oldaxial_disper$tab)
write.csv(oldaxial_disper_df, "oldaxial_disper_location.csv")

#New axial tests
phylo_newaxial_dist = distance(phylo_norm_newaxial, method = "bray")
adonis_data_newaxial = as(sample_data(phylo_norm_newaxial), "data.frame")
phylo_newaxial_adonis = adonis2(phylo_newaxial_dist ~ location_name, data = adonis_data_newaxial)
write.csv(phylo_newaxial_adonis, "phylo_newaxial_adonis_location.csv")
newaxial_disper = permutest(betadisper(phylo_newaxial_dist, adonis_data_newaxial$sample_type))
newaxial_disper_df = as.data.frame(newaxial_disper$tab)
write.csv(newaxial_disper_df, "newaxial_disper_location.csv")

#Gorda tests
phylo_gorda_dist = distance(phylo_norm_gorda, method = "bray")
adonis_data_gorda = as(sample_data(phylo_norm_gorda), "data.frame")
phylo_gorda_adonis = adonis2(phylo_gorda_dist ~ location_name, data = adonis_data_gorda)
write.csv(phylo_gorda_adonis, "phylo_gorda_adonis_location.csv")
gorda_disper = permutest(betadisper(phylo_gorda_dist, adonis_data_gorda$sample_type))
gorda_disper_df = as.data.frame(gorda_disper$tab)
write.csv(gorda_disper_df, "gorda_disper_location.csv")

#MCR tests
phylo_mcr_dist = distance(phylo_norm_mcr, method = "bray")
adonis_data_mcr = as(sample_data(phylo_norm_mcr), "data.frame")
phylo_mcr_adonis = adonis2(phylo_mcr_dist ~ location_name, data = adonis_data_mcr)
write.csv(phylo_mcr_adonis, "phylo_mcr_adonis_location.csv")
mcr_disper = permutest(betadisper(phylo_mcr_dist, adonis_data_mcr$sample_type))
mcr_disper_df = as.data.frame(mcr_disper$tab)
write.csv(mcr_disper_df, "mcr_disper_location.csv")

#### Incubation sample type Taxbar ####
merge_top10_phylo_incu <- function(phylo_norm_incu_metazoa, top=9){
  transformed <- transform_sample_counts(phylo_norm_incu_metazoa, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_norm_incu = subset_samples(phylo_notnorm_taxfix, sample_type == "Incubation")
phylo_norm_incu_metazoa = subset_taxa(phylo_norm_incu, Taxon7 == "Metazoa")

phylo_norm_incu_metazoa_glom = tax_glom(phylo_norm_incu_metazoa, taxrank = "Taxon14")
phylo_norm_incu_metazoa_glom_merged = merge_top10_phylo_incu(phylo_norm_incu_metazoa_glom)
phylo_norm_incu_metazoa_glom_merged_df = psmelt(phylo_norm_incu_metazoa_glom_merged)

phylo_norm_incu_metazoa_glom_merged_df$Taxon14 <- phylo_norm_incu_metazoa_glom_merged_df$Taxon14 %>%
  replace_na('Others')

phylo_norm_incu_metazoa_glom_merged_df_agr = aggregate(Abundance~Sample+Region+Taxon14, data=phylo_norm_incu_metazoa_glom_merged_df, FUN=mean) 

phylo_norm_incu_metazoa_glom_merged_df_agr$Taxon14 <- factor(
  phylo_norm_incu_metazoa_glom_merged_df_agr$Taxon14,
  levels = c("Crustacea", "Cydippida Taxon12", "Eleutheria claparedei", "Gastropoda", "Gnathostomata",
             "Pantachogon haeckeli", "Pleurobrachia pileus Taxon13", "Polychaeta", "Siphonophorae Taxon13", "Others")
)

#Plot
ggplot(phylo_norm_incu_metazoa_glom_merged_df_agr, aes(x = Sample, y = Abundance, fill = Taxon14)) +
  ggtitle("Incubation Taxbar by Region") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top10_merged) +
  facet_wrap(. ~ Region, scales = "free") +
  labs(x = "Sample", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # rotate x labels 90° and center on ticks
    axis.ticks.x = element_line(color = "black"),  # restore x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Special Gorda Ordination ####
phylo_norm_gorda_ord = ordinate(phylo_norm_gorda, "PCoA", "bray")

plot_ordination(phylo_norm_gorda, phylo_norm_gorda_ord, type = "samples", color = "location_name", shape = "sample_type", title = "Gorda Ridge by Location") +
  geom_point(size = 3) +
  scale_color_manual(values = colors_top12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

#### Crustacea Taxbar ####
colors_top21 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#000075",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "black", "gray")

merge_top10_phylo_crustacea <- function(phylo_notnorm_taxfix_crustacea_glom, top=9){
  transformed <- transform_sample_counts(phylo_notnorm_taxfix_crustacea_glom, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

#Prepare dataframe
phylo_notnorm_taxfix_crustacea = subset_taxa(phylo_notnorm_taxfix, Taxon14 == "Crustacea")
phylo_notnorm_taxfix_crustacea_glom = tax_glom(phylo_notnorm_taxfix_crustacea, taxrank = "Taxon16")
phylo_notnorm_taxfix_crustacea_glom_merged = merge_top10_phylo_crustacea(phylo_notnorm_taxfix_crustacea_glom)
phylo_notnorm_taxfix_crustacea_glom_merged_df = psmelt(phylo_notnorm_taxfix_crustacea_glom_merged)

phylo_notnorm_taxfix_crustacea_glom_merged_df$Taxon22 <- phylo_notnorm_taxfix_crustacea_glom_merged_df$Taxon22 %>%
  replace_na('Others')

phylo_notnorm_taxfix_crustacea_glom_merged_df = aggregate(Abundance~Region+sample_type+Taxon22, data=phylo_notnorm_taxfix_crustacea_glom_merged_df, FUN=mean) 

phylo_notnorm_taxfix_crustacea_glom_merged_df$Region <- factor(
  phylo_notnorm_taxfix_crustacea_glom_merged_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

#Plot
ggplot(phylo_notnorm_taxfix_crustacea_glom_merged_df, aes(x = sample_type, y = Abundance, fill = Taxon22)) +
  ggtitle("Crustacea Taxbar by Sample Type") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top21) +
  facet_nested(. ~ Region + sample_type, scales = "free") +
  labs(x = "Sample Type", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Crustacea Coarse Taxbar ####
colors_top3 <- c("#E6194B", "#46F0F0", "#33A02C")

#Prepare dataframe
phylo_notnorm_taxfix_crustacea = subset_taxa(phylo_notnorm_taxfix, Taxon14 == "Crustacea")
phylo_notnorm_taxfix_crustacea_glom = tax_glom(phylo_notnorm_taxfix_crustacea, taxrank = "Taxon15")
phylo_notnorm_taxfix_crustacea_glom_coarse = transform_sample_counts(phylo_notnorm_taxfix_crustacea_glom, function(x) x/sum(x))
phylo_notnorm_taxfix_crustacea_glom_coarse_df = psmelt(phylo_notnorm_taxfix_crustacea_glom_coarse)

phylo_notnorm_taxfix_crustacea_glom_coarse_df = aggregate(Abundance~Region+sample_type+Taxon15, data=phylo_notnorm_taxfix_crustacea_glom_coarse_df, FUN=mean) 

phylo_notnorm_taxfix_crustacea_glom_coarse_df$Region <- factor(
  phylo_notnorm_taxfix_crustacea_glom_coarse_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

#Plot
ggplot(phylo_notnorm_taxfix_crustacea_glom_coarse_df, aes(x = sample_type, y = Abundance, fill = Taxon15)) +
  ggtitle("Crustacea Taxbar by Sample Type") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top3) +
  facet_nested(. ~ Region + sample_type, scales = "free") +
  labs(x = "Sample Type", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Eumalacostraca Taxbar ####
#Prepare dataframe
phylo_notnorm_taxfix_eumala = subset_taxa(phylo_notnorm_taxfix, Taxon16 == "Eumalacostraca")
phylo_notnorm_taxfix_eumala_glom = tax_glom(phylo_notnorm_taxfix_eumala, taxrank = "Taxon18")
phylo_notnorm_taxfix_eumala_glom_coarse_df = psmelt(phylo_notnorm_taxfix_eumala_glom)

phylo_notnorm_taxfix_eumala_glom_coarse_df = aggregate(Abundance~Region+sample_type+Taxon18, data=phylo_notnorm_taxfix_eumala_glom_coarse_df, FUN=mean) 

phylo_notnorm_taxfix_eumala_glom_coarse_df$Region <- factor(
  phylo_notnorm_taxfix_eumala_glom_coarse_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

#Plot
ggplot(phylo_notnorm_taxfix_eumala_glom_coarse_df, aes(x = sample_type, y = Abundance, fill = Taxon18)) +
  ggtitle("eumala Taxbar by Sample Type") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  facet_nested(. ~ Region + sample_type, scales = "free") +
  labs(x = "Sample Type", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Crustacea Read Count ####
colors_top3_matched <- c("Malacostraca" = "#E6194B", 
                         "Maxillopoda" = "#46F0F0", 
                         "Ostracoda" = "#33A02C")

#Prepare dataframe
crustacea_raw = subset_taxa(phylo_notnorm, Taxon14 == "Crustacea")
crustacea_raw_df = psmelt(crustacea_raw)

crustacea_genus_summary_raw = crustacea_raw_df %>%
  group_by(Taxon15) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(crustacea_genus_summary_raw, "CrustaceaSummary.csv", row.names = FALSE)

crustacea_genus_summary_raw$Taxon15 <- factor(crustacea_genus_summary_raw$Taxon15, levels = crustacea_genus_summary_raw$Taxon15[order(-crustacea_genus_summary_raw$total_reads)])

#Plot
ggplot(crustacea_genus_summary_raw, aes(x = Taxon15, y = total_reads, fill = Taxon15)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top3_matched) +
  xlab("Crustacea Class") +
  ylab("Total Reads") +
  ggtitle("Crustacea Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Gastropoda Coarse Taxbar ####
colors_top3 <- c("#E6194B", "#46F0F0", "#33A02C")

#Prepare dataframe
phylo_notnorm_taxfix_gastropoda = subset_taxa(phylo_notnorm_taxfix, Taxon14 == "Gastropoda")
phylo_notnorm_taxfix_gastropoda_glom = tax_glom(phylo_notnorm_taxfix_gastropoda, taxrank = "Taxon16")
phylo_notnorm_taxfix_gastropoda_glom_coarse = transform_sample_counts(phylo_notnorm_taxfix_gastropoda_glom, function(x) x/sum(x))
phylo_notnorm_taxfix_gastropoda_glom_coarse_df = psmelt(phylo_notnorm_taxfix_gastropoda_glom_coarse)

phylo_notnorm_taxfix_gastropoda_glom_coarse_df = aggregate(Abundance~Region+sample_type+Taxon16, data=phylo_notnorm_taxfix_gastropoda_glom_coarse_df, FUN=mean) 

phylo_notnorm_taxfix_gastropoda_glom_coarse_df$Region <- factor(
  phylo_notnorm_taxfix_gastropoda_glom_coarse_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

#Plot
ggplot(phylo_notnorm_taxfix_gastropoda_glom_coarse_df, aes(x = sample_type, y = Abundance, fill = Taxon16)) +
  ggtitle("Gastropoda Taxbar by Sample Type") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top12) +
  facet_nested(. ~ Region + sample_type, scales = "free") +
  labs(x = "Sample Type", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", size = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", size = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

#### Gastropoda Read Count ####
#Prepare dataframe
gastropoda_raw = subset_taxa(phylo_notnorm_taxfix, Taxon14 == "Gastropoda")
gastropoda_raw_df = psmelt(gastropoda_raw)

gastropoda_genus_summary_raw = gastropoda_raw_df %>%
  group_by(Taxon16) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(gastropoda_genus_summary_raw, "GastropodaSummary.csv", row.names = FALSE)

gastropoda_genus_summary_raw$Taxon16 <- factor(gastropoda_genus_summary_raw$Taxon16, levels = gastropoda_genus_summary_raw$Taxon16[order(-gastropoda_genus_summary_raw$total_reads)])

#Plot
ggplot(gastropoda_genus_summary_raw, aes(x = Taxon16, y = total_reads, fill = Taxon16)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("Gastropoda Species") +
  ylab("Total Reads") +
  ggtitle("Gastropoda Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Metazoa Read Count ####
#Prepare dataframe
metazoa_raw = subset_taxa(phylo_notnorm_taxfix, Taxon7 == "Metazoa")
metazoa_raw_df = psmelt(metazoa_raw)

metazoa_genus_summary_raw = metazoa_raw_df %>%
  group_by(Taxon14) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(metazoa_genus_summary_raw, "MetazoaSummary.csv", row.names = FALSE)

metazoa_genus_summary_raw_condensed = metazoa_genus_summary_raw %>%
  mutate(Taxon14 = if_else(row_number() > 19, "Others", Taxon14)) %>%
  group_by(Taxon14) %>%
  summarise(total_reads = sum(total_reads)) %>%
  ungroup()

metazoa_genus_summary_raw_condensed$Taxon14 <- factor(metazoa_genus_summary_raw_condensed$Taxon14, levels = metazoa_genus_summary_raw_condensed$Taxon14[order(-metazoa_genus_summary_raw_condensed$total_reads)])

#Plot
ggplot(metazoa_genus_summary_raw_condensed, aes(x = Taxon14, y = total_reads, fill = Taxon14)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("Gastropoda Species") +
  ylab("Total Reads") +
  ggtitle("Metazoa Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Polychaeta Read Count #####
#Prepare dataframe
polychaeta_raw = subset_taxa(phylo_notnorm, Taxon14 == "Polychaeta")
polychaeta_raw_df = psmelt(polychaeta_raw)

polychaeta_genus_summary_raw = polychaeta_raw_df %>%
  group_by(Taxon16) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(polychaeta_genus_summary_raw, "PolychaetaSummary.csv", row.names = FALSE)

polychaeta_genus_summary_raw$Taxon16 <- factor(polychaeta_genus_summary_raw$Taxon16, levels = polychaeta_genus_summary_raw$Taxon16[order(-polychaeta_genus_summary_raw$total_reads)])

#Plot
ggplot(polychaeta_genus_summary_raw, aes(x = Taxon16, y = total_reads, fill = Taxon16)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("Polychaeta Taxa") +
  ylab("Total Reads") +
  ggtitle("Polychaeta Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Terebellida Read Count ####
#Prepare dataframe
terebellida_raw = subset_taxa(phylo_notnorm, Taxon16 == "Terebellida")
terebellida_raw_df = psmelt(terebellida_raw)

terebellida_genus_summary_raw = terebellida_raw_df %>%
  group_by(Taxon17) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(terebellida_genus_summary_raw, "TerebellidaSummary.csv", row.names = FALSE)

terebellida_genus_summary_raw$Taxon17 <- factor(terebellida_genus_summary_raw$Taxon17, levels = terebellida_genus_summary_raw$Taxon17[order(-terebellida_genus_summary_raw$total_reads)])

#Plot
ggplot(terebellida_genus_summary_raw, aes(x = Taxon17, y = total_reads, fill = Taxon17)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("Terebellida Taxa") +
  ylab("Total Reads") +
  ggtitle("terebellida Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Chelicerata Read Count ####
#Prepare dataframe
chelicerata_raw = subset_taxa(phylo_notnorm, Taxon14 == "Chelicerata")
chelicerata_raw_df = psmelt(chelicerata_raw)

chelicerata_genus_summary_raw = chelicerata_raw_df %>%
  group_by(Taxon17) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(chelicerata_genus_summary_raw, "ChelicerataSummary.csv", row.names = FALSE)

chelicerata_genus_summary_raw$Taxon17 <- factor(chelicerata_genus_summary_raw$Taxon17, levels = chelicerata_genus_summary_raw$Taxon17[order(-chelicerata_genus_summary_raw$total_reads)])

#Plot
ggplot(chelicerata_genus_summary_raw, aes(x = Taxon17, y = total_reads, fill = Taxon17)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("chelicerata Taxa") +
  ylab("Total Reads") +
  ggtitle("Chelicerata Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Rhabditophora Read Count ####
#Prepare dataframe
rhabditophora_raw = subset_taxa(phylo_notnorm_taxfix, Taxon14 == "Rhabditophora")
rhabditophora_raw_df = psmelt(rhabditophora_raw)

rhabditophora_genus_summary_raw = rhabditophora_raw_df %>%
  group_by(Taxon17) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(rhabditophora_genus_summary_raw, "RhabditophoraSummary.csv", row.names = FALSE)

rhabditophora_genus_summary_raw$Taxon16 <- factor(rhabditophora_genus_summary_raw$Taxon16, levels = rhabditophora_genus_summary_raw$Taxon16[order(-rhabditophora_genus_summary_raw$total_reads)])

#Plot
ggplot(rhabditophora_genus_summary_raw, aes(x = Taxon16, y = total_reads, fill = Taxon16)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("rhabditophora Taxa") +
  ylab("Total Reads") +
  ggtitle("Rhabditophora Read Count2") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Enteropneusta Read Count ####
#Prepare dataframe
enteropneusta_raw = subset_taxa(phylo_notnorm, Taxon14 == "Enteropneusta")
enteropneusta_raw_df = psmelt(enteropneusta_raw)

enteropneusta_genus_summary_raw = enteropneusta_raw_df %>%
  group_by(Taxon16) %>%
  summarize(
    total_reads = sum(Abundance),
    asv_count = n_distinct(OTU)
  ) %>%
  arrange(desc(total_reads))

write.csv(enteropneusta_genus_summary_raw, "EnteropneustaSummary.csv", row.names = FALSE)

enteropneusta_genus_summary_raw$Taxon16 <- factor(enteropneusta_genus_summary_raw$Taxon16, levels = enteropneusta_genus_summary_raw$Taxon16[order(-enteropneusta_genus_summary_raw$total_reads)])

#Plot
ggplot(enteropneusta_genus_summary_raw, aes(x = Taxon16, y = total_reads, fill = Taxon16)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("Enteropneusta Taxa") +
  ylab("Total Reads") +
  ggtitle("Enteropneusta Read Count") +
  geom_text(aes(label = round(total_reads, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### Meiofauna vs Macrofauna ####
meiovsmacro = read.csv("MeiovsMacro.csv")
meiovsmacro_agr_readcount = aggregate(ReadCount~Classification, data=meiovsmacro, FUN=sum)
meiovsmacro_agr_readcount$Classification <- factor(meiovsmacro_agr_readcount$Classification, levels = meiovsmacro_agr_readcount$Classification[order(-meiovsmacro_agr_readcount$ReadCount)])

#Bar plot of meiofauna vs macrofauna and other groups readcount
ggplot(meiovsmacro_agr_readcount, aes(x = Classification, y = ReadCount, fill = Classification)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20) +
  xlab("Classification") +
  ylab("Total Reads") +
  ggtitle("Meio vs Macro Read Count") +
  geom_text(aes(label = round(ReadCount, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

meiovsmacro_agr_asvcount = aggregate(ASVCount~Classification, data=meiovsmacro, FUN=sum)
meiovsmacro_agr_asvcount$Classification <- factor(meiovsmacro_agr_asvcount$Classification, levels = meiovsmacro_agr_asvcount$Classification[order(-meiovsmacro_agr_asvcount$ASVCount)])

colors_top20_mod <- c("#46f0f0", "#9a6324", "#1F78B4", "#000075", "#aaffc3",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "gray")

#Barplot for ASV count
ggplot(meiovsmacro_agr_asvcount, aes(x = Classification, y = ASVCount, fill = Classification)) +
  geom_bar(stat = "identity", color = "black") +  # Border around each bar
  scale_fill_manual(values = colors_top20_mod) +
  xlab("Classification") +
  ylab("Total ASVs") +
  ggtitle("Meio vs Macro ASV Count") +
  geom_text(aes(label = round(ASVCount, 1)), vjust = -0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Frame around plot
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  # Bold x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  )

#### New Phyloseq ####

old_tax_table = as.data.frame(tax_table(metazoa_raw))
write.csv(old_tax_table, "old_tax_table.csv", row.names = TRUE)
#Add a new column for classification
new_tax_table = read.csv("old_tax_table.csv")

#Assign classification into the new tax table
for (i in 1:942){
  FakePhylum = new_tax_table$Taxon14[i]
  if (FakePhylum %in% meiovsmacro$Taxon14) {
    meiovsmacro_sub = subset(meiovsmacro, Taxon14 == FakePhylum)
    PreciseTaxon = meiovsmacro_sub$PreciseTaxon[1]
      for (x in 1:nrow(meiovsmacro_sub)) {
        if (new_tax_table[i, PreciseTaxon] == meiovsmacro_sub[x, "Identity"]) {
          new_tax_table$Classification[i] = meiovsmacro_sub$Classification[x]
          new_tax_table$Phylum[i] = meiovsmacro_sub$Phylum[x]
        }
      }
  }
  else {
    new_tax_table$Classification[i] = "Indeterminate"
    new_tax_table$Phylum[i] = "Indeterminate"
  }
}

#Build new phyloseq object from tax table with classifications
row.names(new_tax_table) = new_tax_table$Sample
new_tax_table$Sample = NULL
new_tax_matrix = as.matrix(new_tax_table)
old_otu_table = otu_table(metazoa_raw)
old_meta_table = sample_data(metazoa_raw)

phylo_class <- phyloseq(old_otu_table, tax_table(new_tax_matrix), old_meta_table)

#### Phylo Class Taxbar ####
#Prepare dataframe
phylo_class_glom = tax_glom(phylo_class, taxrank = "Classification")
phylo_class_glom_norm = transform_sample_counts(phylo_class_glom, function(x) x/sum(x))

phylo_class_glom_norm_df = psmelt(phylo_class_glom_norm)
phylo_class_glom_norm_df_agr = aggregate(Abundance~Sample+Region+sample_type+Classification, data=phylo_class_glom_norm_df, FUN=mean)

phylo_class_glom_norm_df$Region <- factor(
  phylo_class_glom_norm_df$Region,
  levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders")
)

#Plot
ggplot(phylo_class_glom_norm_df, aes(x = Sample, y = Abundance, fill = Classification)) +
  ggtitle("Classified Taxbar by Sample Type") +
  geom_bar(stat = "identity", color = "black") +  # black borders around bars
  scale_fill_manual(values = colors_top10_merged) +
  facet_nested(. ~Region+sample_type, scales = "free") +
  labs(x = "Sample Type", y = "Abundance") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # centered, bold title
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold", size = 7),  # slightly smaller bold facet text
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # black frame around each facet
    strip.background = element_rect(color = "black", fill = "white", linewidth = 1),  # black box around facet labels
    panel.grid = element_line(color = "grey70", linewidth = 0.2),  # grey gridlines
    legend.background = element_blank(),  # remove legend box
    legend.key = element_blank()  # remove legend key borders
  )

plot_bar(phylo_class_norm, fill = "Classification")

#### Tree Setup ####
#Import tree
tree <- read_qza("rooted-18S-tree.qza")$data
phy_tree(phylo_notnorm_taxfix) <- tree

ggtree(phy_tree(phylo_notnorm_taxfix)) +
  geom_tiplab(size = 2)

#### Nematode Tree ####
#Subset to only nematodes
nema_tree = subset_taxa(phylo_class, Taxon14 == "Nematoda")

#Prepare dataframe
tax_tree <- as.data.frame(tax_table(nema_tree)) %>%
  rownames_to_column("ASV") %>%
  select(ASV, Family = Taxon21, Genus = Taxon22, Classification) %>%
  mutate(tree_label = ifelse(is.na(Genus), ASV, paste0(Genus, " (", ASV, ")"))) %>%
  select(ASV, tree_label, Family, Classification)

#Plot
ggtree(phy_tree(nema_tree)) %<+% tax_tree +
  geom_tiplab(aes(label = tree_label), size = 3) +
  geom_fruit(geom=geom_star,
             mapping=aes(fill=Family, starshape=Classification, size = 2),
             position="identity",starstroke=0.1) +
  scale_fill_manual(values = colors_top20) +
  theme(
    legend.position = c(0.01, 0.99),          # top left corner inside plot
    legend.justification = c("left", "top"),  # anchor the legend box to bottom-right
    legend.background = element_rect(
      fill = "white",         # white background inside legend box
      color = "black",        # black border around legend box
      size = 0.8,             # border thickness
      linetype = "solid"
    ),
    legend.key = element_rect(fill = "white")  # make legend keys background white too
  )

#### Classification Tree ####
#Build phyloseq
phylo_class <- phyloseq(old_otu_table, tax_table(new_tax_matrix), old_meta_table, tree)

#Prepare dataframe
tax_tree_class <- as.data.frame(tax_table(phylo_class)) %>%
  rownames_to_column("ASV") %>%
  select(ASV, Phylum = Taxon14, Genus = Taxon22, Classification) %>%
  mutate(tree_label = ifelse(is.na(Genus), ASV, paste0(Genus, " (", ASV, ")"))) %>%
  select(ASV, Phylum, tree_label, Classification)

#Plot
ggtree(phy_tree(phylo_class)) %<+% tax_tree_class +
  geom_tiplab(aes(label = tree_label), size = 3) +
  geom_fruit(geom=geom_star,
             mapping=aes(fill=Phylum, starshape=Classification),
             position="identity",starstroke=0.1)+
  theme(legend.position = "none")

#### FASTA to CSV ####
fasta_file <- "ref-seqs-deepsea-DB_052025.fasta"
seqs <- readDNAStringSet(fasta_file)

#Turn fasta file into a dataframe
asv_df <- data.frame(
  ASV_ID = names(seqs),
  Sequence = as.character(seqs),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(SequenceCount = nchar(Sequence))

highlight_id <- "4b7bf6d9767643e6adee58fb9fa0bfcd"

asv_df <- asv_df %>%
  mutate(highlight = ifelse(ASV_ID == highlight_id, "highlight", "other"))

#Plot with a highlight taxon
ggplot(asv_df, aes(x = "", y = SequenceCount)) +
  geom_boxplot(fill = "gray90", outlier.shape = NA) +
  geom_jitter(aes(color = highlight), width = 0.1, alpha = 0.8, size = 2.5) +
  scale_color_manual(values = c("highlight" = "yellow", "other" = "black")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  ylab("Sequence Length (bp)") +
  ggtitle("Distribution of Sequence Lengths with Highlighted ASV")

target_ids <- c("4b7bf6d9767643e6adee58fb9fa0bfcd",
                "0038478be7fb4f097ce93a5e9341af2a",
                "ce991f459b179a1514f559e790b24a9a")

asv_subset <- asv_df %>%
  filter(ASV_ID %in% target_ids)

write.csv(asv_subset, "subset_asvs.csv", row.names = FALSE)

#### Oddball Read Count ####
#Prepare dataframe
panta_raw = subset_taxa(phylo_notnorm, Taxon14 == "Pantachogon haeckeli")
panta_raw_df = psmelt(panta_raw)

panta_genus_summary_raw = panta_raw_df %>%
  group_by(OTU) %>%  # OTU is the ASV ID column in psmelt()
  summarise(TotalReads = sum(Abundance)) %>%
  arrange(desc(TotalReads))

write.csv(panta_genus_summary_raw, "pantaSummary.csv", row.names = FALSE)

#### Tree Without Oddball ####

colors_phyla <- c("#9a6324", "#46f0f0", "#1F78B4", "#ffe119", "#000075",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "gray", "black")

colors_region <- c("#9A6324", "#46F0F0", "#1F78B4", "#AAFFC3", "#E6194B")

colors_sample_type <- c("#33A02C", "#4363D8", "#008080", "#FB9A99", "#f032e6", "#bcf60c")

#Filter phyloseq to remove oddball and indeterminates
phylo_class <- phyloseq(old_otu_table, tax_table(new_tax_matrix), old_meta_table, tree)
phylo_class_filtered = subset_taxa(phylo_class, taxa_names(phylo_class) != "4b7bf6d9767643e6adee58fb9fa0bfcd")
phylo_class_filtered2 = subset_taxa(phylo_class, Phylum != "Indeterminate")

#Filter phyloseq to remove any taxa with less than 600 reads
asv_totals <- taxa_sums(phylo_class_filtered2)
keep_taxa <- names(asv_totals[asv_totals >= 600])
phylo_class_filtered3 <- prune_taxa(keep_taxa, phylo_class_filtered2)

#Modify the asv dataframe so that it reflects the filters on the phyloseq object
class_asv <- as.data.frame(otu_table(phylo_class_filtered3))
class_asv$ASV <- rownames(class_asv)
class_meta <- as(sample_data(phylo_class_filtered3), "data.frame")
class_meta$Sample <- rownames(class_meta)
asv_long <- pivot_longer(class_asv, 
                         cols = -ASV, 
                         names_to = "Sample", 
                         values_to = "Abundance")
asv_long <- merge(asv_long, class_meta, by = "Sample")
asv_long$Presence <- asv_long$Abundance > 0 #Unusesd

#Prepare dataframe
tax_tree_class <- as.data.frame(tax_table(phylo_class_filtered3)) %>%
  rownames_to_column("ASV")

tax_tree_class = merge(tax_tree_class, asv_long, by = "ASV") %>%
  select(ASV, Phylum, Classification, Region, Abundance)

phylum_levels <- setdiff(unique(tax_tree_class$Phylum), "Indeterminate")
phylum_levels <- c(phylum_levels, "Indeterminate") #Make sure indeterminate is listed at the end
tax_tree_class$Phylum <- factor(tax_tree_class$Phylum, levels = phylum_levels)

classification_levels <- c(
  "Meiofauna",
  "Macrofauna",
  "Temporary",
  "Jellyfish",
  "Chordates",
  "Parasites",
  "Indeterminate"
)

tax_tree_class$Classification <- factor(
  tax_tree_class$Classification,
  levels = classification_levels
)

tax_tree_class_cleaned <- tax_tree_class %>%
  group_by(ASV, Region) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    tax_tree_class %>% distinct(ASV, Region, Phylum, Classification),
    by = c("ASV", "Region")
  )

tax_tree_class_wide <- tax_tree_class_cleaned %>%
  mutate(Region = paste0(Region, " Abundance")) %>%  # Rename regions to include "Abundance" after them
  pivot_wider(
    names_from = Region,
    values_from = Abundance,
    values_fill = 0  #Fill missing values with 0
  )

#Scale abundances so that they form an event gradient as alpha values (0-1)
tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(Abundance = rowSums(select(., contains("Abundance")), na.rm = TRUE))

tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(Abundance_Scaled = scales::rescale(Abundance, to = c(0.1, 5)))

tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(`Region Abundance` = scales::rescale(log10(`Old Axial Abundance` + 1), to = c(0, 1)))

tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(`New Axial Alpha` = scales::rescale(log10(`New Axial Abundance` + 1), to = c(0, 1)))

tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(`Gorda Ridge Alpha` = scales::rescale(log10(`Gorda Ridge Abundance` + 1), to = c(0, 1)))

tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(`MCR Alpha` = scales::rescale(log10(`MCR Abundance` + 1), to = c(0, 1)))

tax_tree_class_wide <- tax_tree_class_wide %>%
  mutate(`Siders Alpha` = scales::rescale(log10(`Siders Abundance` + 1), to = c(0, 1)))

#Plot
ggtree(phy_tree(phylo_class_filtered3), layout = "circular", size = 0.15, open.angle = 50, branch.length = "none") %<+% tax_tree_class_wide +
  geom_fruit(
    geom = geom_star,
    mapping = aes(fill = Phylum, starshape = Classification),
    position = "identity",
    offset = 0.025,
    starstroke = 0.1,
    size = 2.5  # increase the star size
  ) +
  scale_fill_manual(values = colors_phyla) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(alpha = `Region Abundance`),   # control darkness by Abundance
    fill = "red",                        # fixed fill color (not mapped → no legend)
    color = "black",                     # black frame around tiles
    width = 0.8,
    height = 0.8,
    offset = 0.03
  ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(alpha = `New Axial Alpha`),   # control darkness by Abundance
    fill = "orange",                        # fixed fill color (not mapped → no legend)
    color = "black",                     # black frame around tiles
    width = 0.8,
    height = 0.8,
    offset = 0.031
  ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(alpha = `Gorda Ridge Alpha`),   # control darkness by Abundance
    fill = "green",                        # fixed fill color (not mapped → no legend)
    color = "black",                     # black frame around tiles
    width = 0.8,
    height = 0.8,
    offset = 0.032
  ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(alpha = `MCR Alpha`),   # control darkness by Abundance
    fill = "blue",                        # fixed fill color (not mapped → no legend)
    color = "black",                     # black frame around tiles
    width = 0.8,
    height = 0.8,
    offset = 0.033
  ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(alpha = `Siders Alpha`),   # control darkness by Abundance
    fill = "violet",                        # fixed fill color (not mapped → no legend)
    color = "black",                     # black frame around tiles
    width = 0.8,
    height = 0.8,
    offset = 0.034
  ) +
  scale_alpha_continuous(range = c(0, 1)) +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = ASV, x = Abundance_Scaled),
    fill = "black",
    stat = "identity",
    orientation = "y",
    width = 0.8,
    offset = 0.035
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",  # stack legends vertically on the right
    plot.margin = margin(10, 100, 10, 100),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      size = 0.8,
      linetype = "solid"
    ),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(hjust = 0.5),
  )

legend_df <- data.frame(
  Region = factor(
    c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders", "Total"),
    levels = c("Old Axial", "New Axial", "Gorda Ridge", "MCR", "Siders", "Total")
  ),
  x = 1,
  y = 1:6
)

#Separate legend figure to edit in photoshop with original ggtree
ggplot(legend_df, aes(x = x, y = y, fill = Region)) +
  geom_tile(color = "black", width = 1, height = 0.6) +  # rectangles with black border
  scale_fill_manual(
    name = "Region",
    values = c(
      "Old Axial" = "red",
      "New Axial" = "orange",
      "Gorda Ridge" = "green",
      "MCR" = "blue",
      "Siders" = "violet",
      "Total" = "black"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",  # stack legends vertically on the right
    plot.margin = margin(10, 100, 10, 100),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      size = 0.8,
      linetype = "solid"
    ),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(hjust = 0.5),
  )
