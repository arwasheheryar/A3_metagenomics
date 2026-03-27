# ============================================================
# BINF6110 Assignment 3 - Shotgun Metagenomics Analysis
# Samples: De Filippis et al. 2019 (SRP126540)
# 3 vegan (SRR8146990, SRR8146961, SRR8146989)
# 3 omnivore (SRR8146956, SRR8146969, SRR8146975)
# ============================================================


# ============================================================
# BINF6110 Assignment 3 - Shotgun Metagenomics Analysis
# Samples: De Filippis et al. 2019 (SRP126540)
# 3 vegan:     SRR8146990, SRR8146961, SRR8146989
# 3 omnivore:  SRR8146956, SRR8146969, SRR8146975
# ============================================================

# ---- 0. Install packages (run this section ONCE, then comment it out) ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq",       update = FALSE, ask = FALSE)
BiocManager::install("biomformat",     update = FALSE, ask = FALSE)
BiocManager::install("MicrobiomeStat", update = FALSE, ask = FALSE)

install.packages("ggplot2",      repos = "https://cran.rstudio.com", quiet = TRUE)
install.packages("vegan",        repos = "https://cran.rstudio.com", quiet = TRUE)
install.packages("dplyr",        repos = "https://cran.rstudio.com", quiet = TRUE)
install.packages("tidyr",        repos = "https://cran.rstudio.com", quiet = TRUE)
install.packages("RColorBrewer", repos = "https://cran.rstudio.com", quiet = TRUE)

# ---- 1. Load libraries ----
library(phyloseq)
library(biomformat)
library(MicrobiomeStat)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ---- 2. Import BIOM table ----
# Find your biom file first — run this to confirm the path:
list.files(path.expand("~"), pattern = "table.biom", 
           recursive = TRUE, full.names = TRUE)

# Then set the correct path (adjust if the search above shows a different location):
biom_file <- "~/assignment3/results/table.biom"

ps <- import_biom(biom_file)

# ---- 3. Verify import ----
ps
tax_table(ps)[1:5, ]
sample_names(ps)
otu_table(ps)[1:5, ]

# ---- 3. Verify import ----
ps
tax_table(ps)[1:5, ]
sample_names(ps)
otu_table(ps)[1:5, ]

BiocManager::install("MicrobiomeStat")
library(MicrobiomeStat)

# ---- 2. Import BIOM table ----
biom_file <- "~/A3_metagenomics/results/table.biom"
biom_data <- read_biom(biom_file)
ps <- import_biom(biom_file)
ps



# ---- 3. Add metadata ----
# Manually define sample metadata: diet group for each SRR ID
sample_metadata <- data.frame(
  SampleID = c("SRR8146990", "SRR8146961", "SRR8146989",
               "SRR8146956", "SRR8146969", "SRR8146975"),
  Diet     = c("Vegan", "Vegan", "Vegan",
               "Omnivore", "Omnivore", "Omnivore"),
  row.names = c("SRR8146990", "SRR8146961", "SRR8146989",
                "SRR8146956", "SRR8146969", "SRR8146975")
)

sample_data(ps) <- sample_data(sample_metadata)

# Verify your object looks right
ps
sample_data(ps)
tax_table(ps)


# ---- 3b. Rename taxonomic ranks (if needed) ----
# kraken-biom uses "Rank1", "Rank2" etc. — rename them
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class",
                             "Order", "Family", "Genus", "Species")


# ============================================================
# FIGURE 1: Taxonomic Abundance Bar Chart (Phylum level)
# ============================================================

# Convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Agglomerate at Phylum level
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)

# Melt into a long dataframe for ggplot
df_phylum <- psmelt(ps_phylum)

# Keep only top 10 phyla by mean relative abundance across samples
top_phyla <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

df_phylum_top <- df_phylum %>%
  mutate(Phylum = ifelse(Phylum %in% top_phyla, Phylum, "Other"))

# Set a colour palette
phylum_colors <- c(brewer.pal(10, "Set3"), "grey70")

# Plot
fig1 <- ggplot(df_phylum_top, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Diet, scales = "free_x") +
  scale_fill_manual(values = phylum_colors) +
  labs(
    title = "Relative Abundance by Phylum",
    x = "Sample",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

fig1
ggsave("~/A3_metagenomics/results/fig1_phylum_abundance.pdf",
       fig1, width = 10, height = 6)


# ============================================================
# FIGURE 2: Alpha Diversity (Observed Richness + Shannon)
# ============================================================

# Calculate alpha diversity metrics
# Note: Bracken/Kraken data retains singletons, so Chao1 is appropriate here
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Chao1"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- merge(alpha_div, sample_metadata, by.x = "SampleID", by.y = "row.names")

# Pivot to long format for faceted plot
alpha_long <- alpha_div %>%
  pivot_longer(cols = c("Observed", "Shannon"),
               names_to = "Metric", values_to = "Value")

# Set factor levels so Vegan comes first
alpha_long$Diet <- factor(alpha_long$Diet, levels = c("Vegan", "Omnivore"))

fig2 <- ggplot(alpha_long, aes(x = Diet, y = Value, colour = Diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, size = 3) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_colour_manual(values = c("Vegan" = "#4DAF4A", "Omnivore" = "#E41A1C")) +
  labs(
    title = "Alpha Diversity by Diet Group",
    x = "Diet",
    y = "Value"
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none")

fig2
ggsave("~/A3_metagenomics/results/fig2_alpha_diversity.pdf",
       fig2, width = 7, height = 5)

# Statistical test (Wilcoxon, given n=3 per group)
wilcox.test(Observed ~ Diet, data = alpha_div)
wilcox.test(Shannon  ~ Diet, data = alpha_div)


# ============================================================
# FIGURE 3: Beta Diversity — Bray-Curtis PCoA + PERMANOVA
# ============================================================

# Ordinate with Bray-Curtis dissimilarity
ord_bray <- ordinate(ps_rel, method = "PCoA", distance = "bray")

# Extract % variance explained by each axis
eig_vals <- ord_bray$values$Eigenvalues
pct_var  <- round(100 * eig_vals / sum(eig_vals[eig_vals > 0]), 1)

fig3 <- plot_ordination(ps_rel, ord_bray, color = "Diet") +
  geom_point(size = 5) +
  scale_colour_manual(values = c("Vegan" = "#4DAF4A", "Omnivore" = "#E41A1C")) +
  labs(
    title = "Beta Diversity: Bray-Curtis PCoA",
    x = paste0("PC1 [", pct_var[1], "%]"),
    y = paste0("PC2 [", pct_var[2], "%]"),
    colour = "Diet"
  ) +
  theme_bw() +
  theme(legend.title = element_text(face = "bold"))

fig3
ggsave("~/A3_metagenomics/results/fig3_beta_pcoa.pdf",
       fig3, width = 7, height = 5)

# PERMANOVA
metadata_df <- as(sample_data(ps_rel), "data.frame")
bray_dist   <- phyloseq::distance(ps_rel, method = "bray")

set.seed(42)
perm_result <- adonis2(bray_dist ~ Diet, data = metadata_df, permutations = 999)
print(perm_result)


# ============================================================
# FIGURE 4: Differential Abundance — ANCOMBC2
# ============================================================

# Run ANCOMBC2 at the Genus level
# prv_cut = 0.5 means a genus must appear in at least 50% of samples
ancom_out <- ancombc2(
  data         = ps,
  tax_level    = "Genus",
  fix_formula  = "Diet",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens  = TRUE,
  prv_cut      = 0.5,       # remove very rare taxa
  lib_cut      = 1000,
  s0_perc      = 0.05,
  group        = "Diet",
  struc_zero   = TRUE,
  neg_lb       = TRUE
)

# Inspect results
head(ancom_out$res)

# Subset to significant results (q < 0.05)
ancom_sig <- subset(ancom_out$res, q_DietVegan < 0.05)
ancom_sig

# Check structural zeroes
ancom_out$zero_ind

# Build the full result dataframe for plotting
ancom_res <- ancom_out$res %>%
  mutate(significant = q_DietVegan < 0.05)

# LFC lollipop plot — all genera, coloured by significance
fig4 <- ggplot(ancom_res,
               aes(x = lfc_DietVegan,
                   y = reorder(taxon, lfc_DietVegan),
                   colour = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lfc_DietVegan - se_DietVegan,
                    xmax = lfc_DietVegan + se_DietVegan),
                width = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#E41A1C"),
                      labels = c("ns", "q < 0.05")) +
  labs(
    title = "Differential Abundance: Vegan vs. Omnivore (ANCOMBC2)",
    x     = "Log Fold Change (Vegan vs. Omnivore)",
    y     = "Genus",
    colour = "Significance"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7))

fig4
ggsave("~/A3_metagenomics/results/fig4_differential_abundance.pdf",
       fig4, width = 9, height = 7)