#!/usr/bin/env Rscript

library(jsonlite)
library(dplyr)
library(ggplot2)
library(gggenes)
library(tidyr)
library(optparse)
library(scales)

option_list = list(
  make_option(c("-i", "--input"), type="character", help="selected_representatives.json"),
  make_option(c("-o", "--output"), type="character", default="operon_synteny.svg")
)

opt = parse_args(OptionParser(option_list=option_list))

data <- fromJSON(opt$input, simplifyDataFrame = FALSE)

df_list <- lapply(data, function(rep) {

  genes <- rep$genes
  if (length(genes) == 0) return(NULL)

  ribo_start <- rep$ribo_start
  ribo_strand <- rep$ribo_strand

  if (ribo_strand == 1) {
    genes <- genes[order(sapply(genes, function(g) g$start))]
  } else {
    genes <- genes[order(sapply(genes, function(g) g$start), decreasing = TRUE)]
  }

  data.frame(
    genome = gsub(".gbff", "", rep$genome),
    organism = rep$organism,
    operon = rep$record,
    product = sapply(genes, function(g) g$product),
    start = sapply(genes, function(g) g$start),
    end = sapply(genes, function(g) g$end),
    strand = sapply(genes, function(g) g$strand),
    ribo_start = ribo_start,
    ribo_strand = ribo_strand
  )
})

df <- bind_rows(df_list)

# Normalise product strings
df$product <- trimws(df$product)
df$product <- gsub("\\s+", " ", df$product)
df$product_key <- tolower(df$product)

# Product collapsing map
product_map <- c(
  "atp-binding cassette domain-containing protein" = "ABC ATP-binding",
  "abc transporter atp-binding protein" = "ABC ATP-binding",
  
  "abc transporter permease" = "ABC permease",
  "apc family permease" = "APC family permease",
  
  "aliphatic sulfonate abc transporter substrate-binding protein" = "ABC substrate-binding",
  "putative urea abc transporter substrate-binding protein" = "ABC substrate-binding",
  
  "agmatinase" = "Putative gdmH",
  "agmatinase family protein" = "Putative gdmH",
  
  "hydrogenase maturation nickel metallochaperone hypa" = "HypA",
  "hydrogenase nickel incorporation protein hypb" = "HypB",
  
  "asma family protein" = "AsmA family protein",
  "asmA family protein" = "AsmA family protein",
  "dmt family transporter" = "DMT family transporter",
  
  "hypothetical protein" = "Hypothetical protein",
  "pas domain s-box protein" = "PAS domain S-box protein"
)

# Apply product short-name mapping
df$product <- recode(df$product_key, !!!product_map, .default = df$product_key)

# Recode display labels
df$product <- recode(df$product,
  "hypothetical protein" = "Hypothetical protein",
  "pas domain s-box protein" = "PAS domain S-box protein",
  .default = df$product
)

df$product <- factor(
  df$product,
  levels = c(
    "Putative gdmH",
    "ABC ATP-binding",
    "ABC permease",
    "ABC substrate-binding",
    "Hypothetical protein",
    "HypA",
    "HypB",
    "APC family permease",
    "AsmA family protein",
    "PAS domain S-box protein",
    "DMT family transporter"
    
  )
)

# ---------------------------
# Fix operon order
# ---------------------------

desired_organism_order <- c(
  "Bradyrhizobium ontarionense",
  "Roseibium aggregatum IAM 12614",
  "Sagittula stellata E-37",
  "Tateyamaria pelophila",
  "Synechocystis sp. PCC 6803",
  "Candidatus Nitrospira inopinata",
  "Oceanobacillus limi",
  "Streptomyces oceani"
)

genome_order <- df %>%
  distinct(genome, organism) %>%
  mutate(
    organism = factor(
      organism,
      levels = rev(desired_organism_order)
    )
  ) %>%
  arrange(organism) %>%
  pull(genome)

df$genome <- factor(df$genome, levels = genome_order)

# ---------------------------
# Position calculator
# ---------------------------
df <- df %>%
  group_by(genome, operon) %>%
  mutate(
    xmin0 = pmin(start, end) - ribo_start,
    xmax0 = pmax(start, end) - ribo_start,
    xmin = ifelse(ribo_strand == 1, xmin0, -xmax0),
    xmax = ifelse(ribo_strand == 1, xmax0, -xmin0)
  ) %>%
  ungroup() %>%
  mutate(
    xmin = pmin(xmin, xmax),
    xmax = pmax(xmin, xmax)
  )

line_df <- df %>%
  group_by(genome, operon) %>%
  summarise(
    xmin = min(xmin, na.rm = TRUE),
    xmax = max(xmax, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------
# Colours
# ---------------------------
col_map <- c(
  # ABC and APC transporter system
  "ABC ATP-binding"      = "#08519C",
  "ABC permease"         = "#3182BD",
  "ABC substrate-binding"= "#9ECAE1",
  "APC family permease"  = "#E3E1CA",

  # gdmH
  "Putative gdmH"        = "#E69F00",

  # Hydrogenase maturation system
  "HypA"                 = "#238B45",
  "HypB"                 = "#74C476",

  # Other proteins
  "AsmA family protein"  = "#984EA3",
  "DMT family transporter" = "#4D4D4D",

  "Hypothetical protein" = "#BDBDBD",
  "PAS domain S-box protein" = "#E41A1C"
)

missing_cols <- setdiff(unique(df$product), names(col_map))
if (length(missing_cols) > 0) {
  warning("Unmapped products: ", paste(missing_cols, collapse = ", "))
  col_map <- c(col_map, setNames(rep("#BDBDBD", length(missing_cols)), missing_cols))
}

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(df, aes(
  xmin = xmin,
  xmax = xmax,
  y = genome,
  fill = product
)) +

  geom_segment(
    data = line_df,
    aes(x = xmin, xend = xmax, y = genome, yend = genome),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.3
  ) +

  geom_gene_arrow() +

  geom_vline(xintercept = 0, linetype = "dashed") +

  labs(
    x = "Distance from riboswitch (bp)",
    y = NULL
  ) +

  scale_y_discrete(
    labels = function(x) {
      org_map <- setNames(df$organism, df$genome)
      org_map[x]
    }
  ) +

  scale_fill_manual(values = col_map, drop = FALSE, na.value = "#BDBDBD") +

  theme_bw() +

  theme(
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    legend.title = element_blank()
  )

ggsave(opt$output, p, width = 10, height = 6, device = svglite::svglite)