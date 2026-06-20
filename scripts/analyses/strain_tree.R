#!/usr/bin/env Rscript

# Libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)  
  library(dplyr)
  library(Biostrings)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(phangorn)
  library(farver)
})

# Command-line arguments
option_list <- list(
  make_option(c("-t", "--tree"), type = "character"),
  make_option(c("-a", "--alignment"), type = "character"),
  make_option(c("-s", "--added_seqs"), type = "character"),
  make_option(c("-n", "--nif"), type = "character"),
  make_option(c("--my_gene_presence"), type = "character", default = NULL),
  make_option(c("--hug_nif"), type = "character", default = NULL),
  make_option(c("--hug_gdmh"), type = "character", default = NULL),
  make_option(c("--taxonomy"), type = "character"),
  make_option(c("-o", "--out"), type = "character"),
  
  make_option(c("-l", "--tax_level"),
              type = "character",
              default = "family"),

  make_option(c("-v", "--overlap_only"),
              action = "store_true",
              default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))

tree_file             <- opt$tree
aln_file              <- opt$alignment
added_seqs_file       <- opt$added_seqs
nif_file              <- opt$nif
my_gene_presence_file <- opt$my_gene_presence
hug_nif_file          <- opt$hug_nif
hug_gdmh_file         <- opt$hug_gdmh
taxonomy_file         <- opt$taxonomy
out_file              <- opt$out
tax_level             <- opt$tax_level
overlap_only          <- opt$overlap_only


# Read tree


tree <- read.tree(tree_file)
tree <- midpoint(tree)


# Remove problematic taxa


bad_taxa <- c(
  "Bacteria_Pseudomonadota_Alphaproteobacteria_Rhodobacterales_Paracoccaceae_Rhodovulum_Rhodovulum_sp._PH10",
  "Bacteria_Pseudomonadota_Alphaproteobacteria_Hyphomicrobiales_Reyranellaceae_Reyranella_Reyranella_massiliensis_521"
)

tree <- drop.tip(tree, bad_taxa)

seqs <- readAAStringSet(aln_file)
headers <- names(seqs)


# Build metadata from taxonomy table


meta <- read.delim(
  taxonomy_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(meta) <- c(
  "Genome",
  "domain",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

meta$label <- paste(
  meta$domain,
  meta$phylum,
  meta$class,
  meta$order,
  meta$family,
  meta$genus,
  gsub("[ /]", "_", meta$species),
  sep = "_"
)

meta <- meta[, c(
  "label",
  "Genome",
  "domain",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)]

# Tree filtering

tree <- drop.tip(tree, setdiff(tree$tip.label, meta$label))
meta <- meta %>% dplyr::filter(label %in% tree$tip.label)

# Added sequence annotation

fasta_lines <- readLines(added_seqs_file)

added_headers <- sub(
  "^>",
  "",
  fasta_lines[grepl("^>", fasta_lines)]
)

meta <- meta %>%
  dplyr::mutate(added = label %in% added_headers)


# NIF data loading

my_nif <- read.delim(nif_file, stringsAsFactors = FALSE)

if (!is.null(hug_nif_file) && file.exists(hug_nif_file)) {

  hug_nif <- read.delim(hug_nif_file, stringsAsFactors = FALSE)

  nif <- dplyr::bind_rows(my_nif, hug_nif)

} else {
  nif <- my_nif
}

meta$Genome <- trimws(as.character(meta$Genome))
nif$Genome <- trimws(as.character(nif$Genome))

stopifnot("Genome" %in% colnames(meta))
stopifnot("Genome" %in% colnames(nif))


# NIF gene processing

nif_genes <- c(
  "NIFB_presence",
  "NIFD_presence",
  "NIFE_presence",
  "NIFH_presence",
  "NIFK_presence",
  "NIFN_presence"
)

nif_genes <- nif_genes[nif_genes %in% colnames(nif)]

nif[nif_genes] <- lapply(nif[nif_genes], as.numeric)


# Expand multi genome entries in NIF table

nif <- nif %>%
  tidyr::separate_rows(Genome, sep = ";") %>%
  dplyr::mutate(Genome = trimws(Genome))

# Collapse duplicate genome entries in NIF table

nif <- nif %>%
  dplyr::group_by(Genome) %>%
  dplyr::summarise(
    dplyr::across(all_of(nif_genes), \(x) max(x, na.rm = TRUE)),
    .groups = "drop"
  )


# Join NIF data

meta <- meta %>%
  dplyr::left_join(
    nif %>% dplyr::select(Genome, all_of(nif_genes)),
    by = "Genome"
  )

meta$nif_any <- rowSums(meta[, nif_genes, drop = FALSE], na.rm = TRUE) > 0


# GDMH processing

my_gdmh_taxa <- unique(my_nif$Genome)

meta$gdmh_hug_present <- "absent"

if (!is.null(hug_gdmh_file) && file.exists(hug_gdmh_file)) {

  gdmh_hug <- read.delim(hug_gdmh_file, stringsAsFactors = FALSE)

  gdmh_hug$GDMH_presence <- as.numeric(gdmh_hug$GDMH_presence)

  gdmh_hug <- gdmh_hug %>%
    tidyr::separate_rows(Genome, sep = ";") %>%
    dplyr::mutate(Genome = trimws(Genome))

  gdmh_hug <- gdmh_hug %>%
    dplyr::group_by(Genome) %>%
    dplyr::summarise(
      GDMH_presence = max(GDMH_presence, na.rm = TRUE),
      .groups = "drop"
    )

  meta <- meta %>%
    dplyr::left_join(
      gdmh_hug %>% dplyr::select(Genome, GDMH_presence),
      by = "Genome"
    )

  meta$gdmh_hug_present <- ifelse(
    !is.na(meta$GDMH_presence) & meta$GDMH_presence == 1,
    "present",
    "absent"
  )
}

meta$gdmh_present <- ifelse(
  meta$Genome %in% my_gdmh_taxa |
    meta$gdmh_hug_present == "present",
  "present",
  "absent"
)

meta$gdmh_present <- factor(meta$gdmh_present, levels = c("absent", "present"))


# Read gene presence/absence of strains from the lab

meta$manual_override <- FALSE

if (!is.null(my_gene_presence_file) &&
    file.exists(my_gene_presence_file)) {

  my_gene_presence <- read.delim(
    my_gene_presence_file,
    stringsAsFactors = FALSE
  )

  my_gene_presence$Species <- trimws(my_gene_presence$Species)

  meta$species <- trimws(meta$species)

  nif_cols <- intersect(
    nif_genes,
    colnames(my_gene_presence)
  )

  idx <- match(meta$species, my_gene_presence$Species)

  matched <- !is.na(idx)

  # mark strains whose annotations are overridden
  meta$manual_override[matched] <- TRUE

  for (col in nif_cols) {

    meta[[col]][matched] <-
      my_gene_presence[[col]][idx[matched]]

  }

  meta$nif_any <-
    rowSums(meta[, nif_genes, drop = FALSE], na.rm = TRUE) > 0

  if ("GDMH_presence" %in% colnames(my_gene_presence)) {

    meta$gdmh_present[matched] <- ifelse(
      my_gene_presence$GDMH_presence[idx[matched]] == 1,
      "present",
      "absent"
    )

    meta$gdmh_present <- factor(
      meta$gdmh_present,
      levels = c("absent", "present")
    )
  }
}

meta$plot_shape <- dplyr::case_when(
  meta$manual_override ~ "override",
  meta$gdmh_present == "present" ~ "present",
  TRUE ~ "absent"
)

meta$plot_shape <- factor(
  meta$plot_shape,
  levels = c("absent", "present", "override")
)


# NIF long format for plotting

nif_long <- meta %>%
  dplyr::select(label, Genome, all_of(nif_genes)) %>%
  tidyr::pivot_longer(
    cols = all_of(nif_genes),
    names_to = "gene",
    values_to = "present"
  ) %>%
  dplyr::mutate(
    present = ifelse(present == 1, "present", "absent")
  )


# Select taxa to highlight

if (overlap_only) {

  highlight_taxa <- meta %>%
    dplyr::filter(gdmh_present == "present" | nif_any == TRUE) %>%
    dplyr::pull(.data[[tax_level]]) %>%
    unique()

} else {

  top_n <- 20

  highlight_taxa <- meta %>%
    dplyr::count(.data[[tax_level]], sort = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(1)
}

# Create plotting column
meta <- meta %>%
  mutate(tax_plot = ifelse(.data[[tax_level]] %in% highlight_taxa,
                           .data[[tax_level]],
                           "Other"))


# Colours

is_grey <- function(hex) {
  rgb <- col2rgb(hex)
  r <- rgb[1]; g <- rgb[2]; b <- rgb[3]
  max_range <- max(r, g, b) - min(r, g, b)   # low = grey-ish
  mean_val  <- mean(c(r, g, b))               # high = light
  max_range < 30 & mean_val > 140             # tunable thresholds
}

n_cols <- length(unique(meta$tax_plot))

if (tax_level == "order") {

  order_palette <- c(
    "#7a8c8c",
    "#d08a8a",
    "#66a366",
    "#009cc0",
    "#69b2af",
    "#FB3C3C",
    "#ffcc80",
    "#9191da",
    "#ff7a7a",
    "#ffffc6",
    "#b8b8ff",
    "#ff5cf7",
    "#eee3ee",
    "#ff99ff",
    "#8ec5ff",
    "#fff34f"
  )

  n_cols <- length(unique(meta$tax_plot))
  
  all_colors <- setNames(
    order_palette[seq_len(n_cols)],
    sort(unique(meta$tax_plot))
  )

} else {

  n_cols <- length(unique(meta$tax_plot))

  all_colors <- hcl.colors(n_cols, palette = "Set 2")

}

palette <- colorRampPalette(order_palette)(
  length(highlight_taxa)
)

tax_colors <- setNames(palette, highlight_taxa)
tax_colors <- c(tax_colors, Other = "grey80")

meta$gdmh_present <- factor(meta$gdmh_present,
                            levels = c("absent", "present"))

# Node support

tree$sh_alrt <- as.numeric(sapply(strsplit(tree$node.label, "/"), `[`, 1))
tree$bootstrap <- as.numeric(sapply(strsplit(tree$node.label, "/"), `[`, 2))

# Keep bootstrap in node.label for ggtree compatibility
tree$node.label <- tree$bootstrap

# Plot

p <- suppressWarnings(
  suppressMessages(
    ggtree(tree, layout="circular", size=0.2) %<+% meta +
      xlim_tree(1.8) +
      
      # Branch colouring by bootstrap support
      geom_tree(aes(color = as.numeric(label))) +
      
      scale_color_gradient2(
        low = "gainsboro",
        mid = "darkgrey",
        high = "black",
        midpoint = 70,
        name = "Bootstrap support\n"
      ) +
      
      # Tip symbols with taxonomy fill
      geom_tippoint(
        aes(
          shape = plot_shape,
          fill = tax_plot,
          alpha = gdmh_present
        ),
        size = 3,
        stroke = 0.4,
        color = "black"
      ) +
      
      scale_shape_manual(
        name = "Tip symbol",
        values = c(
          "absent"   = 21,
          "present"  = 24,
          "override" = 23
        ),
        labels = c(
          "absent"   = expression(italic(gdmH)^"-"),
          "present"  = expression(italic(gdmH)^"+"),
          "override" = "Selected strain"
        ),
        guide = guide_legend(
          override.aes = list(
            fill = "white",
            size = 5,
            alpha = 1,
            color = "black"
          )
        )
      ) +
      
      scale_alpha_manual(
        values = c(
          "absent" = 0.6,
          "present" = 1
        ),
        guide = "none"
      ) +
      
      scale_fill_manual(
        name = paste("Taxonomic", tax_level),
        values = tax_colors,
        breaks = highlight_taxa,
        guide = guide_legend(
          override.aes = list(shape = 21, size = 6, color = "black", alpha = 1)
        )
      ) +

      # NIF presence/absence
      ggnewscale::new_scale_fill() +
      
      geom_fruit(
        data = nif_long %>% filter(present == "present"),
        geom = geom_point,
        mapping = aes(y = label, x = gene, fill = gene),
        shape = 21,
        size = 1.5,
        stroke = 0.3,
        color = "black",
        offset = 0.02,
        pwidth = 0.10
      ) +
      
      scale_fill_manual(
        values = setNames(
          RColorBrewer::brewer.pal(6, "Dark2"),
          c("NIFB_presence","NIFD_presence","NIFE_presence",
            "NIFH_presence","NIFK_presence","NIFN_presence")
        ),
        name = expression(italic(nif)~" presence"),
        labels = c(
          expression(italic(nifB)),
          expression(italic(nifD)),
          expression(italic(nifE)),
          expression(italic(nifH)),
          expression(italic(nifK)),
          expression(italic(nifN))
        ),
        guide = guide_legend(
          override.aes = list(
            shape = 21,
            size = 4,
            stroke = 1,
            color = "black"
          )
        )
      ) +
      
      # Tip labels
      ggnewscale::new_scale_fill() +
      
      geom_tiplab(
        aes(label = species, fill = tax_plot),
        geom = "label",
        color = "black",
        label.padding = unit(0.2, "lines"),
        label.r = unit(0.15, "lines"),
        size = 2,
        linewidth = 0.3,
        align = TRUE,
        offset = 0.20,
        show.legend = FALSE
      ) +
      
      scale_fill_manual(
        values = tax_colors,
        guide = "none"
      ) +
      
      theme(
        legend.position = "right",
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
      )
  )
)



# Save
ext <- tools::file_ext(out_file)

if (ext == "pdf") {
  ggsave(
    out_file,
    plot = p,
    width = 16,
    height = 16,
    device = cairo_pdf
  )
} else if (ext == "svg") {
  ggsave(
    out_file,
    plot = p,
    width = 16,
    height = 16,
    device = svglite::svglite
  )
} else {
  stop("Unsupported file type: use .pdf or .svg")
}

cat(paste("\nTree saved as"), out_file, "\n")