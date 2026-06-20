#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(ggtree)
  library(ggrepel)
  library(dplyr)
  library(readr)
})

# Command-line arguments

option_list <- list(
  make_option(c("-t", "--tree"), type = "character", help = "Tree file in Newick format"),
  make_option(c("-m", "--mapping"), type = "character", help = "Mapping file (tab-delimited)"),
  make_option(c("-l", "--level"), type = "character", help = "Taxonomic level to highlight"),
  make_option(c("-e", "--highlight_euks"),
              action = "store_true",
              default = FALSE,
              help = "Highlight Eukaryota clades"),
  make_option(c("-o", "--out"), type = "character", help = "Output file (.pdf or .svg)")
)

opt <- parse_args(OptionParser(option_list = option_list))

tree_file <- opt$tree
mapping_file <- opt$mapping
tax_level <- opt$level
out_file <- opt$out
highlight_euks <- opt$highlight_euks

# Load tree + midpoint root

tree <- read.tree(tree_file)
tree <- midpoint(tree)

# Load mapping file

mapping <- read_tsv(mapping_file, show_col_types = FALSE)

tax_levels <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
if (!tax_level %in% tax_levels) {
  stop("Invalid taxonomic level. Choose one of: ", paste(tax_levels, collapse = ", "))
}

# Split Taxonomic lineage and extract ranks

mapping <- mapping %>%
  mutate(TaxSplit = strsplit(`Taxonomic lineage`, ",\\s*"))

# helper function: safely extract first match
get_rank <- function(x, pattern) {
  hit <- x[grepl(pattern, x, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  return(hit[1])
}

get_clade <- function(x, n) {
  clades <- x[grepl("\\(clade\\)", x, ignore.case = TRUE)]
  if (length(clades) >= n) {
    return(gsub(" \\(.*\\)$", "", clades[n]))
  }
  return(NA_character_)
}

mapping$superkingdom <- sapply(mapping$TaxSplit, function(x) {
  val <- get_rank(x, "domain")
  if (!is.na(val)) return(val)
  get_clade(x, 1)
})

mapping$kingdom <- sapply(mapping$TaxSplit, function(x) {
  val <- get_rank(x, "kingdom")
  if (!is.na(val)) return(val)
  get_clade(x, 1)
})

mapping$phylum <- sapply(mapping$TaxSplit, function(x) {
  val <- get_rank(x, "phylum")
  if (!is.na(val)) return(val)
  get_clade(x, 2)
})

mapping$class <- sapply(mapping$TaxSplit, function(x) {
  val <- get_rank(x, "class")
  if (!is.na(val)) return(val)
  get_clade(x, 3)
})

mapping$order <- sapply(mapping$TaxSplit, function(x) {
  val <- get_rank(x, "order")
  if (!is.na(val)) return(val)
  get_clade(x, 3)
})

mapping$family <- sapply(mapping$TaxSplit, function(x) {
  get_rank(x, "family")
})

mapping$genus <- sapply(mapping$TaxSplit, function(x) {
  get_rank(x, "genus")
})

mapping$species <- sapply(mapping$TaxSplit, function(x) {
  get_rank(x, "species")
})

# Clean names

mapping <- mapping %>%
  mutate(across(
    c(superkingdom, kingdom, phylum, class, order, family, genus, species),
    ~ gsub(" \\(.*\\)$", "", .)
  ))

# Prepare tip labels with taxonomic info

tip_info <- mapping %>%
  select(
    From,
    Organism,
    superkingdom,
    all_of(tax_level)
  ) %>%
  rename(
    tip = From,
    organism = Organism,
    tax = !!tax_level
  )

matched_tips <- tip_info$tip %in% tree$tip.label
cat("Matched tips:", sum(matched_tips), "out of", length(tree$tip.label), "\n")

na_tips <- tip_info$tip[is.na(tip_info)]
cat("Number of tips with missing", tax_level, ":", length(na_tips), "\n")

tip_info <- tip_info %>%
  filter(!is.na(tax))

# Drop NA tips/branches
tree2 <- ape::drop.tip(tree, setdiff(tree$tip.label, tip_info$tip))

# Extract bootstrap / SH-aLRT from node labels

tree2$sh_alrt <- as.numeric(sapply(strsplit(tree2$node.label, "/"), `[`, 1))
tree2$bootstrap <- as.numeric(sapply(strsplit(tree2$node.label, "/"), `[`, 2))

# Colours

tax_palette <- c(
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

# Organism label selection

label_tips <- tip_info %>%
  filter(
    grepl("Roseibium marinum", organism) |
    grepl("Roseibium aggregatum", organism) |
    grepl("Mameliella alba", organism) |
    grepl("Phaeodactylum tricornutum", organism) |
    grepl("Synechocystis sp.", organism) |
    grepl("Sinorhizobium fredii", organism) |
    grepl("Oceanobacillus", organism) |
    grepl("Candidatus Nitrospira inopinata", organism) |
    grepl("Streptomyces", organism)
  ) %>%
  pull(tip)

# Bootstrap alignment

boot_df <- data.frame(
  node = (Ntip(tree2) + 1):(Ntip(tree2) + Nnode(tree2)),
  bootstrap = tree2$bootstrap
)

# Merge with tree data and plot

p <- suppressWarnings(
  suppressMessages({

    p <- ggtree(tree2, layout = "equal_angle", size = 0.3) %<+% tip_info
    
    p$data <- dplyr::left_join(p$data, boot_df, by = "node")
    
    p <- p +
      geom_nodepoint(
        aes(color = bootstrap, subset = bootstrap >= 70),
        size = 0.8,
        na.rm = TRUE
      ) +
      scale_color_viridis_c(
        limits = c(0, 100),
        name = "Bootstrap"
      )

    if (highlight_euks) {

      euk_tips <- tip_info %>%
        filter(superkingdom == "Eukaryota") %>%
        pull(tip)

      euk_tips <- intersect(euk_tips, tree2$tip.label)

      cat("Eukaryotic tips in tree:", length(euk_tips), "\n")

      if (length(euk_tips) > 1) {

        euk_node <- ape::getMRCA(tree2, euk_tips)

        cat("Eukaryota MRCA node:", euk_node, "\n")

        p <- p +
          geom_hilight(
            node = euk_node,
            fill = "#ffd700",
            alpha = 0.20
          )
      }
    }

    label_data <- p$data %>%
      filter(label %in% label_tips) %>%
      mutate(label = organism)

    p <- p +
      geom_tippoint(
        aes(fill = tax),
        shape = 21,
        color = "black",
        size = 2
      ) +
      geom_tippoint(
        data = label_data,
        aes(fill = tax),
        shape = 24,
        color = "black",
        size = 2,
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = setNames(
          rep(
            tax_palette,
            length.out = length(unique(na.omit(tip_info$tax)))
          ),
          sort(unique(na.omit(tip_info$tax)))
        ),
        name = paste0("Taxonomic ", tax_level)
      ) +
      theme(legend.position = "right")

    p <- p +
      ggrepel::geom_text_repel(
        data = label_data,
        aes(
          x = x,
          y = y,
          label = label
        ),
        size = 2.2,
        force = 2.5,
        force_pull = 0.5,
        max.overlaps = Inf,
        box.padding = 0.3,
        point.padding = 0,
        min.segment.length = 0,
        segment.color = "grey40",
        segment.alpha = 0.8,
        segment.size = 0.4,
        direction = "both",
        seed = 1
      )

    p

  })
)

# Save

ext <- tools::file_ext(out_file)


if (ext == "pdf") {
    ggsave(out_file, p, width = 8, height = 8, device = cairo_pdf)
} else if (ext == "svg") {
    ggsave(out_file, p, width = 8, height = 8, device = svglite::svglite)
} else {
    stop("Unsupported file type: use .pdf or .svg")
}

cat("\nTaxonomy gene tree saved as", out_file, "\n")

# Protein-name plotting data

mapping <- mapping %>%
  mutate(
    protein_group = case_when(

      # SpeB
      grepl("^SpeB", `Protein names`, ignore.case = TRUE) ~ "SpeB",

      # Proclavaminate amidinohydrolase
      grepl("Proclavaminate amidinohydrolase",
            `Protein names`,
            ignore.case = TRUE) ~ "Proclavaminate amidinohydrolase",

      # Arginase/agmatinase/formiminoglutamase family
      grepl(
        "Arginase/agmatinase/formim",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Arginase/agmatinase/formiminoglutamase",

      grepl(
        "Arginase family hydrolase",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Arginase/agmatinase/formiminoglutamase",

      # Agmatinase
      grepl(
        "^Agmatinase",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Agmatinase",

      grepl(
        "Probable agmatinase",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Agmatinase",

      grepl(
        "Agmatinase family protein",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Agmatinase",

      # Agmatine ureohydrolase
      grepl(
        "Agmatine ureohydrolase",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Agmatine ureohydrolase",

      # Arginase family
      grepl(
        "^Arginase family",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Arginase family",

      grepl(
        "Arginase/agmatinase family protein",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Arginase family",
      
      # Shorten GdmH name     
      grepl(
        "^Guanidine hydrolase",
        `Protein names`,
        ignore.case = TRUE
      ) ~ "Guanidine hydrolase (EC 3.5.3.-)",

      # everything else
      TRUE ~ `Protein names`
    )
  )
  
protein_info <- mapping %>%
  select(
    From,
    Organism,
    protein_group
  ) %>%
  rename(
    tip = From,
    organism = Organism,
    protein = protein_group
  ) %>%
  filter(!is.na(protein))

protein_info <- protein_info %>%
  filter(tip %in% tree2$tip.label)

# Protein-name coloured tree

p_protein <- ggtree(tree2, layout = "equal_angle", size = 0.3) %<+% protein_info

p_protein$data <- dplyr::left_join(
  p_protein$data,
  boot_df,
  by = "node"
)

p_protein <- p_protein +
  geom_nodepoint(
    aes(color = bootstrap, subset = bootstrap >= 70),
    size = 0.8,
    na.rm = TRUE
  ) +
  scale_color_viridis_c(
    limits = c(0, 100),
    name = "Bootstrap support"
  )
   
protein_levels <- sort(unique(protein_info$protein))

protein_cols <- setNames(
  rep(
    tax_palette,
    length.out = length(protein_levels)
  ),
  protein_levels
)

label_data_protein <- p_protein$data %>%
  filter(label %in% label_tips) %>%
  mutate(label = organism)

p_protein <- p_protein +
  geom_tippoint(
    aes(fill = protein),
    shape = 21,
    colour = "black",
    size = 2
  ) +
  geom_tippoint(
    data = label_data_protein,
    aes(fill = protein),
    shape = 24,
    colour = "black",
    size = 2,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = protein_cols,
    name = "Protein family/name"
  ) +
  theme(
    legend.position = "right"
  )
  
protein_out <- sub("\\.(pdf|svg)$", "_proteins.\\1", out_file)

if (ext == "pdf") {
  ggsave(protein_out, p_protein, width = 8, height = 8, device = cairo_pdf)
} else if (ext == "svg") {
  ggsave(protein_out, p_protein, width = 8, height = 8, device = svglite::svglite)
} else {
  stop("Unsupported file type: use .pdf or .svg")
}

cat("\nProtein family/name tree saved as", protein_out, "\n")