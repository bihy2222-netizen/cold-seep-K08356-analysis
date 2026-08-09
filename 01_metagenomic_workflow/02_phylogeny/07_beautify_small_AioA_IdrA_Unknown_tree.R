#!/usr/bin/env Rscript

# Beautify the focused AioA/IdrA/Unknown K08356 small tree.
#
# Default input matches the strict group-derep small tree on the cold_seep
# server. The script always writes iTOL annotation files and, when ggtree is
# installed, also writes publication-ready PDF/PNG figures.

default_tree <- paste0(
  "/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/",
  "K08356_four_group_annotation/strict_kofam_best/",
  "tree_small_AioA_IdrA_Unknown_group_derep49_refs136/",
  "AioA_IdrA_Unknown_group_derep49_plus_refs136.trimmed.fasta.treefile"
)

default_target <- "NS|R2111_N500_0-10_bin13-k141_4404301_3"
default_sister <- "IS|SY368YW-8-12_bin16-k141_314641_4"

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript 07_beautify_small_AioA_IdrA_Unknown_tree.R [options]\n\n",
    "Options:\n",
    "  --tree FILE        IQ-TREE .treefile or .contree [server small tree]\n",
    "  --outdir DIR       Output directory [dirname(tree)/beautified_small_tree]\n",
    "  --layout NAME      rectangular or circular [rectangular]\n",
    "  --label-mode MODE  key, mag, all, or none [key]\n",
    "  --target LABEL     Tip label to highlight as the new NS sequence\n",
    "  --sister LABEL     Sister/near-identical tip label to highlight\n",
    "  --help             Show this help\n\n",
    "Outputs:\n",
    "  tip_metadata.tsv\n",
    "  small_tree_beautified.pdf/png, if ggtree is installed\n",
    "  itol_habitat_colorstrip.txt\n",
    "  itol_family_colorstrip.txt\n",
    "  itol_key_symbols.txt\n",
    "  itol_display_labels.txt\n",
    sep = ""
  )
}

parse_args <- function(args) {
  opts <- list(
    tree = default_tree,
    outdir = NA_character_,
    layout = "rectangular",
    label_mode = "key",
    target = default_target,
    sister = default_sister
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--help") {
      usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    if (!name %in% names(opts)) {
      stop("Unknown option: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for option: ", key)
    }
    opts[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  if (is.na(opts$outdir) || opts$outdir == "") {
    opts$outdir <- file.path(dirname(opts$tree), "beautified_small_tree")
  }
  opts
}

need_pkg <- function(pkg, required = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && required) {
    stop(
      "Required R package not installed: ", pkg, "\n",
      "Install it first, or run on an environment that already has it."
    )
  }
  ok
}

read_tip_labels_base <- function(tree_file) {
  newick <- paste(readLines(tree_file, warn = FALSE), collapse = "")
  hits <- gregexpr("[\\(,]([^\\(\\),:;]+):", newick, perl = TRUE)[[1]]
  if (length(hits) == 1 && hits[[1]] == -1) {
    stop("Could not parse tip labels from Newick file: ", tree_file)
  }
  raw <- regmatches(newick, list(hits))[[1]]
  labels <- sub("^[\\(,]", "", raw)
  labels <- sub(":$", "", labels)
  labels <- gsub("^'|'$", "", labels)
  unique(labels)
}

first_match <- function(x, pattern, replacement) {
  out <- rep(NA_character_, length(x))
  hit <- grepl(pattern, x, perl = TRUE)
  out[hit] <- sub(pattern, replacement, x[hit], perl = TRUE)
  out
}

classify_family <- function(label, is_mag) {
  x <- tolower(label)
  out <- rep("Other reference", length(label))
  out[grepl("unknown|unk", x)] <- "Unknown reference"
  out[grepl("idra|iria|idr", x)] <- "IdrA/IriA reference"
  out[grepl("aioa", x)] <- "AioA reference"
  out[grepl("arxa", x)] <- "ArxA reference"
  out[grepl("arra", x)] <- "ArrA reference"
  out[grepl("aor|for|gapor|dmsor", x)] <- "Other DMSOR reference"
  out[is_mag] <- "Group derep MAG K08356"
  out
}

short_label <- function(label) {
  x <- sub("^(IS|AS|ES|NS)\\|", "", label, perl = TRUE)
  x <- sub("_bin", " bin", x, fixed = TRUE)
  x
}

write_itol_colorstrip <- function(path, dataset_label, meta, value_col, colors) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(
    "DATASET_COLORSTRIP\n",
    "SEPARATOR TAB\n",
    "DATASET_LABEL\t", dataset_label, "\n",
    "COLOR\t#333333\n",
    "STRIP_WIDTH\t25\n",
    "MARGIN\t5\n",
    "BORDER_WIDTH\t1\n",
    "BORDER_COLOR\t#FFFFFF\n",
    "LEGEND_TITLE\t", dataset_label, "\n",
    "LEGEND_SHAPES\t", paste(rep(1, length(colors)), collapse = "\t"), "\n",
    "LEGEND_COLORS\t", paste(unname(colors), collapse = "\t"), "\n",
    "LEGEND_LABELS\t", paste(names(colors), collapse = "\t"), "\n",
    "DATA\n",
    sep = "",
    file = con
  )
  for (i in seq_len(nrow(meta))) {
    value <- meta[[value_col]][i]
    if (!is.na(value) && value %in% names(colors)) {
      cat(meta$label[i], colors[[value]], value, sep = "\t", file = con)
      cat("\n", file = con)
    }
  }
}

write_itol_symbols <- function(path, meta, target, sister) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(
    "DATASET_SYMBOL\n",
    "SEPARATOR TAB\n",
    "DATASET_LABEL\tKey sequences\n",
    "COLOR\t#222222\n",
    "LEGEND_TITLE\tKey sequences\n",
    "LEGEND_SHAPES\t2\t2\n",
    "LEGEND_COLORS\t#D62728\t#FF9F1C\n",
    "LEGEND_LABELS\tNew NS K08356\tSister IS K08356\n",
    "DATA\n",
    sep = "",
    file = con
  )
  if (target %in% meta$label) {
    cat(target, "2", "18", "#D62728", "1", "1", "New_NS_K08356", sep = "\t", file = con)
    cat("\n", file = con)
  }
  if (sister %in% meta$label) {
    cat(sister, "2", "16", "#FF9F1C", "1", "1", "Sister_IS_K08356", sep = "\t", file = con)
    cat("\n", file = con)
  }
}

write_itol_labels <- function(path, meta) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  cat("LABELS\nSEPARATOR TAB\nDATA\n", file = con)
  for (i in seq_len(nrow(meta))) {
    cat(meta$label[i], meta$display_label[i], sep = "\t", file = con)
    cat("\n", file = con)
  }
}

plot_with_ggtree <- function(tree, meta, outdir, layout, label_mode) {
  has_plot <- all(vapply(c("ggplot2", "ggtree"), need_pkg, logical(1), required = FALSE))
  if (!has_plot) {
    message("ggtree/ggplot2 not installed; skipped PDF/PNG tree plotting.")
    return(invisible(FALSE))
  }

  library(ggplot2)
  library(ggtree)

  habitat_colors <- c(
    IS = "#1B9E77",
    AS = "#D95F02",
    ES = "#7570B3",
    NS = "#E7298A",
    Reference = "#9E9E9E"
  )
  family_colors <- c(
    "Group derep MAG K08356" = "#2A9D8F",
    "AioA reference" = "#1F78B4",
    "IdrA/IriA reference" = "#C65D2E",
    "Unknown reference" = "#6A4C93",
    "ArxA reference" = "#7F7F7F",
    "ArrA reference" = "#8C564B",
    "Other DMSOR reference" = "#4D4D4D",
    "Other reference" = "#BDBDBD"
  )

  p <- ggtree(tree, layout = layout, size = 0.28, color = "#404040") %<+% meta
  p <- p +
    geom_tippoint(aes(color = family, shape = habitat), size = 1.8, alpha = 0.95) +
    scale_color_manual(values = family_colors, drop = FALSE) +
    scale_shape_manual(values = c(IS = 16, AS = 17, ES = 15, NS = 18, Reference = 1), drop = FALSE) +
    theme_tree2() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(8, 8, 8, 8)
    )

  key_labels <- meta$label[meta$is_key]
  label_tips <- switch(
    label_mode,
    none = character(0),
    key = key_labels,
    mag = meta$label[meta$is_mag | meta$is_key],
    all = meta$label,
    stop("label-mode must be key, mag, all, or none")
  )

  if (length(label_tips) > 0) {
    label_meta <- meta[meta$label %in% label_tips, , drop = FALSE]
    p <- p +
      geom_tiplab(
        data = label_meta,
        aes(label = display_label),
        size = if (label_mode == "all") 1.5 else 2.1,
        align = FALSE,
        linesize = 0.15,
        offset = 0.004
      )
  }

  pdf_file <- file.path(outdir, "small_tree_beautified.pdf")
  png_file <- file.path(outdir, "small_tree_beautified.png")
  width <- if (layout == "circular") 9 else 12
  height <- if (layout == "circular") 9 else 14
  ggsave(pdf_file, p, width = width, height = height, units = "in", limitsize = FALSE)
  ggsave(png_file, p, width = width, height = height, units = "in", dpi = 400, limitsize = FALSE)
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
  invisible(TRUE)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!file.exists(opts$tree)) {
    stop("Tree file not found: ", opts$tree)
  }
  if (!opts$layout %in% c("rectangular", "circular")) {
    stop("--layout must be rectangular or circular")
  }

  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
  labels <- read_tip_labels_base(opts$tree)
  habitat <- first_match(labels, "^(IS|AS|ES|NS)\\|", "\\1")
  is_mag <- !is.na(habitat)
  habitat[!is_mag] <- "Reference"
  family <- classify_family(labels, is_mag)
  bin_id <- first_match(labels, "^(?:IS|AS|ES|NS)\\|(.+_bin[0-9]+)-.+$", "\\1")
  sample_id <- ifelse(is.na(bin_id), NA_character_, sub("_bin[0-9]+$", "", bin_id))

  meta <- data.frame(
    label = labels,
    habitat = habitat,
    family = family,
    is_mag = is_mag,
    bin_id = bin_id,
    sample_id = sample_id,
    is_target = labels == opts$target,
    is_sister = labels == opts$sister,
    stringsAsFactors = FALSE
  )
  meta$is_key <- meta$is_target | meta$is_sister
  meta$display_label <- ifelse(meta$is_key, short_label(meta$label), meta$label)

  habitat_colors <- c(
    IS = "#1B9E77",
    AS = "#D95F02",
    ES = "#7570B3",
    NS = "#E7298A",
    Reference = "#9E9E9E"
  )
  family_colors <- c(
    "Group derep MAG K08356" = "#2A9D8F",
    "AioA reference" = "#1F78B4",
    "IdrA/IriA reference" = "#C65D2E",
    "Unknown reference" = "#6A4C93",
    "ArxA reference" = "#7F7F7F",
    "ArrA reference" = "#8C564B",
    "Other DMSOR reference" = "#4D4D4D",
    "Other reference" = "#BDBDBD"
  )

  metadata_file <- file.path(opts$outdir, "tip_metadata.tsv")
  write.table(meta, metadata_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write_itol_colorstrip(
    file.path(opts$outdir, "itol_habitat_colorstrip.txt"),
    "Habitat",
    meta,
    "habitat",
    habitat_colors
  )
  write_itol_colorstrip(
    file.path(opts$outdir, "itol_family_colorstrip.txt"),
    "Family",
    meta,
    "family",
    family_colors
  )
  write_itol_symbols(file.path(opts$outdir, "itol_key_symbols.txt"), meta, opts$target, opts$sister)
  write_itol_labels(file.path(opts$outdir, "itol_display_labels.txt"), meta)

  if (need_pkg("ape", required = FALSE)) {
    tree <- ape::read.tree(opts$tree)
    plot_with_ggtree(tree, meta, opts$outdir, opts$layout, opts$label_mode)
  } else {
    message("ape not installed; wrote iTOL annotation files only.")
  }

  cat("Tree tips:", length(labels), "\n")
  cat("MAG tips:", sum(meta$is_mag), "\n")
  cat("Reference tips:", sum(!meta$is_mag), "\n")
  cat("Output directory:", opts$outdir, "\n")
  cat("Wrote:", metadata_file, "\n")
}

main()
