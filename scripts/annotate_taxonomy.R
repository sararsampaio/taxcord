#!/usr/bin/env Rscript
# Annotate NCBI BLAST hits with a full taxonomic lineage using taxonomizr.
#
# Usage:
#   Rscript annotate_taxonomy.R <blast_input.tsv> <taxonomy_output.txt> [accessionTaxa.sql]
#
# The input is a tab-delimited BLAST table that includes an `staxids` column.
# The output is pipe-delimited and gains one Taxonomy.<rank> column per rank,
# ready for `metatax condense`.

suppressPackageStartupMessages(library(taxonomizr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: annotate_taxonomy.R <blast_input> <taxonomy_output> [accessionTaxa.sql]")
}

blast_file <- args[1]
output_file <- args[2]
sql_db <- if (length(args) >= 3) args[3] else "accessionTaxa.sql"

# Build the local taxonomy database on first use. This downloads the NCBI
# accession-to-taxid and names/nodes dumps and is several gigabytes.
if (!file.exists(sql_db)) {
  prepareDatabase(sql_db)
}

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

blast <- read.delim(blast_file, sep = "\t", header = TRUE)
blast$Taxonomy <- getTaxonomy(blast$staxids, sql_db, desiredTaxa = ranks)

write.table(blast, output_file, sep = "|", quote = FALSE)
