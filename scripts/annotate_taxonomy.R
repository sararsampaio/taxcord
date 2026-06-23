#!/usr/bin/env Rscript
# Annotate NCBI BLAST hits with a full taxonomic lineage using taxonomizr.
#
# Usage:
#   Rscript annotate_taxonomy.R <blast_input.tsv> <taxonomy_output.txt> [accessionTaxa.sql] [--force]
#
# The input is a tab-delimited BLAST table that includes an `staxids` column.
# The output is pipe-delimited and gains one Taxonomy.<rank> column per rank,
# ready for `metatax condense`.
#
# The accessionTaxa.sql argument is optional and may live in any folder. When
# omitted, the database is found in (or downloaded to) a `data/` folder next to
# this script. Pass --force to re-download even if the database already exists.

# --- Install required packages only if missing -----------------------------
if (!requireNamespace("taxonomizr", quietly = TRUE)) {
  install.packages("taxonomizr", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(taxonomizr))

# --- Locate this script so we can resolve a default data/ folder ------------
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    normalizePath(getwd())
  }
})

# --- Find or download the taxonomy database, return its final path ----------
# Rules:
#   * A user-supplied path that points to an existing .sql file is used as-is,
#     wherever it lives.
#   * Otherwise the database is sought in (or downloaded to) <script_dir>/data.
#   * An existing database is never re-downloaded or overwritten unless force.
resolve_sql_db <- function(user_path, script_dir, force = FALSE) {
  if (!is.null(user_path) && nzchar(user_path)) {
    if (!grepl("\\.sql$", user_path, ignore.case = TRUE)) {
      stop("Invalid SQL database path (must end in .sql): ", user_path)
    }
    sql_db <- normalizePath(user_path, mustWork = FALSE)
    target_dir <- dirname(sql_db)
  } else {
    target_dir <- file.path(script_dir, "data")
    sql_db <- file.path(target_dir, "accessionTaxa.sql")
  }

  # Reuse an existing database unless the user forces a refresh.
  if (file.exists(sql_db) && !force) {
    if (file.access(sql_db, mode = 4) != 0) {
      stop("SQL database exists but is not readable (check permissions): ", sql_db)
    }
    return(sql_db)
  }

  # Need to download/build: ensure the target directory exists and is writable.
  if (!dir.exists(target_dir)) {
    if (!dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create directory for the SQL database (check permissions): ", target_dir)
    }
  }
  if (file.access(target_dir, mode = 2) != 0) {
    stop("No write permission for the SQL database directory: ", target_dir)
  }

  # Build the local taxonomy database. This downloads the NCBI accession-to-taxid
  # and names/nodes dumps and is several gigabytes.
  tryCatch(
    prepareDatabase(sql_db),
    error = function(e) {
      stop("Failed to download/build the NCBI taxonomy database at ",
           sql_db, ": ", conditionMessage(e))
    }
  )
  if (!file.exists(sql_db)) {
    stop("Taxonomy database download did not produce a file at: ", sql_db)
  }
  normalizePath(sql_db)
}

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
args <- args[args != "--force"]
if (length(args) < 2) {
  stop("Usage: annotate_taxonomy.R <blast_input> <taxonomy_output> [accessionTaxa.sql] [--force]")
}

blast_file <- args[1]
output_file <- args[2]
user_sql <- if (length(args) >= 3) args[3] else NULL
sql_db <- resolve_sql_db(user_sql, script_dir, force = force)

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

blast <- read.delim(blast_file, sep = "\t", header = TRUE)
blast$Taxonomy <- getTaxonomy(blast$staxids, sql_db, desiredTaxa = ranks)

write.table(blast, output_file, sep = "|", quote = FALSE)
