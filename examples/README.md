# examples/

Small, synthetic input files for trying the pipeline end to end. **They contain
no real study data** — the taxa are common, publicly-known freshwater organisms
chosen only to exercise the code, and the OTU IDs are made up.

You only need to fabricate the two **branch entry points**; every downstream
file is produced by the tools.

| File | Feeds | Format |
|------|-------|--------|
| [`example_blast.txt`](example_blast.txt) | `taxcord condense` (NCBI branch) | pipe-delimited annotated BLAST hits (many per OTU) |
| [`example_bold.csv`](example_bold.csv) | `taxcord bold-prep` (BOLD branch) | BOLDigger3-style table (one row per OTU) |

The OTU IDs match across both files so the two branches line up at `merge`.

## What each OTU demonstrates

- **OTU0001** — both branches agree fully → resolves to **species** (*Chironomus riparius*).
- **OTU0002** — NCBI hits disagree on species but share the genus; BOLD has genus only → resolves to **genus** (*Baetis*).
- **OTU0003** — low NCBI identity (~88 %) → NCBI resolves only to **order** (Amphipoda), while BOLD has the full species, so `merge` fills the finer ranks from the single available source.
- **OTU0004** — both agree on the species *Potamopyrgus antipodarum* but disagree on **family** (NCBI `Hydrobiidae`, BOLD `Tateidae` — a real reclassification). Plain `merge` drops the species; `merge --gbif-backbone` resolves it against GBIF (accepted family `Tateidae`) and keeps it.
- **OTU0005** — low identity / no BOLD match → little support; dropped at `filter`.

## Run it

```bash
# NCBI branch
taxcord condense    examples/example_blast.txt  ncbi_condensed.txt
taxcord occurrences ncbi_condensed.txt          ncbi_ip.txt
taxcord filter      ncbi_ip.txt                 ncbi.csv

# BOLD branch
taxcord bold-prep   examples/example_bold.csv   bold_condensed.txt
taxcord occurrences bold_condensed.txt          bold_ip.txt
taxcord filter      bold_ip.txt                 bold.csv

# reconcile the two branches (--gbif-backbone keeps OTU0004's species)
taxcord merge       ncbi.csv  bold.csv  consensus.csv --gbif-backbone
```

The `occurrences` step needs network access (it queries GBIF and BOLD); the
counts you get will reflect whatever those databases currently hold for these
taxa.

The `merge --gbif-backbone` step also writes `consensus.reconciliation.tsv`
next to the output, logging the OTU0004 family fix (see the [main README](../README.md#--gbif-backbone-optional-needs-network)).
