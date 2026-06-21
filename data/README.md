# data/

Holds your input and intermediate files. Everything here is git-ignored except
this README, so inputs, outputs, and the NCBI reference database stay local.

Suggested layout:

```
data/
├── blast/          # raw BLAST results (*.blast)
├── taxonomy/       # annotated and condensed tables
└── reference/      # accessionTaxa.sql and NCBI dumps
```

The test fixtures in [`tests/fixtures/`](../tests/fixtures/) are committed
separately so the suite runs without private data.
