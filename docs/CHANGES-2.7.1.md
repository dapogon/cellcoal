# CellCoal 2.7.1 — what changed

A small release: one behaviour change and the documentation.

## Behaviour

- `-d` (single-nucleotide deletions of the sites) is now an error when it is used with a chromosome in blocks (`--blocks`), with a reference (`--fasta`) or in the reads mode (`--bam`, `--fastq`). With a chromosome the deletions are small indels: the message points to `--indel-rate` (with `--indel-ins`, `--indel-len`, `--indel-max`). `-d` keeps working for isolated sites, as before, and the golden outputs do not change.

## Manual

- New section 5.4.9: in which order ADO, the gain along the genome and the coverage act, and how to simulate MDA and PTA (comparison table).
- `examples/config.mda.yaml` and `examples/config.pta.yaml`: illustrative parameter sets (the numbers are not calibrated).
- New figures: main flow of CellCoal 2, tumor regions with colonizations (monoclonal and polyclonal seeding), a chromosome in blocks with regional rates, and indels in the reads.
- Every figure has a number and a caption (the three figures of the quick start are now Figs. 1-3, so the other ones are Figs. 4-14); captions are shown in a smaller shaded block.
- The PDF has page numbers.
- Softer wording on colonization: a metastasis is assumed to be founded once, later seedings are not modeled.
