![Code Ocean Logo](CO_logos/CO_logo_135x72.png)
# Read Quality Control with Base-R Plots

Quality control of paired-end sequencing reads. The capsule runs `fastp` over the
attached FASTQ files, then draws the QC results as figures using base R graphics
only and writes a flat metrics table.

## What it does

1. **Discovery** — searches every folder under `../data` for a read 1 FASTQ by
   pattern (`*_R1*` or `*_1`, gzipped or plain) and pairs it with its read 2 mate,
   so neither folder names nor file names are hard-coded.
2. **Quality control and trimming** — `fastp` with adapter auto-detection for
   paired-end reads:

   | Setting | Value |
   | --- | --- |
   | `--qualified_quality_phred` | 15 |
   | `--n_base_limit` | 5 |
   | `--length_required` | 30 bp |
   | `--detect_adapter_for_pe` | enabled |

3. **Figures and metrics** — `qc_plots.R` reads the `fastp` JSON report and draws
   six figures with `plot`, `lines`, `barplot`, `legend` and `text`.

## Files

| File | Purpose |
| --- | --- |
| `run` | Bash driver: locates reads, runs `fastp`, calls the R script |
| `qc_plots.R` | Reads the JSON report, draws the figures, writes the metrics table |

## Outputs

Written to `../results`:

| Path | Contents |
| --- | --- |
| `figures/fig1_quality_per_cycle.png` | Mean Phred score per cycle, read 1 and read 2, before vs after filtering |
| `figures/fig2_base_content.png` | Base composition (A/C/G/T/N) per cycle, both reads |
| `figures/fig3_gc_per_cycle.png` | GC content per cycle, with the overall level marked |
| `figures/fig4_filtering_and_quality.png` | Read fate during filtering, and the >=Q20 / >=Q30 base fractions before vs after |
| `figures/fig5_duplication.png` | Number of distinct sequences by how often each was observed |
| `figures/fig6_insert_size.png` | Insert size distribution, or a statement when no pair yields one |
| `reports/<sample>_qc_metrics.tsv` | Flat `metric<TAB>value` table of all reported QC values |
| `reports/<sample>_fastp.json` | Machine-readable QC report (input to `qc_plots.R`) |
| `reports/<sample>_fastp.html` | Interactive `fastp` report |
| `reports/<sample>_fastp.log` | `fastp` run log |
| `trimmed/<sample>_trimmed_R[12].fastq.gz` | Reads surviving filtering |

## Running it

```shell
/code/run
```

The script takes no arguments. It searches every folder under `../data` for a
read 1 FASTQ (`*_R1*` or `*_1`, gzipped or plain), pairs it with its read 2 mate,
and names the outputs after it. It exits non-zero if no reads are found, or if the
matching read 2 file is missing.

## Environment

**Starter Environment**
```
R 4.0 or newer
```

**R (CRAN)**
```
jsonlite=2.0.0
```

### Functional

**apt-get**
```
fastp=0.20.1+dfsg-1
```

### Testing (Causes a build error)
**R (CRAN)**
```
fastp=42
```


## Notes

- `qc_plots.R` depends on `jsonlite` only to read the report; every plot is drawn
  with base R graphics.
- The insert-size figure reports "not determinable" when `fastp` cannot estimate
  an insert for any read pair, which happens when read 1 and read 2 do not behave
  as mates.
- The metrics table reports post-filter read counts from the `fastp` report; it
  does not recount the trimmed FASTQ files.



![Code Ocean](CO_logos/CO_logo_68x36.png)
