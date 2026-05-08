# VariantPrioritiser.jl

Julia rewrite of the `ProcessParentsAlamut` / `ProcessParentsAlamutR04` workflow for `r04` GRCh38 data.

Current process:

- reads VEP-annotated small-variant `vcf` / `vcf.gz`
- auto-discovers `r04_manta/*.SV.vcf.gz` alongside the batch when present
- uses `r04_metrics`, `OMIM`, `PanelApp`, imprinted gene, contamination, homozygosity, and relatedness inputs
- writes an interactive HTML report or TSV
- supports trio, singleton, one-known-parent, and manually specified affected/unaffected family structures

## Layout

- `config/defaults.toml`: config file thresholds and resource loacations
- `src/Config.jl`: config reader
- `src/CLI.jl`: command-line parsing
- `src/VCFLoader.jl`: small-variant VEP loader
- `src/MantaLoader.jl`: `Manta` SV loader
- `src/FamilyInference.jl`: family inference and relatedness parsing
- `src/Filtering.jl`: prioritisation rules
- `src/Report.jl`: HTML / TSV report generation
- `bin/prioritise.jl`: entry point

## Typical usage

From the `julia/` directory:

```bash
julia --project=. bin/prioritise.jl \
  --html \
  --pipeline-prefix r04 \
  ../testbatchtrio/r04_vep/WG0501-WG0502-WG0503_vep_annotated.vcf.gz \
  WG0502 WG0503 WG0501
```

Positional sample arguments keep the historical ordering:

- trio: `MOTHER FATHER PROBAND`
- singleton: `PROBAND`

If the path does not contain `rNN`, pass `--pipeline-prefix r04`.

## Manual family specification

For non-standard pedigrees, specify the family directly instead of relying on `control` / relatedness inference:

```bash
julia --project=. bin/prioritise.jl \
  --html \
  --pipeline-prefix r04 \
  --parent1 MOTHER \
  --parent2 FATHER \
  --affected PROBAND,AFFECTED_RELATIVE \
  --unaffected-list UNAFFECTED_RELATIVE \
  batch/r04_vep/family_vep_annotated.vcf.gz
```

Available manual family options:

- `--parent1 SAMPLE`
- `--parent2 SAMPLE`
- `--affected sample1,sample2,...`
- `--unaffected-list sample3,sample4,...`

These override `control`-based inference.

## Supported analysis modes

### Trio

Possible outputs:

- `Homozygous Variants`
- `Compound Heterozygous Variants`
- `De Novo Variants`
- `Variants In Imprinted Genes`
- `Manta Deletions`
- `Manta Duplications`
- `Manta Insertions`

Notes:

- de novo SNVs/indels are driven by intersection with `r04_denovocnn`
- proband is shown first in the report and genotype columns

### Singleton

Possible outputs:

- `Homozygous Variants`
- `Variants Which Could Be Compound Heterozygous`
- `Variants In Imprinted Genes`
- `Manta Deletions`
- `Manta Insertions`

Notes:

- no de novo box
- no phased compound-het claim; this box means `2+` coding or canonical splice variants in the same gene
- singleton `Manta` is intentionally stricter and excludes duplications

### One known parent plus proband

Behavior is similar to trio, except:

- no de novo box
- possible compound-het calls are based on one heterozygous variant present in the known parent and one absent from the known parent in the same gene

### Other affected / unaffected family structures

These are handled as cosegregation analyses:

- variants must be present in all affected samples
- unaffected samples are used as exclusions
- report output includes a `Segregating Variants` box

## Report features

HTML output includes:

- sample cards with sample comments and QC
- relatedness / parent-child validation
- contamination and homozygosity summaries
- OMIM / PanelApp annotations
- imprinted-gene flagging
- ROH overlap marking for prioritised homozygous variants
- section-level `Copy Table` buttons
- persistent `Assessed` checkboxes using browser `localStorage`
- full-row highlight when marked as assessed

Use `--html` for HTML output and `--tsv` for tabular output.
