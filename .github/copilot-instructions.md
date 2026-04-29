# TraceQC AI Coding Instructions

## Project Overview
TraceQC is an R package for quality control of CRISPR lineage tracing sequence data. It provides pipelines for bulk DNA and single-cell RNA-seq data, integrating Python-based sequence alignment with R-based mutation analysis and visualization.

## Architecture
- **R Layer**: Core functions in `R/` directory wrap Python scripts and handle data processing/visualization
- **Python Layer**: Alignment scripts in `inst/py/` use Biopython's pairwise2 for global sequence alignment
- **Data Flow**: FASTQ/BAM → alignment → mutation detection → filtering → visualization/output

## Key Workflows
- **Alignment**: Use `sequence_alignment()` for bulk or `sequence_alignment_for_10x()` for 10x single-cell data
- **Thresholding**: Run `sequence_permutation()` to determine alignment score cutoffs via permutation testing
- **Mutation Calling**: Apply `seq_to_character()` to identify insertions/deletions/substitutions
- **Filtering**: For single-cell, use `filter_mutations_per_UMI()` and `filter_mutations_per_cell()` with frequency thresholds
- **Output**: Format results with `format_mutation_df()` for standardized tabular output

## Reference File Format
Reference files contain:
- Line 1: Full construct sequence (e.g., `ATGGACTATCATATGCTTACCG...`)
- Subsequent lines: Region definitions (e.g., `target 24 140`, `spacer 88 107`)

Parse with `parse_ref_file()` which returns list with `refseq` and `regions` data.frame.

## Data Conventions
- Alignment output: Columns include `target_seq`, `target_ref`, `score`, `CB` (cell barcode), `UB` (UMI)
- Mutation data.frames: `type` (insertion/deletion/substitution), `start`, `length`, `mutate_to`, `count`
- Use `group_by(target_seq, target_ref) %>% summarise(count=n(), score=max(score))` to deduplicate identical sequences

## Dependencies & Setup
- R package built with roxygen2; install with `devtools::install_github("LiuzLab/TraceQC")`
- Requires Python packages: `pip install biopython pandas tqdm pysam`
- Specify `python_path` if not "python3"; alignment uses multiprocessing with `ncores`

## Common Patterns
- Access package data: `system.file("extdata/test_data/...", package="TraceQC")`
- Visualization: Pass parsed reference to plotting functions (e.g., `plot_construct(ref)`)
- Error handling: Check alignment scores with permutation-derived cutoffs to filter noise
- Performance: Group identical sequences before mutation calling for efficiency

## Examples
- Load reference: `ref <- parse_ref_file("path/to/ref.txt")`
- Align bulk: `sequence_alignment(input.fastq, ref_file, output.txt)`
- Call mutations: `mutations <- seq_to_character(aligned_reads, use_CPM=FALSE, alignment_score_cutoff=threshold)`
- Visualize: `plot_deletion_hotspot(mutations, ref)`</content>
<parameter name="filePath">/Users/zd/Documents/GitHub/TraceQC/.github/copilot-instructions.md