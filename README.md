# RNA editing QTL pipeline

<div align="justify">

A Snakemake pipeline for mapping cis RNA editing QTLs (edQTLs) in each tissue. It quantifies A-to-I RNA editing at known sites using RNA-seq BAM files, builds and normalizes a tissue-specific editing-ratio matrix, and maps cis-edQTLs with [tensorQTL](https://github.com/broadinstitute/tensorqtl), using genotype PCs, phenotype (editing-ratio) PCs, sex, and age as covariates.

The pipeline is based on Li et al. (2022), “RNA editing underlies genetic risk of common inflammatory diseases,” and is also implemented in [Pantry](https://github.com/PejLab/Pantry) (PAN-TRanscriptomic phenotYping), a framework for multimodal transcriptomic phenotyping and QTL mapping.

Use this repository to run the RNA editing pipeline independently. To analyze RNA editing alongside six additional RNA phenotypes, use Pantry.

It was written for GTEx v8 (BAM files aligned with STAR to the hg38 reference genome, and WGS genotypes), but should work for any cohort with BAM files, a VCF, and sample annotation.

## Installation

```sh
mamba env create -f envs/environment.yml
mamba activate rna_editing_env
```

The environment installs a CPU build of PyTorch. tensorQTL runs on CPU, but a CUDA GPU is much faster, mainly for the permutation step. To use a GPU, run the tensorQTL rules in an environment with a CUDA build of PyTorch and check with `python -c "import torch; print(torch.cuda.is_available())"`.

## Input

Paths are set in `config/config.yaml`. Input data is not included in the repository.

- RNA-seq BAMs, one per sample, in `data/bam/` (named `{sample_id}.{bam_extension}`)
- Genotypes as a VCF, with subject IDs as sample names
- Reference genome FASTA (the same build used for alignment)
- A BED file of known editing sites (a filtered hg38 list is included in `src/`; full hg38 and hg19 lists are available from [REDIportal](https://rediportal.cloud.ba.infn.it/atlas/index.html))
- A sample annotation file (`sample_id`, `subject_id`, `tissue_cln`)
- A subject covariates file (`subject_id`, `sex`, `age`)

The annotation and covariates files in `src/` are small examples of the format.

## Configuration

All settings are in `config/config.yaml`. Mapping and base quality cutoffs are at the top of `scripts/query_editing_level.pl`.

There are two normalization options (`normalization`):

- `quantile_int` (default): quantile normalization of samples, then a rank-based inverse normal transform of each site across samples, as in Pantry and the GTEx eQTL pipeline. Every site ends up normally distributed, so a few extreme samples can't drive an association.
- `zscore_qqnorm`: z-score of each site, then quantile normalization of each sample to a normal distribution, as in LeafCutter and Li et al. (2022).

## Running

Run from the repository root:

```sh
snakemake -j 16 -n                                  # dry run
snakemake -j 16 --printshellcmds --latency-wait 30  # run
```

To limit how many tensorQTL jobs run at the same time (useful on a single GPU), add `--resources tensorqtl_slots=1`. On CPU, leave it out so tissues can run in parallel.

## Output

Results go to `results/{tissue}/`:

- `{tissue}.cis_edQTLs.txt.gz`: top variant per site from permutations, with q-values
- `{tissue}.cis_independent_edQTLs.txt.gz`: conditionally independent edQTLs
- `{tissue}.cis_edqtl_signif.txt.gz`: all significant site-variant pairs
- `{tissue}.cis_edqtl_all_pvals.txt.gz`: nominal p-values for all tested pairs
- `nominal/`: full nominal statistics, one parquet file per chromosome

Editing levels and matrices are in `intermediate/`, genotype files and covariates in `data/genotype/{tissue}/`, and logs in `logs/`. Site IDs are `chrN_position`.

## Credits and references

This pipeline builds on code and methods from:

- Qin Li, Michael J. Gloudemans, Jonathan M. Geisinger, Boming Fan, François Aguet, Tao Sun, Gokul Ramaswami, Yang I. Li, Jin-Biao Ma, Jonathan K. Pritchard, Stephen B. Montgomery, Jin Billy Li. RNA editing underlies genetic risk of common inflammatory diseases. *Nature* 608, 569–577 (2022). https://doi.org/10.1038/s41586-022-05052-x
- Daniel Munro, Nava Ehsan, Seyed Mehdi Esmaeili-Fard, Alexander Gusev, Abraham A. Palmer, Pejman Mohammadi. Multimodal analysis of RNA sequencing data powers discovery of complex trait genetics. *Nature Communications* 15, 10387 (2024). https://doi.org/10.1038/s41467-024-54840-8

## License

MIT, see [LICENSE](LICENSE).

</div>
