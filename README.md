<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-issentinel_logo_dark.png">
    <img alt="nf-core/issentinel" src="docs/images/nf-core-issentinel_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/issentinel)
[![GitHub Actions CI Status](https://github.com/nf-core/issentinel/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/issentinel/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/issentinel/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/issentinel/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/issentinel/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/issentinel)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23issentinel-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/issentinel)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/issentinel** is a bioinformatics pipeline designed to detect new Insertion Sequence (IS) element insertions in bacterial genomes using high-throughput sequencing data (paired-end FastQ reads).

The pipeline uses an engineered reference genome containing a "decoy contig" of the insertion sequence of interest, alongside masking of existing reference insertion sequence locations. This setup ensures that reads spanning new insertion junctions map with high mapping quality (MAPQ) to both the genomic target and the decoy sequence. Structural variant calling is subsequently performed using Delly to identify these junctions as translocations or insertions, and candidate insertions are annotated, compiled, and filtered.

### Pipeline Workflow
1. **Reference Engineering (`PREPARE_REFERENCE`)**: Masks existing insertion sequence coordinates in the reference genome and appends the target IS sequence as a decoy contig.
2. **Reference Indexing (`INDEX_REFERENCE`)**: Indexes the engineered reference genome using `bwa index` and `samtools faidx`.
3. **Read Alignment (`BWA_ALIGN`)**: Aligns paired-end reads to the engineered decoy reference using BWA-MEM.
4. **Structural Variant Calling (`DELLY_CALL`)**: Calls translocations (TRA/BND) and insertions (INS) using Delly.
5. **Candidate Extraction (`EXTRACT_IS_CANDIDATES`)**: Parses the VCF output, filtering for insertions supported by paired-end and split-read evidence, and annotating them with features like insertion coordinates, orientation, double-strand break (DSB) junctions, allele frequency (VAF), and heteropopulation status.
6. **Result Compilation (`COMPILE_RESULTS`)**: Aggregates candidate insertions across all samples into a single master summary table.
7. **Region of Interest Filtering (`FILTER_REGION`)**: (Optional) Filters candidates within specific target coordinate windows (e.g., drug resistance loci).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:
```csv
sample,fastq_1,fastq_2
ERR13259960,reads/ERR13259960_1.fastq.gz,reads/ERR13259960_2.fastq.gz
ERR13260062,reads/ERR13260062_1.fastq.gz,reads/ERR13260062_2.fastq.gz
```

Each row represents a sample with its corresponding paired-end FastQ reads.

Now, you can run the pipeline using the command below. Note that all parameters starting from `--run_pipeline` are optional (they are shown here explicitly, but omitting them will fall back to the default settings for MTB H37Rv, IS6110, and the mmpL5-Rv0678 region):

```bash
nextflow run nf-core/issentinel \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --run_pipeline true \
   --run_compile true \
   --run_filter true \
   --decoy_name 'IS6110' \
   --ref_genome 'reference/H37RV.fna' \
   --gbk_file 'reference/H37RV_REF.gb' \
   --is_fasta 'reference/IS6110.fasta' \
   --filter_name 'mmpL5-Rv0678' \
   --filter_chrom 'NC_000962.3' \
   --filter_start 775586 \
   --filter_end 779487
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/issentinel/usage) and the [parameter documentation](https://nf-co.re/issentinel/parameters).

## Pipeline Parameters

`nf-core/issentinel` is designed to be bacterial-agnostic and insertion sequence (IS)-agnostic. While it defaults to analyzing the *Mycobacterium tuberculosis* (H37Rv) genome with the `IS6110` insertion sequence and filtering for the `mmpL5-Rv0678` genomic region, you can fully customize the reference files, decoy names, and genomic coordinates.

### Reference and IS Options
* `--decoy_name`: Name of the insertion sequence decoy contig (default: `'IS6110'`).
* `--ref_genome`: Path to the raw FASTA reference file (default: `"${projectDir}/reference/H37RV.fna"`).
* `--gbk_file`: Path to the GenBank reference record used for candidate annotation (default: `"${projectDir}/reference/H37RV_REF.gb"`).
* `--is_fasta`: Path to the FASTA sequence of the Insertion Sequence element to mask and decoy (default: `"${projectDir}/reference/IS6110.fasta"`).

### Run Execution Toggles
* `--run_pipeline`: Run BWA alignment, SV calling with Delly, and candidate extraction (default: `true`).
* `--run_compile`: Compile all sample candidates into a single master summary table (default: `true`).
* `--run_filter`: Filter the compiled master table for variants falling within a specific genomic Region of Interest (default: `false`).

### Genomic Region Filtering
When `--run_filter` is set to `true`, you can define custom target coordinate windows:
* `--filter_name`: Custom name of the target region/genes (default: `'mmpL5-Rv0678'`).
* `--filter_chrom`: Chromosome/contig name matching the reference FASTA (default: `'NC_000962.3'`).
* `--filter_start`: 1-based start position of the target region window (default: `775586`).
* `--filter_end`: 1-based end position of the target region window (default: `779487`).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/issentinel/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/issentinel/output).

## Credits

nf-core/issentinel was originally written by Fernando Falat.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#issentinel` channel](https://nfcore.slack.com/channels/issentinel) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations


If you use nf-core/issentinel for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
