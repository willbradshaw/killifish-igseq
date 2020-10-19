# Antibody repertoire sequencing reveals systemic and mucosal immunosenescence in the short-lived turquoise killifish
[William J. Bradshaw](https://orcid.org/0000-0002-7345-3866), Michael L. Poeschla, Aleksandra Placzek, and [Dario Riccardo Valenzano](https://orcid.org/0000-0002-8761-8289)

DOI: https://doi.org/10.1101/2020.08.21.261248 (bioRxiv preprint)

## Abstract

Aging individuals exhibit a pervasive decline in adaptive immune function, with important implications for health and lifespan. Previous studies have found a pervasive loss of immune-repertoire diversity in human peripheral blood; however, little is known about repertoire aging in other immune compartments, or in species other than humans. Here, we perform the first study of immune-repertoire aging in an emerging model of vertebrate aging, the African turquoise killifish (Nothobranchius furzeri). Despite their extremely short lifespans, these killifish exhibit complex and individualised heavy-chain repertoires, with a generative process capable of producing millions of productive receptor sequences. Whole-body killifish repertoires decline rapidly in within-individual diversity with age, while between-individual variability increases. Large, expanded B-cell clones exhibit far greater diversity loss with age than small clones, suggesting an important difference in the age-sensitivity of different B-cell populations. Compared to the whole body, the immune repertoires of isolated intestinal samples exhibit much more dramatic age-related phenotypes, apparently due to an elevated prevalence of age-sensitive expanded clones. Our results highlight the importance of organ-specific dynamics in adaptive immunosenescence.

## Usage

### Pre-processing

The pre-processing pipeline takes raw IgSeq Illumina reads and carries out (among other things) quality filtering, UMI clustering, assignment of VDJ identities, and clonotyping, returning a Change-O database of unique IgH sequences from each sample. For more details, see Supplementary Note 6 from [the preprint][preprint].

The pre-processing pipeline must be run separately on each dataset to be analysed. Three such datasets were generated for this dataset: the pilot dataset (Fig S1 in [the preprint][preprint]), the ageing dataset (Fig. 2A in the preprint) and the gut dataset (Fig. 4A in the preprint). All three datasets are available via [NCBI][].

[preprint]: https://doi.org/10.1101/2020.08.21.261248
[NCBI]: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA662612

To execute the pre-processing pipeline:

- Clone this repo on your local system.
- Download and extract the raw data from [NCBI][].
- Create a separate run directory for each dataset.
- Copy the appropriate config file from `preprocessing/configs_run` (`config_pilot.yaml` for the pilot dataset, `config_ageing.yaml` for the ageing dataset, etc) to the corresponding run directory, renaming it to `config_preprocess.yaml`.
- Edit `config_preprocess.yaml` to reflect the location of the appropriate raw reads files (in the `samples` entry).
- Edit the `auxiliary_file`, `primers` and `tsa` entries of `config_preprocess.yaml` to reflect the location of these files (in the `preprocess/source/input_files` directory in this repo) on your system.
- Edit the `base_path_preprocess` entry in `config_preprocess.yaml` to reflect the absolute path to `preprocessing/source/snakefiles/base` on your system.
- Copy `preprocessing/Snakefile` to each run directory.
- [Install conda][conda] and [Snakemake][snake] on your system.
- Run Snakemake from within the appropriate run directory:
```snakemake --use-conda --cores```

NB: The preprocessing pipeline is by far the most computationally intensive component of the analysis process. It is recommended to run it on some kind of high-performance computing cluster. Depending on your cluster architecture, you may wish to use a different Snakemake command; see the [Snakemake website][cluster] for more details.

[cluster]: https://snakemake.readthedocs.io/en/v5.1.4/executable.html#cluster-execution
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[snake]: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

### Analysis

The analysis pipeline takes the sequence databases output by the pre-processing pipeline, and infers diversity spectra and (using [IGoR][igor]) generative models. Unlike the pre-processing pipeline, it can be run on multiple datasets simultaneously.

[igor]: https://github.com/qmarcou/IGoR

To generate the outputs required to generate the figures in [our preprint][preprint], the analysis pipeline can be run on the pre-processed data files included in this repo, or on the outputs of your own pre-processing runs. To do this:

- Clone this repo on your local system.
- Navigate to the `analysis/` directory.
- If using your own pre-processed data files, edit the `data_changeo` entry in `config_analysis.yaml` to reflect the paths to the output sequence databases for each dataset (within each preprocessing run directory, `outfiles/preprocess/changeo/output/seqs-all.tab`).
- If not already installed, install Conda and Snakemake on your system as described above.
- Run Snakemake from within the `analysis/` directory:
```snakemake --use-conda --cores```

### Figure generation

The pipeline to generate the figures used in [our preprint][preprint] can be run on the processed data files included in this repo, or on the corresponding output files of the pre-processing and analysis pipelines. To do this:

- Clone this repo on your local system.
- Navigate to the `figures/` directory.
- Edit `config_figures.yaml` to include paths to the required processed data files; by default it will use the files included in this repository.
- If not already installed, install Conda and Snakemake on your system as described above.
- Run Snakemake from within the `figures/` directory:
```snakemake --use-conda --cores```
