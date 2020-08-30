# Antibody repertoire sequencing reveals systemic and mucosal immunosenescence in the short-lived turquoise killifish
[William J. Bradshaw](https://orcid.org/0000-0002-7345-3866), Michael L. Poeschla, Aleksandra Placzek, and [Dario Riccardo Valenzano](https://orcid.org/0000-0002-8761-8289)

DOI: https://doi.org/10.1101/2020.08.21.261248 (bioRxiv preprint)

## Abstract

Aging individuals exhibit a pervasive decline in adaptive immune function, with important implications for health and lifespan. Previous studies have found a pervasive loss of immune-repertoire diversity in human peripheral blood; however, little is known about repertoire aging in other immune compartments, or in species other than humans. Here, we perform the first study of immune-repertoire aging in an emerging model of vertebrate aging, the African turquoise killifish (Nothobranchius furzeri). Despite their extremely short lifespans, these killifish exhibit complex and individualised heavy-chain repertoires, with a generative process capable of producing millions of productive receptor sequences. Whole-body killifish repertoires decline rapidly in within-individual diversity with age, while between-individual variability increases. Large, expanded B-cell clones exhibit far greater diversity loss with age than small clones, suggesting an important difference in the age-sensitivity of different B-cell populations. Compared to the whole body, the immune repertoires of isolated intestinal samples exhibit much more dramatic age-related phenotypes, apparently due to an elevated prevalence of age-sensitive expanded clones. Our results highlight the importance of organ-specific dynamics in adaptive immunosenescence. 

## Usage

Code for running the pre-processing and analysis pipelines will be made available in the near future, along with links to the raw data on SRA.

Processed data for figure generation is provided in this repository. To generate figures from this data:

- [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on your system.
- Navigate to the `figures` directory.
- Edit `config_figures.yaml` to include paths to the required processed data files; by default it will use the files included in this repository.
- To the model with all available cores, run:
```snakemake --use-conda --cores```
- To specify a number of cores, run:
```snakemake --use-conda --cores <n_cores>```
