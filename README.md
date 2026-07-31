# P_Malariae_Transcriptomics
Source code and methodological record for the transcriptomic analysis of *Plasmodium malariae* in the mosquito *Anopheles Gambiae*

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

## Repository Scope & Intent

The scripts provided here are intended to act as an audit trail. They document the exact commands, parameters, R packages, and filtering thresholds utilized during our study. **They are not intended to be a turnkey, software pipeline**, as they were executed in a specific high-performance computing (HPC) environment and rely on specific local configurations.

Researchers canuse these scripts as a reference to understand our methodology in detail, or to adapt the code snippets for their own analytical workflows.


## Project Description

This repository contains the computational pipelines, custom scripts, and bioinformatics workflows developed for our study. It serves as a transparent methodological record of all the analysis made.

For a complete overview of the study's background, methodology, and findings, please refer to our full publication. 

**Abstract:**
>   To survive and multiply in different host niches, malaria parasites require sets of proteins that are expressed during the life cycle in a timely manner. Plasmodium malariae is a human malaria parasite that is present in most malaria-endemic regions. This species remains enigmatic due to the inability to maintain it in long-term culture, the low parasitaemia observed in human infections, and its frequent co-infection with other malaria species.
> 
>  We investigate the transcriptome of P. malariae during transmission in the vector by analyzing RNA-seq data from midgut and salivary gland parasite stages obtained in two experimental mosquito infections with field P. malariae isolates. Our analysis identified 3,699 expressed genes, of which 263 are developmentally regulated and 309 are P. malariae specific, including genes with unknown functions and without orthologs in other Plasmodium. We detected unique expression patterns of the ApiAP2 family of transcription factors, many of which appear to function as master regulators during sporogony. Several members of multigene families, such as PIR, PHIST, and fam-l, were found to be expressed, and new families potentially showing transcriptional heterogeneity were described.
> 
>Our findings indicate that the P. malariae transcriptome is largely conserved; however, a distinct set of expressed genes lacking detectable orthologs may contribute to differences in transmission and virulence. Given the scarcity of genomic and transcriptomic data for this pathogen and considering the global distribution of P. malariae and its potential role in sustaining malaria transmission, this study has direct implications for malaria control and elimination efforts.


## Code Archive Breakdown

The analysis workflow is documented across the following modular scripts:
*   **`01_XXXXX.sh`**: Records the raw data processing (KrakenUniq, SortMeRNA), quantification (CoCo), and the Differential Expression Analysis (DESeq2) comparing midgut and salivary gland stages.
*   **`02_XXXXX.sh`**: Details the competitive mapping approach (HISAT2), quality filtering, duplicate marking (Picard), and read extraction.
*   **`03_XXXXXs.R`**: Shows the Gene Ontology (GO) enrichment analysis (using `TopGO`) and the integration of `PANNZER2` functional predictions for unannotated genes.
*   **`04_XXXXX.sh`**: Contains the logic for analyzing ApiAP2 transcription factor expression and the exact `HOMER` commands used for promoter motif discovery.
*   **`05_XXXXX.R`**: Demonstrates the Reciprocal Best Hit (RBH) ortholog detection across *Plasmodium* species using the `orthologr` package.
*   **`06_XXXXX.sh`**: Documents the identification of Hypervariable Genes (HVGs) using genomic locations, BLAST, and integration with Malaria Cell Atlas (MCA) single-cell data.
*   **`07_XXXXX.R`**: Details the comparative analysis between our mosquito-stage data and previously published blood-stage datasets.

## Data Availability

To adhere to code repository best practices, **no raw data, intermediate BAM files, or generated tables are hosted here.** 

*   Raw RNA-Seq `.fastq` files associated with this study are deposited in GEO under accession GSE231944 .
*   Reference genomes were obtained from [PlasmoDB](https://plasmodb.org/) and [VectorBase](https://vectorbase.org/).


## Citation

If you use the code, pipelines, or methodological approaches provided in this repository for your research, please cite our publication:

**Ruiz, J. L., Díaz-Terenti, B., Lefèvre, T., Yerbanga, R. S., Otto, T. D., & Gómez-Díaz, E.** (2026). *The dynamic transcriptome of Plasmodium malariae in the mosquito*. DOI: [https://doi.org/10.1101/2025.05.05.652270]

**BibTeX:**
```bibtex
@article{ruiz_2026_pmalariae,
  title = {The dynamic transcriptome of Plasmodium malariae in the mosquito},
  author = {Ruiz, J. L. and D{\'i}az-Terenti, B. and Lef{\`e}vre, T. and Yerbanga, R. S. and Otto, T. D. and G{\'o}mez-D{\'i}az, E.},
  year = {2026},
  doi = {[https://doi.org/10.1101/2025.05.05.652270]}
}



