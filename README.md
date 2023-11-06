# SARS-CoV-2 Origin Study

## Description
This bioinformatics project is focused on tracing the origins of SARS-CoV-2, the virus known for causing COVID-19, by comparing its genetic material with that of coronaviruses found in animals. We use a technique called the Levenshtein Distance algorithm to measure the number of differences between the amino acid sequences of two different viral genomes. This serves as a proxy for determining the minimum number of mutations that must have occurred within the viral RNA to transform a coronavirus common within other animals into SARS-CoV-2, and thus, the most likely original animal host.

## Results
Our research suggests that pangolins or bats are the most likely sources of SARS-CoV-2, with slightly more evidence pointing towards the former. For a detailed exploration of our study, I have made the [full report](https://andylabs.org/static/media/SARS-CoV-2-Origin-Study.pdf) available on my site.

## Data Visualizations

### Minimum Mutation Analysis Table
![Minimum Mutation Analysis Table](https://user-images.githubusercontent.com/91595477/205288133-45c5e0f7-2816-46c5-86ff-dc9dc40b8e71.png)
*Figure 1: An Excel snapshot capturing the minimum number of mutations required for various coronaviruses to match the reference SARS-CoV-2 strain. The dataset includes five samples from each animal species tested, along with MERS and SARS for comparison, and five other SARS-CoV-2 variants.*

### Comparative Mutation Frequency Graph
![Comparative Mutation Frequency Graph](https://user-images.githubusercontent.com/91595477/205288069-233a7c10-4ac0-4d14-ba0e-d741635f4b66.png)
*Figure 2: A bar graph illustrating the comparative analysis of mutation frequencies. The graph compares the number of mutations across different coronaviruses, including samples from each animal species tested, relative to the reference SARS-CoV-2 strain.*

## Repository Contents

- `Main.py`: This script runs the core analysis. It facilitates the process of RNA translation into amino acids, aligns the genetic sequences, and calculates their differences.
- `FASTAs/`: This directory contains FASTA files, which are formatted sequences of RNA genomes from various animal coronaviruses. All FASTA files have been found and extracted from the [National Center for Biotechnology Information (NCBI) database](https://www.ncbi.nlm.nih.gov/nuccore).

## Installation and Setup
1. **Download or clone the repository**: Clone or download this repository to your local machine to get started.
   
   ```bash
    git clone https://github.com/AndyAnderson8/SARS-CoV-2-Origin-Study.git
    cd SARS-CoV-2-Origin-Study
    ```

2. **Running analysis**: Simply launch the `Project.py` file.
   
    ```bash
    python Project.py
    ```

## License
[MIT](https://github.com/AndyAnderson8/SARS-CoV-2-Origin-Study/blob/main/LICENSE.txt)
