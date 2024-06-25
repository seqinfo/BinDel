# BinDel: software tool for detecting clinically significant microdeletions in low-coverage WGS-based NIPT samples
BinDel is distributed under the Attribution-NonCommercial-ShareAlike 4.0 International ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)) license.

Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
All rights reserved, unauthorised usage and distribution are prohibited.
Contact: priitpaluoja@gmail.com/priit.palta@gmail.com



Here we present the BinDel, a novel region-aware microdeletion detection software package developed to infer clinically relevant microdeletion risk in low-coverage whole-genome sequencing NIPT data. 

Our [paper](https://doi.org/10.1101/2022.09.20.22280152) describes the BinDel algorithm and how it was tested. We quantified the impact of sequencing coverage, fetal DNA fraction, and region length on microdeletion risk detection accuracy. We also estimated BinDel accuracy on known microdeletion samples and clinically validated aneuploidy samples. 


## Installation
#### Docker: [![Docker Repository on Quay](https://quay.io/repository/priitpaluoja/bindel/status "Docker Repository on Quay")](https://quay.io/repository/priitpaluoja/bindel)
<details><summary>Installation on Ubuntu 22.04</summary>
<p>
 
The following is tested with [ubuntu-22.04.1-live-server-amd64](https://releases.ubuntu.com/22.04/).

#### Install R as shown in [DigitalOcean](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04). [From DigitalOcean](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04):
```bash
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list
sudo apt update
sudo apt install --no-install-recommends r-base
```
#### Install BinDel dependencies and [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)
```bash
sudo apt -y install r-cran-devtools r-bioc-biostrings r-cran-dplyr r-bioc-genomicalignments r-bioc-genomicranges r-cran-ggplot2  r-bioc-iranges r-cran-magrittr r-cran-purrr r-cran-readr r-bioc-rsamtools r-cran-stringr  r-cran-tidyr git r-bioc-bsgenome  libcairo2-dev libxt-dev
```
#### Install BSgenome.Hsapiens.UCSC.hg38 and BinDel
```R
sudo -i R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
devtools::install_github("seqinfo/BinDel", upgrade = "never")
```
</p>
</details>


<details><summary>Installation on Windows 10</summary>
<p>

1. Install [R](https://cran.r-project.org/bin/windows/base/).
2. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
3. Install devtools and BinDel in R:
```R
# In R:
install.packages("devtools") 
devtools::install_github("seqinfo/BinDel")
```
</p>
</details>


## Alignment/Mapping
BinDel requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Usage
### Reference creation
BinDel requires the creation of a reference set file. A reference set file is a file that contains known euploid NIPT samples. The read counts of these samples are used to compare the sample of interest with the healthy reference group. The creation of the reference file requires known euploid NIPT samples in `.bam` format and
a [coordinates file](example/locations.info.tsv) describing bin lengths and regions of interest. 

A reference file with at least 200 samples is recommended to analyse regular WGS NIPT. The reference sample set size is not constrained, but PCA-based normalisation yields more accurate analysis when an increased number of euploid fetus samples are used as a reference.

```R
# In R:
BinDel::write_reference(c("sample1.bam", "sample2.bam"), "coordinates.tsv", "reference.gz")
```

In addition to generating a single reference file in tab-separated text format, it's possible to create individual reference files per sample and subsequently merge them. This approach can be advantageous in scenarios involving high-performance computing (HPC) workflows or when memory constraints are a concern.

For instance, one can create a reference with the first sample as follows:
```R
# In R
BinDel::write_reference(c("sample1.bam"), "coordinates.tsv", "reference_with_header.gz", anonymise=F)
```
For subsequent sample(s), the reference can be generated without including the column names/header:
```R
# In R
BinDel::write_reference(c("sample2.bam"), "coordinates.tsv", "reference_no_header.gz", col_names = F, anonymise=F)
```
Finally, to consolidate the reference files, they can be merged using the following command (in bash):

```bash
cat reference_with_header reference_no_header.gz > reference.gz
```

### Running BinDel
```R
# In R:
BinDel::infer_normality("sample.bam", "reference.gz")
```
Note that `infer_normality` with small reference set size and high cumulative PCA settings (regional, and `cumulative_variance_per_region = c(50)`) may lead to `NA` microdeletion risk scores (omitted from the output) or the inability to call microdeletion risk accurately.

Based on the [preprint](https://www.medrxiv.org/content/10.1101/2022.09.20.22280152v2.full-text), the recommended analysis approach would be adopting a microdeletion-region-specific cut-off threshold when NIPT samples featuring microdeletions are available for calibration. This calibration should encompass parameters such as `cumulative_variance`, `cumulative_variance_per_region`, and the cut-off threshold, tailored to align with the application’s specific requirements, ensuring a balance between sensitivity and specificity that aligns with the intended analysis purposes. In cases where a region-specific approach is not feasible, we advocate the adoption of PCA95% (`cumulative_variance = 95.0`) in conjunction with a microdeletion cut-off threshold set at 90%. 

### Interpreting output
Results location: The output file name and location are determined by the function `infer_normality` named parameter `output_file_path` (presumes that the folder exists).

The `md_prob` is a probability of microdeletion risk with a risk score between 0-100, where 100 represents the highest risk. In non-standard NIPT analysis (exploratory or novel method development), one might benefit from using `dist_md` Mahalanobis distance (set `output_intermediate_scores = T`), which is not capped at 100 but requires distribution analysis to determine where to set high microdeletion cut-off. Please note that the risk score interpretation may vary between different NIPT WGS protocols.
