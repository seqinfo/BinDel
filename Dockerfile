FROM ubuntu:24.04
CMD ["bash"]
ENV DEBIAN_FRONTEND noninteractive
RUN  apt-get update \
  && apt-get install -y --no-install-recommends \
    r-base \
    r-cran-devtools \
    r-bioc-biostrings \
    r-cran-dplyr \
    r-bioc-genomicalignments \
    r-bioc-genomicranges \
    r-bioc-iranges \
    r-cran-magrittr \
    r-cran-purrr \
    r-cran-readr \
    r-bioc-rsamtools \
    git \
    r-bioc-bsgenome \
  && R -e 'options(timeout = 15000);if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager");BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=F);devtools::install_github("seqinfo/BinDel", upgrade = "never")' \
  && apt remove -y r-cran-devtools git \
  && apt-get autoremove -y \
  && apt-get clean -y
