FROM ubuntu:22.04
CMD ["bash"]
ENV DEBIAN_FRONTEND noninteractive
RUN  apt-get update \
  && apt-get install -y wget gpg software-properties-common \
  && rm -rf /var/lib/apt/lists/* \
  && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
  && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
  && apt-get install -y --no-install-recommends \
    libxt-dev \
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
    r-cran-stringr \
    r-cran-tidyr \
    git \
    r-bioc-bsgenome \
  && R -e 'options(timeout = 5000);if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager");BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=F);devtools::install_github("seqinfo/BinDel", upgrade = "never")' \
  && apt remove -y wget gpg software-properties-common r-cran-devtools git \
  && apt-get autoremove -y \
  && apt-get clean -y
