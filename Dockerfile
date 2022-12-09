FROM ubuntu:22.04
CMD ["bash"]
ENV DEBIAN_FRONTEND noninteractive
RUN  apt-get update \
  && apt-get install -y wget gpg software-properties-common\
  && rm -rf /var/lib/apt/lists/*
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get install -y --no-install-recommends r-base
RUN apt-get install -y r-cran-devtools r-bioc-biostrings r-cran-dplyr r-bioc-genomicalignments r-bioc-genomicranges r-cran-ggplot2  r-bioc-iranges r-cran-magrittr r-cran-purrr r-cran-readr r-bioc-rsamtools r-cran-stringr  r-cran-tidyr git r-bioc-bsgenome 
RUN R -e  'options(timeout = 5000);if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager");BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=F)#;devtools::install_github("seqinfo/BinDel", upgrade = "never")'
RUN R -e  'devtools::install_github("seqinfo/BinDel", upgrade = "never")'
