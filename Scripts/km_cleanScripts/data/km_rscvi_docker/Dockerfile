FROM lsteuernagel/r_scvi_docker_v2:v9

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libbz2-dev \
    libfftw3-dev \
    libgsl-dev \
    libprotobuf-dev \
    protobuf-compiler \
    libtbb-dev \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

# Create and export custom R lib path
ENV R_LIBS_SITE=/usr/local/lib/R_custom
RUN mkdir -p $R_LIBS_SITE

# Install Matrix and remotes explicitly with lib path override
RUN Rscript -e ".libPaths(Sys.getenv('R_LIBS_SITE')); install.packages(c('Matrix', 'remotes'), repos='https://cran.rstudio.com')" \
&& Rscript -e ".libPaths(Sys.getenv('R_LIBS_SITE')); install.packages('sctransform', dependencies = TRUE)" \
&& Rscript -e ".libPaths(Sys.getenv('R_LIBS_SITE')); library(remotes); remotes::install_github('mojaveazure/seurat-disk')" \
&& Rscript -e ".libPaths(Sys.getenv('R_LIBS_SITE')); library(remotes); remotes::install_github('lsteuernagel/mapscvi')" \
&& Rscript -e "library(mapscvi)"
