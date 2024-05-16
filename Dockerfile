# Use R Base
FROM rocker/tidyverse:4

# Set metadata for the image
LABEL maintainer="Felipe González Casabianca <minigonche@gmail.com>"
LABEL version="1.0"
LABEL description="Docker image for the Nextflow Community Network Analyazer Pipeline"

RUN R -e "install.packages(c('missMDA','readxl','tidyr','dplyr', 'igraph','optparse','vegan'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('GraceYoon/SPRING')"





