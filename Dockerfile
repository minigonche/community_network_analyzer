# Use R Base
FROM rocker/tidyverse:4

# Set metadata for the image
LABEL maintainer="Felipe González Casabianca <minigonche@gmail.com>"
LABEL version="1.0"
LABEL description="Docker image for the Nextflow Community Network Analyazer Pipeline"

RUN apt-get update && apt-get install -y cmake

RUN R -e "install.packages(c('missMDA'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('optparse'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('readxl','tidyr','dplyr', 'igraph','vegan'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('GraceYoon/SPRING')"





