FROM ubuntu

## Installing dependencies

RUN apt-get update
RUN apt-get -y install r-base
RUN apt-get -y install aptitude libcurl4-openssl-dev  libxml2-dev libxml2-dev
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('randomForest')"

## Make sure all files that your execution scripts need are in the Docker image

COPY Dry_Run_SC3.R /Dry_Run_SC3.R
COPY score_sc3.sh /score_sc3.sh
COPY model_storage /model_storage