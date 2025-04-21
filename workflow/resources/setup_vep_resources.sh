#!/bin/bash

set -euo pipefail

#
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Constants for URLs
HUMAN_ANCESTOR_FA="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa"
PHYLOCSF_GERP="https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz"

# Downloader selection
if command -v wget &> /dev/null; then
  downloader="wget"
elif command -v curl &> /dev/null; then
  downloader="curl"
else
  echo "Neither wget nor curl is installed. Aborting."
  exit 1
fi

# Download VCF using the chosen downloader (wget or curl)
mkdir -p plugin_resources/loftee && cd $_

if [ "$downloader" = "wget" ]; then
  ${downloader} \
	--timestamping \
	--append-output=./wget_loftee_resources.log \
	--show-progress \
	${HUMAN_ANCESTOR_FA}.gz
  ${downloader} \
	--timestamping \
	--append-output=./wget_loftee_resources.log \
	--show-progress \
	${HUMAN_ANCESTOR_FA}.gz.fai
  ${downloader} \
	--timestamping \
	--append-output=./wget_loftee_resources.log \
	--show-progress \
	${HUMAN_ANCESTOR_FA}.gz.gzi
  ${downloader} \
	--timestamping \
	--append-output=./wget_loftee_resources.log \
	--show-progress \
	${PHYLOCSF_GERP}
else
  ${downloader} -vL -o human_ancestor.fa.gz \
    ${HUMAN_ANCESTOR_FA}.gz 2>&1 | tee -a curl_loftee_resources.log
  ${downloader} -vL -o human_ancestor.fa.gz.fai \
    ${HUMAN_ANCESTOR_FA}.gz.fai 2>&1 | tee -a curl_loftee_resources.log
  ${downloader} -vL -o human_ancestor.fa.gz.gzi \
	${HUMAN_ANCESTOR_FA}.gz.gzi 2>&1 | tee -a curl_loftee_resources.log
  ${downloader} -vL -o phylocsf_gerp.sql.gz \
	${PHYLOCSF_GERP} 2>&1 | tee -a curl_loftee_resources.log
fi

# Unzip the downloaded phyloCSF GERP file
gunzip -f phylocsf_gerp.sql.gz
# Check if the file was downloaded successfully
if [ ! -f phylocsf_gerp.sql ]; then
  echo "Error: phylocsf_gerp.sql was not downloaded successfully."
  exit 1
fi

# Set up MaxEntScan resources
cd ${SCRIPT_DIR}
# mkdir -p plugin_resources/maxentscan && cd $_




