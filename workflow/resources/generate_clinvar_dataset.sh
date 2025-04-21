#!/bin/bash

# Check if assembly argument is provided
if [ "$1" = "" ]
then
    echo -e "No argument about assembly version. Please add GRCh37 or GRCh38. \n"
    echo -e "Usage: "
    echo -e "        $0 [GRCh37|GRCh38] \n"
    exit 1
fi

set -euo pipefail

assembly="$1"
output_fn="clinvar_${assembly}.germline.nocoflicted.bcf.gz"
VCF="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly}/clinvar.vcf.gz"
MD5="${VCF}.md5"
TBI="${VCF}.tbi"

# Function to download files
function download_files() {
  echo "Downloading ..."

  # Check if wget is installed, if not, use curl
  if command -v wget &> /dev/null; then
    downloader="wget"
  elif command -v curl &> /dev/null; then
    downloader="curl"
  else
    echo "Neither wget nor curl is installed. Aborting."
    exit 1
  fi

  # Download VCF file
  if [ "$downloader" = "wget" ]; then
    ${downloader} --timestamping --append-output=./wget_clnvcf.log --show-progress ${VCF}
  else
    ${downloader} -vL -o clinvar.vcf.gz ${VCF} 2>&1 | tee -a curl_clnvcf.log
  fi
  
  md5sum clinvar.vcf.gz > downloaded_vcf.md5
  sleep 2

  # Download MD5 file
  if [ "$downloader" = "wget" ]; then
    ${downloader} --timestamping --append-output=./wget_clnvcf.log --show-progress ${MD5}
  else
    ${downloader} -vL -o clinvar.vcf.md5 ${MD5} 2>&1 | tee -a curl_clnvcf.log
  fi
  sleep 2

  # Download TBI file
  if [ "$downloader" = "wget" ]; then
    ${downloader} --timestamping --append-output=./wget_clnvcf.log --show-progress ${TBI}
  else
    ${downloader} -vL -o clinvar.vcf.gz.tbi ${TBI} 2>&1 | tee -a curl_clnvcf.log
  fi
}

# Function to check the MD5 checksum of the downloaded VCF file
function check_md5() {
  echo "Checking MD5 ..."
  if md5sum --check downloaded_vcf.md5 | grep -Ev 'OK$' >/dev/null; then
    echo "The downloaded file may be corrupted."
    exit 2
  else
    echo "md5sum check -> OK!"
  fi
}

# Function to filter and sort the VCF file
function filter_sort_bcf() {
  echo "Filtering ..."
  bcftools view \
    --output-type u \
    --include \
      'INFO/ORIGIN=="1"
        & (INFO/CLNSIG=="Pathogenic"
            | INFO/CLNSIG=="Likely_pathogenic"
            | INFO/CLNSIG=="Pathogenic/Likely_pathogenic"
            | INFO/CLNSIG=="Benign" 
            | INFO/CLNSIG=="Likely_benign" 
            | INFO/CLNSIG=="Benign/Likely_benign") 
        & INFO/CLNREVSTAT!="no_assertion_criteria_provided" 
        & INFO/CLNREVSTAT!="no_classification_provided" 
        & INFO/CLNREVSTAT!="no_classification_for_the_individual_variant"' \
    clinvar.vcf.gz \
  | bcftools sort \
      --output-type b \
      --output ${output_fn} \
      --write-index
}

# Function to move files to an output directory
function moving_files() {
  echo "Moving files ..."
  local datetime=$(date +"%Y%m%d-%I%M%S")
  output_dir="Filtered_BCF_${assembly}_${datetime}"
  mkdir -p ${output_dir}
  mv ${output_fn}* ./${output_dir}
  mv wget_clnvcf.log ./${output_dir}
  # Remove downloaded files and MD5
  rm -rf clinvar.vcf.gz* *.md5
}

# Main function
function main() {
  download_files
  check_md5
  filter_sort_bcf
  moving_files
  echo -e "Completed!"
}

# Run the main function
main

