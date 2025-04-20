# NAR_2025
A scoring framework for splicing single‑nucleotide variants (SNVs) with a **PriorityScore** (range –10 to 14). Values ≥1 are considered “screening‑positive.” This repository underpins our Nucleic Acids Research (NAR) 2025 submission.

## Description



## Requirements
- **Docker** (tested with version 27.2.1)  
- **Nextflow** (tested with version 24.10.5.5935)
- **bcftools** (tested with version bcftools 1.21, using htslib 1.21)

### For set up
- **wget**
- **git**

We tested this workflow script on Intel Mac (MacOS version 15.3.2) with Inel CPU and Ubuntu (version 24.04.2).


## Installation and try our framework
### STEP 1. Clone this repository
```bash
git clone https://github.com/ysut/NAR_2025.git
```

### STEP 2. Setup docker images
This step takes time.
```bash
# Move to docker directry in this repository
cd NAR_2025/workflow/docker

# Run setup shell script for build docker images
./build_docker_images.sh
```

### STEP 3. Download VEP resources using docker images built above
You can change the `vep_data` directry path where you like.
Please run the command below and follow the comments displayed.
```bash
docker run --rm -v ${HOME}/vep_data:/data ps_vep:113.4 INSTALL.pl
```

### STEP 3. Download LOFTEE resources using a custom script

### STEP 4. Download and filter the ClinVar's variant data
```bash
cd resources && ./do
```

### Let's try a framework
```bash
nextflow main.nf --input_vcf ../examples/example.vcf.gz --output_dir ../examples --resources resources
```
