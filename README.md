# NAR_2025
Splicing single‑nucleotide variants (SNVs) with a **PriorityScore** (range –10 to 14). Values ≥1 are considered “screening‑positive.” This repository underpins our Nucleic Acids Research (NAR) 2025 submission.

## Description



## Requirements
- **Docker** (tested with version 27.2.1)  
- **Nextflow** (tested with version 24.10.5.5935)
- **wget


## Installation and try our framework

```bash
git clone https://github.com/ysut/NAR_2025.git
```

#### STEP 1. Setup docker images
This step takes time.
```bash
NAR_2025/workflow/docker/build_docker_images.sh
```

#### STEP 2. Download VEP resources using docker images built above
You can change the `vep_data` directry path where you like.
Please run the command below and follow the comments displayed.
```bash
docker run --rm -v vep_data:/data ps_vep:113.4 INSTALL.pl
```

#### STEP 3. Download LOFTEE resources using a custom script.


#### STEP 4. 


#### Try framework
```bash
nextflow main.nf --input_vcf ../examples/example.vcf.gz --output_dir ../examples
```