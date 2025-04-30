# NAR_2025
A scoring framework for splicing single‑nucleotide variants (SNVs) with a **PriorityScore** (range –10 to 14). Values ≥1 are considered “screening‑positive.” This repository underpins our Nucleic Acids Research (NAR) 2025 submission.

## Requirements
- **Docker** (tested with version 27.2.1)  
- **Nextflow** (tested with version 24.10.5.5935)
- **bcftools** (tested with version bcftools 1.21, using htslib 1.21)

### For set up
- **wget** or **curl**
- **git**

We tested this workflow script on Intel Mac (MacOS version 15.3.2) with Inel CPU and Ubuntu (version 24.04.2).

## Installation and try our framework
### STEP 1. Clone this repository
```bash
git clone https://github.com/ysut/NAR_2025.git
wd=$(pwd)
wf="${wd}/NAR_2025/workflow"           # Set the path to workflow dictory
rc="${wd}/NAR_2025/workflow/resources"  # Set the path to resources directory 
```

### STEP 2. Setup docker images
This step takes time.
```bash
# Move to docker directry in this repository
cd ${wf}/docker

# Run setup shell script for build docker images
./build_docker_images.sh
```

### STEP 3. Download VEP resources using docker images built above
Run the following command and follow the on‑screen instructions.  
Then, please select a number for `homo_sapiens_merged_vep_113_GRCh37.tar.gz`. 
```bash
docker run -it --rm -v ${HOME}/vep_data:/data ps_vep:113.4 INSTALL.pl
```
💡 The path of `${HOME}/vep_data` is an example.  
You can change the `vep_data` directory path to any location you prefer.


### STEP 4. Download LOFTEE resources
```bash
cd ${rc}
./setup_vep_resources.sh
```

### STEP 5. Download and filter the ClinVar's variant data
```bash
cd ${rc}
./generate_clinvar_dataset.sh
```

### STEP 6. Edit the nextflow.config file
Please edit five paths, i.e., `reference`, `annotation_gtf`, `vep_data`, `vep_plugin_resources`, and `ps_resources`.
A `nextflow.config` file is located in `NAR_2025/workflow/scripts`.

### Let's try a framework
```bash
nextflow ${wf}/scripts/main.nf --input_vcf ${wf}/examples/example.vcf.gz --output_dir ${wf}/examples
