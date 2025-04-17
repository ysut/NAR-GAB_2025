# NAR_2025
Splicing single‑nucleotide variants (SNVs) with a **PriorityScore** (range –10 to 14). Values ≥1 are considered “screening‑positive.” This repository underpins our Nucleic Acids Research (NAR) 2025 submission.

## Features

- **In‑place VCF annotation**  
  Reads an existing VCF, looks up per‑variant PriorityScores in a pandas DataFrame, and writes them into a new `INFO=PriorityScore` field without dropping any original columns or samples.

- **High performance**  
  - Uses a one‑time build of an in‑memory dict for O(1) lookups  
  - Optional parallel DataFrame processing via `pandarallel`

- **VCF compatibility**  
  - Adds a proper `##INFO=<ID=PriorityScore,Number=1,Type=Integer,…>` header  
  - Omits the tag entirely for missing scores (downstream `vcfanno` can `fill(0)`)

- **Extensible workflow**  
  - Dynamic handling of new contigs (two‑pass header collection)  
  - Easily integrated into larger Nextflow or Snakemake pipelines

## Installation

```bash
git clone https://github.com/YourUser/NAR_2025.git
cd NAR_2025
pip install -r requirements.txt