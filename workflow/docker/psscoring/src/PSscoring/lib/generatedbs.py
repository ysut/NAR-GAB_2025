import os
import sys
import requests

from pathlib2 import Path
import gffutils
import gffutils.pybedtools_integration
from absl import app
from absl import flags
from absl import logging


FLAGS = flags.FLAGS
flags.DEFINE_string(
    'output_dir', '.', 'Path to output directory', short_name='o')
flags.DEFINE_string(
    'release', '43', 'GENCODE release version (e.g., 43)', short_name='r')
flags.DEFINE_string(
    'assembly', 'GRCh37', 'Assembly version (GRCh37 or GRCh38)', short_name='a')


def download_gencode_files(release: str, assembly: str, output_dir: str) -> Path:
    """
    Download GENCODE GTF file from the GENCODE website.
    """
    # For using this script by itself, set the base URL to the GENCODE FTP site with "https"
    # BASE_URL = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}"
    # For using this script in the Docker container, set the base URL to the GENCODE FTP site with "http"
    BASE_URL = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}"

    if assembly == 'GRCh37':
        gtf_fn, gff_fn = f"gencode.v{release}lift37.annotation.gtf.gz", f"gencode.v{release}lift37.annotation.gff3.gz"
        gtf_url, gff_url = f"{BASE_URL}/GRCh37_mapping/{gtf_fn}", f"{BASE_URL}/GRCh37_mapping/{gff_fn}"
    elif assembly == 'GRCh38':
        gtf_fn, gff_fn = f"gencode.v{release}.annotation.gtf.gz", f"gencode.v{release}.annotation.gff3.gz"
        gtf_url, gff_url = f"{BASE_URL}/{gtf_fn}", f"{BASE_URL}/{gff_fn}"
    else:
        raise ValueError("Assembly must be either 'GRCh37' or 'GRCh38'.")

    # Check if the GTF file already exists
    gtf_path, gff_path = Path(f"{output_dir}/{gtf_fn}"), Path(f"{output_dir}/{gff_fn}")
    
    if gtf_path.exists():
        logging.info(f"GTF file already exists: {gtf_path}")
    else:
        logging.info(f"Downloading GTF file from {gtf_url}")
        response = requests.get(gtf_url, stream=True)
        if response.status_code == 200:
            with open(f"{output_dir}/{gtf_fn}", 'wb') as f:
                f.write(response.content)
            logging.info(f"GTF file downloaded to {output_dir}")
        else:
            raise f"Failed to download GTF file: {response.status_code}"
    
    if gff_path.exists():
        logging.info(f"GFF file already exists: {gff_path}")
    else:
        logging.info(f"Downloading GFF3 file from {gff_url}")
        response = requests.get(gff_url, stream=True)
        if response.status_code == 200:
            with open(f"{output_dir}/{gff_fn}", 'wb') as f:
                f.write(response.content)
            logging.info(f"GFF3 file downloaded to {output_dir}")
        else:
            raise f"Failed to download GFF3 file: {response.status_code}"

    return gtf_path

def generate_intoron_gtf(db: gffutils.FeatureDB, output: str) -> None:
    introns = db.create_introns(exon_featuretype='exon', 
                                new_featuretype='intron', 
                                merge_attributes=True, 
                                numeric_sort=True)
    pybed = gffutils.pybedtools_integration.to_bedtool(introns)
    pybed.saveas(output)
    
    return None

def main(argv):
    del argv  # Unused.
    gtf_path: Path = download_gencode_files(FLAGS.release, FLAGS.assembly, FLAGS.output_dir)
    
    # Set the output file names
    gtf_base_name: str = gtf_path.name.rstrip('.gtf.gz')
    db_anno_gencode = f"{FLAGS.output_dir}/{gtf_base_name}.gtf.db"
    intron_gtf = f"{FLAGS.output_dir}/{gtf_base_name}.intron.gtf.gz"
    db_anno_intron = f"{FLAGS.output_dir}/{gtf_base_name}.intron.gtf.db"

    if not os.path.exists(db_anno_gencode):
        db = gffutils.create_db(str(gtf_path), db_anno_gencode,
                                disable_infer_genes=True, 
                                disable_infer_transcripts=True,
                                keep_order=True)
    
    if not os.path.exists(db_anno_intron):
        # Create intron information file as GTF
        db = gffutils.FeatureDB(db_anno_gencode)
        generate_intoron_gtf(db, intron_gtf)

        # Create intron DB from above GTF
        gffutils.create_db(intron_gtf, db_anno_intron,
                        disable_infer_genes=True, 
                        disable_infer_transcripts=True, 
                        keep_order=True,
                        merge_strategy="merge")

if __name__ == '__main__':
    app.run(main)