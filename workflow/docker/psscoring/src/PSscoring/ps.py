#!/usr/bin/env python

import subprocess
import glob
import os
import re

import numpy as np
import pandas as pd
from pathlib2 import Path
from pandarallel import pandarallel
import gffutils
import pysam
from logging import getLogger, config
import yaml

from lib import posparser, splaiparser, predeffect, anno_clinvar
from lib.preprocess import parse_vcf
from lib.scoring import Scoring
from lib.vcfwriter import write_vcf


#===============================================================================
# Functions 
#===============================================================================
def set_gtf_db(db_list: list, gff_list: list) -> tuple:
    if len(db_list) == 0 or len(gff_list) == 0:
        subprocess.run(
            ["python", "/opt/psscoring/lib/generatedbs.py", "--output_dir", FLAGS.resources, 
            "--release", FLAGS.release, "--assembly", FLAGS.assembly], 
            shell=False, check=True, stdin=subprocess.DEVNULL)

    if FLAGS.assembly == 'GRCh37':
        db_anno_gencode = f"{FLAGS.resources}/gencode.v{FLAGS.release}lift37.annotation.gtf.db"
        db_anno_intron = f"{FLAGS.resources}/gencode.v{FLAGS.release}lift37.annotation.intron.gtf.db"
        gencode_gff = f"{FLAGS.resources}/gencode.v{FLAGS.release}lift37.annotation.gff3.gz"
    elif FLAGS.assembly == 'GRCh38':
        db_anno_gencode = f"{FLAGS.resources}/gencode.v{FLAGS.release}.annotation.gtf.db"
        db_anno_intron = f"{FLAGS.resources}/gencode.v{FLAGS.release}.annotation.intron.gtf.db"
        gencode_gff = f"{FLAGS.resources}/gencode.v{FLAGS.release}.annotation.gff3.gz"
    else:
        raise ValueError("Assembly must be either 'GRCh37' or 'GRCh38'.")
    
    return gffutils.FeatureDB(db_anno_gencode), gffutils.FeatureDB(db_anno_intron), gencode_gff

def setup_logging(output_vcf: str, verbose: bool):
    config_path = '/opt/psscoring/logging.yaml'
    with open(config_path, 'r') as f:
        log_cfg = yaml.safe_load(f)

    out_dir = os.path.dirname(os.path.abspath(output_vcf))
    base = os.path.splitext(os.path.basename(output_vcf))[0]
    log_file = os.path.join(out_dir, f"{base}.log")
    log_cfg['handlers']['file']['filename'] = log_file

    if verbose:
        log_cfg['loggers']['__main__']['handlers'] = ['console', 'file']
        log_cfg['loggers']['__main__']['level']    = 'DEBUG'
        log_cfg['handlers']['console']['level']     = 'DEBUG'
    else:
        log_cfg['loggers']['__main__']['handlers'] = ['file']

    config.dictConfig(log_cfg)


def map_and_calc_score(row, score_map: dict) -> int:
    """
    PriortiyScore is the sum of the "clinvar_screening", "insilico_screening", and "recalibrated_splai"
    """
    if row['insilico_screening'] == "Not available":
        return np.nan

    return int(score_map[row['recalibrated_splai']]) + int(score_map[row['insilico_screening']]) + int(score_map[row['clinvar_screening']])

#===============================================================================
# Arugments parser using absl-py 
#===============================================================================
from absl import app
from absl import flags

FLAGS = flags.FLAGS
flags.DEFINE_string(
    'input', None, 'Path to input VCF file', short_name='i')
flags.DEFINE_string(
    'output', None, 'Path to output VCF file', short_name='o')
flags.DEFINE_string(
    'resources', None, 'Path to resources directory', short_name='r')
flags.DEFINE_string(
    'release', '43', 'Release version (e.g., 43)', short_name='R')
flags.DEFINE_string(
    'assembly', 'GRCh37', 'Assembly version (GRCh37 or GRCh38)', short_name='a')
flags.DEFINE_boolean(
    'raw_tsv', False, 'Output raw TSV file')
flags.DEFINE_boolean(
    'verbose', False, 'Verbose logging')

flags.DEFINE_integer(
    'n_workers', 2, 'Number of workers for parallel processing in pandas')
flags.DEFINE_float(
    'min_score_aldl', 0.02, 'Minimum SpliceAI score for AL or DL')
flags.DEFINE_float(
    'max_score_aldl', 0.2, 'Maximum SpliceAI score for AL or DL')
flags.DEFINE_float(
    'min_score_agdg', 0.01, 'Minimum SpliceAI score for AG or DG')
flags.DEFINE_float(
    'max_score_agdg', 0.05, 'Maximum SpliceAI score for AG or DG')
flags.DEFINE_integer(
    'min_gain_exon_len', 25, 'Minimum length of gained exon')
flags.DEFINE_integer(
    'max_gain_exon_len', 500, 'Maximum length of gained exon')
flags.DEFINE_float(
    'activation_score_ag', 0.2, 'Activation score for AG')
flags.DEFINE_float(
    'activation_score_dg', 0.2, 'Activation score for DG')


#===============================================================================
# Main function
#===============================================================================
def main(argv):
    del argv  # Unused.
    setup_logging(FLAGS.output, FLAGS.verbose)
    logger = getLogger(__name__)
    os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 
    pandarallel.initialize(nb_workers=FLAGS.n_workers, 
                           progress_bar=False, verbose=1, use_memory_fs=False
                           ) 

    # Display the input arguments
    logger.info(f"""
                Input args
                ----------------
                Input VCF    : {FLAGS.input}
                Output VCF   : {FLAGS.output}
                Resources dir: {FLAGS.resources}
                """)

    thresholds_SpliceAI_parser: dict = {
        'TH_min_sALDL': FLAGS.min_score_aldl, 
        'TH_max_sALDL': FLAGS.max_score_aldl, 
        'TH_min_sAGDG': FLAGS.min_score_agdg, 
        'TH_max_sAGDG': FLAGS.max_score_agdg,
        'TH_min_GExon': FLAGS.min_gain_exon_len, 
        'TH_max_GExon': FLAGS.max_gain_exon_len,
        'TH_sAG': FLAGS.activation_score_ag, 
        'TH_sDG': FLAGS.activation_score_dg
    }

    # raw_vcf: str = FLAGS.input
    fp = Path(FLAGS.input)
    fp_stem, fp_dir = fp.stem, fp.parent

    ## eLoF genes list (only HGNC IDs) 
    elof_path = "/opt/psscoring/eLoF_genes.tsv"
    elofs = pd.read_table(
        elof_path, usecols=['HGNC_ID'], sep='\t')
    elofs_hgnc_ids_with_prefix = elofs['HGNC_ID'].unique().tolist()
    elofs_hgnc_ids = [re.sub('HGNC:', '', hgnc) for hgnc in elofs_hgnc_ids_with_prefix]
    
    # Find gencode GTF file databases (gencode.*.annotation.gtf.db) in resources directory.
    db_list = glob.glob(f"{FLAGS.resources}/gencode.*.annotation.gtf.db")
    gff_list = glob.glob(f"{FLAGS.resources}/gencode.*.annotation.gff3.gz")
    db, db_intron, gencode_gff = set_gtf_db(db_list, gff_list)
    
    # Generate gff3 index file using pysam
    tbi_path = f"{gencode_gff}.tbi"
    if not os.path.exists(tbi_path):
        logger.info("Re-compressing and sorting GFF3 for BGZF+Tabix...")
        sorted_bgz = f"{gencode_gff}.sorted.gz"
        cmd = (
            f"gunzip -c {gencode_gff} | "
            f"sort -k1,1 -k4,4n | "
            f"bgzip -c > {sorted_bgz}"
        )
        subprocess.run(cmd, shell=True, check=True)
        os.replace(sorted_bgz, gencode_gff)

        logger.info("Indexing sorted BGZF-compressed GFF3 with tabix...")
        subprocess.run(
            ["tabix", "-f", "-p", "gff", gencode_gff],
            check=True
        )
        logger.info("Tabix index created successfully.")

    # Find a processed ClinVar bcf file in resources directory.
    clinvar_file_list = glob.glob(f"{FLAGS.resources}/Filtered_BCF_{FLAGS.assembly}_*-*/clinvar_{FLAGS.assembly}.germline.nocoflicted.bcf.gz")
    clinvar_file_index = glob.glob(f"{FLAGS.resources}/Filtered_BCF_{FLAGS.assembly}_*-*/clinvar_{FLAGS.assembly}.germline.nocoflicted.bcf.gz.*i")
    if len(clinvar_file_list) == 0 or len(clinvar_file_index) == 0:
        raise FileNotFoundError(
            f"Cannot find the processed ClinVar bcf file in {FLAGS.resources} directory. "
            f"Please check the directory and try again."
            f"You can generate it using the 'ss_generate_clinvar_dataset.sh' script.")
    else:
        clinvar_file = clinvar_file_list[0]
        clinvar_file_index = clinvar_file_index[0]
        logger.debug("ClinVar bcf file: %s", clinvar_file)
        logger.debug("ClinVar bcf index file: %s", clinvar_file_index)
    
    # Find CCRs file in resources directory.
    ccrs_auto_file_list = glob.glob(f"{FLAGS.resources}/ccrs.autosomes.*.bed.gz")
    ccrs_x_file_list = glob.glob(f"{FLAGS.resources}/ccrs.xchrom.*.bed.gz")
    if len(ccrs_auto_file_list) == 0 or len(ccrs_x_file_list) == 0:
        subprocess.run(['/opt/psscoring/dlccrs.sh', FLAGS.resources], 
                       shell=False, check=True, stdin=subprocess.DEVNULL)
        ccrs_auto = glob.glob(f"{FLAGS.resources}/ccrs.autosomes.*.bed.gz")[0]
        ccrs_x = glob.glob(f"{FLAGS.resources}/ccrs.xchrom.*.bed.gz")[0]
    else:
        ccrs_auto = ccrs_auto_file_list[0]
        ccrs_x = ccrs_x_file_list[0]

    ## Convert to pandas DataFrame from a input VCF file
    df = parse_vcf(raw_vcf=FLAGS.input, db=db)

    logger.info('Calculate the distance to the nearest splice site in intron variant...')
    df['IntronDist'] = df.apply(
        posparser.signed_distance_to_exon_boundary, 
        db=db, db_intron=db_intron, axis=1)

    logger.info('Classify "Canonical" splice site or "Non-canonical" splice site...')
    df = posparser.classifying_canonical(df)

    df['Ex_or_Int'] = np.where(
        df['IntronDist'] == "[Warning] Invalid ENST ID", "[Warning] Invalid ENST ID",
        np.where(df['IntronDist'].isnull(), 'Exonic', 'Intronic'))

    tbx_anno = pysam.TabixFile(gencode_gff)
    df['exon_loc'] = df.apply(
        posparser.calc_exon_loc, tabixfile=tbx_anno, enstcolname='ENST', axis=1)
    df = pd.concat([df, df['exon_loc'].str.split(':', expand=True)], axis=1)
    df.rename(columns={0: 'ex_up_dist', 1: 'ex_down_dist'}, inplace=True)
    df.drop(columns=['exon_loc'], inplace=True)

    #2-2. Select minimum distance from upstream distance and downstream distance
    df['exon_pos'] = df.parallel_apply(posparser.select_exon_pos, axis=1)
    #2-3. Relative exon location
    df['prc_exon_loc'] = df.parallel_apply(posparser.calc_prc_exon_loc, axis=1)

    #2-4. Decision exonic splice sites (1 nt in acceptor site or 3 nts on Donor site)
    df['exon_splice_site'] = df.parallel_apply(posparser.extract_splicing_region, axis=1)

    #3.   Additional Splicing information
    logger.info('Annotating splicing information...')
    #3-1. Annotate splicing type ('Exonic Acceptor' etc.)
    df['SpliceType'] = df.parallel_apply(posparser.select_donor_acceptor, axis=1)

    #5.   Annotate ClinVar varaints interpretations
    logger.info('Annotating ClinVar varaints interpretations...')
    # clinvar_file = '../../../clinvar/Filtered_BCF_GRCh37_20241211-044124/clinvar_GRCh37.germline.nocoflicted.bcf.gz'
    cln_bcf = pysam.VariantFile(clinvar_file)
    df['clinvar_same_pos'] = df.apply(
        anno_clinvar.anno_same_pos_vars, cln_bcf=cln_bcf, axis=1)
    df['clinvar_same_motif'] = df.apply(
        anno_clinvar.anno_same_motif_vars, cln_bcf=cln_bcf, axis=1)
    df['same_motif_clinsigs'] = df['clinvar_same_motif'].parallel_apply(
        anno_clinvar.extract_same_motif_clinsigs)

    logger.info('Parsing SpliceAI results...')
    logger.info('Annotating Exon/Intron position information...')
    df['ExInt_INFO'] = df.apply(
        splaiparser.calc_exint_info, db=db, db_intron=db_intron, axis=1)

    #6-3. Predict splicing effects
    df['Pseudoexon'] = df.apply(
        splaiparser.pseudoexon_activation,
        thresholds=thresholds_SpliceAI_parser, 
        db_intron=db_intron,
        axis=1)

    df['Part_IntRet'] = df.parallel_apply(
        splaiparser.partial_intron_retention,
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    df['Part_ExDel'] = df.parallel_apply(
        splaiparser.partial_exon_deletion,
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    df['Exon_skipping'] = df.parallel_apply(
        splaiparser.exon_skipping, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)
                                            
    df['Int_Retention'] = df.parallel_apply(
        splaiparser.intron_retention, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    df['multiexs'] = df.parallel_apply(
        splaiparser.multi_exon_skipping, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    #7.   Annotate aberrant splicing size (bp)
    logger.info('Annotating aberrant splicing size (bp)...')
    #7-1. Annotate size of 
    df['Size_Part_ExDel'] = df.parallel_apply(
        splaiparser.anno_partial_exon_del_size, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    #7-3. Annotate size of partial intron retention
    df['Size_Part_IntRet'] = df.parallel_apply(
        splaiparser.anno_partial_intron_retention_size, 
        thresholds=thresholds_SpliceAI_parser,
        axis=1)

    #7-2. Annotate size of pseudoexon
    df['Size_pseudoexon'] = df.parallel_apply(
        splaiparser.anno_gained_exon_size, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    #7-4. Annotate size of intron retention
    df['Size_IntRet'] = df.parallel_apply(
        splaiparser.anno_intron_retention_size, 
        thresholds=thresholds_SpliceAI_parser,
        axis=1)

    #7-5. Annotate size of exon skipping
    df['Size_skipped_exon'] = df.parallel_apply(
        splaiparser.anno_skipped_exon_size, 
        thresholds=thresholds_SpliceAI_parser,
        axis=1)

    df['variant_id'] = df['CHROM'].astype(str) + '-' \
        + df['POS'].astype(str) + '-' + df['REF'] + '-' + df['ALT']

    #8.   Evaluate splicing effects
    logger.info('Predicting CDS change...')
    #8-1. Predict CDS change
    df['CDS_Length'] = df.apply(predeffect.calc_cds_len, db=db, axis=1)
    df['is_10%_truncation'] = df.apply(predeffect.calc_cds_len_shorten, axis=1)

    #8-2. Determine if the gene is included in eLoFs genes
    df['is_eLoF'] = df.parallel_apply(
        predeffect.elofs_judge, elofs_hgnc_ids=elofs_hgnc_ids, axis=1
        )

    #8-3. Determine causing NMD or not
    df['is_NMD_at_Canon'] = df.parallel_apply(predeffect.nmd_judge, axis=1)

    cannot_predict: str = 'Cannot predict splicing event'
    df['Size_Part_ExDel'] = df['Size_Part_ExDel'].replace(cannot_predict, np.nan)
    df['Size_Part_IntRet'] = df['Size_Part_IntRet'].replace(cannot_predict, np.nan)
    df['Size_pseudoexon'] = df['Size_pseudoexon'].replace(cannot_predict, np.nan)
    df['Size_IntRet'] = df['Size_IntRet'].replace(cannot_predict, np.nan)
    df['Size_skipped_exon'] = df['Size_skipped_exon'].replace(cannot_predict, np.nan)

    df['is_Frameshift_Part_ExDel'] = df['Size_Part_ExDel'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_Part_IntRet'] = df['Size_Part_IntRet'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_pseudoexon'] = df['Size_pseudoexon'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_IntRet'] = df['Size_IntRet'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_skipped_exon'] = df['Size_skipped_exon'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift'] = df[['is_Frameshift_Part_ExDel', 
                            'is_Frameshift_Part_IntRet', 
                            'is_Frameshift_pseudoexon', 
                            'is_Frameshift_IntRet', 
                            'is_Frameshift_skipped_exon'
                            ]].any(axis=1)

    #9.   CCRs
    logger.info('Setting up CCRs info...')
    #9-1. Annotate truncated regions 
    df['skipped_region'] = df.parallel_apply(
        splaiparser.anno_skipped_regions, axis=1)
    df['deleted_region'] = df.parallel_apply(
        splaiparser.anno_deleted_regions, 
        thresholds=thresholds_SpliceAI_parser, axis=1)

    #9-2. Intersect with CCRs
    logger.info('Annotating CCRs score')
    df = predeffect.anno_ccr_score(df, autoccr=ccrs_auto, xccr=ccrs_x)

    # Extract data with SymbolSource == 'HGNC'
    df = df[df['SymbolSource'] == 'HGNC']

    logger.info('Scoring...')
    scoring = Scoring()
    df['insilico_screening'] = df.parallel_apply(scoring.insilico_screening, axis=1)
    df['clinvar_screening'] = df.parallel_apply(scoring.clinvar_screening, axis=1)
    df['recalibrated_splai'] = df.parallel_apply(scoring.recal_scores_in_canon, axis=1)

    solution = {'s1': 9.0, 's2': 6.0, 's3': 0.0, 's4': -5.0, 
                's5': -3.0, 's6': 0.0, 's7': 2.0, 's8': 3.0, 's9': 2.0,
                's10': 4.0, 's11': 2.0, 's12': -1.0, 's13': 0.0, 's14': 1.0, 
                's15': -5.0, 's0': 0.0}

    df['PriorityScore'] = df.parallel_apply(map_and_calc_score, args=(solution,), axis=1)
    df = df[['CHROM', 'POS', 'REF', 'ALT', 'PriorityScore']]

    logger.info('Writing VCF file...')
    write_vcf(df, FLAGS.input, FLAGS.output)

    if FLAGS.raw_tsv:
        logger.info('Saving raw TSV file...')
        df.to_csv(
            f"{fp_dir}/{fp_stem}.raw.tsv", index=False, sep='\t')

    print("Done!")

if __name__ == '__main__':
    app.run(main)
