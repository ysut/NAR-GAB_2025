import re

import gffutils
import pandas as pd
from cyvcf2 import VCF

from . import posparser
# from .deco import print_filtering_count


def parse_vcf(raw_vcf: str, db: gffutils.interface.FeatureDB) -> pd.DataFrame:
    """Parse VCF file and extract relevant information
    Args:
        raw_vcf (str): Path to the VCF file
        db (str): Path to the GTF database
    Returns:
        pd.DataFrame: DataFrame containing parsed VCF information
    """

    header = VCF(raw_vcf).header_iter()
    for h in header:
        try:
            h['ID']
        except KeyError:
            continue
        else:
            if h['ID'] == 'CSQ':
                vep_cols_list = h['Description'].split('Format: ')[1].rstrip('"').split('|')
            elif h['ID'] == 'SpliceAI':
                splai_cols_list = h['Description'].split('Format: ')[1].rstrip('"').split('|')
            else:
                pass

    cols = [
        'CHROM', 'POS', 'REF', 'ALT', 'GeneSymbol', 'SymbolSource', 'HGNC_ID', 
        'ENST', 'HGVSc', 'Consequence', 'EXON', 'INTRON', 'Strand',
        'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 
        'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'maxsplai', 'loftee' 
    ]
    vepidx: dict = {col: i for i, col in enumerate(vep_cols_list)}
    splaidx: dict = {col: i for i, col in enumerate(splai_cols_list)}
    df = pd.DataFrame(columns=cols)

    for v in VCF(raw_vcf):
        vep: list = v.INFO.get('CSQ').split('|')

        # Get HGVSc from VEP
        try:
            hgvsc = re.search('(?<=:).*',vep[vepidx['HGVSc']])[0]
        except TypeError:
            hgvsc = "NA"

        # Get SpliceAI scores
        if v.INFO.get('SpliceAI'):
            splai: list = v.INFO.get('SpliceAI').split(',')[0].split('|')
        else:
            splai = ['NA'] * len(splai_cols_list)

        # Convert strand to +/- 
        strand = lambda s: '+' if s == '1' else '-'

        # Get max SpliceAI scores
        ds_ag: float = splai[splaidx['DS_AG']]
        ds_al: float = splai[splaidx['DS_AL']]
        ds_dg: float = splai[splaidx['DS_DG']]
        ds_dl: float = splai[splaidx['DS_DL']]
        if splai[splaidx['DP_AG']] == 'NA':
            maxsplai: str = "NA"
        maxsplai: float = max(ds_ag, ds_al, ds_dg, ds_dl)
    
        # Add df row
        df = pd.concat(
            [df, pd.DataFrame(
                [
                    [
                        v.CHROM, v.POS, v.REF, v.ALT[0], 
                        vep[vepidx['SYMBOL']], vep[vepidx['SYMBOL_SOURCE']], 
                        vep[vepidx['HGNC_ID']], vep[vepidx['Feature']], hgvsc, 
                        vep[vepidx['Consequence']], vep[vepidx['EXON']], 
                        vep[vepidx['INTRON']], strand(vep[vepidx['STRAND']]), 
                        ds_ag, ds_al, ds_dg, ds_dl, 
                        splai[splaidx['DP_AG']], splai[splaidx['DP_AL']], 
                        splai[splaidx['DP_DG']], splai[splaidx['DP_DL']],
                        maxsplai, vep[vepidx['LoF']], 
                    ]
                ],
                columns=cols
                )
            ], ignore_index=True)
        
    df.drop_duplicates(inplace=True)

    # Annotate full ENST IDs with GTF database
    df['ENST_Full'] = df.apply(posparser.fetch_enst_full, db=db, axis=1)
    df = df.fillna({'loftee': 'NANANANANNA'})

    return df
