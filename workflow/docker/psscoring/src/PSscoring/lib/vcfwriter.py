from cyvcf2 import VCF, Writer
import pandas as pd
from pandas import Int64Dtype


def write_vcf(df: pd.DataFrame, raw_vcf: str, output_vcf: str) -> None:
    """
    Write a VCF file with the priority score for pathogenic splicing SNVs.
    Args:
        df (pd.DataFrame): DataFrame containing the priority scores.
        raw_vcf (str): Path to the input VCF file.
        output_vcf (str): Path to the output VCF file.
    """

    # Check if the input DataFrame is empty
    if df.empty:
        raise ValueError("The input DataFrame is empty. Please provide a valid DataFrame.")
    # Check if the required columns are present in the DataFrame
    required_columns = ["CHROM", "POS", "REF", "ALT", "PriorityScore"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"The input DataFrame is missing the following required columns: {', '.join(missing_columns)}")
    # Check if the input VCF file exists
    if not isinstance(raw_vcf, str):
        raise ValueError("The input VCF file path must be a string.")

    # cast to int for the priority score
    df["PriorityScore"] = df["PriorityScore"].astype(Int64Dtype())

    mapping = {
        (row.CHROM, row.POS, row.REF, row.ALT): int(row.PriorityScore)
        for row in df.itertuples(index=False)
        if not pd.isna(row.PriorityScore)
    }

    vcf_in  = VCF(raw_vcf)
    vcf_out = Writer(output_vcf, vcf_in)
    vcf_out.add_to_header(
        '##INFO=<ID=PriorityScore,Number=1,Type=Integer,'
        'Description="Priority score for pathogenic splicing SNVs '
        '(range -10 to 14; values ≥1 are screening-positive)">'
    )

    vcf_in.add_info_to_header({
        "ID": "PriorityScore", "Number": "1", "Type": "Integer",
        "Description": "Priority score for pathogenic splicing SNVs (range -10 to 14; values ≥1 are considered screening-positive)"
    })
    vcf_out.write_header()


    for var in vcf_in:
        key = (var.CHROM, var.POS, var.REF, var.ALT[0])
        if key in mapping:
            var.INFO["PriorityScore"] = mapping[key]
        vcf_out.write_record(var)
    
    vcf_out.close()
    vcf_in.close()