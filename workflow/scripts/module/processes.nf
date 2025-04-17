process SPLICEAI {
    input:
    tuple path(input_vcf), path(input_tbi), 
          path(reference_fasta), path(annotation_gtf)

    output:
    tuple path('*.splai.vcf'), path(reference_fasta)

    script:
    """
    /bin/bash -c " \\
    source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate spliceai && \\
    spliceai \\
      -I ${input_vcf} \\
      -O ${input_vcf.baseName}.splai.vcf \\
      -R ${reference_fasta} \\
      -A ${annotation_gtf} \\
      -D 4999 \\
      -M 0
    "
    """
}

process VEP {
    input: 
    tuple path(input_vcf), path(reference_fasta)
    
    output:
    path('*.splai.vep.vcf')

    script:
    """
    /opt/vep/src/ensembl-vep/vep \\
      --dir_cache /data \\
      --cache \\
      --offline \\
      --no_stats \\
      --merged \\
      --gencode_basic \\
      --variant_class \\
      --canonical \\
      --symbol \\
      --numbers \\
      --vcf \\
      --pick_allele \\
      --force_overwrite \\
      --use_given_ref \\
      --assembly ${params.assembly} \\
      --fasta ${reference_fasta} \\
      --plugin MaxEntScan,/plugin_resources/maxentscan/fordownload \\
      --plugin LoF,loftee_path:/vep/loftee,human_ancestor_fa:/plugin_resources/loftee/human_ancestor.fa.gz,conservation_file:/plugin_resources/loftee/phylocsf_gerp.sql \\
      -i ${input_vcf} \\
      -o ${input_vcf.baseName}.vep.vcf
    """
}

process PS {
    publishDir "${params.out_root}", mode: 'copy', overwrite: true

    input:
    path(input_vcf)

    output:
    path "*.psscored.vcf"

    script:
    """
    bash -c "
      source /opt/conda/etc/profile.d/conda.sh && \\
      conda activate psscoring && \\
      /opt/psscoring/ps.py \\
        --input ${input_vcf} \\
        --output ${input_vcf.baseName}.psscored.vcf \\
        --resources /ps_resources \\
        --verbose
    "
    """
}