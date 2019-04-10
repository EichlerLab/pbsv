"""
Call variants from alignments with pbsv.
"""

#
# Definitions
#

def _pbsv_call_get_chrom_list():
    """
    Get a list of chromosomes from the reference file.
    """

    # Check for FAI
    fai_file = '{}.fai'.format(config['reference'])

    if not os.path.isfile(fai_file):
        raise RuntimeError('Cannot find FAI for reference: {}'.format(fai_file))

    # Get a list of chromosomes
    with open(fai_file, 'r') as in_file:
        chrom_list = [line.split('\t', 1)[0].strip() for line in in_file]

    chrom_list = [chrom for chrom in chrom_list if chrom]  # Remove empty elements

    if not chrom_list:
        raise RuntimeError('Found 0 chromosomes in FAI file: {}'.format(fai_file))

    # Return list
    return chrom_list

def _pbsv_get_all_svsig_files(wildcards):
    """
    Get a list of all svsig files. Searches wildcards for "sample" and "chrom". For each of these
    wildcards that is missing, each possible value is filled in and one file is returned for all combinations
    of these. For wildcards that are defined, only svsig files matching that wildcard are returned.

    For example, to get svsig files for all movies on chromosome "chr1", then only `wildcards.chrom` should be defined;
    `wildcards.sample` would not be defined, and all possible values for the sample would automatically
    be filled in.

    :param wildcards: Snakemake wildcards.

    :return: A list of svsig.gz files.
    """

    # Input file pattern
    svsig_file_pattern = 'call/svsig/{sample}/discover_{chrom}.svsig.gz'

    # Get subread table
    df_subread = pd.read_table('init/subread_table.tab', header=0)

    # Get list of chromosomes
    chrom = wildcards.get('chrom', None)

    if chrom:
        chrom_list = [chrom]
    else:
        chrom_list = _pbsv_call_get_chrom_list()

    # Get for a single sample if in wildcards, otherwise, get a list of samples
    sample = wildcards.get('sample', None)

    if sample:
        sample_list = [sample]
    else:
        sample_list = list(SAMPLE_TABLE['SAMPLE'])

    # Init output file list
    svsig_file_list = list()

    # Process each sample
    for sample in sample_list:

        # Add to svsig_file_list
        for chrom in chrom_list:
            svsig_file_list.append(svsig_file_pattern.format(
                **{
                    'sample': sample,
                    'chrom': chrom
                }
            ))

    # Return
    return svsig_file_list


#############
### Rules ###
#############


#
# Per sample callset
#

# pbsv_call_sample_final_bnd
#
# Make final BND output file.
#
# Note: awk command is a workaround for version 1.2 and may not be needed for newer versions.
# See: https://github.com/PacificBiosciences/pbbioconda/issues/60
rule pbsv_call_sample_final_bnd:
    input:
        vcf='temp/call/unmerged_vcf/sample_{sample}/bnd/variants_all.vcf'
    output:
        calls='pbsv_sample_{sample}_bnd.vcf.gz'
    shell:
        """awk -vOFS="\\t" '($1 !~ /^#/) {{gsub(",", ";", $7)}} {{print}}' {input.vcf} | """
        """bcftools view -O z -o {output.calls}; """
        """tabix {output.calls}"""

# pbsv_call_sample_final_vcf
#
# Merge and sort SV calls.
#
# Note: sed and awk command is a workaround for version 1.2 and may not be needed for newer versions.
# See: https://github.com/PacificBiosciences/pbbioconda/issues/60
rule pbsv_call_sample_final_vcf:
    input:
        vcf=expand('temp/call/unmerged_vcf/sample_{{sample}}/sv/variants_{chrom}.vcf', chrom=_pbsv_call_get_chrom_list())
    output:
        calls='pbsv_sample_{sample}_sv.vcf.gz'
    shell:
        """bcftools concat -O v {input.vcf} | """
        """awk -vOFS="\\t" '($1 !~ /^#/) {{gsub(",", ";", $7)}} {{print}}' | """
        """sed '12i##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise SV breakpoint">' | """
        """bcftools sort -O z -o {output.calls}; """
        """tabix {output.calls}"""

# pbsv_call_sample_unmerged_bnd
#
# Make initial SV calls (INS, DEL, and INV) for one chromosome.
rule pbsv_call_sample_unmerged_bnd:
    input:
        ref_fa=config['reference'],
        svsig=_pbsv_get_all_svsig_files
    output:
        vcf=temp('temp/call/unmerged_vcf/sample_{sample}/bnd/variants_all.vcf')
    shell:
        """pbsv call --types BND -j 4 {input.ref_fa} {input.svsig} {output.vcf}"""

# pbsv_call_sample_unmerged_sv
#
# Make initial SV calls (INS, DEL, and INV) for one chromosome.
rule pbsv_call_sample_unmerged_sv:
    input:
        ref_fa=config['reference'],
        svsig=_pbsv_get_all_svsig_files
    output:
        vcf=temp('temp/call/unmerged_vcf/sample_{sample}/sv/variants_{chrom}.vcf')
    shell:
        """pbsv call --types INS,DEL,INV -j 4 {input.ref_fa} {input.svsig} {output.vcf}"""


#
# Joint callset
#

# pbsv_call_joint_final_bnd
#
# Make final BND output file.
#
# Note: awk command is a workaround for version 1.2 and may not be needed for newer versions.
# See: https://github.com/PacificBiosciences/pbbioconda/issues/60
rule pbsv_call_joint_final_bnd:
    input:
        vcf='temp/call/unmerged_vcf/joint_all/bnd/variants_all.vcf'
    output:
        calls='pbsv_joint_all_bnd.vcf.gz'
    shell:
        """awk -vOFS="\\t" '($1 !~ /^#/) {{gsub(",", ";", $7)}} {{print}}' {input.vcf} | """
        """bcftools view -O z -o {output.calls}; """
        """tabix {output.calls}"""

# pbsv_call_joint_merge_sv
#
# Merge and sort SV calls.
#
# Note: sed and awk command is a workaround for version 1.2 and may not be needed for newer versions.
# See: https://github.com/PacificBiosciences/pbbioconda/issues/60
rule pbsv_call_joint_merge_sv:
    input:
        vcf=expand('temp/call/unmerged_vcf/joint_all/sv/variants_{chrom}.vcf', chrom=_pbsv_call_get_chrom_list())
    output:
        calls='pbsv_joint_all_sv.vcf.gz'
    shell:
        """bcftools concat -O v {input.vcf} | """
        """awk -vOFS="\\t" '($1 !~ /^#/) {{gsub(",", ";", $7)}} {{print}}' | """
        """sed '12i##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise SV breakpoint">' | """
        """bcftools sort -O z -o {output.calls}; """
        """tabix {output.calls}"""

# pbsv_call_joint_unmerged_bnd
#
# Make initial SV calls (INS, DEL, and INV) for one chromosome.
rule pbsv_call_joint_unmerged_bnd:
    input:
        ref_fa=config['reference'],
        svsig=_pbsv_get_all_svsig_files
    output:
        vcf=temp('temp/call/unmerged_vcf/joint_all/bnd/variants_all.vcf')
    shell:
        """pbsv call --types BND -j 4 {input.ref_fa} {input.svsig} {output.vcf}"""

# pbsv_call_joint_unmerged_sv
#
# Make initial SV calls (INS, DEL, and INV) for one chromosome.
rule pbsv_call_joint_unmerged_sv:
    input:
        ref_fa=config['reference'],
        svsig=_pbsv_get_all_svsig_files
    output:
        vcf=temp('temp/call/unmerged_vcf/joint_all/sv/variants_{chrom}.vcf')
    shell:
        """pbsv call --types INS,DEL,INV -j 4 {input.ref_fa} {input.svsig} {output.vcf}"""


#
# Discover
#

# pbsv_call_discover
#
# Discover SV signals in alignments.
rule pbsv_call_discover:
    input:
        bam='align/{sample}/mapped_reads.bam'
    output:
        svsig='call/svsig/{sample}/discover_{chrom}.svsig.gz'
    benchmark:
        'call/svsig/{sample}/bm/discover_{chrom}.tab'
    run:

        # Get tandem repeat bed, if configured
        tandem_bed = config.get('tandem_bed', None)

        if tandem_bed is not None:
            if not os.path.isfile(tandem_bed):
                raise RuntimeError(
                    'Cannot find configured tandem repeat BED file ("tandem_bed" in config.json): {}'.format(tandem_bed)
                )

            tandem_opt = '--tandem-repeats {} '.format(tandem_bed)
        else:
            tandem_opt = ' '

        # Run PBSV discover
        shell(
            """pbsv discover {tandem_opt}--region {wildcards.chrom} {input.bam} {output.svsig}"""
        )
