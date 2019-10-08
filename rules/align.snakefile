"""
Align cells for PBSV.
"""

def _pbsv_align_get_bams_per_sample(wildcards):
    """
    Get all BAMs for one sample.
    """

    # Set BAM file pattern
    bam_pattern = 'temp/align/{SAMPLE}/movie/movie_{MOVIE}.bam'

    # Get data table
    cell_table = pd.read_csv('init/cell_table.tab', sep='\t', header=0)

    cell_table = cell_table.loc[cell_table['SAMPLE'] == wildcards.sample]

    if cell_table.shape[0] == 0:
        raise RuntimeError('No input for sample in cell table: {}'.format(wildcards.sample))

    # Return a list of BAM files
    return [bam_pattern.format(**row) for index, row in cell_table.iterrows()]

def _pbsv_align_pbmm2_params(wildcards):
    """
    Get datatype-specific parameters for pbmm2 (CCS or subreads).

    :param wildcards: Rule wildcards.

    :return: PBSV parameters.
    """

    if wildcards.sample not in SAMPLE_TABLE.index:
        raise RuntimeError('Missing sample table entry for: {sample}'.format(**wildcards))

    sample_type = SAMPLE_TABLE.loc[wildcards.sample].squeeze()['TYPE'].lower()

    if sample_type == 'ccs':
        return '--preset CCS'
    elif sample_type == 'subreads':
        return '--median-filter'
    else:
        raise RuntimeError('Unrecognized sequence data type for sample {}: {}'.format(wildcards.sample, sample_type))




# pbsv_align_merge_bam
#
# Merge bams per sample (collapse movies).
rule pbsv_align_merge_bam:
    input:
        bam=_pbsv_align_get_bams_per_sample
    output:
        bam='align/{sample}/mapped_reads.bam',
        bai='align/{sample}/mapped_reads.bam.bai'
    shell:
        """samtools merge -O BAM -@ 6 {output.bam} {input.bam}; """
        """samtools index {output.bam}"""

# pbsv_align_get_bam
#
# Split samples into table files.
rule pbsv_align_get_bam:
    input:
        ref_fa=config['reference'],
        tab='init/{sample}/cell_table.tab'
    output:
        bam=temp('temp/align/{sample}/movie/movie_{movie}.bam'),
        bai=temp('temp/align/{sample}/movie/movie_{movie}.bam.bai')
    benchmark:
        'align/{sample}/bm/align/movie_{movie}.tab'
    run:

        # Read cell table
        df = pd.read_csv(input.tab, sep='\t', header=0, index_col='MOVIE', squeeze=True)

        try:
            movie = np.int32(wildcards.movie)
        except:
            raise RuntimeError('Movie wildcard is not an integer: {movie}'.format(**wildcards))

        if movie not in df.index:
            raise RuntimeError('Movie {} is not in the sample cell table: {}'.format(wildcards.movie, input.tab))

        df_movie = df.loc[movie].squeeze()

        # Set parameter dictionary
        param_dict = {
            'SAMPLE': wildcards.sample,
            'MOVIE': movie,
            'INPUT_FILE': df_movie['DATA'],
            'REF_FA': input.ref_fa,
            'OUTPUT_FILE': output.bam,
            'PARAMS': _pbsv_align_pbmm2_params(wildcards)
        }

        # Report temp directory
        print('Aligning: {sample} movie {movie}'.format(**wildcards))

        # Align
        shell((
            """pbmm2 align {INPUT_FILE} {REF_FA} {OUTPUT_FILE} --sort --sample '{SAMPLE}' {PARAMS} -j 4 -J 2; """
            """samtools index {OUTPUT_FILE}"""
        ).format(**param_dict))
