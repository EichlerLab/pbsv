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
    cell_table = pd.read_csv('init/cell_table.tsv.gz', sep='\t', header=0)

    cell_table = cell_table.loc[cell_table['SAMPLE'] == wildcards.sample]

    if cell_table.shape[0] == 0:
        raise RuntimeError('No input for sample in cell table: {}'.format(wildcards.sample))

    # Return a list of BAM files
    return [bam_pattern.format(**row) for index, row in cell_table.iterrows()]


def _pbsv_align_pbmm2_get_type(wildcards):
    """
    Get sample type.

    :param wildcards: Rule wildcards.

    :return: Sample type (e.g. "ccs" or "subreads").
    """

    if wildcards.sample not in SAMPLE_TABLE.index:
        raise RuntimeError('Missing sample table entry for: {sample}'.format(**wildcards))

    return SAMPLE_TABLE.loc[wildcards.sample].squeeze()['TYPE'].lower()


def _pbsv_align_pbmm2_align_params(wildcards):
    """
    Get datatype-specific parameters for pbmm2 align (CCS or subreads).

    :param wildcards: Rule wildcards.

    :return: PBSV parameters.
    """

    sample_type = _pbsv_align_pbmm2_get_type(wildcards)

    if sample_type == 'ccs':
        return '--preset CCS -L 0.1 -c 0'
    elif sample_type == 'subreads':
        return '--preset SUBREAD --median-filter'
    else:
        raise RuntimeError('Unrecognized sequence data type for sample {}: {}'.format(wildcards.sample, sample_type))


def _pbsv_align_pbmm2_index_params(wildcards):
    """
    Get datatype-specific parameters for pbmm2 index (CCS or subreads).

    :param wildcards: Rule wildcards.

    :return: PBSV parameters.
    """

    sample_type = wildcards.type

    if sample_type == 'ccs':
        return '--preset CCS'
    elif sample_type == 'subreads':
        return '--preset SUBREAD'
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
        mmi=lambda wildcards: 'data/ref_{type}.mmi'.format(type=_pbsv_align_pbmm2_get_type(wildcards)),
        tsv='init/{sample}/cell_table.tsv.gz'
    output:
        bam=temp('temp/align/{sample}/movie/movie_{movie}.bam'),
        bai=temp('temp/align/{sample}/movie/movie_{movie}.bam.bai')
    benchmark:
        'align/{sample}/bm/align/movie_{movie}.tab'
    run:

        # Read cell table
        df = pd.read_csv(input.tsv, sep='\t', header=0, index_col='MOVIE', squeeze=True)

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
            'REF_MMI': input.mmi,
            'OUTPUT_FILE': output.bam,
            'PARAMS': _pbsv_align_pbmm2_align_params(wildcards)
        }

        # Report temp directory
        print('Aligning: {sample} movie {movie}'.format(**wildcards))

        # Align
        shell((
            """pbmm2 align --sort --sample '{SAMPLE}' {PARAMS} -j 6 -J 2 {REF_MMI} {INPUT_FILE} {OUTPUT_FILE}; """
            """samtools index {OUTPUT_FILE}"""
        ).format(**param_dict))

# pbsv_align_index_reference
#
# Index the reference. Generates an MMI for mapping with pbmm2.
rule pbsv_align_index_reference:
    input:
        ref_fa='data/ref.fa'
    output:
        mmi='data/ref_{type}.mmi'
    run:

        # Index (generate MMI)
        # Set parameter dictionary
        param_dict = {
            'INPUT_FILE': input.ref_fa,
            'OUTPUT_FILE': output.mmi,
            'PARAMS': _pbsv_align_pbmm2_index_params(wildcards)
        }

        # Index
        shell((
            """pbmm2 index {PARAMS} {INPUT_FILE} {OUTPUT_FILE}"""
        ).format(**param_dict))


rule pbsv_align_get_ref_fa:
    input:
        ref_fa=config['reference']
    output:
        ref_fa='data/ref.fa'
    run:

        # If compressed, uncompress
        ref_tok = (input.ref_fa.lower().split('.'))

        if ref_tok[-1] == 'gz':

            if len(ref_tok) < 3:
                raise RuntimeError(f'Unrecognized reference format: Ends with gz, but no extension before it indicating the file type: {input.ref_fa}')

            if ref_tok[-2] not in {'fa', 'fasta', 'fn'}:
                raise RuntimeError(f'Unrecognized reference format: Ends with gz, but unrecognized FASTA extension before it: {input.ref_fa}')

            shell(
                """zcat {input.ref_fa} > {output.ref_fa}"""
            )

            fa_in_name = output.ref_fa

        elif ref_tok[-1] in {'fa', 'fasta', 'fn'}:

            shell(
                """ln -s $(readlink -f {input.ref_fa}) {output.ref_fa}"""
            )

        else:
            raise RuntimeError(f'Unrecognized reference format: No recognized FASTA extension (gz or uncompressed): {input.ref_fa}')
