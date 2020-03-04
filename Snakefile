"""
Pipeline for running PBSV.
"""

import os

#
# Init
#

# Initialize, import, and read configs
include: 'rules/include.snakefile'

# Set LD_LIBRARY_PATH

if 'ld_path' in config:
    os.environ['LD_LIBRARY_PATH'] = config['ld_path']

#
# Rules
#

def _get_all_list(wildcards):
    df = pd.read_csv('samples.tab')

# pbsv_all
#
# Run whole pipeline
rule pbsv_all:
    input:
        calls_sv=expand('pbsv_sample_{sample}_sv.vcf.gz', sample=SAMPLE_TABLE['SAMPLE']),
        calls_bnd=expand('pbsv_sample_{sample}_bnd.vcf.gz', sample=SAMPLE_TABLE['SAMPLE']),
        calls_dup=expand('pbsv_sample_{sample}_dup.vcf.gz', sample=SAMPLE_TABLE['SAMPLE'])


include: 'rules/align.snakefile'
include: 'rules/call.snakefile'
