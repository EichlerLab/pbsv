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

# pbsv_all
#
# Run whole pipeline
rule pbsv_all:
    input:
        calls_sv='pbsv_sv.vcf.gz',
        calls_bnd='pbsv_bnd.vcf.gz'


include: 'rules/align.snakefile'
include: 'rules/call.snakefile'
