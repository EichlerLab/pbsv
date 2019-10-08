"""
Setup pipeline. Should be included in all snakefiles.
"""

#
# Import packages
#

import datetime
import gzip
import os
import sys
import tempfile

import pandas as pd
import numpy as np


#
# Read config
#

CONFIG_FILE_NAME = 'config.json'

configfile: CONFIG_FILE_NAME


#
# Set locations
#

RULES_DIR = os.path.dirname(workflow.snakefile)

if os.path.basename(workflow.snakefile) == 'Snakefile':
    RULES_DIR = os.path.join(RULES_DIR, 'rules')

PIPELINE_DIR = os.path.dirname(RULES_DIR)

WORKING_DIR = os.getcwd()

sys.path.append(PIPELINE_DIR)


#
# Import pipeline packages
#

# Import from a local package from PIPELINE_DIR here. System path is setup to find them.



#
# Init
#

# Bash "strict mode"
shell.prefix('set -euo pipefail; ')

# Set paths from configuration parameters (if given)
LD_LIBRARY_PATH = config.get('ld_path', None)
PATH = config.get('path', None)

if LD_LIBRARY_PATH is not None:
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH

if PATH is not None:
    os.environ['PATH'] = PATH

# Save environment
PROCESS_ENV = os.environ.copy()


#
# Read sample table
#

# Definitions
SAMPLE_TABLE_FILE = os.path.join(WORKING_DIR, 'samples.tab')
SAMPLE_TABLE_COLUMNS = ['SAMPLE', 'FOFN', 'TYPE']

# Check for file
if not os.path.exists(SAMPLE_TABLE_FILE):
    raise RuntimeError('Missing sample table: {}'.format(SAMPLE_TABLE_FILE))

# Read
SAMPLE_TABLE = pd.read_csv(SAMPLE_TABLE_FILE, sep='\t', header=0)

# Check for missing columns
missing_cols = [col for col in SAMPLE_TABLE_COLUMNS if col not in SAMPLE_TABLE.columns]

if missing_cols:
    raise RuntimeError('Missing sample table column(s) "{}": {}'.format(', '.join(missing_cols), SAMPLE_TABLE_FILE))

del(missing_cols)

# Read
SAMPLE_TABLE.set_index('SAMPLE', inplace=True, drop=False)
SAMPLE_TABLE = SAMPLE_TABLE.loc[:, SAMPLE_TABLE_COLUMNS]
SAMPLE_TABLE.set_index('SAMPLE', inplace=True, drop=False)

# Check for duplicate sample names
if len(set(SAMPLE_TABLE['SAMPLE'])) != SAMPLE_TABLE.shape[0]:
    raise RuntimeError('Found duplicate sample names in sample table: {}'.format(SAMPLE_TABLE_FILE))


#
# Set TEMPDIR
#

TEMP_DIR = config.get('tempdir', None)

if TEMP_DIR is None or TEMP_DIR == '' or TEMP_DIR == 'None':
    TEMP_DIR = tempfile.gettempdir()
else:
    TEMP_DIR = os.path.abspath(TEMP_DIR)

if os.path.isdir(TEMP_DIR) and os.path.samefile(TEMP_DIR, '.'):
    # Defaults to local directory if a temp cannot be found. This
    # should not occur on real systems, but don't clutter the working
    # directory if it does.
    TEMP_DIR = os.path.join(TEMP_DIR, 'temp')

#
# Set include flag
#

INCLUDE_SNAKEFILE = True
