# PBSV pipeline

## About

Runs the PacBio pbsv tool (https://github.com/PacificBiosciences/pbsv) to call structural variants. Pipeline is up to
date as of PBSV 2.2.0.

## Setup

First, change to a working directory where all output data will be saved. This can be a clean directory with nothing in
it.

### Configuration

Create a configuration file called `config.json`. Add the following text adjusting the path to the reference FASTA:

```
{
  "reference": "/net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa",
  "tandem_bed": "/net/eichler/vol27/projects/structural_variation/nobackups/resources/svpop/data/anno/trf/trf_regions_200_0.bed",
}
```

"tandem_bed" points to a bed file of tandem repeat regions. In this example, the UCSC TRF track was merged (bedtools merge)
concatenating with records that touch or are within 200 bp.

An optional parameter "min_svlen" sets the minimum SV/Indel length. The default value in PBSV is currently "20", and this
can be adjusted as low as "10" (< 10 causes PBSV to fail and return an error message).


### Define input samples

Create a `samples.tab` file with two tab-delimited columns:

* SAMPLE: Name of the sample
* FOFN: A path to an FOFN file that points to all input BAM files
* TYPE: "ccs" if the reads are CCS, and "subreads" if the reads are subreads.

The column headings (SAMPLE and FOFN) must be the first line. The FOFN file is a list of absolute paths to each input
BAM, one per line.

### Set environment

Define a variable that gives the full path to the pbsv pipeline code, which is the directory that contains `Snakefile`
and this `README.md` file. The pipeline itself does not use the variable, but commands in this README will.

Example:
`PIPELINE_DIR=/net/eichler/vol27/projects/structural_variation/nobackups/pipelines/pbsv/201910`

This section assumes the pbsv pipeline is not in the working directory, which is the recommended usage. That means the
current directory is where `samples.tab` is and where all output files will go, but the pipeline code is in
a different directory so that it does not need to be copied each time it is run. If the pipeline code is in the same
directory as the output files, you can set `PIPELINE_DIR` to '.' or you can ignore it completely and adjust the commands
below.

### Load prerequisites

A set of pacbio tools, including pbsv, minimap2, bcftools, samtools, and a python environment with Pandas installed must
be available before running the pipeline.

Pipeline was last tested with these Eichler lab modules:
```
module load htslib/1.9
module load bcftools/1.9
module load samtools/1.9
module load pbconda/201910
module load miniconda/4.5.12
```

If you are running on CentOS 7, replace `htslib/1.9` with `htslib/1.9-20`.


## Run

Once `samples.tab` has been defined, the prerequisite modules are loaded, and `PIPELINE_DIR` is set to point to the
pipeline directory where `Snakefile` exists, then the pipeline is ready to run.

The pipeline is executed in two steps. The first is a quick initialization process that finds paths to all input and
and saves them in the `results/init` directory. After this step, changes to `samples.tab` will be ignored, and the
pipeline will read the information in `results/init`. The second step then aligns and runs pbsv on the input samples.

### Run init

`snakemake -s ${PIPELINE_DIR}/Snakefile.init`

### Run pipeline

`snakemake -s ${PIPELINE_DIR}/Snakefile`


### Run distributed

To distribute jobs over the cluster, make sure DRMAA_LIBRARY_PATH is set in the environment (see below). Use the command
below. You may want to modify the number of concurrent jobs (-j).

`mkdir -p log; snakemake -s ${PIPELINE_DIR}/Snakefile -j 30 -k --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -e ./log -o ./log -pe serial {cluster.cpu} -l mfree={cluster.mem} -l h_rt={cluster.rt} -w n -S /bin/bash" -w 60 -u ${PIPELINE_DIR}/cluster.eichler.json --config ld_path=${LD_LIBRARY_PATH}`

If DRMAA_LIBRARY_PATH is not set in your environment, run this before distributing jobs:
`export DRMAA_LIBRARY_PATH=/opt/uge/lib/lx-amd64/libdrmaa.so.1.0`
