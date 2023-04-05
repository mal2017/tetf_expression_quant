# TETF_EXPRESSION_QUANT: processing of DGRP expression data for the Ellison lab TETF project

## dependencies

- `snakemake` (tested with v7.21.0 on Ubuntu 22 LTS and 7.18.2 on CentOS 7)
- `mamba` (tested with version 1.0.0 on Ubuntu 22 LTS and 0.9.1 on CentOS 7). `conda` alone should work just fine.
- `singularity` (tested with v3.1.0-1 on CentOS 7) or `apptainer` (tested with v1.1.3 on Ubuntu 22 LTS)

## usage

```
snakemake --use-conda --use-singularity --cores <threads>
```

A test run (just DGRP_100 and DGRP_101) can be performed by adding `--config RUN_TYPE=test --SALMON_REPLICATES=1`

## Description

This workflow qunatifies expression using salmon. It depends on several other `tetf_*` repositories developed
for the project and is intended to be run as a submodule/subworkflow of the parent repo: `tetf`.
