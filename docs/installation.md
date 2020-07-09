# Installation
<!-- MarkdownTOC -->

- [Nextflow](#nextflow)
- [Pipeline software dependencies](#pipeline-software-dependencies)
- [Pipeline code](#pipeline-code)

<!-- /MarkdownTOC -->

## Nextflow

The `czbiohub/sc2-msspe-bioinfo` pipeline uses Nextflow, which must be present on the system where the pipeline is launched. See [nextflow.io](https://www.nextflow.io/) for the latest installation instructions.

You can also install Nextflow using Bioconda:

```
conda install -c bioconda nextflow
```

## Pipeline software dependencies

There is a docker container [`czbiohub/sc2-msspe`](https://hub.docker.com/repository/docker/czbiohub/sc2-msspe) available on DockerHub.

Alternatively, pipeline dependencies can be installed into a conda environment `sc2-msspe`:

```
conda env create -f environment.yaml
```

## Pipeline code

Nextflow will automatically fetch the pipeline code from GitHub if `czbiohub/sc2-msspe-bioinfo` is specified as the pipeline name.

-------------

Documentation adapted from [`nf-core`](https://nf-co.re/usage/installation).