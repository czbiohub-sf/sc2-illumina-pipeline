FROM nfcore/base:dev
LABEL authors="Samantha Hao, Jack Kamm" \
      description="Docker image containing all software requirements for the nf-core/msspe pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-msspe-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-msspe-1.0dev > nf-core-msspe-1.0dev.yml
