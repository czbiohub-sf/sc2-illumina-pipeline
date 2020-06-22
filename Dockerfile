FROM nfcore/base:1.9
LABEL authors="Jack Kamm and Samantha Hao" \
      description="Docker image containing all software requirements for the nf-core/msspe pipeline"

# Install the conda environment
COPY environment.yaml /
RUN conda env create -f /environment.yaml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/sc2-msspe/bin/:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-msspe-1.0dev > nf-core-msspe-1.0dev.yaml

# Install fastv
wget http://opengene.org/fastv/fastv && chmod a+x ./fastv && mv fastv /bin

# Install baltic
RUN pip install git+git://github.com/evogytis/baltic.git
