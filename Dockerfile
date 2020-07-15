FROM nfcore/base:1.9
LABEL authors="Jack Kamm and Samantha Hao" \
      description="Docker image containing all software requirements for the nf-core/msspe pipeline"

RUN apt-get update && \
	apt-get install -y unzip autoconf build-essential curl && \
	rm -rf /var/lib/apt/lists/*

RUN PERL_MM_USE_DEFAULT=1 cpan install Inline Inline::C

RUN wget https://raw.githubusercontent.com/nawrockie/vadr/master/vadr-install.sh && \
	sh vadr-install.sh linux && \
	rm vadr-install.sh

RUN mkdir -p /data/vadr-models-corona && \
	cd /data && \
	mkdir -p vadr-models-corona && \
	wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/coronaviridae/CURRENT/vadr-models-corona-1.1-1.tar.gz -O vadr-models-corona.tar.gz && \
	tar xvf vadr-models-corona.tar.gz -C vadr-models-corona --strip-components 1 && \
	rm vadr-models-corona.tar.gz

# Install the conda environment
COPY environment.yaml /
RUN conda env create -f /environment.yaml && conda clean -ay

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/sc2-msspe/bin/:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-msspe-1.0dev > nf-core-msspe-1.0dev.yaml

# Install baltic
RUN pip install git+git://github.com/evogytis/baltic.git
