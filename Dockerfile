FROM continuumio/miniconda3:latest
LABEL authors="Jack Kamm and Samantha Hao" \
      description="Docker image containing all requirements for sars-cov-2 MSSPE pipeline"

COPY environment.yaml /
RUN /opt/conda/bin/conda env create -f /environment.yaml && /opt/conda/bin/conda clean -a
ENV PATH /opt/conda/envs/sc2-msspe/bin:$PATH
