FROM nfcore/base
MAINTAINER Qi Zhao <zhaoqi@sysucc.org.cn>
LABEL authors="zhaoqi@sysucc.org.cn" \
    description="Docker image containing all requirements for the sysucc/rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-lncpipe-1.0dev/bin:$PATH


# Install DAtools 

RUN curl http://omicsbio.info/pub/DATOOLS/DAtools_v2.7.4.jar  && \
    chmod 777 DAtools_v2.7.4.jar  && \
    mv DAtools_v2.7.4.jar /opt/datools.jar
