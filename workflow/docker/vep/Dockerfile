ARG VEP_RELEASE=113.4
ARG SAMTOOLS_VERSION=1.21
ARG LOFTEE_VERSION=1.0.3
ARG OPT_SRC=/opt/src

FROM ubuntu:22.04 AS samtools_builder
ARG SAMTOOLS_VERSION
WORKDIR /opt/tmp
USER root
RUN set -ex \
    # Install dependencies
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        autoconf \
        automake \
        bzip2 \ 
        curl \
        make \
        gcc \
        perl \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses5-dev \
        libdeflate-dev \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sLk --output /opt/tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -jxvf /opt/tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /opt/tmp \
    && cd /opt/tmp/samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/opt/tools/samtools-${SAMTOOLS_VERSION}  \
    && make \
    && make install 
    # /opt/tools/samtools-${SAMTOOLS_VERSION}/bin/*

FROM ubuntu:24.10 AS loftee_builder
ARG LOFTEE_VERSION
WORKDIR /opt/tmp
# USER root
RUN set -ex \
    # Install dependencies
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        curl \
        unzip \
    # libncurses5-dev \
    > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sLk --output /opt/tmp/master.zip \
        https://github.com/konradjk/loftee/archive/refs/tags/v${LOFTEE_VERSION}.zip \
    && unzip /opt/tmp/master.zip -d /opt/tmp \
    # Rename loftee folder
    && ls -l /opt/tmp \
    && rm -rf /opt/tmp/loftee-${LOFTEE_VERSION}/LICENSE /opt/tmp/loftee-${LOFTEE_VERSION}/README.md \
    && mv /opt/tmp/loftee-${LOFTEE_VERSION}/ /opt/tmp/loftee/ \
    && ls -l /opt/tmp/loftee
    
#3. Finalize
FROM ensemblorg/ensembl-vep:release_${VEP_RELEASE} AS vep
ARG SAMTOOLS_VERSION
ARG LOFTEE_VERSION

ENV PATH=${PATH}:/opt/samtools-${SAMTOOLS_VERSION}/bin
ENV PERL5LIB=${PERL5LIB_TMP}:/vep/loftee

COPY --chown=vep:vep --from=samtools_builder /opt/tools/samtools-${SAMTOOLS_VERSION}/* /opt/samtools-${SAMTOOLS_VERSION}
COPY --chown=vep:vep --from=loftee_builder /opt/tmp/loftee /vep/loftee

# RUN ln -s $OPT/.vep /data
USER vep
WORKDIR /data