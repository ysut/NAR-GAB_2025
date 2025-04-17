FROM debian:bullseye-slim

#1. Set build arguments
# miniforge
ARG MINIFORGE_NAME=Miniforge3
ARG MINIFORGE_VERSION=23.1.0-4
ARG ENV_NAME=spliceai

#2. Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV CONDA_DIR=/opt/conda 
ENV PATH=${CONDA_DIR}/bin:${PATH}

WORKDIR /usr/tools/tmp

# RUN
RUN set -ex \
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        bzip2 \
        ca-certificates \
        git \
        tini \
        wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && wget --no-hsts --quiet \
        https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-$(uname -m).sh \
        -O /tmp/miniforge.sh \
    && /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniforge.sh \
    && conda clean --tarballs --index-cache --packages --yes \
    && find ${CONDA_DIR} -follow -type f -name '*.a' -delete \
    && find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete \
    && conda clean --force-pkgs-dirs --all --yes \
    && echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc \
    && echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

# Install conda packages using environment.yml
COPY ./environment.yaml ./
RUN set -ex \
    && conda update --name base --channel conda-forge conda \
    && conda env create --file environment.yaml \
    && conda clean --force-pkgs-dirs --all --yes \
    && rm -rf /usr/tools/tmp

# Set working directory and entrypoint
WORKDIR /
COPY ./entrypoint.sh /usr/tools/entrypoint.sh
RUN chmod +x /usr/tools/entrypoint.sh

ENTRYPOINT ["/usr/tools/entrypoint.sh"]
# CMD ["spliceai", "-h"]