FROM condaforge/miniforge3:24.7.1-2

ENV CONDA_DEFAULT_ENV="psscoring"
ENV PATH="/opt/psscoring:$PATH"

WORKDIR /opt

RUN groupadd -r psscoring && \
    useradd -r -g psscoring -m -d /home/psscoring psscoring

COPY ./environment.yaml ./
RUN set -ex \
    && conda update --name base --channel conda-forge conda \
    && conda env create --file environment.yaml \
    && conda clean --force-pkgs-dirs --all --yes \
    && rm ./environment.yaml \
    && echo "conda activate ${CONDA_DEFAULT_ENV}" >> ~/.bashrc

# Copy the scoring script
COPY src/PSscoring /opt/psscoring
RUN set -ex \
    && chmod +x /opt/psscoring/ps.py

WORKDIR /home

CMD ["/bin/bash"]