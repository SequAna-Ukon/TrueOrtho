FROM mambaorg/micromamba:1.5.8

USER root

# 1. Install system utilities and OpenBLAS runtime libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    wget \
    curl \
    git \
    libopenblas0-serial \
    libopenblas-dev \
 && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# 2. Copy environment file into container
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# 3. Build environment and clean up cache
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# 4. Set environment paths & ensure dynamic libraries in /opt/conda/lib are found
ENV PATH="/opt/conda/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/lib:/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}"

WORKDIR /workspace
