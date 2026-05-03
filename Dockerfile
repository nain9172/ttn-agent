FROM nvidia/cuda:12.8.0-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget curl git && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge (supports aarch64)
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH=/opt/conda/bin:$PATH

WORKDIR /app

COPY env.yaml .

RUN conda env create -f env.yaml && conda clean -afy

SHELL ["conda", "run", "-n", "gm3", "/bin/bash", "-c"]

ENV CONDA_DEFAULT_ENV=gm3
ENV PATH=/opt/conda/envs/gm3/bin:$PATH

COPY . .

CMD ["bash"]
