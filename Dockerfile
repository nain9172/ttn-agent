# TTN agent：網站服務（web）與 GPU 分析（worker）共用的映像。
#
# 映像只放程式與 conda 環境；資料、模型、結果一律放在 docker volume（/vol），
# 授權資料（HGMD、cardiodb）絕不能進映像。用法見 docs/DEPLOY_DOCKER.md。
#
# 目標平台：DGX Spark（aarch64、GB10、CUDA 13）。environment.yml 與 wheels/ 的 flash-attn 都是 aarch64 版。
FROM nvidia/cuda:13.0.0-devel-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget git tini procps && \
    rm -rf /var/lib/apt/lists/*

# Miniforge（conda-forge，支援 aarch64）
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && rm /tmp/miniforge.sh

WORKDIR /app

# conda 環境（與 setup_env.sh 相同：environment.yml → gm4）。這一層最久（30–60 分鐘），放在程式碼之前以利快取
COPY environment.yml .
RUN /opt/conda/bin/conda env create -n gm4 -f environment.yml && /opt/conda/bin/conda clean -afy
ENV PATH=/opt/conda/envs/gm4/bin:/opt/conda/bin:$PATH \
    CONDA_DEFAULT_ENV=gm4

# flash-attn：PyPI 沒有 aarch64 版，用 repo 裡預先編好的 wheel；
# 再移除 evo2 依賴帶進來、只有空殼的 transformer-engine（有它 Evo2 反而無法載入）
COPY wheels/ /tmp/wheels/
RUN pip install --no-cache-dir /tmp/wheels/flash_attn-*-linux_aarch64.whl && rm -rf /tmp/wheels && \
    if pip show transformer_engine >/dev/null 2>&1 && ! python -c "import transformer_engine.pytorch" 2>/dev/null; then \
        pip uninstall -y transformer_engine transformer-engine-cu12 transformer-engine-cu13 || true; fi

# 執行期需要、但 conda / pip 不會裝的系統函式庫（vLLM 載入 Gemma 4 時會用到 OpenCV，它需要 libGL）。
# 刻意放在 conda 環境之後：之後要補新的系統函式庫只會重建這一層之後的部分，不必重建整個 conda 環境
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgl1 libglib2.0-0 && \
    rm -rf /var/lib/apt/lists/*

# 建置時就檢查 pipeline 需要的套件都載入得了（缺系統函式庫會在這裡失敗，而不是等到分析時才發現）。
# 建置時沒有 GPU，這裡只做 import，不初始化 CUDA
# 在 /tmp 執行：Evo2 載入時會在工作目錄建立 activations_debug.log，留在 /app 會變成 root 擁有、一般帳號寫不進去
RUN cd /tmp && python -c "import cv2, torch, transformers, flash_attn, evo2, vllm; from vllm import LLM; print('imports OK')" && \
    rm -f /tmp/activations_debug.log

# 程式碼。.dockerignore 排除了 data/、outputs/、site/、research/ 與所有資料檔；
# data/ 只會帶進 git 追蹤的 transcript_interval.xlsx，移到 seed/，啟動時才複製進 volume
COPY . .
# 程式用的 data/、outputs/、site/、llm_logs/ 在建置時就指向 volume，執行時不需要 root 權限
RUN mkdir -p /app/seed && if [ -d /app/data ]; then mv /app/data /app/seed/data; fi && \
    for d in data outputs site llm_logs; do rm -rf "/app/$d"; ln -s "/vol/$d" "/app/$d"; done && \
    chmod +x /app/docker/entrypoint.sh /app/setup_and_run.sh /app/setup_env.sh /app/run_batch_inprocess.sh && \
    chmod -R a+rX /app && \
    # 有些套件會在工作目錄寫 log（例：Evo2 的 activations_debug.log）；以一般帳號執行時也要能建立新檔
    chmod 1777 /app

# 模型權重、編譯快取與網站都在 volume；寫在映像裡，docker compose exec 進來的指令也看得到
ENV TTN_VOLUME=/vol \
    TTN_SERVER_MODE=1 \
    TTN_PORT=8080 \
    TTN_SITE_DIR=/vol/site \
    HF_HOME=/vol/models/huggingface \
    XDG_CACHE_HOME=/vol/cache \
    TRITON_CACHE_DIR=/vol/cache/triton \
    VLLM_CACHE_ROOT=/vol/cache/vllm \
    MPLCONFIGDIR=/vol/cache/matplotlib \
    TORCHINDUCTOR_CACHE_DIR=/vol/cache/torchinductor
EXPOSE 8080
VOLUME ["/vol"]

ENTRYPOINT ["/usr/bin/tini", "--", "/app/docker/entrypoint.sh"]
CMD ["web"]
