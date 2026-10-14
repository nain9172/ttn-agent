#!/usr/bin/env python3
"""釋放檔案快取，讓 vLLM 的顯存預檢看得到足夠的可用記憶體。

DGX Spark（GB10）的 CPU 與 GPU 共用記憶體，torch.cuda.mem_get_info 會把作業系統的
page cache 算成已使用。讀過模型權重、PDF、報告之後，page cache 動輒數十 GB，Gemma 4
（需要 0.85 × 總量）就會被預檢擋下。

不需要 root：對自己讀得到的檔案呼叫 posix_fadvise(DONTNEED)，請核心丟掉它們的快取頁
（只影響乾淨的快取，不會遺失資料）。可以免密碼 sudo 時再整機 drop_caches。

用法：python scripts/release_page_cache.py [額外目錄 ...]
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def advise(paths) -> int:
    n = 0
    for base in paths:
        if not base.exists():
            continue
        for dirpath, _, files in os.walk(base):
            for name in files:
                try:
                    fd = os.open(os.path.join(dirpath, name), os.O_RDONLY)
                except OSError:
                    continue
                try:
                    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                    n += 1
                except OSError:
                    pass
                finally:
                    os.close(fd)
    return n


def gpu_free() -> str:
    try:
        import torch

        free, total = torch.cuda.mem_get_info(0)
        return f"{free / 2**30:.1f} / {total / 2**30:.1f} GiB"
    except Exception:  # noqa: BLE001
        return "?"


def main() -> int:
    before = gpu_free()
    hf_home = Path(os.environ.get("HF_HOME", Path.home() / ".cache" / "huggingface"))
    targets = [hf_home, ROOT, Path(sys.prefix)] + [Path(p) for p in sys.argv[1:]]
    n = advise(targets)
    dropped = subprocess.run(["sudo", "-n", "sh", "-c", "sync; echo 3 > /proc/sys/vm/drop_caches"],
                             capture_output=True).returncode == 0
    print(f"    已釋放 {n:,} 個檔案的快取{'，並清除整機 page cache' if dropped else ''}；"
          f"GPU 可用記憶體 {before} → {gpu_free()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
