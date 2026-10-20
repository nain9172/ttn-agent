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


def meminfo_gib(key: str) -> float:
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith(key + ":"):
                return int(line.split()[1]) / 2**20
    return 0.0


def need_gib() -> float:
    """接下來要載入的模型需要多少 GPU 可用記憶體：NEED_GIB（setup_and_run.sh 會設），否則 Gemma 的需求。"""
    if os.environ.get("NEED_GIB"):
        return float(os.environ["NEED_GIB"])
    sys.path.insert(0, str(ROOT))
    import config
    import torch

    return torch.cuda.mem_get_info(0)[1] / 2**30 * float(config.LLM_GPU_MEMORY_UTILIZATION)


def reclaim(margin_gib: float = 6.0) -> float:
    """可用記憶體不夠、但差額能靠清檔案快取補回來時，暫時配置匿名記憶體逼核心回收快取，再全部釋放。

    GB10 的 GPU 與 CPU 共用記憶體，檔案快取會被算成 GPU「已使用」；vLLM 載入前只看 mem_get_info，
    不會等核心回收快取就放棄。別人讀過的檔案（或容器看不到的檔案）無法用 fadvise 清，又沒有 sudo
    可以 drop_caches 時，這是唯一不需權限的作法：配置量約等於模型需要的量（核心只在空閒記憶體
    用完時才回收快取），並隨時保留 margin_gib 的可用記憶體，不會把別的程式擠到 OOM。
    被其他程式佔住（快取補不回差額）時不做。回傳實際配置過的 GiB。
    """
    try:
        import torch

        free = torch.cuda.mem_get_info(0)[0] / 2**30
        need = need_gib()
    except Exception:  # noqa: BLE001
        return 0.0
    deficit = need - free + 1.0
    if deficit <= 0 or meminfo_gib("Cached") < deficit:
        return 0.0
    target = free + deficit
    chunks, got = [], 0.0
    try:
        while got + 1 <= target and meminfo_gib("MemAvailable") - 1 > margin_gib:
            chunks.append(b"\x01" * 2**30)  # 寫入每一頁，才會真的佔用實體記憶體
            got += 1
    except MemoryError:
        pass
    finally:
        chunks.clear()
    return got


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
    # docker：data/、outputs/ 是指向 volume 的 symlink，os.walk 不會跟進去，直接把 volume 加進來
    if os.environ.get("TTN_VOLUME"):
        targets.append(Path(os.environ["TTN_VOLUME"]))
    n = advise(targets)
    try:
        dropped = subprocess.run(["sudo", "-n", "sh", "-c", "sync; echo 3 > /proc/sys/vm/drop_caches"],
                                 capture_output=True).returncode == 0
    except OSError:  # 容器裡沒有 sudo
        dropped = False
    pushed = 0.0 if dropped else reclaim()
    print(f"    已釋放 {n:,} 個檔案的快取{'，並清除整機 page cache' if dropped else ''}"
          f"{f'，並以 {pushed:.0f} GiB 暫時配置逼出其餘快取' if pushed else ''}；"
          f"GPU 可用記憶體 {before} → {gpu_free()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
