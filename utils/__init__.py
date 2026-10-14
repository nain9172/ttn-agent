"""TTN Variant AI Agent - utility modules.

各模組直接 import（例如 `from utils.evo2_predictor import Evo2Predictor`）。
這裡刻意不預先載入任何子模組：否則只想用 utils.pubmind_lookup 的輕量腳本
（如 build_site.py）也會被迫載入 matplotlib、Bio 等重套件。
"""

__version__ = "1.0.0"
