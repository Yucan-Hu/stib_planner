"""
STIB Planner — 拟杆菌 STIB 基因编辑方案生成器（Demo-0）。

纯 Python 标准库。无第三方依赖。

快速开始:
    python -m stib.cli selftest
    python -m stib.cli fetch NC_004663.1
    python -m stib.cli plan --genome data/NC_004663.1.fasta \\
        --gene tdk --motif ATTCGGAAGAAGACGTGGTATCA --out plan_tdk.md
"""

__version__ = "0.0.1"

__all__ = ["knowledge", "sequence", "planner", "render", "cli"]
