#!/usr/bin/env python3
"""
自包含 demo —— 不需要下载真实基因组即可看到完整产出。

原理:
  构造一个合成的迷你基因组（~3 kb），把 Table S2 里 sgRNA26 的
  "PAM + spacer" 按正确构型埋进去，然后跑完整规划流程。

  这样做的价值: 你可以**手工核对**工具算出的 PAM 坐标与预测插入位点，
  因为底物是我们自己造的。

用法:
  python scripts/demo_offline.py            # 打印方案 Markdown
  python scripts/demo_offline.py --out demo_plan.md
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stib.knowledge import KnowledgeBase                      # noqa: E402
from stib.planner import make_plan                            # noqa: E402
from stib.render import render_markdown, render_json          # noqa: E402
from stib.sequence import (                                   # noqa: E402
    GenomeIndex, find_pams, INSERT_WINDOW, SPACER_LEN,
)

# Table S2 已发表序列
SG26_SPACER = "ATTCGGAAGAAGACGTGGTATCA"
SG26_PAM = "GTT"

#: sgRNA26 的 PAM 在合成基因组中的起点（可手工核对）
PAM_AT = 1000

#: 目标基因区间（demo 用，模拟一个 600 bp 的基因）
GENE_START = 700
GENE_END = 1300


def build_synthetic_genome() -> str:
    """
    构造合成基因组。

    构型:  ... target | GTT | ATTCGGAAGAAGACGTGGTATCA | ...
            ^              ^
            PAM_AT=1000    spacer 起点 = PAM_AT + 3 = 1003
    """
    # 用去重复的填充序列，避免种子区唯一性检查被误判
    def filler(n: int, seed: int = 1) -> str:
        out = []
        i = seed
        bases = "ACGT"
        while len(out) < n:
            # 简单 LCG，产生可复现但不重复的序列
            i = (i * 1103515245 + 12345) & 0x7FFFFFFF
            b = bases[(i >> 16) & 3]
            # 避免 4 连相同碱基
            if len(out) >= 4 and out[-1] == out[-2] == out[-3] == out[-4] == b:
                b = "ACGT"[("ACGT".index(b) + 1) % 4]
            out.append(b)
        return "".join(out)

    left = filler(PAM_AT, seed=7)
    right = filler(3000 - PAM_AT - 3 - SPACER_LEN, seed=99)
    genome = left + SG26_PAM + SG26_SPACER + right
    assert len(genome) == 3000, len(genome)
    return genome


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", help="写出 Markdown 到该路径")
    ap.add_argument("--json", help="写出 JSON 到该路径")
    args = ap.parse_args()

    genome = build_synthetic_genome()
    g = GenomeIndex(genome, "synthetic_demo_3kb")

    print("=" * 72)
    print("合成基因组自检（可手工核对）")
    print("=" * 72)
    print(f"  长度              : {len(g):,} bp")
    print(f"  埋入的 PAM 位置   : {PAM_AT} (期望 {PAM_AT})")
    print(f"  埋入的 PAM 序列   : {genome[PAM_AT:PAM_AT + 3]} (期望 {SG26_PAM})")
    print(f"  埋入的 spacer     : {genome[PAM_AT + 3:PAM_AT + 3 + SPACER_LEN]}")
    print(f"                   (期望 {SG26_SPACER})")
    print(f"  spacer 长度       : {len(SG26_SPACER)} (期望 {SPACER_LEN})")

    # 扫描全基因组，确认只埋了一处
    pams = list(find_pams(genome))
    gtt = [p for p in pams if p[1] == "GTT"]
    print(f"  全基因组 GTT 数   : {len(gtt)}")
    print(f"  GTT 位置          : {[p for p, _ in gtt]}")

    ok = (genome[PAM_AT:PAM_AT + 3] == SG26_PAM
          and genome[PAM_AT + 3:PAM_AT + 3 + SPACER_LEN] == SG26_SPACER)
    print(f"\n  构型核对: {'PASS' if ok else 'FAIL'}")

    # --- 跑规划 ---
    kb = KnowledgeBase()
    plan = make_plan(
        genome=g,
        gene="demo_gene",
        start=GENE_START,
        end=GENE_END,
        species="Bacteroides thetaiotaomicron",
        top_n=4,
        kb=kb,
    )

    print("\n" + "=" * 72)
    print("方案摘要")
    print("=" * 72)
    t = plan.target
    print(f"  候选扫描数       : {t['n_candidates_scanned']}")
    print(f"  推荐 spacer      : {t['selected_spacer']}")
    print(f"  推荐 PAM         : {t['selected_pam']} @ {t['selected_pam_start']}")
    print(f"  预测插入位点     : {t['predicted_insert_site']}")
    print(f"  氯霉素浓度       : {t['cm_conc']}")
    print(f"  高等级警示       : {sum(1 for w in plan.warnings if w.severity == 'high')}")
    print(f"  风险条目         : {len(plan.risks)}")
    print(f"  判定点           : {len(plan.decisions)}")

    # 核对：推荐的那条 spacer 是否就是埋进去的 sgRNA26
    hit = any(s.spacer == SG26_SPACER for s in plan.sgrnas)
    print(f"\n  推荐集合含 sgRNA26: {'是' if hit else '否（排序把它排到 top_n 之外，属正常）'}")

    md = render_markdown(plan)
    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(md)
        print(f"\n  方案已写出: {args.out}")
    else:
        print("\n" + "=" * 72)
        print("完整方案 Markdown")
        print("=" * 72 + "\n")
        print(md)

    if args.json:
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_json(plan))
        print(f"  JSON 已写出: {args.json}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
