#!/usr/bin/env python3
"""
对已取回的基因组验证工具的核心逻辑。

三类检查:
  1. **已发表锚点** — Table S2 中属于本物种的 sgRNA 必须能在基因组中找到，
     且其 **spacer 起点前 3 nt** 等于该 sgRNA 的 PAM
  2. **PAM 扫描一致性** — 用工具自己的 scan_sgrnas 独立复算，
     确认能重新发现同一位点
  3. **插入窗口** — 计算预测插入位点，人工核对

⚠️ PAM 位置约定（ShCAST / Type V-K）:
     5'-... [PAM 3nt] [protospacer 23nt] ...-3'
   即 **PAM 在 spacer 上游**。因此核对方式是「spacer 起点前 3 nt == PAM」。
   （v0.1 的脚本查成了下游，导致误报，已修正。）

用法:
  python scripts/verify_against_paper.py --genome data/NC_004663.1.fasta
  python scripts/verify_against_paper.py --genome ... --dump 200
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stib.sequence import (                                   # noqa: E402
    GenomeIndex, read_fasta, scan_sgrnas, revcomp, dedupe_sgrnas,
    INSERT_WINDOW, SPACER_LEN,
)

#: Table S2 中属于 *B. thetaiotaomicron* VPI-5482 的锚点。
#: 只列本物种的 —— sgRNA60/65/76/78/79/80 属于其他物种，不在本基因组中。
BT_ANCHORS: list[tuple[str, str, str, str, str]] = [
    # (名称, spacer, 期望 PAM, spacer_location, target_site)
    ("sgRNA26", "ATTCGGAAGAAGACGTGGTATCA", "GTT", "tdk", "Bt01"),
    ("sgRNA39", "TGGTTTAATCGGGATGACGGCTT", "GTT", "bt1754", "Bt02"),
    ("sgRNA47", "GTTCGTGTTTGCCGTGTATCATT", "GTT", "neutral site 5", "Bt03"),
    ("sgRNA31", "ATATGGGTCCCGATGTATGCCAG", "GTT", "lpxF", "Bt04"),
    ("sgRNA44", "CAGATCGGTGCTGGGAGTGGTTC", "GTT", "neutral site 44", "Bt05"),
    ("sgRNA48", "ATGGTGGCGTTTCCATCACAAGA", "GTT", "BT3804", "Bt06"),
    ("sgRNA50", "TCATACTCTTTGATATTGCTATG", "GTC", "BT0234", "Bt07"),
    ("sgRNA52", "CTTCAATGGAATACATTCGCCTT", "GTT", "lpxF", "Bt08"),
    ("sgRNA55", "TTATCAGGTTGATGCGCAAAGAT", "GTT", "neutral site 55", "Bt09"),
    ("sgRNA57", "TATAACGTTTGCGCTTATATGCC", "GTT", "BT1438", "Bt10"),
]


def pam_before(g: GenomeIndex, sp_start: int) -> str:
    """返回 spacer 起点**上游** 3 nt（即 PAM 位置）。"""
    return g.seq[sp_start - 3:sp_start] if sp_start >= 3 else ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--dump", type=int, default=0,
                    help="打印每个锚点前后 N bp 上下文")
    args = ap.parse_args()

    seqs = read_fasta(args.genome)
    if not seqs:
        raise SystemExit("ERROR: 未读到序列")
    sid = next(iter(seqs))
    g = GenomeIndex(seqs[sid], sid)
    print(f"基因组: {g.id}  {len(g):,} bp\n")

    fails = 0

    # ---------- 检查 1: 已发表锚点（PAM 在 spacer 上游）----------
    print("[1] 已发表 sgRNA 锚点（Table S2，仅本物种）")
    print("    PAM 约定: 5'-[PAM]-[spacer 23nt]-3'（PAM 在 spacer 上游）\n")
    for name, spacer, exp_pam, loc, site in BT_ANCHORS:
        if len(spacer) != SPACER_LEN:
            print(f"  [WARN] {name}: spacer 长度 {len(spacer)} != {SPACER_LEN}")

        hits_fwd = g.find_all(spacer)
        hits_rev = g.find_all(revcomp(spacer))

        if not hits_fwd and not hits_rev:
            print(f"  [FAIL] {name:<8} {site:<6} {loc:<16} 未找到")
            fails += 1
            continue

        pam_ok = False
        detail: list[str] = []

        for h in hits_fwd:
            prev = pam_before(g, h)
            ok = (prev == exp_pam)
            pam_ok = pam_ok or ok
            detail.append(f"+{h:,} PAM={prev}{'OK' if ok else 'X'}")

        for h in hits_rev:
            # 负链命中：反向互补序列位于正链 [h, h+23)。
            # 该链上的 PAM 在正链 h+23 处，取其反向互补即为该链 PAM。
            after = g.seq[h + len(spacer): h + len(spacer) + 3]
            prev_rc = revcomp(after) if len(after) == 3 else "?"
            ok = (prev_rc == exp_pam)
            pam_ok = pam_ok or ok
            detail.append(f"(rc)@{h:,} PAM={prev_rc}{'OK' if ok else 'X'}")

        if not pam_ok:
            fails += 1
        flag = "PASS" if pam_ok else "FAIL"
        print(f"  [{flag}] {name:<8} {site:<6} {loc:<16} "
              f"命中 {len(hits_fwd) + len(hits_rev)} 处: {'; '.join(detail[:3])}")

        if args.dump and hits_fwd:
            h = hits_fwd[0]
            s = max(0, h - args.dump)
            e = min(len(g), h + len(spacer) + args.dump)
            print(f"         上下文 [{s:,}-{e:,}]:")
            print(f"           {g.seq[s:h]}|{spacer}|{g.seq[h+len(spacer):e]}")

    # ---------- 检查 2: 工具自身的扫描一致性 ----------
    print("\n[2] 工具自身扫描一致性（scan_sgrnas 应能重新发现锚点）")
    sp26 = "ATTCGGAAGAAGACGTGGTATCA"
    h26 = g.find_all(sp26)
    if h26:
        anchor = h26[0]
        cands = dedupe_sgrnas(scan_sgrnas(g.seq, anchor, anchor + len(sp26), flank=50))
        found = [c for c in cands if c.spacer == sp26]
        if found:
            c = found[0]
            ok = (c.pam == "GTT") and (c.pam_start == anchor - 3)
            print(f"  [{'PASS' if ok else 'FAIL'}] sgRNA26 被重新发现: "
                  f"PAM={c.pam} @ {c.pam_start:,}（期望 GTT @ {anchor - 3:,}）")
            if not ok:
                fails += 1
        else:
            print(f"  [FAIL] scan_sgrnas 未能发现 sgRNA26（扫描 {len(cands)} 条候选）")
            fails += 1
    else:
        print("  [SKIP] sgRNA26 未找到")

    # ---------- 检查 3: 插入窗口 ----------
    print("\n[3] 插入窗口预测")
    if h26:
        pam_start = h26[0] - 3
        for lo, hi in (INSERT_WINDOW, (56, 63)):
            mid = pam_start + 3 + (lo + hi) // 2
            print(f"  PAM 起点 {pam_start:,} → 窗口 {lo}-{hi} bp → 预测插入位点 ≈ {mid:,}")
        print("  两个窗口来自论文不同章节（Figure 2 vs Figure 4），需实测确认")

    print(f"\n{'=' * 62}")
    print(f"结果: {'ALL PASS' if fails == 0 else f'{fails} 项未通过'}")
    print(f"{'=' * 62}")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
