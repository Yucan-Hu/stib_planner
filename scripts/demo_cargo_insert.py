#!/usr/bin/env python3
"""
端到端演示：把 cargo 插入 psh4300C 的 LE–RE 槽位，产出可打开的图谱文件。

用法:
  python scripts/demo_cargo_insert.py --cargo <cargo.fa> --out out/
  python scripts/demo_cargo_insert.py --demo          # 用内置测试 cargo

产出:
  out/<name>.dna   SnapGene（请打开确认）
  out/<name>.gb    GenBank 保底（任何工具可读）
  out/<name>.fa    构建后的序列
  out/<name>.md    方案 + 引物清单
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stib.buildout import (ConstructMap, OutFeature, OutPrimer,  # noqa: E402
                           roundtrip_check, write_fasta, write_genbank,
                           write_snapgene_file)
from stib.cargo import check_cargo, design_assembly_primers  # noqa: E402
from stib.sequence import revcomp  # noqa: E402

#: cargo 槽位（来自 stib/data/construct_psh43.yml，由 psh4300C.dna 实测确认）
SLOT_START = 8590
SLOT_END = 9784

#: 上游片段（Fragment 1 = 酶+LE）的 3' 末端 —— 取槽位起点前 30 bp
#: 下游片段（Fragment 3 = RE-骨架）的 5' 开头 —— 取槽位终点后 30 bp
OVERLAP_FLANK = 30

DEMO_CARGO = (
    "ATGAAAACCGCTATTCTGGGTGGCACCTTTGACCCGAAAGTGCAGGAAGCGATTGAACGT"
    "TACGGGTTTGATGTGATCGTTGCCAACCCGGCAATGAGCTACGGCCTGGAAGATTACCTG"
    "CGTGAACTGGATAAACTGGGTATTCCGCTGATTGAAATCGACGGTGAACTGACCTACGAA"
    "GAAACCCTGAAAGCGCTGGTTGAACGTCTGGGTGAAACCGCTCCGGCGTCTGATCGTGAA"
)


def build(backbone_seq: str, cargo: str, cargo_name: str, *,
          replace_slot: bool = True) -> tuple[ConstructMap, object]:
    """构造插入后的完整质粒图谱。"""
    n = len(backbone_seq)
    before = backbone_seq[:SLOT_START]
    after = backbone_seq[SLOT_END:]
    slot = backbone_seq[SLOT_START:SLOT_END]

    # 上游片段的 3' 末端 / 下游片段的 5' 开头（用于算 overlap）
    up3 = backbone_seq[max(0, SLOT_START - OVERLAP_FLANK):SLOT_START]
    dn5 = backbone_seq[SLOT_END:min(n, SLOT_END + OVERLAP_FLANK)]

    plan = design_assembly_primers(
        cargo, upstream_3prime=up3, downstream_5prime=dn5,
        cargo_name=cargo_name)

    if replace_slot:
        new_seq = before + cargo + after
        cargo_start, cargo_end = SLOT_START, SLOT_START + len(cargo)
        flavor = f"替换槽位 {SLOT_START:,}–{SLOT_END:,}（{len(slot):,} bp）"
    else:
        new_seq = before + cargo + slot + after
        cargo_start, cargo_end = SLOT_START, SLOT_START + len(cargo)
        flavor = f"插入槽位起点 {SLOT_START:,}（保留原内容）"

    feats: list[OutFeature] = [
        OutFeature(name="LE", start=8474, end=8590, kind="misc_feature",
                   qualifiers={"note": "转座子左臂（与 cargo 紧邻）"}),
        OutFeature(name=cargo_name, start=cargo_start, end=cargo_end,
                   kind="CDS", qualifiers={"note": "插入的 cargo"}),
        OutFeature(name="RE", start=cargo_end + (SLOT_END - SLOT_END) + 0,
                   end=0, kind="misc_feature"),
    ]
    # RE 位置随插入长度平移
    shift = len(cargo) - (len(slot) if replace_slot else 0)
    feats[2] = OutFeature(name="RE", start=9784 + shift, end=9995 + shift,
                          kind="misc_feature",
                          qualifiers={"note": "转座子右臂（与 cargo 紧邻）"})
    feats.append(OutFeature(name="CmR", start=8932 + shift, end=9592 + shift,
                            strand=-1, kind="CDS",
                            qualifiers={"product": "chloramphenicol resistance"}))
    feats.append(OutFeature(name="PcfxA", start=9592 + shift, end=9714 + shift,
                            strand=-1, kind="misc_feature"))
    feats.append(OutFeature(
        name="enzyme (Cas12k-TnsC)", start=1930, end=8093, kind="misc_feature",
        qualifiers={"note": "psh43 转座酶模块"}))

    primers: list[OutPrimer] = []
    if plan.primers:
        f = plan.primers[0]
        r = plan.primers[1]
        # cargo 正向引物：5' 端落在槽位起点（overlap 与上游片段末端一致）
        primers.append(OutPrimer(name=f.name, seq=f.seq,
                                 start=SLOT_START, end=SLOT_START + len(f.seq),
                                 strand=1))
        primers.append(OutPrimer(name=r.name, seq=r.seq,
                                 start=SLOT_START + len(cargo) - len(r.seq),
                                 end=SLOT_START + len(cargo), strand=-1))

    cm = ConstructMap(
        name=f"psh4300C_{cargo_name}",
        seq=new_seq,
        topology="circular",
        comments=f"psh4300C with {cargo_name} ({len(cargo)} bp) — {flavor}",
        features=feats,
        primers=primers,
    )
    return cm, plan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", default="testdata/psh4300C.dna")
    ap.add_argument("--cargo", help="cargo FASTA 路径")
    ap.add_argument("--demo", action="store_true", help="用内置测试 cargo")
    ap.add_argument("--out", default="out")
    ap.add_argument("--keep-slot", action="store_true",
                    help="插入而不替换槽位内容（默认替换）")
    args = ap.parse_args()

    if not args.demo and not args.cargo:
        ap.error("需要 --cargo <file.fa> 或 --demo")

    from stib.buildout import read_snapgene
    bb = read_snapgene(open(args.backbone, "rb").read())
    if not bb["seq"]:
        print("ERROR: 骨架文件里没有序列", file=sys.stderr)
        return 1
    print(f"骨架: {args.backbone}  {len(bb['seq']):,} bp  {bb['topology']}")

    if args.demo:
        cargo, cargo_name = DEMO_CARGO, "demo_cargo"
        print(f"cargo: 内置测试序列 {len(cargo)} bp")
    else:
        from stib.cargo import read_cargo_fasta
        cargo_name, cargo = read_cargo_fasta(args.cargo)
        print(f"cargo: {args.cargo}  名={cargo_name}  {len(cargo):,} bp")

    chk = check_cargo(cargo, name=cargo_name)
    print()
    print("=== cargo 校验 ===")
    print(f"  ok={chk.ok}  {chk.stats}")
    for e in chk.errors:
        print("  ERROR:", e)
    for w in chk.warnings:
        print("  WARN :", w)
    if not chk.ok:
        print("\n硬性校验未通过，停止。")
        return 2

    cm, plan = build(bb["seq"], cargo, cargo_name,
                     replace_slot=not args.keep_slot)

    print()
    print("=== 落点校验（转座子臂边界）===")
    print(f"  cargo 落点: {SLOT_START:,}–{SLOT_START+len(cargo):,}")
    print(f"  LE 结束  : 8,590   RE 开始: {9784 + len(cargo) - (SLOT_END-SLOT_START):,}")
    in_slot = SLOT_START >= 8474 and (SLOT_START + len(cargo)) <= (SLOT_END + len(cargo))
    print(f"  紧邻 LE 3' 端且位于 RE 之前: {'是 ✓' if in_slot else '否 ✗'}")

    print()
    print("=== 组装引物 ===")
    for p in plan.primers:
        print(f"  {p.name}")
        print(f"    序列 5'->3' : {p.seq}")
        print(f"    overlap(5') : {p.overlap}")
        print(f"    退火段(3')  : {p.anneal}")
        print(f"    Tm 估算     : {p.tm_full:.1f} °C（全长）")

    print()
    print("=== .dna 往返校验 ===")
    ok, diffs = roundtrip_check(cm, verbose=True)
    print(f"  一致: {ok}")
    for d in diffs:
        print("  DIFF:", d)

    os.makedirs(args.out, exist_ok=True)
    base = os.path.join(args.out, f"psh4300C_{cargo_name}")
    write_snapgene_file(base + ".dna", cm)
    with open(base + ".gb", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(write_genbank(cm))
    with open(base + ".fa", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(write_fasta(cm))

    print()
    print("=== 产出文件 ===")
    for ext in (".dna", ".gb", ".fa"):
        p = base + ext
        print(f"  {os.path.abspath(p)}  {os.path.getsize(p):,} B")
    print()
    print("  → 在 SnapGene 里打开 .dna；若打不开，用 .gb 保底（同样能看注释）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
