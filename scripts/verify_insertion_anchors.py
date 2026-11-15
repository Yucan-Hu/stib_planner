"""
任意插入位点（锚序列 / 坐标）的行为回归。

覆盖：
  A. 用两侧锚序列能精确复现默认位点 8,866–8,888
  B. 位点落在 LE–RE 之外 → **硬拒**（转座必败，不能只告警）
  C. 位点覆盖保留元件（rrnB 终止子等）→ **告警并列出被删元件**
  D. 锚序列太短 / 多处匹配 → 报错且列出候选位置
  E. 只给一侧锚 → 报错并说明还缺什么
  F. 纯插入（start == end）允许
  G. 坐标越界 → 报错
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from stib.buildout import read_snapgene
from stib.cargo import read_cargo_fasta
from stib.cargo_assembly import (build_construct_map, resolve_insertion_anchors,
                                 validate_insertion_region)
from stib.construct import features_from_snapgene, load_spec

BACKBONE = "testdata/psh4300C.dna"
SPEC = "stib/data/construct_psh43.yml"
CARGO = "testdata/demo_cargo.fa"

rows: list[tuple[bool, str, str]] = []


def check(ok: bool, name: str, detail: str = "") -> None:
    rows.append((ok, name, detail))


rec = read_snapgene(open(BACKBONE, "rb").read())
seq = rec["seq"].upper()
spec = load_spec(SPEC)
feats = features_from_snapgene(rec)
cname, cargo = read_cargo_fasta(CARGO)

LE_END, RE_START = spec.le.core_end, spec.re.core_start
print("=" * 74)
print("任意插入位点回归（锚序列 / 坐标）")
print("=" * 74)
print(f"  载体 {len(seq):,} bp；LE 3' 端 {LE_END:,}；RE 5' 端 {RE_START:,}")

# ── A. 锚序列复现默认位点 ──
print("\n[A] 用 SOP 的两条引物当锚序列")
a, b = resolve_insertion_anchors(
    seq, "ACATGCAACTGTCTAGAATGCGAG", "GTGGCCTATCAAGATTACCCG")
check((a, b) == (8866, 8888), "锚序列复现 8,866–8,888", f"得到 {a:,}–{b:,}")
print(f"    左锚 3' 端 = {a:,}；右锚 5' 端 = {b:,}")

cm, plan, warns = build_construct_map(
    spec, seq, feats, cargo=cargo, cargo_name=cname,
    anchor_left="ACATGCAACTGTCTAGAATGCGAG",
    anchor_right="GTGGCCTATCAAGATTACCCG",
    cargo_strand=-1)
sch = plan.scheme
check(sch is not None and sch["edit_start"] == 8866 and sch["edit_end"] == 8888,
      "build_construct_map 采用锚序列给出的位点")
check(sch is not None and all(r["verdict"] == "OK" for r in sch["specificity"]),
      "该位点下 6 条引物全部唯一")
check(len(cm.seq) == len(seq) - 22 + len(cargo), "组装总长 = 载体 - 22 + cargo",
      f"{len(cm.seq):,}")

# ── B. LE–RE 之外 → 硬拒 ──
print("\n[B] 位点在 LE–RE 之外（应是硬拒）")
outside_l = seq[300:320]
outside_r = seq[330:350]
try:
    build_construct_map(spec, seq, feats, cargo=cargo, cargo_name=cname,
                        anchor_left=outside_l, anchor_right=outside_r)
    check(False, "LE–RE 之外被拒绝", "竟然没有报错")
except ValueError as exc:
    msg = str(exc)
    check("不在 LE 与 RE 之间" in msg or "LE" in msg, "LE–RE 之外被拒绝",
          msg.splitlines()[0][:70])
    print(f"    拒绝理由：{msg.splitlines()[0][:90]}")

# ── C. 覆盖保留元件 → 告警并列出 ──
print("\n[C] 位点覆盖 rrnB T1 终止子（应告警并列出元件）")
# 让区间落在 8780–8800，横跨 rrnB T1（8759–8846）
errs, ws = validate_insertion_region(spec, feats, 8780, 8800)
joined = " ".join(ws)
check(not errs, "区间本身合法（在 LE–RE 内）", "；".join(errs))
check("删掉" in joined and ("rrnB" in joined or "terminator" in joined),
      "告警列出了会被删掉的元件", joined[:80] if joined else "（无告警）")
print(f"    告警：{joined[:100]}")

# ── D. 锚太短 / 多处匹配 ──
print("\n[D] 锚序列质量问题")
try:
    resolve_insertion_anchors(seq, "ACGT", "GTGGCCTATCAAGATTACCCG")
    check(False, "过短锚被拒绝")
except ValueError as exc:
    check("太短" in str(exc), "过短锚被拒绝", str(exc).splitlines()[0][:60])

# 找一个真的重复的 12-mer
dup12 = None
for i in range(0, len(seq) - 12):
    probe = seq[i:i + 12]
    cnt = seq.count(probe)
    if cnt >= 2:
        dup12 = probe
        break
if dup12:
    try:
        resolve_insertion_anchors(seq, dup12, "GTGGCCTATCAAGATTACCCG")
        check(False, "多处匹配的锚被拒绝", dup12)
    except ValueError as exc:
        check("处匹配" in str(exc), "多处匹配的锚被拒绝",
              f"{dup12} → {str(exc).splitlines()[0][:60]}")
        print(f"    重复 12-mer {dup12}：{str(exc).splitlines()[0][:80]}")

# ── E. 只给一侧 ──
print("\n[E] 只给一侧锚")
try:
    build_construct_map(spec, seq, feats, cargo=cargo, cargo_name=cname,
                        anchor_left="ACATGCAACTGTCTAGAATGCGAG")
    check(False, "只给一侧锚时被拒绝")
except ValueError as exc:
    check("还需要" in str(exc), "只给一侧锚时被拒绝", str(exc).splitlines()[0][:60])

# ── F. 纯插入 ──
print("\n[F] 纯插入（start == end）")
errs, ws = validate_insertion_region(spec, feats, 8888, 8888)
check(not errs, "纯插入被允许")
check(any("长度为 0" in w for w in ws), "纯插入有提示")
print(f"    提示：{[w for w in ws if '长度为 0' in w][:1]}")

# ── G. 坐标越界 / 越出 LE–RE ──
print("\n[G] 坐标输入")
try:
    build_construct_map(spec, seq, feats, cargo=cargo, cargo_name=cname,
                        replace_region=(100, 50))
    check(False, "反向坐标被拒绝")
except ValueError as exc:
    check("越界" in str(exc) or "不在" in str(exc), "反向坐标被拒绝",
          str(exc).splitlines()[0][:60])

# 在 LE–RE 内但会切掉 CmR：允许，但必须**告警**（不能硬拒）
cm_, plan_, warns_ = build_construct_map(
    spec, seq, feats, cargo=cargo, cargo_name=cname,
    replace_region=(9000, 9500))
joined_ = " ".join(warns_)
check("CmR" in joined_ and ("删掉" in joined_ or "移除" in joined_),
      "坐标切掉 CmR 时给出告警", joined_[:70])
check(not any(f.name == "CmR" for f in cm_.features),
      "切掉 CmR 后图谱里确实没有 CmR")
# 而 LE–RE 之外的坐标必须硬拒
try:
    build_construct_map(spec, seq, feats, cargo=cargo, cargo_name=cname,
                        replace_region=(300, 330))
    check(False, "坐标落在 LE 之前被硬拒")
except ValueError as exc:
    check("LE" in str(exc), "坐标落在 LE 之前被硬拒", str(exc).splitlines()[0][:60])

# ── 汇总 ──
print("\n" + "=" * 74)
width = max(len(n) for _, n, _ in rows)
fails = sum(1 for ok, _, _ in rows if not ok)
for ok, n, d in rows:
    print(f"  [{'PASS' if ok else 'FAIL'}] {n:<{width}}  {d}")
print("=" * 74)
print(f"  {len(rows) - fails}/{len(rows)} passed")
print("RESULT:", "ALL PASS" if fails == 0 else f"{fails} FAILED")
sys.exit(0 if fails == 0 else 1)
