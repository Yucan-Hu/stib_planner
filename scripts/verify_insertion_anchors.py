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

# ── H. 默认路线的 CmR 结论必须唯一且正确（v0.0.6 修复的 P0）──
print("\n[H] 不给 --site 时的默认路线（不能自相矛盾）")
cm_d, plan_d, warns_d = build_construct_map(
    spec, seq, feats, cargo=cargo, cargo_name=cname)
joined_d = " ".join(warns_d)
check(sch is not None, "默认调用能生成方案")
check(plan_d.scheme is not None and plan_d.scheme["edit_start"] == 8866,
      "默认位点是 after_selection（8,866–8,888）",
      f"实际 {plan_d.scheme['edit_start'] if plan_d.scheme else '?'}")
check(any(f.name == "CmR" for f in cm_d.features),
      "默认路线**保留** CmR")
say_keep = "CmR 完好" in joined_d or "保留氯霉素筛选" in joined_d
say_cut = "移除 CmR" in joined_d or "失去氯霉素筛选" in joined_d
check(say_keep and not say_cut, "默认路线只说「保留」、不说「移除」",
      f"保留={say_keep} 移除={say_cut}")
print(f"    CmR 相关结论：{'保留' if say_keep else '—'}；自相矛盾={say_keep and say_cut}")

# 显式选 after_le 时必须说「移除」，且不能再出现「完好」
cm_l, plan_l, warns_l = build_construct_map(
    spec, seq, feats, cargo=cargo, cargo_name=cname, site_id="after_le")
joined_l = " ".join(warns_l)
check(not any(f.name == "CmR" for f in cm_l.features), "after_le 确实移除 CmR")
check("移除 CmR" in joined_l, "after_le 明确报告会移除 CmR")
check("CmR 完好" not in joined_l, "after_le 不再出现「CmR 完好」",
      "仍含「CmR 完好」" if "CmR 完好" in joined_l else "")

# ── I. 锚序列归一化与方向（v0.0.6 修的两个真 bug）──
print("\n[I] 锚序列的空白与方向")
L_A = "ACATGCAACTGTCTAGAATGCGAG"
R_A = "GTGGCCTATCAAGATTACCCG"
plus_left = seq[8866 - 24:8866]
check(resolve_insertion_anchors(seq, plus_left, R_A) == (8866, 8888),
      "+链左锚（无空格）→ 8,866")
check(resolve_insertion_anchors(seq, plus_left[:12] + " " + plus_left[12:], R_A)
      == (8866, 8888),
      "内部空格不再造成边界偏移", "早期版本会 +1 bp")
check(resolve_insertion_anchors(seq, L_A, "GTGGCCTAT\nCAAGATTACCCG") == (8866, 8888),
      "换行能被识别")
_d = {}
resolve_insertion_anchors(seq, L_A, R_A, details=_d)
check(_d["left"]["strand"] == "-" and _d["right"]["strand"] == "+",
      "details 回报了每条锚的匹配链",
      f"{_d['left']['strand']}/{_d['right']['strand']}")
check(_d["left"]["boundary"] == 8866 and _d["right"]["boundary"] == 8888,
      "details 回报的边界与返回值一致")

# ── J. 整段替换走两片段方案（after_le）──
print("\n[J] after_le 的全段替换（两片段）")
cm_fe, plan_fe, warns_fe = build_construct_map(
    spec, seq, feats, cargo=cargo, cargo_name=cname, site_id="after_le")
sch_fe = plan_fe.scheme
check(sch_fe is not None and sch_fe.get("kind") == "2fragment",
      "after_le 采用两片段方案",
      f"kind={sch_fe.get('kind') if sch_fe else None}")
if sch_fe:
    n_pr = sum(len(fr["primers"]) for fr in sch_fe["fragments"])
    check(n_pr == 4, "共 4 条引物（载体 2 + cargo 2）", f"实际 {n_pr}")
    check(len(sch_fe["checks"]) == 2, "只有 2 处同源区", f"实际 {len(sch_fe['checks'])}")
    check(all(c["ok"] for c in sch_fe["checks"]), "两处同源区序列一致")
    check(all(c.get("arm_sites", 1) == 1 for c in sch_fe["checks"]),
          "两处同源臂在骨架上唯一")
    bad = [r for r in sch_fe["specificity"] if r["verdict"] != "OK"]
    check(not bad, "4 条引物的 3' 端全部唯一",
          "；".join(f"{r['name']}={r['anchor_sites']}" for r in bad) or "全部 = 1")
    check(sch_fe["total_len"] == len(seq) - (9784 - 8590) + len(cargo),
          "组装总长 = 载体 - 被替换段 + cargo", f"{sch_fe['total_len']:,}")
    used = [pr["name"] for fr in sch_fe["fragments"] for pr in fr["primers"]
            if pr.get("from_sop")]
    check(not used, "两片段方案不占用 GA725 / Cargo1 / GA778", f"实际 {used}")
    # AT 富集的接点必须被明确告知
    joined = " ".join(sch_fe["warnings"])
    check("放宽" in joined, "AT 富集接点的放宽已写明")
    check("骨架间隔" in joined, "给出了「留骨架间隔」的替代参数")
    for fr in sch_fe["fragments"]:
        for pr in fr["primers"]:
            print(f"    {pr['dir']} {pr['name']:<12} L={len(pr['seq']):<3} "
                  f"Tm={pr['tm']:.1f}")
check(not any(f.name == "CmR" for f in cm_fe.features), "after_le 确实移除了 CmR")

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
