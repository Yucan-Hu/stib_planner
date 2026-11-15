"""
引物特异性回归测试 —— 锁住 2026-11 发现的真实缺陷。

背景（实测）：
  psh4300C 的 8,909–8,935 与 9,711–9,737 是一段**正向重复**。
  早期版本把插入位点右边界取在 8,909，于是：
    - F3 正向引物（Frag3_F）在 20–26 nt 全部有 **2 个**结合位点
    - J2 同源臂也重复 2 处
  后果：PCR 会从两处起始，较短的产物删掉 RE —— 而 Gibson 与菌落 PCR 都正常，
        直到接合才发现零转化子。

正确边界 = GA778 的 5' 端（8,888），它天然唯一，J2 臂也唯一。

本脚本断言：
  A. 正确边界下，6 条引物的 3' 端锚定区在各自模板上**全部唯一**
  B. 三处同源臂在骨架上**全部唯一**
  C. F3 正向引物直接来自 SOP 的 GA778（不是推导）
  D. 旧边界（8,909）下**不会**静默输出双位点引物 —— 要么自动换边界，
     要么报错；总之任何输出都必须满足 3' 端唯一
  E. 重复区本身被正确识别（作为事实记录）
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from stib.buildout import read_snapgene
from stib.cargo import (PRIMER_ANCHOR_LEN, build_3fragment_scheme,
                        read_cargo_fasta, _sites)

BACKBONE = "testdata/psh4300C.dna"
CARGO = "testdata/demo_cargo.fa"
GA725_START = 1904

rows: list[tuple[bool, str, str]] = []


def check(ok: bool, name: str, detail: str = "") -> None:
    rows.append((ok, name, detail))


rec = read_snapgene(open(BACKBONE, "rb").read())
seq = rec["seq"].upper()
name, cargo = read_cargo_fasta(CARGO)

# 构建规范里的已知引物（与真实运行一致）
known = {
    "GA725": {"seq": "GACCAGTGTCTTGTGTCTCTTC", "anneal_pos": 1904},
    "Cargo1": {"seq": "ACATGCAACTGTCTAGAATGCGAG", "anneal_pos": 8866},
    "GA778": {"seq": "GTGGCCTATcaagattacccg", "anneal_pos": 8888},
}

print("=" * 72)
print("引物特异性回归（psh4300C + demo_cargo）")
print("=" * 72)

# ── E. 重复区事实 ──
print("\n[0] 重复区事实核对")
look = seq[8909:8935]
dup = _sites(look, seq)
check(len(dup) == 2 and dup[0] == 8909,
      "8,909 处 26-mer 确有 2 个位点（8,909 与 9,711）", f"位点={dup}")
print(f"    26-mer {look} → 位点 {dup}")

# ── A/B/C. 正确边界 8,866–8,888 ──
print("\n[1] 正确边界 edit_end = 8,888（GA778 的 5' 端）")
sch = build_3fragment_scheme(
    seq, 8866, 8888, cargo,
    ga725_start=GA725_START, overlap_len=22,
    cargo_name=name, known=known)

spec = sch["specificity"]
bad = [r for r in spec if r["verdict"] != "OK"]
check(len(spec) == 6, "共 6 条引物", f"实际 {len(spec)}")
check(not bad, "6 条引物的 3' 端锚定区全部唯一",
      "；".join(f"{r['name']}={r['anchor_sites']}" for r in bad) or "全部 = 1")
for r in spec:
    print(f"    {r['name']:<14}{r['template']:<9}len={r['length']:<3} "
          f"3'端={r['anchor_sites']}  全长={r['full_sites']}  {r['verdict']}")

arms = sch["checks"]
check(all(c["ok"] for c in arms), "三处同源区序列一致")
check(all(c.get("arm_sites", 1) == 1 for c in arms), "三处同源臂在骨架上均唯一",
      "；".join(f"{c['junction'].split()[0]}={c.get('arm_sites')}" for c in arms
               if c.get("arm_sites", 1) != 1) or "全部 = 1")
for c in arms:
    print(f"    {c['junction']:<40} {c['len']} bp  臂位点 {c.get('arm_sites')}")

f3 = sch["fragments"][2]["primers"][0]
check(f3["name"] == "GA778" and f3.get("from_sop") is True,
      "F3 正向引物直接用 SOP 的 GA778", f"实际 {f3['name']} from_sop={f3.get('from_sop')}")
f1 = sch["fragments"][0]
check(f1["primers"][0]["name"] == "GA725" and f1["primers"][0].get("from_sop"),
      "F1 正向引物用 SOP 的 GA725")
check(f1["primers"][1]["name"] == "Cargo1" and f1["primers"][1].get("from_sop"),
      "F1 反向引物用 SOP 的 Cargo1")
check(sch["edit_start"] == 8866 and sch["edit_end"] == 8888,
      "插入位点为 8,866–8,888", f"{sch['edit_start']}–{sch['edit_end']}")

# ── D. 旧边界不应静默输出双位点引物 ──
print("\n[2] 旧边界 edit_end = 8,909 的行为（必须不静默输出双位点引物）")
try:
    old = build_3fragment_scheme(
        seq, 8866, 8909, cargo,
        ga725_start=GA725_START, overlap_len=22,
        cargo_name=name, known=known)
    bad_old = [r for r in old["specificity"] if r["verdict"] != "OK"]
    check(not bad_old, "旧边界下输出的引物仍然全部唯一",
          "；".join(f"{r['name']}={r['anchor_sites']}" for r in bad_old) or "全部 = 1")
    check(old["edit_end"] == 8909, "边界未被改动（8,909 本身就可行）",
          f"实际 {old['edit_end']}")
    # 8,909 处 J2 同源臂确实与 9,715 重复 —— 必须**显式告警**（不静默通过）
    j2 = [c for c in old["checks"] if c["junction"].startswith("J2")][0]
    check(j2.get("arm_sites", 1) > 1, "旧边界的 J2 同源臂确实重复",
          f"位点数 {j2.get('arm_sites')}")
    check(any("同源区在骨架上出现" in w for w in old["warnings"]),
          "重复的同源臂已写入警告", "警告中应含「同源区在骨架上出现」")
    for r in old["specificity"]:
        print(f"    {r['name']:<14}3'端={r['anchor_sites']}  {r['verdict']}")
except ValueError as exc:
    check(True, "旧边界被明确拒绝（而不是输出错引物）",
          str(exc).splitlines()[0])
    print(f"    已拒绝：{str(exc).splitlines()[0]}")

# ── 汇总 ──
print("\n" + "=" * 72)
width = max(len(n) for _, n, _ in rows)
fails = 0
for ok, n, d in rows:
    if not ok:
        fails += 1
    mark = "PASS" if ok else "FAIL"
    print(f"  [{mark}] {n:<{width}}  {d}")
print("=" * 72)
print(f"  {len(rows) - fails}/{len(rows)} passed")
print("RESULT:", "ALL PASS" if fails == 0 else f"{fails} FAILED")
sys.exit(0 if fails == 0 else 1)
