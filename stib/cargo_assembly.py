"""
STIB Planner — cargo 插入的完整组装方案（v0.7）。

═══ 两种策略（课题组 v0.7 澄清）═══

**策略 A：keep_cmr（默认）—— 保留氯霉素筛选**

    LE ──[ 新 cargo ]── rrnB 终止子 ── PcfxA ── CmR ── RE ── 骨架
         ↑ 新 cargo 插在这里

  做法（SOP 的两对引物）：
     片段1 = GA725 / CARGO1   （酶 … LE … 至 cargo 5' 端）
     片段2 = cargo 本身
     片段3 = GA778 / CARGO2   （cargo 3' 端 … RE … PcfxA … CmR … 骨架）
  → CmR 留在片段3 里，**仍可用氯霉素筛选**。

**策略 B：replace_selection —— 连 CmR 表达框一起换掉**

  把 LE 与 RE 之间（8,590–9,784）全部替换为新 cargo。
  → **失去氯霉素筛选**，必须另设筛选手段（报告会显式告警）。

═══ 引物设计规则 ═══

cargo 引物带 5' 同源臂（Gibson overlap）：
  正向 = [上游末端 overlap] + cargo 前 N bp
  反向 = [下游开头 overlap 的反向互补] + cargo 后 N bp 的反向互补

⚠️ 反向引物上加的是**反向互补** —— 手工设计时最容易搞反的一步。

═══ 输出必须保留骨架全部注释 ═══

构造最终图谱时，从载体文件取**全部**特征（sgRNA 盒、酶、启动子、LE/RE、
CmR、AmpR、ori…），按插入位置**整体平移**，而不是只放新建的几个特征。
跨越插入点的元件会显式告警，不静默处理。
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

from .cargo import AssemblyPlan, check_cargo, design_assembly_primers
from .sequence import revcomp
from .construct import (ConstructSpec, features_from_snapgene,
                        features_to_construct_map, shift_features)

#: 策略名
STRATEGY_KEEP_CMR = "keep_cmr"
STRATEGY_REPLACE_SELECTION = "replace_selection"


@dataclass
class InsertionResult:
    """一次 cargo 插入的完整结果。"""

    strategy: str
    cargo_name: str
    cargo_len: int
    #: 被替换/插入的载体区间
    edit_start: int
    edit_end: int
    new_seq: str
    features: list[dict]
    primers: list[dict]
    checks: Any = None
    warnings: list[str] = field(default_factory=list)
    steps: list[str] = field(default_factory=list)
    selection_note: str = ""

    @property
    def delta(self) -> int:
        return len(self.new_seq) - 0      # 由调用方算，见 build()


def build_construct_map(spec: ConstructSpec, backbone_seq: str,
                        backbone_features: list[dict], *,
                        cargo: str, cargo_name: str,
                        strategy: str = STRATEGY_KEEP_CMR,
                        overlap_len: int = 22,
                        flank: int = 30,
                        replace_region: tuple[int, int] | None = None,
                        site_id: str | None = None,
                        cargo_strand: int | None = None):
    """
    构造插入后的图谱。

    参数:
      site_id       插入位点 id（见 spec.cargo_slot.sites）；None → 用策略默认
      cargo_strand  cargo 应处的链（+1 / -1）。None → 用位点自身建议值。
                    **若 cargo 需表达，应设为与 CmR 同向（本载体为 -1）。**

    返回 (ConstructMap, AssemblyPlan, warnings)。
    """
    n = len(backbone_seq)
    cargo = cargo.upper().strip()

    le_end = spec.le.core_end
    re_start = spec.re.core_start
    cm = spec.cmr_region()
    warnings: list[str] = []

    site = spec.site_named(site_id) if site_id else None
    if site_id and site is None:
        avail = "、".join(s.id for s in spec.cargo_slot.sites)
        raise ValueError(f"未知插入位点 {site_id!r}；可选：{avail}")
    strand = cargo_strand if cargo_strand is not None else (
        site.cargo_strand if site else 1)

    # ---- 决定编辑区间 ----
    if replace_region is not None:
        edit_start, edit_end = int(replace_region[0]), int(replace_region[1])
        strategy = "custom"
    elif site is not None:
        edit_start, edit_end = site.region
        strategy = f"site:{site.id}"
        warnings.append(f"📍 **插入位点**：{site.label}（{site.id}）"
                        f" —— 区间 {edit_start:,}–{edit_end:,}")
        if site.note:
            warnings.append(f"    位点说明：{' '.join(site.note.split())}")
    elif strategy == STRATEGY_REPLACE_SELECTION:
        edit_start, edit_end = le_end, re_start
        warnings.append(
            "⚠️ **已选择替换筛选盒**：LE 与 RE 之间的全部内容（含 PcfxA-CmR）"
            "都会被新 cargo 取代。**此后无法再用氯霉素筛选**，"
            "须另设筛选手段（如报告基因、营养缺陷互补）。")
    else:
        default_site = spec.site_named("after_le")
        if default_site is not None:
            site = default_site
            edit_start, edit_end = default_site.region
            strategy = f"site:{default_site.id}"
        else:
            edit_start, edit_end = spec.gap_between_le_and_selection()
        sel = f"{cm[0]:,}–{cm[1]:,}" if cm else "?"
        warnings.append(
            f"✅ **保留氯霉素筛选**：新 cargo 插在 LE 之后（{le_end:,}）、"
            f"筛选盒 PcfxA-CmR（{sel}）之前。CmR 不受影响。")

    # ---- 区间内已有元件：只作**提示**，不当作长度上限 ----
    #
    # ⚠️ v0.9 更正：早期版本把这段的终止子算作「占位」，进而报出
    #    「最大连续空档 91 bp / cargo 装不下」，让课题组以为长度受限。
    #    这是理解偏差：
    #      - 本工具的插入是**替换** edit 区间，cargo 长度不受区间长度限制；
    #      - 区间里原有的元件会被这一操作**一并移除**，这才是要说的事。
    #    所以现在只报告「本操作会动到哪些元件」，不再给长度上限。
    blockers = spec.blockers_in(edit_start, edit_end)
    if blockers:
        bl_txt = "；".join(
            f"{b['name']}（{b['start']:,}–{b['end']:,}，{b['end'] - b['start']} bp）"
            for b in blockers)
        warnings.append(
            f"ℹ️ **本次操作会移除区间内已有的元件**：{bl_txt}。"
            "若这些元件（转录终止子/筛选盒）需要保留，请改用其它位点。")
        if any("terminator" in b["name"].lower() for b in blockers):
            warnings.append(
                "    其中含**转录终止子** —— 移除后可能发生读通（read-through），"
                "影响下游表达，请评估。")

    # CmR 是否会落在被替换区间里
    cmr_el = next((e for e in spec.elements if e.get("name") == "CmR"), None)
    if cmr_el is not None:
        cut = not (edit_end <= int(cmr_el["start"]) or edit_start >= int(cmr_el["end"]))
        if cut:
            warnings.append(
                "🔴 **本次操作会移除 CmR**（氯霉素抗性基因，位于 "
                f"{cmr_el['start']:,}–{cmr_el['end']:,}）—— "
                "**无法再用氯霉素筛选**，须另设筛选手段。请确认这是你要的。")
        else:
            warnings.append("✅ CmR 完好，仍可用氯霉素筛选。")

    if not (0 <= edit_start < edit_end <= n):
        raise ValueError(f"编辑区间非法: [{edit_start}, {edit_end}) vs 长度 {n}")

    # ---- cargo 方向 ----
    # ⚠️ 若 cargo 需表达且要与 CmR 同向，插入序列取**反向互补**，
    #    这样它在载体上沿与 CmR 相同的方向读取。
    if strand < 0:
        inserted = revcomp(cargo)
        warnings.append(
            "➡️ **cargo 方向：负链（与 CmR 同向）** —— 插入的是所给序列的"
            "**反向互补**，使其转录方向与 CmR 一致。请确认符合你的表达意图。")
    else:
        inserted = cargo
        warnings.append("➡️ **cargo 方向：正链**。若该 cargo 需要表达且希望"
                        "与 CmR 同向，应改用 `--cargo-strand -1`。")

    # ---- 组装引物的 overlap 依据 ----
    up3 = backbone_seq[max(0, edit_start - flank):edit_start]
    dn5 = backbone_seq[edit_end:min(n, edit_end + flank)]
    if len(up3) < 15 or len(dn5) < 15:
        warnings.append(
            f"⚠️ 编辑区间两侧可用序列偏短（左 {len(up3)} bp / 右 {len(dn5)} bp），"
            "overlap 可能不够特异。")

    plan = design_assembly_primers(
        inserted, upstream_3prime=up3, downstream_5prime=dn5,
        cargo_name=cargo_name, overlap_len=overlap_len)
    warnings.extend(plan.warnings)

    # ---- 拼接新序列 ----
    new_seq = backbone_seq[:edit_start] + inserted + backbone_seq[edit_end:]
    delta = len(inserted) - (edit_end - edit_start)

    # ---- 平移骨架注释（保留全部）----
    shifted, sh_warn = shift_features(
        backbone_features, edit_start=edit_start, edit_end=edit_end, delta=delta)
    warnings.extend(sh_warn)

    # ---- 加 cargo 自身的特征 ----
    from .buildout import ConstructMap, OutFeature, OutPrimer
    feats = features_to_construct_map(shifted, len(new_seq))
    feats.append(OutFeature(
        name=cargo_name, start=edit_start, end=edit_start + len(inserted),
        strand=strand, kind="CDS",
        qualifiers={"note": (f"插入的 cargo（原序列 {len(cargo)} bp，"
                             f"{'负链' if strand < 0 else '正链'}）"
                             + ("；已取反向互补以获得与 CmR 同向"
                                if strand < 0 else ""))}))

    # ---- 引物：定位到图谱坐标 ----
    primers: list[OutPrimer] = []
    if plan.primers:
        f, r = plan.primers[0], plan.primers[1]
        primers.append(OutPrimer(
            name=f.name, seq=f.seq, start=edit_start,
            end=min(edit_start + len(f.seq), len(new_seq)), strand=strand))
        r_start = max(0, edit_start + len(inserted) - len(r.seq))
        primers.append(OutPrimer(
            name=r.name, seq=r.seq, start=r_start,
            end=min(edit_start + len(inserted), len(new_seq)), strand=-strand))

    cm = ConstructMap(
        name=f"{os.path.splitext(os.path.basename(spec.meta.get('backbone_file', 'backbone')))[0]}_{cargo_name}",
        seq=new_seq,
        topology=spec.meta.get("topology", "circular"),
        comments=(f"{cargo_name}（原 {len(cargo)} bp，strand {strand:+d}）插入 "
                  f"{edit_start:,}–{edit_end:,}（净变化 {delta:+,} bp）；{strategy}"),
        features=feats,
        primers=primers,
    )
    return cm, plan, warnings



def cargo_report(spec: ConstructSpec, cm, plan: AssemblyPlan,
                 warnings: list[str], strategy: str) -> str:
    """生成人类可读的组装方案（Markdown）。"""
    L: list[str] = []
    L.append(f"# cargo 组装方案 — {cm.name}\n")
    L.append(f"**载体**：{spec.meta.get('name', '?')}"
             f"（{spec.meta.get('backbone_file', '?')}，"
             f"{spec.meta.get('backbone_length', '?'):,} bp，"
             f"{spec.meta.get('topology', '?')}）\n")
    L.append(f"**cargo**：{plan.cargo_name}（{plan.cargo_len:,} bp，"
             f"GC {plan.cargo_gc:.0%}）\n")
    L.append(f"**策略**：`{strategy}`\n")
    L.append("\n## 1. 校验与告警\n")
    for w in warnings:
        L.append(f"- {w}\n")
    if plan.checks and plan.checks.errors:
        L.append("\n**硬性错误**（不出引物）：\n")
        for e in plan.checks.errors:
            L.append(f"- ❌ {e}\n")

    L.append("\n## 2. 组装引物（可直接下单）\n")
    L.append("| 名称 | 序列 (5'→3') | 长度 | overlap (5') | 退火段 (3') | Tm 估算 |\n")
    L.append("|---|---|---|---|---|---|\n")
    for p in plan.primers:
        L.append(f"| {p.name} | `{p.seq}` | {p.length} | `{p.overlap}` | "
                 f"`{p.anneal}` | {p.tm_full:.1f} °C |\n")
    L.append("\n> Tm 为 **估算值**（Wallace 修正式），首次使用须梯度 PCR 确认。\n")
    L.append("> 反向引物 5' 端是下游 overlap 的**反向互补** —— 手工设计最易搞反。\n")

    L.append("\n## 3. 片段与加量（来自 SOP）\n")
    L.append("| 片段 | 内容 | 长度 | 建议加量 |\n|---|---|---|---|\n")
    for fr in plan.fragments:
        ln = f"{fr['len']:,}" if fr.get("len") else "（见 SOP）"
        L.append(f"| {fr['name']} | {fr['desc']} | {ln} | {fr['pmol']} pmol |\n")
    L.append(f"\n> pmol = weight(ng) / (length(bp) × 0.65)\n")

    L.append("\n## 4. 流程\n")
    for s in plan.steps:
        L.append(f"{s}\n")

    L.append("\n## 5. 图谱内容\n")
    L.append(f"- 总长 {len(cm.seq):,} bp；特征 {len(cm.features)} 个；"
             f"引物 {len(cm.primers)} 条\n")
    L.append("\n| 元件 | 区间 | 链 |\n|---|---|---|\n")
    for f in sorted(cm.features, key=lambda x: x.start):
        L.append(f"| {f.name} | {f.start:,}–{f.end:,} | "
                 f"{'+' if f.strand > 0 else '−'} |\n")
    return "".join(L)
