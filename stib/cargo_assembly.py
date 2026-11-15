"""
STIB Planner — cargo 插入的完整组装方案（v0.7）。

═══ 两种策略（课题组 v0.7 澄清）═══

**策略 A：keep_cmr（默认）—— 保留氯霉素筛选**

    LE ──[ 新 cargo ]── rrnB 终止子 ── PcfxA ── CmR ── RE ── 骨架
         ↑ 新 cargo 插在这里

  做法（SOP 的两对引物）：
     片段1 = GA725 + Cargo1   （酶 … LE … 至插入位点**左**边界 8,866）
     片段2 = cargo 本身（两侧带同源臂）
     片段3 = GA778 + Frag3_R  （插入位点**右**边界 8,888 … RE … PcfxA … CmR … 骨架）
  → CmR 留在片段3 里，**仍可用氯霉素筛选**。
  ★ 插入位点的左右边界由 **Cargo1 的 3' 端**与 **GA778 的 5' 端**确定，
    不是任意坐标；取到重复区会让引物出现双结合位点。

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

from .cargo import (AssemblyPlan, build_3fragment_scheme, check_cargo,
                    design_assembly_primers)
from .sequence import revcomp
from .construct import (ConstructSpec, features_from_snapgene,
                        features_to_construct_map, shift_features)

#: 策略名
STRATEGY_KEEP_CMR = "keep_cmr"
STRATEGY_REPLACE_SELECTION = "replace_selection"


# ---------------------------------------------------------------------------
# 任意插入位点：用两侧「锚序列」定位
# ---------------------------------------------------------------------------
#
# 为什么需要这个：规范里的两个位点只是**两个常用落点**，不是全部可能。
# 用户想说"插在这两段序列之间"时，让他给两条紧邻插入点的序列，
# 比让他数坐标可靠得多 —— 换一个载体文件（注释版本不同）也不会错位。
#
# 语义（与 SOP 的两条引物一致）：
#   anchor_left  —— 紧邻插入点**左侧**的序列，其 **3' 端**就是插入点
#                   与 Cargo1 的语义相同（反向引物，3' 端落在左边界）
#   anchor_right —— 紧邻插入点**右侧**的序列，其 **5' 端**就是插入点
#                   与 GA778 的语义相同（正向引物，5' 端落在右边界）


def _match_sites(seq: str, probe: str) -> list[tuple[int, str]]:
    """probe 及其反向互补在 + 链上的全部匹配：[(起始坐标, '+'|'-')]。"""
    probe = probe.upper()
    out: list[tuple[int, str]] = []
    i = seq.find(probe)
    while i >= 0:
        out.append((i, "+"))
        i = seq.find(probe, i + 1)
    rc = revcomp(probe)
    if rc != probe:                      # 回文探针避免重复计
        i = seq.find(rc)
        while i >= 0:
            out.append((i, "-"))
            i = seq.find(rc, i + 1)
    return sorted(out)


def _locate_unique(seq: str, probe: str, label: str) -> tuple[int, str]:
    """把锚序列定位到**唯一**位点；0 处或多处都报错并列出。"""
    probe = probe.upper().replace(" ", "").strip()
    if len(probe) < 12:
        raise ValueError(
            f"{label} 只有 {len(probe)} nt，太短 —— 至少给 12 nt，"
            f"建议 20–30 nt（与引物退火段等长最可靠）。")
    sites = _match_sites(seq, probe)
    if not sites:
        raise ValueError(
            f"{label} 在载体上找不到：{probe}\n"
            f"  确认它来自**这个**载体文件（大小写、正反方向都可能影响）。")
    if len(sites) > 1:
        where = "、".join(f"{p:,}({s}链)" for p, s in sites)
        raise ValueError(
            f"{label} 在载体上有 {len(sites)} 处匹配，无法唯一确定插入点：{where}\n"
            f"  请给更长的序列（多 5–10 nt 通常就能区分），"
            f"或改用坐标 `--replace-region START END`。")
    return sites[0]


def resolve_insertion_anchors(backbone_seq: str,
                              anchor_left: str,
                              anchor_right: str) -> tuple[int, int]:
    """
    用两侧锚序列求出插入区间 ``[start, end)``。

    左锚的 3' 端 = 左边界；右锚的 5' 端 = 右边界。
    两条锚的方向由它们在载体上的匹配链自动判定，用户按任一方向给都行。
    """
    lp, lstrand = _locate_unique(backbone_seq, anchor_left, "左锚序列")
    rp, rstrand = _locate_unique(backbone_seq, anchor_right, "右锚序列")

    # 左锚：3' 端贴着插入点
    left = lp + len(anchor_left.strip()) if lstrand == "+" else lp
    # 右锚：5' 端贴着插入点
    right = rp if rstrand == "+" else rp + len(anchor_right.strip())

    if left > right:
        raise ValueError(
            f"两条锚序列的顺序反了：左锚推出边界 {left:,}，右锚推出 {right:,}。\n"
            f"  左锚应是插入点**上游**（5' 侧）的序列，右锚是**下游**（3' 侧）的序列。")
    return left, right


def validate_insertion_region(
    spec: ConstructSpec,
    backbone_features: list[dict],
    edit_start: int,
    edit_end: int,
) -> tuple[list[str], list[str]]:
    """
    校验自定义插入区间。返回 (hard_errors, warnings)。

    - **硬拒**：区间不在 LE 的 3' 端与 RE 的 5' 端之间 —— 那里的序列是转座子
      末端识别区，插进去必然转座失败，而且失败要到接合才暴露。
    - **告警**：区间覆盖了原本要保留的元件（CmR / PcfxA / rrnB 终止子等），
      逐条列出被删掉的元件名与坐标。
    """
    errors: list[str] = []
    warnings: list[str] = []

    le_end = spec.le.core_end
    re_start = spec.re.core_start

    if not (le_end <= edit_start <= edit_end <= re_start):
        errors.append(
            f"插入区间 {edit_start:,}–{edit_end:,} 不在 LE 与 RE 之间"
            f"（允许范围 {le_end:,}–{re_start:,}）。\n"
            f"  LE/RE 是转座子的末端识别序列，cargo 必须落在两者**之间**。"
            f"插进末端会在接合后表现为零转化子，且此前所有 PCR 都正常。")

    if edit_start == edit_end:
        warnings.append(
            f"ℹ️ 区间长度为 0（纯插入，不删除任何碱基）—— 允许，"
            f"但请确认这正是你要的效果。")

    # 被区间覆盖的元件
    hit: list[str] = []
    for f in backbone_features or []:
        if isinstance(f, dict):
            fs, fe = f.get("start"), f.get("end")
            quals = f.get("qualifiers") or {}
            raw = quals.get("label") or quals.get("gene") or f.get("name")
        else:
            fs, fe = getattr(f, "start", None), getattr(f, "end", None)
            raw = getattr(f, "name", None)
        if fs is None or fe is None:
            continue
        if fs < edit_end and fe > edit_start:
            label = raw[0] if isinstance(raw, (list, tuple)) and raw else (
                raw if isinstance(raw, str) else "未命名")
            if label in ("source", None):
                continue
            hit.append(f"{label}（{fs:,}–{fe:,}）")

    if hit:
        warnings.append(
            "⚠️ 该区间会**删掉**下列原本存在的元件：" + "、".join(hit) +
            "\n  若其中包含 PcfxA-CmR，则**不能再靠氯霉素筛选**，"
            "须确认你有替代筛选手段；若含终止子，注意转录通读风险。")

    return errors, warnings


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
                        cargo_strand: int | None = None,
                        anchor_left: str | None = None,
                        anchor_right: str | None = None):
    """
    构造插入后的图谱。

    参数:
      site_id       插入位点 id（见 spec.cargo_slot.sites）；None → 用策略默认
      replace_region 直接给坐标（高级）；给了就优先于 site_id
      anchor_left   紧邻插入点左侧的序列，其 **3' 端** = 左边界（与 Cargo1 同语义）
      anchor_right  紧邻插入点右侧的序列，其 **5' 端** = 右边界（与 GA778 同语义）
                    —— 锚序列比坐标稳：换注释版本也不会错位；给了就优先于 site_id
      cargo_strand  cargo 应处的链（+1 / -1）。None → 用位点自身建议值。
                    **若 cargo 需表达，应设为与 CmR 同向（本载体为 -1）。**

    优先级：anchor_left/right > replace_region > site_id > 策略默认。

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
    if anchor_left and anchor_right:
        edit_start, edit_end = resolve_insertion_anchors(
            backbone_seq, anchor_left, anchor_right)
        strategy = "anchors"
        errors, warns = validate_insertion_region(
            spec, backbone_features, edit_start, edit_end)
        if errors:
            raise ValueError("\n".join(errors))
        warnings.extend(warns)
        warnings.append(
            f"📍 **插入位点**：由两侧锚序列定位 —— 区间 "
            f"{edit_start:,}–{edit_end:,}（替换 {edit_end - edit_start} bp）\n"
            f"    左锚 3' 端 = {edit_start:,}；右锚 5' 端 = {edit_end:,}")
    elif anchor_left or anchor_right:
        missing = "anchor_right（右锚）" if anchor_left else "anchor_left（左锚）"
        raise ValueError(
            f"只给了一侧锚序列 —— 还需要 {missing}。\n"
            f"  插一个片段需要确定**两个**边界：\n"
            f"    左锚 = 紧邻插入点左侧的序列（其 3' 端是插入点）\n"
            f"    右锚 = 紧邻插入点右侧的序列（其 5' 端是插入点）\n"
            f"  若只想给坐标，用 --replace-region START END。")
    elif replace_region is not None:
        edit_start, edit_end = int(replace_region[0]), int(replace_region[1])
        if not (0 <= edit_start <= edit_end <= n):
            raise ValueError(
                f"坐标越界：{edit_start:,}–{edit_end:,} 不在 0–{n:,} 内。")
        strategy = "custom"
        errors, warns = validate_insertion_region(
            spec, backbone_features, edit_start, edit_end)
        if errors:
            raise ValueError("\n".join(errors))
        warnings.extend(warns)
        warnings.append(
            f"📍 **插入位点**：手动指定坐标 —— 区间 {edit_start:,}–{edit_end:,}"
            f"（替换 {edit_end - edit_start} bp）")
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
        cargo_name=cargo_name, overlap_len=overlap_len,
        vector_primers=getattr(spec, "primers", None))
    warnings.extend(plan.warnings)

    # ---- 三片段 Gibson：完整 6 条引物 + 三处同源区检查 ----
    # ⚠️ 配对：GA725 + Cargo1（Fragment 1）、GA778 + Frag3_R（Fragment 3）。
    #    插入位点两端就是 Cargo1 的 3' 端与 GA778 的 5' 端 —— 不能随手取，
    #    取错会让 F3 被迫另设正向引物，而该区常常不唯一（见 construct_psh43.yml）。
    try:
        ga = (getattr(spec, "primers", None) or {}).get("GA725") or {}
        plan.scheme = build_3fragment_scheme(
            backbone_seq, edit_start, edit_end, inserted,
            ga725_start=int(ga.get("anneal_pos") or 1904),
            overlap_len=overlap_len, cargo_name=cargo_name,
            ga725_seq=str(ga.get("seq") or "").strip() or None,
            known=getattr(spec, "primers", None))
        warnings.extend(plan.scheme["warnings"])
    except Exception as exc:  # 设计失败不应让整份方案崩溃
        warnings.append(
            f"⚠️ 三片段引物设计失败：{exc}\n"
            f"  这只影响**载体片段的扩增引物**，cargo 本身的引物仍然有效。\n"
            f"  常见原因：该边界落在重复区（本载体的 LE/RE、rrnB 终止子、"
            f"CmR 附近都有重复），使 3' 端在模板上不唯一。\n"
            f"  建议：① 改用 `after_selection` 位点（边界由 Cargo1 / GA778 决定，"
            f"已实测唯一）；② 或人工设计该端的引物后自行核对特异性。")

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
    L.append("### 2.1 cargo 引物（2 条）\n")
    L.append("| 名称 | 序列 (5'→3') | 长度 | overlap (5') | 退火段 (3') | Tm 估算 |\n")
    L.append("|---|---|---|---|---|---|\n")
    for p in plan.primers:
        L.append(f"| {p.name} | `{p.seq}` | {p.length} | `{p.overlap}` | "
                 f"`{p.anneal}` | {p.tm_full:.1f} °C |\n")
    L.append("\n> Tm 为 **估算值**（Wallace 修正式），首次使用须梯度 PCR 确认。\n")
    L.append("> 反向引物 5' 端是下游 overlap 的**反向互补** —— 手工设计最易搞反。\n")

    # ---- 2.2 载体扩增引物：三片段组装的另两条，最容易漏订 ----
    sch = getattr(plan, "scheme", None)   # §2.2 与 §2.3 都用，先取出
    if plan.vector_primers:
        n_vp = len(plan.vector_primers)
        # 本方案实际用到哪几条 SOP 引物 —— 任意位点下大部分并不适用
        used: set[str] = set()
        if sch:
            for fr in sch["fragments"]:
                for pr in fr["primers"]:
                    if pr.get("from_sop"):
                        used.add(pr["name"])
        unused = [vp["name"] for vp in plan.vector_primers if vp["name"] not in used]
        L.append(f"\n### 2.2 SOP 记录的载体引物（规范里共 {n_vp} 条）\n")
        if sch:
            L.append(f"> 本方案**实际用到**：{'、'.join(sorted(used)) if used else '（无，全部由序列推导）'}。\n")
            if unused:
                L.append(f"> 未用到（只在该引物自身的端点成立时才适用）："
                         f"{'、'.join(unused)}。\n")
        L.append("\n⚠️ **它们不是一对。** SOP 是三片段组装：\n"
                 "- **GA725 + Cargo1** 扩增 Fragment 1（Cargo1 反向，3' 端落在插入位点左边界）\n"
                 "- **GA778 + Frag3_R** 扩增 Fragment 3（GA778 正向，5' 端落在右边界）\n\n"
                 "**插入位点的左右边界就是由这两条决定的**（默认位点下为 8,866 与 8,888）。"
                 "完整 6 条引物见 §2.3。\n\n")
        L.append("| 名称 | 序列 (5'→3') | 长度 | Tm≈ | 产物 | 用途 |\n")
        L.append("|---|---|---|---|---|---|\n")
        for vp in plan.vector_primers:
            tm = "—" if vp.get("tm") is None else str(vp["tm"])
            ps = f"{vp['product_size']:,} bp" if vp.get("product_size") else "—"
            use = vp.get("description") or vp.get("role") or "—"
            L.append(f"| {vp['name']} | `{vp['seq']}` | {vp['length']} | "
                     f"{tm} | {ps} | {use} |\n")
        src = os.path.basename(getattr(spec, "path", "") or "") or "构建规范"
        L.append(f"\n> 序列、Tm 与产物长度引自构建规范 `{src}`（SOP 表），"
                 "**不是本工具估算**。\n")
        for vp in plan.vector_primers:
            if vp.get("has_lowercase"):
                L.append(f"> 📝 **{vp['name']}** 序列含小写段 —— SOP 用大小写区分"
                         "「同源臂 / 退火段」，但**本载体的 GA778 整条都与模板一致**"
                         "（8,888–8,909），并不存在人工尾巴，下单按大写即可。\n")
            if vp.get("note"):
                L.append(f"> 📝 **{vp['name']}**：{vp['note']}\n")

    # ---- 2.3 三片段引物组（6 条）+ 同源区检查 ----
    if sch:
        L.append("\n### 2.3 三片段 Gibson 引物组（**共 6 条**）与同源区检查\n")
        L.append("⚠️ **GA725 与 GA778 不是一对引物。** SOP 是三片段、两两配对：\n"
                 "Fragment 1 = **GA725 + Cargo1**；Fragment 2 = cargo 本身；"
                 "Fragment 3 = **GA778 + Frag3_R**。\n"
                 "线性重组顺序 **F1 → cargo → F3 → 回到 F1**，因此一共有 **3 处**同源区"
                 "（不是 2 处）：J1 / J2 由 cargo 引物的同源臂提供，**J3 是环化用的，最容易漏。**\n\n")
        L.append(f"产物总长（组装后）：{sch['total_len']:,} bp"
                 f"　｜　插入位点：{sch['edit_start']:,}–{sch['edit_end']:,}"
                 f"（替换 {sch['edit_end'] - sch['edit_start']} bp）\n\n")
        L.append("| 片段 | 引物 | 方向 | 序列 (5'→3') | 长度 | Tm≈ | 片段长 | 来源 |\n")
        L.append("|---|---|---|---|---|---|---|---|\n")
        for fr in sch["fragments"]:
            for pr in fr["primers"]:
                src_tag = "SOP 原引物" if pr.get("from_sop") else "由序列推导"
                L.append(f"| {fr['name']} | **{pr['name']}** | {pr['dir']} | "
                         f"`{pr['seq']}` | {len(pr['seq'])} | {pr['tm']:.1f} | "
                         f"{fr['length']:,} bp | {src_tag} |\n")

        # ---- 引物特异性：每条在自己的模板上有几个结合位点 ----
        spec_rows = sch.get("specificity") or []
        if spec_rows:
            L.append(f"\n**引物特异性检查**（3' 末端 {spec_rows[0]['anchor_len']} nt "
                     "在模板双链上的结合位点数 —— 必须为 1）\n\n")
            L.append("| 引物 | 模板 | 长度 | 3' 端位点数 | 全长位点数 | 位置 | 判定 |\n")
            L.append("|---|---|---|---|---|---|---|\n")
            for row in spec_rows:
                pos = row["anchor_pos"][0] if row["anchor_pos"] else "—"
                mark = "✅" if row["verdict"] == "OK" else f"❌ {row['verdict']}"
                L.append(f"| {row['name']} | {row['template']} | {row['length']} | "
                         f"{row['anchor_sites']} | {row['full_sites']} | {pos} | {mark} |\n")
            bad = [r for r in spec_rows if r["verdict"] != "OK"]
            if bad:
                L.append("\n> 🔴 **有引物不唯一，不要下单。** PCR 会从多个位置起始扩增，"
                         "产物是混合物；最坏情况下会得到缺失关键元件的片段，"
                         "而 Gibson 与菌落 PCR 都还能正常，直到接合才发现零转化子。\n")
            else:
                L.append("\n> ✅ 全部 6 条引物的 3' 端在各自模板上唯一。\n")

        L.append("\n**同源区检查**（逐 junction 比对，必须完全一致才能组装）\n\n")
        L.append("| Junction | 左端序列 | 右端序列 | 长度 | GC | 臂在骨架上 | 结果 |\n")
        L.append("|---|---|---|---|---|---|---|\n")
        for c in sch["checks"]:
            mark = "✅ 一致" if c["ok"] else "❌ **不一致**"
            arm = "唯一" if c.get("arm_sites", 1) == 1 else f"⚠ {c['arm_sites']} 处"
            L.append(f"| {c['junction']} | `{c['a']}` | `{c['b']}` | {c['len']} | "
                     f"{c['gc']:.0%} | {arm} | {mark} |\n")
        ok_all = all(c["ok"] for c in sch["checks"])
        L.append(f"\n> 结论：**{'全部 3 处同源区一致，可按上表下单' if ok_all else '存在不一致，先解决再下单'}**\n")
        if sch.get("overlap_shift"):
            L.append(f"\n> 🔧 **自动调整**：{sch['overlap_shift']}\n")
        f1_len = sch["fragments"][0]["length"]
        f3_len = sch["fragments"][2]["length"]
        L.append(f"\n> 📏 **片段长度对照**：SOP 记 Fragment 1 ≈ 6,986 bp、"
                 f"Fragment 3 ≈ 3,076 bp；本方案 F1 = {f1_len:,} bp、F3 = {f3_len:,} bp。"
                 "F1 差 24 bp 是 Cargo1 的引物长度（其 5' 端到 8,890），"
                 "**属预期差异**。\n")

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
