"""
STIB Planner — cargo 校验与 Gibson 组装引物设计。

═══ 组装拓扑（来自课题组 SOP「039_STIB 载体构建.md」）═══

    Fragment 1  [ STIB 酶 ………… LE ]           GA725 PCR（6986 bp）
    Fragment 2  [        cargo       ]          用户提供的 FASTA
    Fragment 3  [ RE … PcfxA … CmR … 骨架 ]     GA778 PCR（3076 bp）

线性重组顺序：  酶 — LE — [cargo] — RE — PcfxA — CmR — 骨架

═══ ⚠️ 本模块最想防的一个错误 ═══

**把 cargo 放进 LE 或 RE 内部。**

后果：转座子末端识别被破坏 → 转座直接失败。而 PCR 条带、Gibson 转化、
菌落 PCR 全都正常 —— 直到做接合才发现转化子为零，白做两周。

本模块因此把「cargo 落点」作为**硬性校验**，而不是提示。

═══ Gibson overlap 是什么 ═══

两个相邻片段要在重叠区有一模一样的序列。所以：

  上游片段 3' 末端 ……[ overlap ]┐
                                ├─ 同一段序列（同向，非反向互补）
  下游片段 5' 末端 ……[ overlap ]┘

对 cargo 来说：
  cargo 正向引物 = [上游片段末端的 overlap] + [cargo 前 20-25 bp]
  cargo 反向引物 = [下游片段开头的 overlap 的**反向互补**] + [cargo 后 20-25 bp]

⚠️ 反向引物上加的那段是 overlap 的反向互补，因为引物本身朝左。
   这是手工设计时最容易搞反的一步。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

DEFAULT_OVERLAP_LEN = 22
MIN_OVERLAP_LEN = 15
MAX_OVERLAP_LEN = 30

#: 引物退火段建议长度（不含 overlap 尾巴）
PRIMER_ANNEAL_LEN = 22
PRIMER_TM_MIN = 55.0
PRIMER_TM_MAX = 70.0

#: cargo 中禁止出现的位点（会干扰本载体的 BsaI Golden Gate）
FORBIDDEN_MOTIFS: list[dict] = [
    {"motif": "GGTCTC", "name": "BsaI",
     "why": "本载体用 BsaI-HF2 Golden Gate 更换 sgRNA；cargo 内的 BsaI 位点会被切开"},
    {"motif": "GAGACC", "name": "BsaI(rc)",
     "why": "同上（反向互补链）"},
]


# ---------------------------------------------------------------------------
# 结果结构
# ---------------------------------------------------------------------------

@dataclass
class CargoCheck:
    ok: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    stats: dict = field(default_factory=dict)


@dataclass
class AssemblyPrimer:
    """一条组装引物。"""

    name: str
    seq: str
    #: 5' 端的 overlap 尾巴（Gibson 同源区），可能为空
    overlap: str = ""
    #: 3' 端的退火段（与模板互补的部分）
    anneal: str = ""
    tm_anneal: float = 0.0
    tm_overlap: float = 0.0
    tm_full: float = 0.0
    length: int = 0
    note: str = ""

    def as_dict(self) -> dict:
        return {
            "name": self.name, "seq": self.seq, "length": self.length,
            "overlap": self.overlap, "anneal": self.anneal,
            "tm_anneal": round(self.tm_anneal, 1),
            "tm_overlap": round(self.tm_overlap, 1),
            "tm_full": round(self.tm_full, 1),
            "note": self.note,
        }


@dataclass
class AssemblyPlan:
    """一次 cargo 插入的完整组装方案。"""

    cargo_name: str
    cargo_len: int
    cargo_gc: float
    primers: list[AssemblyPrimer] = field(default_factory=list)
    #: 载体侧扩增引物（GA725 / GA778 等）—— 来自构建规范，**同样需要下单**
    vector_primers: list[dict] = field(default_factory=list)
    #: 三片段 Gibson 的完整引物组 + 同源区检查（build_3fragment_scheme 产出）
    scheme: dict | None = None
    checks: CargoCheck | None = None
    fragments: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    steps: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "cargo": {"name": self.cargo_name, "length": self.cargo_len,
                      "gc": round(self.cargo_gc, 3)},
            "checks": ({"ok": self.checks.ok, "errors": self.checks.errors,
                        "warnings": self.checks.warnings, "stats": self.checks.stats}
                       if self.checks else None),
            "primers": [p.as_dict() for p in self.primers],
            "vector_primers": self.vector_primers,
            "fragments": self.fragments,
            "warnings": self.warnings,
            "steps": self.steps,
        }


# ---------------------------------------------------------------------------
# cargo 校验
# ---------------------------------------------------------------------------

def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def _tm_basic(seq: str) -> float:
    """
    简易 Tm（Wallace 规则 + 长度修正）。

    ⚠️ 这是**估算**，不是最近邻热力学模型。首次使用必须梯度 PCR 确认。
    本工具不引入 primer3 依赖以保持零依赖；宁可标清"估算"也不假装精确。
    """
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = _gc(seq)
    # 长度 ≥14 时用 Wallace 修正式
    if len(seq) < 14:
        at = len(seq) - seq.count("G") - seq.count("C")
        return 2 * at + 4 * (seq.count("G") + seq.count("C"))
    return 64.9 + 41.0 * (seq.count("G") + seq.count("C") - 16.4) / len(seq)


def check_cargo(cargo: str, *, name: str = "cargo") -> CargoCheck:
    """
    校验 cargo 序列是否可用于本载体的 Gibson 组装。

    硬性项不通过即 `ok=False` —— 调用方应拒绝出引物。
    """
    errors: list[str] = []
    warnings: list[str] = []
    s = cargo.upper().strip()

    if not s:
        errors.append("cargo 序列为空")
        return CargoCheck(ok=False, errors=errors)

    bad = [c for c in set(s) if c not in "ACGT"]
    if bad:
        errors.append(
            f"cargo 含非 ACGT 字符: {''.join(sorted(bad))}。"
            "若这是简并碱基或含 N，请先确定为具体序列 —— 引物无法带简并碱基设计。")

    if len(s) < 30:
        warnings.append(f"cargo 仅 {len(s)} bp，是否漏了序列？")

    # 禁止位点
    for m in FORBIDDEN_MOTIFS:
        k = m["motif"]
        hits = [i for i in range(len(s) - len(k) + 1) if s.startswith(k, i)]
        if hits:
            errors.append(
                f"cargo 含 **{m['name']} 位点** {k} × {len(hits)}"
                f"（位置 {', '.join(str(h) for h in hits[:5])}）—— {m['why']}。\n"
                "  处理办法：(a) 同义突变去掉该位点（推荐，不改氨基酸）；"
                "(b) 改用其它 Golden Gate 酶并相应更换载体。")

    gc = _gc(s)
    if gc < 0.25:
        warnings.append(f"cargo GC 仅 {gc:.0%}，偏低；PCR 可能困难，建议加 GC buffer 或 DMSO")
    if gc > 0.70:
        warnings.append(f"cargo GC 达 {gc:.0%}，偏高；同样建议优化 PCR 条件")

    # 长的同聚物
    for base in "ACGT":
        run = 1
        mx = 1
        for i in range(1, len(s)):
            run = run + 1 if s[i] == s[i - 1] else 1
            mx = max(mx, run)
        if mx >= 10:
            warnings.append(f"cargo 含 {base}×{mx} 同聚物，测序与 PCR 可能受影响")

    stats = {
        "length": len(s), "gc": round(gc, 3),
        "n_forbidden_sites": sum(
            sum(1 for i in range(len(s) - len(m["motif"]) + 1)
                if s.startswith(m["motif"], i))
            for m in FORBIDDEN_MOTIFS),
    }
    return CargoCheck(ok=not errors, errors=errors, warnings=warnings, stats=stats)


# ---------------------------------------------------------------------------
# Gibson 引物设计
# ---------------------------------------------------------------------------

def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _pick_overlap(seq: str, *, prefer: int = DEFAULT_OVERLAP_LEN,
                  from_end: str = "3prime") -> str:
    """
    从片段末端挑一段 overlap。

    from_end="3prime" → 取该片段最后 prefer bp（供下游片段作同源）
    from_end="5prime" → 取该片段最前 prefer bp
    """
    prefer = max(MIN_OVERLAP_LEN, min(MAX_OVERLAP_LEN, prefer))
    if from_end == "3prime":
        return seq[-prefer:]
    return seq[:prefer]


def normalize_vector_primers(raw: dict | None) -> list[dict]:
    """
    把构建规范（`construct_psh43.yml` 的 `primers:`）规范化成可展示的**载体扩增引物**。

    ⚠️ 存在的理由：Gibson 是**三片段**组装 —— cargo 之外还要用 GA725 / GA778
    从骨架 PCR 扩增 Fragment 1 与 Fragment 3。这两条**同样要下单合成**，
    早期版本只在流程里提名字、不给序列，导致用户漏订 → 构建卡住。

    只做字段归一，**不修改序列**：规范里的大小写（SOP 表用来区分同源臂/退火段）
    原样保留，由展示层显式提示用户核对。
    """
    out: list[dict] = []
    for name, v in (raw or {}).items():
        if not isinstance(v, dict):
            continue
        seq = str(v.get("seq", "")).strip()
        if not seq:
            continue
        ps = v.get("product_size")
        out.append({
            "name": str(name),
            "seq": seq,
            "length": len(seq),
            "tm": v.get("tm"),
            "product_size": int(ps) if isinstance(ps, (int, float)) else None,
            "anneal_pos": v.get("anneal_pos"),
            "template": str(v.get("template", "")).strip(),
            "description": str(v.get("description", "")).strip(),
            "role": str(v.get("role", "")).strip(),
            "note": str(v.get("note", "")).strip(),
            "has_lowercase": any(c.islower() for c in seq),
        })
    return out


def design_assembly_primers(cargo: str, *,
                            upstream_3prime: str,
                            downstream_5prime: str,
                            cargo_name: str = "cargo",
                            overlap_len: int = DEFAULT_OVERLAP_LEN,
                            anneal_len: int = PRIMER_ANNEAL_LEN,
                            prefix_fwd: str = "",
                            prefix_rev: str = "",
                            vector_primers: dict | None = None) -> AssemblyPlan:
    """
    为「上游片段 + cargo + 下游片段」的 Gibson 组装设计 cargo 引物。

    参数:
      cargo               cargo 序列（FASTA 读入的序列，已去换行）
      upstream_3prime     上游片段（酶+LE）的 **3' 末端**序列 —— 用于算 overlap
      downstream_5prime   下游片段（RE+骨架）的 **5' 开头**序列 —— 用于算 overlap
                          ⚠️ 这两个值必须来自实际 PCR 产物（GA725 / GA778），
                             不能靠猜；错了会导致 Gibson 失败
      overlap_len         overlap 长度（SOP / NEB 推荐 15-25）
      anneal_len          退火段长度
      prefix_fwd/rev      可选的额外 5' 尾巴（如酶切位点）；一般留空

    返回 AssemblyPlan。

    ⚠️ 本函数**不做**序列层面的个性化优化（密码子、RBS、连接肽），
    那类需求按用户具体意图单独设计。这里只保证**组装正确**。
    """
    cargo = cargo.upper().strip()
    up_ov = _pick_overlap(upstream_3prime.upper(), prefer=overlap_len,
                          from_end="3prime")
    dn_ov = _pick_overlap(downstream_5prime.upper(), prefer=overlap_len,
                          from_end="5prime")

    checks = check_cargo(cargo, name=cargo_name)

    plan = AssemblyPlan(
        cargo_name=cargo_name, cargo_len=len(cargo), cargo_gc=_gc(cargo),
        checks=checks,
    )

    if not checks.ok:
        plan.warnings.append(
            "cargo 未通过硬性校验，**不产出引物**。请先解决上面的问题。")
        return plan

    # ---- cargo 正向引物：overlap(上游) + cargo 前 anneal_len ----
    a_fwd = cargo[:anneal_len]
    fwd_seq = prefix_fwd + up_ov + a_fwd

    # ---- cargo 反向引物：overlap(下游 rc) + cargo 后 anneal_len 的 rc ----
    a_rev_template = cargo[-anneal_len:]
    a_rev = _revcomp(a_rev_template)          # 引物朝左，用反向互补
    rev_seq = prefix_rev + _revcomp(dn_ov) + a_rev

    plan.primers = [
        AssemblyPrimer(
            name=f"{cargo_name}_F",
            seq=fwd_seq,
            overlap=prefix_fwd + up_ov,
            anneal=a_fwd,
            tm_anneal=_tm_basic(a_fwd),
            tm_overlap=_tm_basic(up_ov),
            tm_full=_tm_basic(fwd_seq),
            length=len(fwd_seq),
            note="5' 端 overlap 与上游片段（酶+LE）的 3' 末端一致",
        ),
        AssemblyPrimer(
            name=f"{cargo_name}_R",
            seq=rev_seq,
            overlap=prefix_rev + _revcomp(dn_ov),
            anneal=a_rev,
            tm_anneal=_tm_basic(a_rev),
            tm_overlap=_tm_basic(_revcomp(dn_ov)),
            tm_full=_tm_basic(rev_seq),
            length=len(rev_seq),
            note="5' 端 overlap 是下游片段（RE+骨架）5' 开头的**反向互补**"
                 "—— 这一步手工设计最容易搞反",
        ),
    ]

    # ---- 自检：overlap 是否真的能对上 ----
    # 组装后拼接序列里，上游末端的 overlap 应与 cargo F 引物带来的 overlap 相同
    if not fwd_seq.startswith(prefix_fwd + up_ov):
        plan.warnings.append("内部错误：正向引物 overlap 前缀不符")
    if upstream_3prime.upper()[-len(up_ov):] != up_ov:
        plan.warnings.append("内部错误：上游 overlap 取段不一致")

    # ---- Tm 提示 ----
    for p in plan.primers:
        if p.tm_anneal < PRIMER_TM_MIN:
            plan.warnings.append(
                f"{p.name} 退火段 Tm 估算 {p.tm_anneal:.1f} °C 偏低（建议 ≥{PRIMER_TM_MIN}）"
                "—— 可把退火段加长；注意这是估算值，首次须梯度 PCR")
        if p.tm_anneal > PRIMER_TM_MAX:
            plan.warnings.append(
                f"{p.name} 退火段 Tm 估算 {p.tm_anneal:.1f} °C 偏高（建议 ≤{PRIMER_TM_MAX}）")

    # ---- 片段清单（含 SOP 的加量建议）----
    plan.fragments = [
        {"name": "Fragment 1", "desc": "STIB 酶 + LE（GA725 PCR）",
         "len": None, "pmol": 0.16,
         "note": "来自 psh4300C 模板；实际长度见 SOP（6986 bp）"},
        {"name": "Fragment 2", "desc": f"cargo（{cargo_name}）",
         "len": len(cargo), "pmol": 0.16,
         "note": "本方案设计的引物用于扩增此片段"},
        {"name": "Fragment 3", "desc": "RE + PcfxA + CmR + 骨架（GA778 PCR）",
         "len": None, "pmol": 0.08,
         "note": "长片段按 SOP 用 0.08 pmol"},
    ]

    # ---- 流程（来自 SOP）----
    plan.steps = [
        "1. 用本方案的两条 cargo 引物 PCR 扩增 cargo（高保真酶，如 Phanta Max）",
        "2. 用 GA725 / GA778 从 psh4300C 扩增 Fragment 1 与 Fragment 3",
        "3. 跑胶确认条带单一；用 DMT 酶消化 1 小时后纯化回收",
        f"4. Gibson 组装：50 °C，三片段 40–60 分钟（SOP 推荐）",
        "5. 加量：Fragment 1 / cargo 各 0.16 pmol，Fragment 3（长片段）0.08 pmol",
        "   pmol = weight(ng) / (length(bp) × 0.65)",
        "6. 转化 S17-1λpir 感受态（**必须含 pir**，R6K 载体才能复制）",
        "7. 涂 100 µg/mL 氨苄青霉素 LB 平板；测序确认 cargo 边界",
    ]

    # ---- 载体侧扩增引物（GA725 / GA778）：与 cargo 引物一样必须下单 ----
    plan.vector_primers = normalize_vector_primers(vector_primers)
    if plan.vector_primers:
        names = " / ".join(str(vp["name"]) for vp in plan.vector_primers)
        plan.steps[1] = (
            f"2. 用 **{names}** 及各自的伙伴引物扩增载体的两个片段"
            " —— 完整 6 条引物与同源检查见 §2.3")

    # ---- 边界提醒（本模块的核心目的）----
    plan.warnings.insert(0,
        "⚠️ **cargo 必须落在 LE 与 RE 之间**，不能进入 LE 或 RE 内部，"
        "也不能插到 LE 之前（酶那一侧）。"
        "放错位置会破坏转座子末端识别 → 转座失败，"
        "而 PCR、Gibson、菌落 PCR 全都正常，直到接合才发现转化子为零。")
    plan.warnings.append(
        "⚠️ cargo 引物的同源臂取自**插入位点两翼**（由本工具计算），"
        "Cargo1 / GA778 只负责扩增载体片段，不参与 cargo 替换。"
        "务必核对报告 §2.3 的同源区比对（三片段 3 处、两片段 2 处）。")

    return plan


# ---------------------------------------------------------------------------
# 三片段 Gibson：完整 6 条引物 + 三处同源区检查
# ---------------------------------------------------------------------------

#: 判定特异性的锚定长度（引物 3' 末端）。
#: PCR 延伸由 3' 端锚定 —— 只要这一段在模板上唯一，即使 5' 端有偶然匹配
#: 也不会从别处起始扩增。
PRIMER_ANCHOR_LEN = 15

#: 允许的引物长度区间。
#: ⚠️ 早期版本硬编码 20–26，这**正是** Frag3_F 出现双结合位点的直接原因：
#:    psh4300C 的 8,909 处在 20–26 这个区间内**每一档都不唯一**，
#:    而 L≥28 立刻变唯一 —— 窗口太窄，再好的打分也救不了。
PRIMER_MIN_LEN = 18
PRIMER_MAX_LEN = 32

#: 同源臂/引物里禁止出现的位点（本载体的 Golden Gate 用 BsaI-HF2）。
_BAD_MOTIFS = ("GGTCTC", "GAGACC")


def _sites(primer: str, template: str, k: int | None = None) -> list[int]:
    """
    引物（或其 3' 末端 k nt）在**双链**模板上的全部物理结合位点，返回 + 链坐标。

    一个物理位点表现为两种形式之一：引物序列出现在 + 链（引物结合 − 链模板），
    或它的反向互补出现在 + 链（引物结合 + 链模板）。两种都计一次。
    """
    probe = (primer if k is None else primer[-k:]).upper()
    if not probe:
        return []
    found: set[int] = set()
    i = template.find(probe)
    while i >= 0:
        found.add(i)
        i = template.find(probe, i + 1)
    rc = _revcomp(probe)
    # 回文探针：+ 链命中已覆盖同一物理位点，不再重复计
    if rc != probe:
        i = template.find(rc)
        while i >= 0:
            found.add(i)
            i = template.find(rc, i + 1)
    return sorted(found)


def _primer_defects(seq: str) -> list[str]:
    """返回引物序列上的硬伤（空列表 = 干净）。"""
    out: list[str] = []
    up = seq.upper()
    for m in _BAD_MOTIFS:
        if m in up:
            out.append(f"含 BsaI 位点 {m}")
    for base in "ACGT":
        if base * 5 in up:
            out.append(f"含 ≥5 连 {base}")
    if up.endswith("T"):
        out.append("3' 端为 T（延伸效率偏低）")
    return out


def _design_at(
    template: str,
    *,
    pos: int,
    direction: str,
    tm_target: float = 60.0,
    tm_range: tuple[float, float] = (56.0, 66.0),
    gc_range: tuple[float, float] = (0.40, 0.60),
    anchor: int = PRIMER_ANCHOR_LEN,
    min_len: int = PRIMER_MIN_LEN,
    max_len: int = PRIMER_MAX_LEN,
    prefer_full_unique: bool = True,
) -> dict | None:
    """
    在固定端点处设计一条**特异性经过校验**的引物。

    端点语义：
      direction='F' —— 引物的 5' 端固定在 `pos`（扩增子左边界），向 3' 方向取长度 L
      direction='R' —— 引物的 5' 端固定在 `pos`（扩增子右边界），3' 端落在 pos−L

    **硬性淘汰**（不是扣分）：
      - 3' 末端 `anchor` nt 在双链模板上必须**唯一**
      - 不得含 BsaI 位点、≥5 连碱基

    排序：先看全长是否也唯一（更稳），再按 |Tm − tm_target| 接近程度。
    返回 dict（含 sites / full_sites / site_pos，供输出特异性表），无解返回 None。
    """
    n = len(template)
    best: tuple | None = None

    for L in range(min_len, max_len + 1):
        if direction == "F":
            sub = template[pos:pos + L]
            primer = sub.upper()
        else:
            sub = template[pos - L:pos]
            primer = _revcomp(sub.upper())

        if len(sub) < L or any(c not in "ACGT" for c in primer):
            continue

        gc = _gc(primer)
        if not (gc_range[0] <= gc <= gc_range[1]):
            continue
        tm = _tm_basic(primer)
        if not (tm_range[0] <= tm <= tm_range[1]):
            continue

        defects = _primer_defects(primer)
        if defects:
            continue

        anchor_sites = _sites(primer, template, anchor)
        if len(anchor_sites) != 1:
            continue                       # ← 硬性淘汰：会从多处起始扩增

        full_sites = _sites(primer, template)
        rank = (
            0 if (prefer_full_unique and len(full_sites) == 1) else 1,
            abs(tm - tm_target),
            L,
        )
        if best is None or rank < best[0]:
            best = (rank, {
                "seq": primer,
                "length": L,
                "tm": tm,
                "gc": gc,
                "anchor_sites": len(anchor_sites),
                "anchor_pos": anchor_sites,
                "full_sites": len(full_sites),
                "full_pos": full_sites,
            })
    return None if best is None else best[1]


def _fwd_at(seq: str, x: int, **kw) -> dict | None:
    """正向引物：5' 端固定在 x（扩增子左边界）。"""
    return _design_at(seq, pos=x, direction="F", **kw)


def _rev_ending(seq: str, y: int, **kw) -> dict | None:
    """反向引物：5' 端固定在 y（扩增子右边界），3' 端落在 y−L。"""
    return _design_at(seq, pos=y, direction="R", **kw)


def _design_with_relax(template: str, *, pos: int, direction: str) -> tuple[dict | None, str | None]:
    """
    先按常规窗口设计；不行再逐级放宽，并**说明放宽了什么**。

    ⚠️ 为什么需要放宽：`after_le` 这类「整段替换」的端点落在转座子末端，
    而 LE/RE 本身 AT 富集（本载体 LE 末端 30 nt 的 GC 只有 20%）。
    端点又被构建方式钉死，不能随便挪 —— 于是只能在"放宽引物参数"与
    "留出骨架间隔"之间选一个。这里选择**放宽 + 明确报告**，
    并在报告里给出「下移几 bp 就能用常规引物」的替代方案。
    """
    got = _design_at(template, pos=pos, direction=direction)
    if got is not None:
        return got, None
    got = _design_at(template, pos=pos, direction=direction,
                     gc_range=(0.22, 0.68), tm_range=(54.0, 70.0),
                     min_len=18, max_len=45)
    if got is not None:
        return got, ("此处 GC 40–60% / Tm 56–66 °C 无解（多半是 AT 富集的转座子末端），"
                     "已放宽到 GC 22–68% / Tm 54–70 °C，并允许加长到 45 nt")
    return None, None


def _junction_alternatives(template: str, *, junction: int, direction: str,
                           span: int = 20, limit: int = 3) -> list[dict]:
    """
    在端点附近找满足**常规**窗口的替代端点。

    代价是留出「骨架间隔」：反向引物下移 → LE 与 cargo 之间多出 N bp；
    正向引物上移 → cargo 与 RE 之间多出 N bp。
    """
    out: list[dict] = []
    if direction == "R":          # 3' 端；下移 = 把端点往右挪
        rng = range(junction + 1, junction + span + 1)
    else:                          # 5' 端；上移 = 把端点往左挪
        rng = range(junction - 1, junction - span - 1, -1)
    for y in rng:
        if y <= 0:
            continue
        got = _design_at(template, pos=y, direction=direction)
        if got is None:
            continue
        out.append({"pos": y, "gap": abs(y - junction),
                    "seq": got["seq"], "length": got["length"],
                    "tm": got["tm"], "gc": got["gc"]})
        if len(out) >= limit:
            break
    return out


def build_2fragment_scheme(backbone: str, edit_start: int, edit_end: int,
                           cargo: str, *,
                           overlap_len: int = DEFAULT_OVERLAP_LEN,
                           cargo_name: str = "cargo",
                           known: dict | None = None) -> dict:
    """
    全段替换（`after_le` / `--replace-selection`）用的**两片段** Gibson。

    为什么要单独一个方案：这条路要把 LE 与 RE 之间整段换掉，两端的接点被
    构建方式**钉死**在 LE 的 3' 端与 RE 的 5' 端，不能像三片段那样把边界
    挪到好用的位置。于是：

        片段 1（载体）：backbone[edit_end:] + backbone[:edit_start]   （跨环线性）
        片段 2（cargo）：cargo 本身，两侧带同源臂

    引物：
      vec_F —— 5' 端在 edit_end（**自由端**，长度可挑）
      vec_R —— 3' 端在 edit_start（**固定端**，必须紧贴左边界）

    ⚠️ 这两条引物会落在转座子末端附近，而 LE/RE 是 AT 富集的。
    常规 GC/Tm 窗口可能无解 —— 此时会放宽并在 `warnings` 里写明，
    同时给出「端点挪 N bp 就能用常规引物」的替代参数（代价是留骨架间隔）。

    返回结构与 `build_3fragment_scheme` 对齐（fragments / checks /
    specificity / warnings / edit_start / edit_end / total_len），
    额外带 `kind: "2fragment"`。
    """
    n = len(backbone)
    warnings: list[str] = []

    vec_len = (n - edit_end) + edit_start
    if vec_len <= 0:
        raise ValueError(
            f"两片段方案不可行：替换区间 {edit_start:,}–{edit_end:,} "
            f"覆盖了整个载体，没有载体片段可扩。")

    # ---------- 载体片段的两条引物 ----------
    vec_f, relax_f = _design_with_relax(backbone, pos=edit_end, direction="F")
    vec_r, relax_r = _design_with_relax(backbone, pos=edit_start, direction="R")
    for label, got, relax, side in (
        ("vec_F", vec_f, relax_f, "右接点（RE 侧）"),
        ("vec_R", vec_r, relax_r, "左接点（LE 侧）"),
    ):
        if got is None:
            alts = _junction_alternatives(
                backbone, junction=(edit_end if side.startswith("右") else edit_start),
                direction=("F" if side.startswith("右") else "R"))
            hint = ""
            if alts:
                hint = ("可用的替代端点：" + "、".join(
                    f"{a['pos']:,}（留 {a['gap']} bp 骨架间隔，L={a['length']} "
                    f"Tm={a['tm']:.1f} GC={a['gc']*100:.0f}%）" for a in alts))
            raise ValueError(
                f"{label} 在 {side} 找不到可用引物（连放宽窗口也无解）。\n"
                f"  该接点被构建方式钉死，不能随意挪动。{hint or ''}")
        warnings.append(
            f"{label}：{side} L={got['length']} nt，Tm≈{got['tm']:.1f} °C，"
            f"GC {got['gc']*100:.0f}%，3' 端位点 "
            f"{len(_sites(got['seq'], backbone, PRIMER_ANCHOR_LEN))} 处")
        if relax:
            warnings.append(f"⚠️ {label}：{relax}")
            alts = _junction_alternatives(
                backbone, junction=(edit_end if side.startswith("右") else edit_start),
                direction=("F" if side.startswith("右") else "R"))
            if alts:
                warnings.append(
                    f"    若不接受加长/AT 富集的引物，可把该端点挪到 "
                    f"{alts[0]['pos']:,}（留 {alts[0]['gap']} bp 骨架间隔）—— "
                    f"那时 L={alts[0]['length']}、Tm={alts[0]['tm']:.1f}、"
                    f"GC {alts[0]['gc']*100:.0f}%，落在常规窗口内。\n"
                    f"    ⚠️ 代价是 LE/RE 与 cargo 之间**多出 {alts[0]['gap']} bp 骨架序列**，"
                    f"该间隔是否影响转座末端识别**未经验证**，需自行评估。")

    # ---------- cargo 片段 ----------
    up_ov = backbone[edit_start - overlap_len:edit_start]
    dn_ov = backbone[edit_end:edit_end + overlap_len]
    cargo_f = up_ov + cargo[:21]
    cargo_r = _revcomp(dn_ov) + _revcomp(cargo[-21:])

    frags = [
        {"name": "Fragment 1（载体）", "desc": "载体线性化（跨环，已去掉被替换段）",
         "span": [edit_end, n], "length": vec_len,
         "primers": [
             {"name": "vec_F", "role": "载体正向，5' 端在右接点（自由端）",
              "dir": "F", **{k: vec_f[k] for k in ("seq", "tm", "gc", "length")},
              "from_sop": False, "pos": edit_end},
             {"name": "vec_R", "role": "载体反向，3' 端在左接点（固定端）",
              "dir": "R", **{k: vec_r[k] for k in ("seq", "tm", "gc", "length")},
              "from_sop": False, "pos": edit_start},
         ]},
        {"name": "Fragment 2（cargo）", "desc": f"cargo（{cargo_name}）+ 两侧同源臂",
         "span": [edit_start, edit_end], "length": len(cargo),
         "primers": [
             {"name": f"{cargo_name}_F", "role": "cargo 正向（含上游同源臂）",
              "dir": "F", "seq": cargo_f, "length": len(cargo_f),
              "tm": _tm_basic(cargo_f), "gc": _gc(cargo_f),
              "from_sop": False, "pos": edit_start},
             {"name": f"{cargo_name}_R", "role": "cargo 反向（含下游同源臂 rc）",
              "dir": "R", "seq": cargo_r, "length": len(cargo_r),
              "tm": _tm_basic(cargo_r), "gc": _gc(cargo_r),
              "from_sop": False, "pos": edit_end},
         ]},
    ]

    checks = [
        {"junction": "J1  vec_R 3' ↔ cargo_F 5'（左接点）",
         "a": backbone[edit_start - overlap_len:edit_start],
         "b": cargo_f[:overlap_len]},
        {"junction": "J2  cargo_R 5' ↔ vec_F 5'（右接点）",
         "a": backbone[edit_end:edit_end + overlap_len],
         "b": _revcomp(cargo_r[:overlap_len])},
    ]
    for c in checks:
        c["ok"] = c["a"].upper() == c["b"].upper()
        c["len"] = len(c["a"])
        c["gc"] = _gc(c["a"])
        arms = _sites(c["a"], backbone)
        c["arm_sites"] = len(arms)
        c["arm_pos"] = arms
        if c["ok"] and c["len"] < MIN_OVERLAP_LEN:
            warnings.append(f"{c['junction']} 同源区仅 {c['len']} bp，低于推荐下限 {MIN_OVERLAP_LEN}")
        if c["ok"] and c["arm_sites"] > 1:
            warnings.append(
                f"⚠️ {c['junction']} 的同源区在骨架上出现 {c['arm_sites']} 处"
                f"（{c['arm_pos']}）—— Gibson 时有两个可退火位置，可能装配错。")
    if not all(c["ok"] for c in checks):
        warnings.append("🔴 **存在不一致的同源区** —— 按此下单 Gibson 会失败。")

    specificity: list[dict] = []
    for fr in frags:
        for pr in fr["primers"]:
            tmpl = cargo if fr["name"].startswith("Fragment 2") else backbone
            row = _spec_row(pr["name"], pr["seq"], tmpl,
                            from_sop=False,
                            template_label="cargo" if tmpl is cargo else "psh4300C")
            specificity.append(row)
            if row["verdict"] != "OK":
                warnings.append(
                    f"🔴 {pr['name']} 的 3' 端 {PRIMER_ANCHOR_LEN} nt 有 "
                    f"{row['anchor_sites']} 个结合位点 —— 会从多处起始扩增，不要下单。")

    return {"kind": "2fragment", "fragments": frags, "checks": checks,
            "specificity": specificity, "warnings": warnings,
            "overlap_shift": None,
            "edit_start": edit_start, "edit_end": edit_end,
            "total_len": vec_len + len(cargo)}


def _spec_row(name: str, seq: str, template: str, *, from_sop: bool,
              template_label: str) -> dict:
    """生成特异性表的一行：这条引物在模板上有几个结合位点。"""
    anchor = _sites(seq, template, PRIMER_ANCHOR_LEN)
    full = _sites(seq, template)
    if len(anchor) == 1:
        verdict = "OK"
    else:
        verdict = f"⚠ {len(anchor)} 个结合位点"
    return {
        "name": name,
        "seq": seq,
        "length": len(seq),
        "template": template_label,
        "anchor_len": PRIMER_ANCHOR_LEN,
        "anchor_sites": len(anchor),
        "anchor_pos": anchor,
        "full_sites": len(full),
        "verdict": verdict,
        "from_sop": from_sop,
    }


def build_3fragment_scheme(backbone: str, edit_start: int, edit_end: int,
                           cargo: str, *, ga725_start: int,
                           overlap_len: int = DEFAULT_OVERLAP_LEN,
                           cargo_name: str = "cargo",
                           ga725_seq: str | None = None,
                           known: dict | None = None) -> dict:
    """
    设计**三片段 Gibson 所需的全部 6 条引物**，并逐 junction 检查同源区。

    片段与配对（SOP 原意）：

        Fragment 1 = GA725(F) + Cargo1(R)     酶 … LE ……→ 插入位点左边界
        Fragment 2 = cargo_F + cargo_R        cargo 本身（两侧带同源臂）
        Fragment 3 = GA778(F) + Frag3_R(R)    插入位点右边界 …… RE … CmR … 骨架

    环状顺序 F1 → F2 → F3 → F1，因此有 **3 处**同源区。

    ★ 端点语义（v0.0.4 更正，重要）：
      `edit_start` = 插入位点左边界 = **Cargo1 的 3' 端**
      `edit_end`   = 插入位点右边界 = **GA778 的 5' 端**
      这两个边界**不是随意取的**：落在别处会让 F3 被迫另设计一条正向引物，
      而该区域常常不唯一 —— psh4300C 的 8,909 就是这种情况
      （8,909–8,935 与 9,711–9,737 是正向重复，20–26 nt 全部双位点）。
      用 GA778 的 5' 端（8,888）则天然唯一，且 J2 同源臂也唯一。

    所有**推导**出来的引物都经过特异性校验：3' 末端 `PRIMER_ANCHOR_LEN` nt
    在模板双链上必须唯一，否则淘汰。端点处找不到唯一引物时，
    J3 同源区长度会自动微移并在 `warnings` / `overlap_shift` 里报告。

    参数:
      ga725_start   Fragment 1 正向引物 5' 端坐标（SOP 记 1,904）
      known         构建规范里的已知引物（GA725 / Cargo1 / GA778），优先采用

    返回 {'fragments', 'checks', 'specificity', 'arms', 'warnings', ...}
    """
    n = len(backbone)
    f1_start = ga725_start
    f3_far = ga725_start + overlap_len
    warnings: list[str] = []
    frags: list[dict] = []

    # ---------- Fragment 1 正向：GA725（SOP 优先）----------
    f1_f: dict | None = None
    g7_raw = ((known or {}).get("GA725") or {}).get("seq")
    for cand_raw in (ga725_seq, g7_raw):
        if not cand_raw:
            continue
        cand = str(cand_raw).upper().replace(" ", "")
        if backbone.upper().startswith(cand, f1_start):
            f1_f = {"seq": cand, "length": len(cand), "tm": _tm_basic(cand),
                    "gc": _gc(cand), "from_sop": True}
            break
    if f1_f is None:
        got = _fwd_at(backbone, f1_start)
        if got is None:
            raise ValueError(
                f"Fragment 1 正向引物在 {f1_start:,} 处找不到**唯一**序列。\n"
                f"  该端点 3' 端 {PRIMER_ANCHOR_LEN} nt 在模板上不唯一，"
                f"或 GC/Tm 不合格。请换一个 ga725_start，或显式提供引物序列。")
        f1_f = {**got, "from_sop": False}

    # ---------- Fragment 1 反向：Cargo1（SOP 优先），3' 端落在 edit_start ----------
    f1_r: dict | None = None
    f1_r_name, f1_r_role = "Frag1_R", "插入位点左边界（反向）"
    c1_raw = ((known or {}).get("Cargo1") or {}).get("seq")
    if c1_raw:
        c1 = str(c1_raw).upper().replace(" ", "")
        if _revcomp(backbone[edit_start:edit_start + len(c1)].upper()) == c1:
            f1_r = {"seq": c1, "length": len(c1), "tm": _tm_basic(c1),
                    "gc": _gc(c1), "from_sop": True}
            f1_r_name, f1_r_role = "Cargo1", "SOP 原引物（Fragment 1 反向）"
    if f1_r is None:
        got = _rev_ending(backbone, edit_start)
        if got is None:
            raise ValueError(
                f"Fragment 1 反向引物在 {edit_start:,} 处找不到**唯一**序列。\n"
                f"  这是插入位点的左边界，通常应由 Cargo1 承担 —— "
                f"检查构建规范里 Cargo1 的序列/坐标是否与载体一致。")
        f1_r = {**got, "from_sop": False}

    frags.append({
        "name": "Fragment 1", "desc": "STIB 酶 + LE（至插入位点左边界）",
        "span": [f1_start, edit_start],
        "length": edit_start - f1_start,
        "primers": [
            {"name": "GA725", "role": "SOP 原引物" if f1_f["from_sop"] else "由序列推导",
             "dir": "F", **{k: f1_f[k] for k in ("seq", "tm", "gc", "length")},
             "from_sop": f1_f["from_sop"], "pos": f1_start},
            {"name": f1_r_name, "role": f1_r_role,
             "dir": "R", **{k: f1_r[k] for k in ("seq", "tm", "gc", "length")},
             "from_sop": f1_r["from_sop"], "pos": edit_start},
        ],
    })

    # ---------- Fragment 2：cargo ± 同源臂 ----------
    up_ov = backbone[edit_start - overlap_len:edit_start]
    dn_ov = backbone[edit_end:edit_end + overlap_len]
    cargo_f = up_ov + cargo[:21]
    cargo_r = _revcomp(dn_ov) + _revcomp(cargo[-21:])
    frags.append({
        "name": "Fragment 2", "desc": f"cargo（{cargo_name}）+ 两侧同源臂",
        "span": [edit_start, edit_end], "length": len(cargo),
        "primers": [
            {"name": f"{cargo_name}_F", "role": "cargo 正向（含上游同源臂）",
             "dir": "F", "seq": cargo_f, "tm": _tm_basic(cargo_f), "gc": _gc(cargo_f),
             "from_sop": False, "pos": edit_start},
            {"name": f"{cargo_name}_R", "role": "cargo 反向（含下游同源臂 rc）",
             "dir": "R", "seq": cargo_r, "tm": _tm_basic(cargo_r), "gc": _gc(cargo_r),
             "from_sop": False, "pos": edit_end},
        ],
    })

    # ---------- Fragment 3 正向：GA778（SOP 优先），5' 端落在 edit_end ----------
    f3_f: dict | None = None
    f3_f_name, f3_f_role = "Frag3_F", "插入位点右边界（正向）"
    g8_raw = ((known or {}).get("GA778") or {}).get("seq")
    if g8_raw:
        g8 = str(g8_raw).upper().replace(" ", "")
        if backbone[edit_end:edit_end + len(g8)].upper() == g8:
            f3_f = {"seq": g8, "length": len(g8), "tm": _tm_basic(g8),
                    "gc": _gc(g8), "from_sop": True}
            f3_f_name, f3_f_role = "GA778", "SOP 原引物（Fragment 3 正向）"
    if f3_f is None:
        got = _fwd_at(backbone, edit_end)
        if got is None:
            raise ValueError(
                f"Fragment 3 正向引物在 {edit_end:,} 处找不到**唯一**序列。\n"
                f"  常见原因：这个边界落在重复区里。建议改用 SOP 的 GA778 端点"
                f"（检查构建规范里 GA778 的 anneal_pos 是否等于 edit_end）。")
        f3_f = {**got, "from_sop": False}

    # ---------- Fragment 3 反向：端点无解时自动微移 J3 同源区 ----------
    f3_r = _rev_ending(backbone, f3_far)
    overlap_shift: str | None = None
    if f3_r is None:
        for cand_end in range(ga725_start + MIN_OVERLAP_LEN,
                              ga725_start + PRIMER_MAX_LEN + 1):
            if cand_end == f3_far:
                continue
            got = _rev_ending(backbone, cand_end)
            if got is not None:
                overlap_shift = (
                    f"J3 同源区由 {overlap_len} bp 自动调整为 "
                    f"{cand_end - ga725_start} bp —— "
                    f"原端点 {f3_far:,} 处找不到唯一的反向引物")
                f3_far = cand_end
                f3_r = got
                break
    if f3_r is None:
        raise ValueError(
            f"Fragment 3 反向引物在 {ga725_start + MIN_OVERLAP_LEN}–"
            f"{ga725_start + PRIMER_MAX_LEN} 范围内都找不到唯一序列。\n"
            f"  该区（GA725 上游）可能整体存在重复。请人工指定 Frag3_R，"
            f"或换一个 ga725_start。")
    # 推导路径的 dict 不带 from_sop —— 这里补齐，避免下游 KeyError
    f3_r.setdefault("from_sop", False)
    if overlap_shift:
        warnings.append(overlap_shift)

    ov3 = f3_far - f1_start            # J3 同源区长度（可能已被微移改动）
    frags.append({
        "name": "Fragment 3", "desc": "插入位点右边界 + RE + PcfxA + CmR + 骨架（环状）",
        "span": [edit_end, n], "length": (n - edit_end) + f3_far,
        "primers": [
            {"name": f3_f_name, "role": f3_f_role,
             "dir": "F", **{k: f3_f[k] for k in ("seq", "tm", "gc", "length")},
             "from_sop": f3_f["from_sop"], "pos": edit_end},
            {"name": "Frag3_R", "role": "与 GA725 的 5' 端共享 J3 同源区",
             "dir": "R", **{k: f3_r[k] for k in ("seq", "tm", "gc", "length")},
             "from_sop": f3_r["from_sop"], "pos": f3_far},
        ],
    })

    # ---------- 三处同源区检查 ----------
    checks = [
        {"junction": "J1  F1 3' ↔ cargo_F 5'",
         "a": backbone[edit_start - overlap_len:edit_start],
         "b": cargo_f[:overlap_len]},
        {"junction": "J2  cargo_R 5' ↔ F3 5'",
         "a": backbone[edit_end:edit_end + overlap_len],
         "b": _revcomp(cargo_r[:overlap_len])},
        {"junction": "J3  F3 3' ↔ F1 5'（环化，最易漏）",
         "a": backbone[f3_far - ov3:f3_far],
         "b": backbone[f1_start:f1_start + ov3]},
    ]
    for c in checks:
        c["ok"] = c["a"].upper() == c["b"].upper()
        c["len"] = len(c["a"])
        c["gc"] = _gc(c["a"])
        # 同源臂在骨架上是否也重复 —— 重复的臂有两个退火去处，是独立风险
        arm_sites = _sites(c["a"], backbone)
        c["arm_sites"] = len(arm_sites)
        c["arm_pos"] = arm_sites

    for c in checks:
        if c["ok"] and ("GGTCTC" in c["a"].upper() or "GAGACC" in c["a"].upper()):
            warnings.append(f"{c['junction']} 的同源区含 BsaI 位点 —— Golden Gate 换 sgRNA 时会切")
        if c["ok"] and c["len"] < MIN_OVERLAP_LEN:
            warnings.append(f"{c['junction']} 同源区仅 {c['len']} bp，低于 NEB 推荐下限 {MIN_OVERLAP_LEN}")
        if c["ok"] and c["arm_sites"] > 1:
            warnings.append(
                f"⚠️ {c['junction']} 的同源区在骨架上出现 {c['arm_sites']} 处"
                f"（{c['arm_pos']}）—— Gibson 时这段臂有两个可退火位置，"
                f"可能装配错。若能在不改变插入位点的前提下避开更好；"
                f"无法避开时请在测序阶段确认该 junction。")
    if not all(c["ok"] for c in checks):
        warnings.append("🔴 **存在不一致的同源区** —— 按此下单 Gibson 会失败。见下表红叉项。")

    # ---------- 特异性表（每条引物在**自己的模板**上有几个结合位点）----------
    specificity: list[dict] = []
    for fr in frags:
        for pr in fr["primers"]:
            tmpl_label = "cargo" if fr["name"] == "Fragment 2" else "psh4300C"
            tmpl = cargo if tmpl_label == "cargo" else backbone
            row = _spec_row(pr["name"], pr["seq"], tmpl,
                            from_sop=pr.get("from_sop", False),
                            template_label=tmpl_label)
            specificity.append(row)
            if row["verdict"] != "OK":
                warnings.append(
                    f"🔴 {pr['name']} 的 3' 端 {PRIMER_ANCHOR_LEN} nt 在"
                    f" {tmpl_label} 上有 {row['anchor_sites']} 个结合位点 —— "
                    f"会从多处起始扩增，不要下单。")

    return {"fragments": frags, "checks": checks, "specificity": specificity,
            "warnings": warnings, "overlap_shift": overlap_shift,
            "f1_start": f1_start, "f3_far": f3_far, "overlap_len_j3": ov3,
            "edit_start": edit_start, "edit_end": edit_end,
            # F1 与 F3 在 J3 处共享 ov3 bp；直接相加会虚高，必须减掉一次。
            "total_len": ((edit_start - f1_start) + len(cargo)
                          + ((n - edit_end) + f3_far) - ov3)}


# ---------------------------------------------------------------------------
# FASTA 读取（cargo 用）
# ---------------------------------------------------------------------------

def read_cargo_fasta(path: str) -> tuple[str, str]:
    """
    读 cargo FASTA，返回 (name, sequence)。

    支持多行序列、小写、含空白。**只接受一条序列** ——
    多条时报错而不是默默取第一条，因为 cargo 取错是灾难性的。
    """
    from .genome import _open_text

    entries: list[tuple[str, list[str]]] = []
    name: str | None = None
    chunks: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, chunks))
                # FASTA 惯例：标题行**第一个 token** 才是标识符，其余是描述。
                # 若整行都当名字，描述会一起进文件名（实测踩到）。
                _hdr = line[1:].strip().split()
                name = _hdr[0] if _hdr else "cargo"
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        entries.append((name, chunks))

    if not entries:
        raise SystemExit(f"ERROR: {path} 里没有 FASTA 序列（应以 '>' 开头）")
    if len(entries) > 1:
        names = "、".join(n for n, _ in entries[:6])
        raise SystemExit(
            f"ERROR: {path} 含 {len(entries)} 条序列（{names}…）。\n"
            "cargo 只能是一条序列 —— 请确认你要插入的是哪一条，"
            "或拆成多个文件分别设计。\n"
            "（不自动取第一条：cargo 取错会导致整个构建白做。）")

    nm, ch = entries[0]
    seq = "".join(ch).upper()
    return nm, seq
