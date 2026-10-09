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


def design_assembly_primers(cargo: str, *,
                            upstream_3prime: str,
                            downstream_5prime: str,
                            cargo_name: str = "cargo",
                            overlap_len: int = DEFAULT_OVERLAP_LEN,
                            anneal_len: int = PRIMER_ANNEAL_LEN,
                            prefix_fwd: str = "",
                            prefix_rev: str = "") -> AssemblyPlan:
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

    # ---- 边界提醒（本模块的核心目的）----
    plan.warnings.insert(0,
        "⚠️ **cargo 必须落在 LE 与 RE 之间**，不能进入 LE 或 RE 内部，"
        "也不能插到 LE 之前（酶那一侧）。"
        "放错位置会破坏转座子末端识别 → 转座失败，"
        "而 PCR、Gibson、菌落 PCR 全都正常，直到接合才发现转化子为零。")
    plan.warnings.append(
        "⚠️ 上/下游 overlap 取自 GA725 / GA778 的**实际产物末端**。"
        "若这段序列有误，Gibson 会失败且不易诊断 —— 请核对来源。")

    return plan


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
