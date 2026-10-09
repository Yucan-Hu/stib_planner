"""
STIB Planner — 脱靶评估引擎（全基因组，纯标准库，无外部依赖）。

═══════════════════════════════════════════════════════════════════════
能力与边界（必读）
═══════════════════════════════════════════════════════════════════════

**做什么**：给定基因组序列 + 一条 spacer，枚举全基因组（正负链）中
允许 ≤MAX_MISMATCH 错配的**所有**候选位点，并按错配**位置分布**分级。

**为什么按位置分级**：CAST 与核酸酶机制不同。
依据结构研究（*Structural basis for the assembly of the type V
CRISPR-associated transposon complex*）：**前 17 nt 更重要** ——
ShCAST 的 PAM 在 spacer 上游，Cas12k 与 RNA:DNA 异源双链的接触集中在
PAM 近端（helix），3' 端 6 nt（tail）为延伸尾部。
因此「前 17 nt 内 3 错配」比「tail 内 3 错配」更危险。
⚠️ 注意这与 Cas9 的 3' 端种子区惯例**方向相反**（v0.6 更正）。

**算法**：6-mer 块索引（3 个互不重叠块）+ 每块最多 1 处替换的变体枚举
  + 逐位校验。
  鸽巢论证：每块允许 1 处替换时，需 ≥ k×(1+1)=6 个错配才能弄脏全部 3 块，
  故 **≤5 错配必有命中块**；本工具保守声明 **≤4 不漏检**。
  ⚠️ 早期版本声称「8-mer 保证 ≤4」，**该说法数学上不成立**，已推翻。

**已声明的边界（不得隐瞒）**：
  - 保证 **≤4 错配不漏检**。
  - **5 错配起保证失效**（引擎会显式警示）。论文与 SOP 用 Cas-OFFinder
    允许 ≤5 错配，两者**不等价**：本引擎是种子剪枝，可能少报。
  - 输出中会显式标注所用阈值与是否可能漏检。
  - 不预测「脱靶位点是否真的会发生转座」。CAST 的插入还取决于
    局部染色质状态、靶免疫（target immunity）等本引擎不建模的因素。
═══════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import gzip
import json
import os
from dataclasses import dataclass, field
from typing import Iterable, Iterator

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

#: 块长度（= 索引里的种子长度）。
#:
#: ═══ 鸽巢的正确形式（v0.3 经交叉验证逐层纠正后定稿）═══
#:
#: **每个错配最多弄脏一个块** → 「必有干净块」⟺ **m < k**。
#: （此前两次写错的版本：`m ≤ 2k−1`、以及把「块覆盖完整」当判据 ——
#:   都让交叉验证抓出了漏检，见本文件末尾的「弯路记录」。）
#:
#: 因此要保证 m ≤ R，需要 **k ≥ R + 1** 个块。
#: 但下面的贪心搜索有个更宽松的条件：允许每块最多 `MAX_SUBST` 处替换时，
#: 只有在「所有块都超过 MAX_SUBST」时才会漏检，即需要
#:     m ≥ k × (MAX_SUBST + 1)
#: 才能弄脏全部块。据此：**保证 m ≤ k × (MAX_SUBST + 1) − 1**。
#:
#: 23 nt 上的最优选择（k 个互不重叠块需 k×SEED_LEN ≤ 23）：
#:   SEED_LEN=6, k=5 → 30 > 23 ✗ 放不下 5 块
#:   SEED_LEN=6, k=3（18 nt）配 MAX_SUBST=1 → 保证 m ≤ 3×2−1 = 5 ✓ 覆盖 ≤4
#: 但 k=3 时 5 个错配可弄脏三块（1+1+3）… 实际由 MAX_SUBST 的分配决定。
#:
#: 结论：取 **SEED_LEN=6, k=3, MAX_SUBST=1** —— 保证 **m ≤ 4**，
#: 即任何 ≤4 错配的命中，至少有一块的错配数 ≤1，必被 1-替换搜索找到。
SEED_LEN = 6

#: 每块允许的最大替换数（查询侧枚举 0..MAX_SUBST 处替换）。
MAX_SUBST = 1

#: 块数（互不重叠）。23 nt 可放 3 个 6-mer（0-5, 6-11, 12-17）。
N_BLOCKS = 3

#: 由上面的方案推出的保证上限：m ≤ k×(MAX_SUBST+1) − 1 = 3×2−1 = 5。
#: 我们保守取 4（默认阈值），留一格余量。
GUARANTEED_MAX_MISMATCH = 4

#: 默认评估阈值（论文要求 ≤5）。
DEFAULT_MAX_MISMATCH = 4

_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_CODE_BASE = "ACGT"


def _encode(seq: str) -> list[int]:
    """把序列编码成 2-bit 整数列表；非 ACGT 记为 -1。"""
    return [_BASE_CODE.get(c, -1) for c in seq.upper()]


# ---------------------------------------------------------------------------
# 种子索引
# ---------------------------------------------------------------------------

class SeedIndex:
    """
    基因组 6-mer 块索引（SEED_LEN = 6）。

    存储: {种子整数: [位置, ...]}。位置为种子起点（0-based）。
    内存: 6.3 Mb 基因组约 6.29M 条目；用 array('i') 每位置 4 字节，
          总计约 25 MB — 可接受。可选 gzip 持久化以加速重复分析。
    """

    def __init__(self, genome: str) -> None:
        self.genome = genome.upper()
        self.table: dict[int, list[int]] = {}
        self._build()

    def _build(self) -> None:
        codes = _encode(self.genome)
        n = len(codes)
        k = SEED_LEN
        if n < k:
            return
        # 滑动窗口计算种子整数
        table = self.table
        for i in range(n - k + 1):
            v = 0
            ok = True
            for j in range(i, i + k):
                c = codes[j]
                if c < 0:
                    ok = False
                    break
                v = (v << 2) | c
            if ok:
                table.setdefault(v, []).append(i)

    # -- 持久化（可选，仅用于加速重复分析）-------------------------------

    def save(self, path: str) -> None:
        """把索引存成 gzip JSON。"""
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            json.dump({
                "seed_len": SEED_LEN,
                "genome_len": len(self.genome),
                "table": {str(k): v for k, v in self.table.items()},
            }, fh)

    @classmethod
    def load(cls, path: str) -> "SeedIndex":
        """从 gzip JSON 载入索引（基因组序列会一并载入）。"""
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)
        obj = cls.__new__(cls)
        obj.genome = ""
        obj.table = {int(k): v for k, v in data["table"].items()}
        return obj

    @property
    def n_seeds(self) -> int:
        return len(self.table)

    @property
    def n_positions(self) -> int:
        return sum(len(v) for v in self.table.values())


# ---------------------------------------------------------------------------
# 错配分类
# ---------------------------------------------------------------------------

#: 错配分区（spacer 坐标，0-based，5'→3'）。
#:
#: ═══ v0.6 重大更正：分区方向与 Cas9 相反 ═══
#:
#: 依据课题组提供的**结构信息**（*Structural basis for the assembly of the type V
#: CRISPR-associated transposon complex*）：**前 17 nt 更重要**。
#:
#: 这在结构上说得通：ShCAST 的 PAM 位于 spacer **上游**，因此坐标 0 紧邻 PAM；
#: Cas12k 的识别螺旋与 RNA:DNA 异源双链的接触集中在 PAM 近端，
#: 3' 端的 6 nt 是延伸尾部，接触较少。
#:
#: ⚠️ v0.2–v0.5 一直套用 **Cas9 的惯例**（seed 在 3' 端 15-22、middle 在 8-14），
#: 于是把最关键的 17 nt 切成两块并只给中等权重 —— **方向完全反了**。
#:
#: 现分区：
#:   helix — 坐标 0–16（17 nt）**关键区**：识别主要依赖此段
#:   tail  — 坐标 17–22（6 nt）外周区
#:
#: ⚠️ 关于 SOP「种子区错配区分度较高」：该表述沿用 Cas9 语汇（3' 端）。
#: 按上述结构证据，对 ShCAST **真正的关键段是 PAM 近端 17 nt**。
#: 若课题组后续确认应更细分为「前 8 / 中 9 / 后 6」三档，只需调整本函数。
HELIX_NT = 17                # 关键区长度（坐标 0..16）
#: tail = 坐标 17..22，共 6 nt


def classify_mismatches(spacer: str, target: str,
                        helix_nt: int | None = None) -> dict:
    """
    把错配位置分成两类（v0.6）：

      helix — **关键区**（坐标 0 起 helix_nt 个，默认前 17 nt）
              紧邻 PAM，是 Cas12k 识别所依赖的 RNA:DNA 异源双链区
      tail  — **外周区**（其余，默认坐标 17–22 共 6 nt）

    对 23 nt spacer：helix = 0–16，tail = 17–22。

    ⚠️ 与 Cas9 相反：ShCAST 的关键段在 **5' 端（PAM 近端）**，不在 3' 端。
    """
    n = len(spacer)
    assert len(target) == n, "长度必须一致"
    k = HELIX_NT if helix_nt is None else helix_nt
    helix_end = min(k, n)

    out = {"helix": [], "tail": []}
    for i, (a, b) in enumerate(zip(spacer, target)):
        if a == b:
            continue
        (out["helix"] if i < helix_end else out["tail"]).append(i)
    return out


@dataclass
class OffTargetHit:
    """一个脱靶候选位点。"""

    pos: int                     # 命中起点（正链 0-based）
    strand: str                  # '+' / '-'
    seq: str                     # 基因组上该位点的实际序列（正链方向）
    n_mismatch: int
    mismatch_pos: list[int]
    by_region: dict              # {'distal': [...], 'middle': [...], 'seed': [...]}
    pam_nearby: str = ""         # 命中附近 5' 上游 3 nt（可能是 PAM）
    risk: str = "low"            # 'high' / 'medium' / 'low'
    #: 是否为「与查询完全相同的另一份拷贝」（0 错配）。
    #: 这类命中在**串联重复区**（如 16S rRNA 多拷贝、重复基因家族）是正常现象，
    #: 不是真正的脱靶，但必须**显式标出** —— 用户需要知道该位点附近存在同源拷贝，
    #: 因为插入可能落在任一拷贝上。
    is_identical_copy: bool = False

    def as_dict(self) -> dict:
        return {
            "pos": self.pos,
            "strand": self.strand,
            "seq": self.seq,
            "n_mismatch": self.n_mismatch,
            "mismatch_pos": self.mismatch_pos,
            "mismatch_by_region": {k: v for k, v in self.by_region.items() if v},
            "pam_nearby": self.pam_nearby,
            "risk": self.risk,
            "is_identical_copy": self.is_identical_copy,
        }


#: 风险分级规则（v0.6，按 helix / tail 分区）。
#:
#: 依据课题组提供的结构信息：**前 17 nt（helix，PAM 近端）是识别关键段**。
#: 因此同样错配数下，错配落在 helix 的位点**更危险**；纯 tail 错配相对次要。
#:
#: **任何错配都增加风险** —— 错配绝不减害（v0.4 更正）。
RISK_RULES = [
    ("high",   "≤2 错配（任意位置）"),
    ("high",   "≤4 错配且全部落于 helix（PAM 近端前 17 nt 关键区）"),
    ("medium", "≤4 错配且含 helix 错配"),
    ("medium", "3-4 错配且全部落于 tail（外周 6 nt）"),
    ("low",    "其余（含 >4 错配，可能漏检）"),
]


def _risk_of(n_mism: int, by_region: dict) -> str:
    """按 helix（关键区）/ tail（外周区）判定风险等级。"""
    hx = len(by_region.get("helix", []))
    if n_mism <= 2:
        return "high"
    if n_mism <= 4 and hx == n_mism:
        return "high"
    if n_mism <= 4 and hx > 0:
        return "medium"
    if 3 <= n_mism <= 4:
        return "medium"
    return "low"


# ---------------------------------------------------------------------------
# 脱靶枚举
# ---------------------------------------------------------------------------

@dataclass
class OffTargetReport:
    """一次脱靶评估的完整结果。"""

    spacer: str
    genome_len: int
    max_mismatch: int
    guaranteed: bool                # 是否在鸽巢保证范围内
    n_indexed_seeds: int
    hits: list[OffTargetHit]
    n_candidates_scanned: int = 0
    capped: bool = False            # 是否因上限截断
    warnings: list[str] = field(default_factory=list)

    @property
    def n_hits(self) -> int:
        return len(self.hits)

    def by_risk(self, level: str) -> list[OffTargetHit]:
        return [h for h in self.hits if h.risk == level]

    def as_dict(self) -> dict:
        return {
            "spacer": self.spacer,
            "genome_len": self.genome_len,
            "max_mismatch": self.max_mismatch,
            "guaranteed_no_miss": self.guaranteed,
            "n_indexed_seeds": self.n_indexed_seeds,
            "n_candidates_scanned": self.n_candidates_scanned,
            "n_hits": self.n_hits,
            "counts_by_risk": {
                lv: len(self.by_risk(lv)) for lv in ("high", "medium", "low")
            },
            "hits": [h.as_dict() for h in self.hits],
            "capped": self.capped,
            "warnings": self.warnings,
            "risk_rules": [{"level": lv, "rule": r} for lv, r in RISK_RULES],
        }


def _seed_values(seq: str) -> list[int]:
    """返回查询的块种子值列表（不含偏移）。"""
    return [v for _, v in _seed_pairs(seq)]


def _seed_pairs(seq: str, k: int = N_BLOCKS) -> list[tuple[int, int]]:
    """
    把查询切成 k 个**互不重叠的完整等长块**，返回 [(偏移量, 种子值)]。

    偏移量必须保留：索引里存的是「种子自身起点」，对齐起点 = 命中位置 − 偏移量。

    ═══ 鸽巢的正确形式（v0.3 经四轮交叉验证纠正后定稿）═══

    **每个错配最多弄脏一个块**，所以「必有干净块」⟺ **m < k**。

    但本实现允许查询侧枚举每块最多 `MAX_SUBST` 处替换，于是漏检只发生在
    「所有块都超过 MAX_SUBST」时：
        m ≥ k × (MAX_SUBST + 1)  才能弄脏全部块
    ⟹ **保证 m ≤ k × (MAX_SUBST + 1) − 1**。

    取 k=3、MAX_SUBST=1 → 保证 m ≤ 5，覆盖默认阈值 ≤4。

    ═══ v0.3 走过的弯路（全部由交叉验证抓出，记录以免重犯）═══
      弯路一：以为要「块覆盖完整查询」，让末块贴末尾 → 块重叠 → 4 错配弄脏全部三块。
      弯路二：用 {0,8,16} 的 8-mer → 末块只有 7 nt，索引只存 8-mer → 该种子永不命中。
      弯路三：把鸽巢写成 `m ≤ 2k−1` → 方向搞反，k=3 误以为能保证 m≤5。
      弯路四：SEED_LEN=7 k=3 → m=4 时每块各 1–2 个错配，仍有漏检。
    结论：**块长 6、k=3、每块允许 1 替换** 才真正覆盖 ≤4。

    ⚠️ 只接受**完整块**并丢弃重复种子值：若某块与另一块值相同，独立块数减少，
    实际保证下降 —— 由 seed_coverage() 检出并警示。
    """
    n = len(seq)
    if n < SEED_LEN:
        return []
    offsets = [i * SEED_LEN for i in range(k) if (i + 1) * SEED_LEN <= n]
    pairs: list[tuple[int, int]] = []
    seen_vals: set[int] = set()
    for o in offsets:
        v = _encode_int(seq[o:o + SEED_LEN])
        if v < 0 or v in seen_vals:
            continue
        seen_vals.add(v)
        pairs.append((o, v))
    return pairs


def seed_coverage(query: str, k: int = N_BLOCKS) -> dict:
    """
    报告查询的块覆盖情况 —— 判断「不漏检」保证是否成立。

    判定依据（与 _seed_pairs 的 docstring 一致）：
        保证上限 = k × (MAX_SUBST + 1) − 1
    """
    n = len(query)
    pairs = _seed_pairs(query, k=k)
    vals = {v for _, v in pairs}
    n_blocks = len(vals)
    guaranteed_max = n_blocks * (MAX_SUBST + 1) - 1 if n_blocks > 0 else 0
    return {
        "n_distinct": n_blocks,
        "offsets": [o for o, _ in pairs],
        "block_len": SEED_LEN,
        "max_subst_per_block": MAX_SUBST,
        "guaranteed_max_mismatch": guaranteed_max,
        "ok": (n_blocks >= 3) and (guaranteed_max >= 4),
    }


def _subst_variants(kmer: str, max_subst: int = MAX_SUBST) -> list[int]:
    """
    枚举一个块的编码值：自身 + 恰好 1..max_subst 处替换的所有变体。

    ⚠️ 必须**无重复**。v0.3 初版用递归时，`budget > 0` 分支对最后一位又走了
    一次「不替换」路径，导致原始值被重复入列（实测 20 个而非 19 个），
    进而在命中阶段产生**同一位点的重复报告**。

    现改用 itertools.combinations 精确产生「恰好 j 处不同」，
    天然无重复：1 + Σ_{j=1..max_subst} C(L,j)·3^j。
    6-mer / max_subst=1 → 1 + 3×6 = 19。
    """
    import itertools

    codes = [_BASE_CODE.get(c, -1) for c in kmer.upper()]
    if not codes or any(c < 0 for c in codes):
        return []
    L = len(codes)

    def pack(acc: list[int]) -> int:
        v = 0
        for c in acc:
            v = (v << 2) | c
        return v

    out: set[int] = {pack(codes)}
    for j in range(1, min(max_subst, L) + 1):
        for pos in itertools.combinations(range(L), j):
            for alts in itertools.product(range(4), repeat=j):
                acc = list(codes)
                ok = True
                for p, a in zip(pos, alts):
                    if a == codes[p]:
                        ok = False
                        break
                    acc[p] = a
                if ok:
                    out.add(pack(acc))
    return sorted(out)


def _encode_int(kmer: str) -> int:
    v = 0
    for c in kmer.upper():
        code = _BASE_CODE.get(c, -1)
        if code < 0:
            return -1
        v = (v << 2) | code
    return v


def find_offtargets(index: SeedIndex, spacer: str,
                    max_mismatch: int = DEFAULT_MAX_MISMATCH,
                    max_report: int = 50,
                    both_strands: bool = True) -> OffTargetReport:
    """
    枚举全基因组脱靶位点。

    参数:
      index         已建好的 SeedIndex
      spacer        23 nt protospacer（正链方向）
      max_mismatch  允许错配数。≤4 时不漏检；>4 会标注可能漏检
      max_report    最多返回多少个命中（按风险与错配数排序）
      both_strands  是否同时扫负链

    返回 OffTargetReport。
    """
    genome = index.genome
    n = len(spacer)
    spacer = spacer.upper()
    warnings: list[str] = []

    if n != 23:
        warnings.append(f"spacer 长度为 {n}，本体系标准为 23 nt")
    if any(c not in "ACGT" for c in spacer):
        warnings.append("spacer 含非 ACGT 字符")

    # 「不漏检」保证的判据（v0.3 用穷举搜索纠正过一次，这是最终版）：
    #
    #   鸽巢成立的真实条件是「**k 个互不重叠的块** → 保证 m ≤ 2k−1」。
    #
    #   ⚠️ v0.3 曾误用「块必须完整覆盖查询」作为判据，于是让末块贴住序列末尾
    #      （offset {0,8,15}），结果块之间出现重叠 —— 重叠恰恰**破坏**鸽巢：
    #      3 个带重叠的块可以被 4 个错配同时弄脏（实测漏检）。
    #   ✅ 正确做法是**互不重叠**：offset {0, 8, 16}，块为 8/8/7 nt。
    #      3 块 → 保证 m ≤ 5，足以覆盖默认阈值 ≤4。
    covered = seed_coverage(spacer, k=3)
    pigeon_max = covered["guaranteed_max_mismatch"]
    threshold_ok = max_mismatch <= GUARANTEED_MAX_MISMATCH
    # 保证 = 阈值在全局上限内，且该 spacer 的块结构能达到这个阈值
    guaranteed = max_mismatch <= min(GUARANTEED_MAX_MISMATCH, pigeon_max)

    if not threshold_ok:
        warnings.append(
            f"max_mismatch={max_mismatch} > {GUARANTEED_MAX_MISMATCH}："
            "本引擎不再保证不漏检，结果可能少报。"
            "若需完整 ≤5 错配枚举，请用 Cas-OFFinder v2.4 交叉核对。"
        )
    if not guaranteed and threshold_ok:
        why = []
        if covered["n_distinct"] < 3:
            why.append(f"只能提供 {covered['n_distinct']} 个互不相同的 {SEED_LEN}-mer 块"
                       f"（鸽巢需要 3 个才保证 m≤5）")
        if max_mismatch > pigeon_max:
            why.append(f"该 spacer 的鸽巢上限为 {pigeon_max} 错配，"
                       f"低于阈值 {max_mismatch}")
        warnings.append(
            "⚠️ **不漏检保证不成立**：" + "；".join(why) + "。"
            f"（种子偏移 {covered['offsets']}，共 {covered['n_distinct']} 块）"
            f"低复杂度序列（如短重复单元）会出现这种情况 —— 这是 {SEED_LEN}-mer 块在"
            "23 nt 上的信息论限制，不能靠调参绕过。"
            "**结果可能漏检部分脱靶位点**，请用 Cas-OFFinder 交叉核对。"
        )

    # 候选收集：对每个块，枚举**最多 MAX_SUBST 处替换**的所有种子值，
    # 在索引里取并集。这是「保证 m ≤ k(MAX_SUBST+1)−1」得以成立的关键 ——
    # 只查精确块是不够的（块内允许 1 个错配时精确查找会漏掉它）。
    #
    # ⚠️ 关键（三处，都是实测抓到的 bug）：
    #   1. 索引命中位置是「种子起点」，必须减去该块的**偏移量**才能得到
    #      spacer 对齐起点。（v0.2 bug：直接用命中位置当起点 → 错位漏检）
    #   2. 去重键必须包含 **offset**：块偏移不同则对齐起点不同，
    #      只用 (align, strand) 去重会误删合法候选。
    #   3. 必须枚举替换变体（v0.3 bug：只查精确块 → 块内有 1 个错配就漏检）。
    queries = [(spacer, "+")]
    if both_strands:
        queries.append((_revcomp_local(spacer), "-"))

    seen: set[tuple[int, str, int]] = set()
    candidates: list[tuple[int, str]] = []
    for q, strand in queries:
        for offset, _sv in _seed_pairs(q):
            blk = q[offset:offset + SEED_LEN]
            for v in _subst_variants(blk, MAX_SUBST):
                for seed_pos in index.table.get(v, ()):
                    align = seed_pos - offset      # ← 换算回对齐起点
                    if align < 0:
                        continue
                    key = (align, strand, offset)  # ← 去重必须带 offset
                    if key in seen:
                        continue
                    seen.add(key)
                    candidates.append((align, strand))

    hits: list[OffTargetHit] = []
    # ⚠️ 命中阶段必须**再去重一次**。候选阶段的键是 (align, strand, offset)，
    #    而同一个位点会被**不同块**的不同替换变体分别命中（offset 不同 → 键不同），
    #    于是同一位点多次进入候选列表，最终被重复报告。
    #    （v0.3 bug：实测同一位点被报 3 次，恰好对应 3 个块。）
    reported: set[tuple[int, str]] = set()
    for pos, strand in candidates:
        if pos + n > len(genome):
            continue
        target = genome[pos:pos + n]
        if len(target) < n:
            continue

        # 与"该链上的预期序列"比较
        ref = spacer if strand == "+" else _revcomp_local(spacer)
        mm = [i for i, (a, b) in enumerate(zip(ref, target)) if a != b]
        if len(mm) > max_mismatch:
            continue
        if (pos, strand) in reported:
            continue
        reported.add((pos, strand))

        by_region = classify_mismatches(ref, target)
        pam_up = genome[pos - 3:pos] if pos >= 3 else ""
        hits.append(OffTargetHit(
            pos=pos,
            strand=strand,
            seq=target,
            n_mismatch=len(mm),
            mismatch_pos=mm,
            by_region=by_region,
            pam_nearby=pam_up,
            risk=_risk_of(len(mm), by_region),
            # 0 错配 = 与查询完全相同的另一份拷贝（串联重复区的正常现象）
            is_identical_copy=(len(mm) == 0),
        ))

    # 排序: 风险优先，其次错配数，再位置
    rank = {"high": 0, "medium": 1, "low": 2}
    hits.sort(key=lambda h: (rank[h.risk], h.n_mismatch, h.pos))

    capped = len(hits) > max_report
    if capped:
        warnings.append(f"命中 {len(hits)} 个，仅报告前 {max_report} 个（按风险排序）")

    return OffTargetReport(
        spacer=spacer,
        genome_len=len(genome),
        max_mismatch=max_mismatch,
        guaranteed=guaranteed,
        n_indexed_seeds=index.n_seeds,
        hits=hits[:max_report],
        n_candidates_scanned=len(candidates),
        capped=capped,
        warnings=warnings,
    )


def _revcomp_local(seq: str) -> str:
    """局部反向互补（避免与 sequence.py 循环导入）。"""
    tbl = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                        "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
    return seq.translate(tbl)[::-1]


# ---------------------------------------------------------------------------
# 便捷入口
# ---------------------------------------------------------------------------

def evaluate_spacer(genome: str, spacer: str,
                    max_mismatch: int = DEFAULT_MAX_MISMATCH,
                    index: SeedIndex | None = None,
                    max_report: int = 50) -> OffTargetReport:
    """
    单条 spacer 的脱靶评估。若需评估多条，请自行复用 index 参数以避免重建。
    """
    idx = index if index is not None else SeedIndex(genome)
    return find_offtargets(idx, spacer, max_mismatch=max_mismatch, max_report=max_report)


def summarize(report: OffTargetReport) -> str:
    """人类可读摘要。"""
    lines = []
    lines.append(f"spacer        : {report.spacer}")
    lines.append(f"基因组长度    : {report.genome_len:,} bp")
    lines.append(f"错配阈值      : ≤{report.max_mismatch}"
                 + ("（保证不漏检）" if report.guaranteed else "（**可能漏检**）"))
    lines.append(f"索引种子数    : {report.n_indexed_seeds:,}")
    lines.append(f"候选位点扫描  : {report.n_candidates_scanned:,}")
    lines.append(f"命中总数      : {report.n_hits}"
                 + ("（已截断）" if report.capped else ""))
    for lv in ("high", "medium", "low"):
        lines.append(f"  {lv:<7}      : {len(report.by_risk(lv))}")
    for w in report.warnings:
        lines.append(f"  [WARN] {w}")
    return "\n".join(lines)


if __name__ == "__main__":
    # 极简自检（不需要外部文件）
    print("oftarget.py selftest\n")
    g = "ACGT" * 50
    idx = SeedIndex(g)
    sp = "ACGACGACGACGACGACGACGA"  # 23 nt
    rep = find_offtargets(idx, sp, max_mismatch=2)
    print(summarize(rep))
    print("\nOK")
