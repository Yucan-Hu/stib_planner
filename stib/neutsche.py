"""
STIB Planner — 中性位点查找与 cargo 插入设计。

═══ 用途 ═══

给**大片段 cargo**（外源功能基因、通路、报告基因）找一个安全的落点，
使插入不破坏宿主功能。这与「基因失活」是**两类不同任务**：

  基因失活（knockout）  目标 = 某个已知基因，插入落在 CDS 中段
  中性位点（knock-in）  目标 = 基因间区，插入**不落在任何基因内部**

本模块处理后者。

═══ ⚠️ 关于「已验证中性位点」的诚实声明 ═══

论文中已验证的中性位点（Bt05 @157,758、Bt09 @2,325,135 等）是**在拟杆菌里
实际做过并被证明可插入**的位点 —— 它们是**经验验证过的**，价值远高于
本模块**从序列特征预测**出来的候选。

因此：

  ✅ 优先推荐论文已验证的位点（见 `KNOWN_NEUTRAL_SITES`）
  ⚠️ 序列特征预测只能**给候选**，不能声称「已验证安全」

本模块的评分只反映可解释的序列/注释特征，**不预测插入效率**。

═══ 评分依据（每条都可审阅）═══

  1. **完全落在基因间区** —— 硬性要求。插进任何 CDS/tRNA/rRNA 都可能破坏功能
  2. **靠近基因间区中心** —— 两侧留出等长余量，减小误伤邻近元件的概率
  3. **远离必需基因** —— 若注释里标了必需性（或缺注释时用保守清单）
  4. **避开可移动元件** —— IS / 转座酶 / 前噬菌体附近易发生重排
  5. **GC 含量与基因组均值接近** —— 极端 GC 区在拟杆菌中常伴随特殊结构
  6. **（可选）同向转录** —— 插在转录单元内部会打断操纵子

以上均为启发式，**权重取整、可审阅**。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# ---------------------------------------------------------------------------
# 论文已验证的中性位点（真实数据，优先于任何预测）
# ---------------------------------------------------------------------------

#: 论文中实际使用并验证过的中性位点。
#: 坐标为 0-based，来自 scripts/verify_against_paper.py 的锚点核对（10/10 PASS）。
KNOWN_NEUTRAL_SITES: list[dict] = [
    {"label": "neutral site 44", "pos": 157_758, "strand": "+",
     "note": "论文 sgRNA44（Bt05）目标；已验证可插入"},
    {"label": "neutral site 55", "pos": 2_325_135, "strand": "+",
     "note": "论文 sgRNA55（Bt09）目标；已验证可插入"},
    {"label": "neutral site 5", "pos": 4_005_688, "strand": "-",
     "note": "论文 sgRNA47（Bt03）目标；已验证可插入（反向链）"},
]


# ---------------------------------------------------------------------------
# 结果结构
# ---------------------------------------------------------------------------

@dataclass
class NeutralSite:
    """一个候选中性位点。"""

    start: int                    # 基因间区起点（0-based）
    end: int                      # 基因间区终点（开区间）
    insert_pos: int               # 建议插入位置（0-based，sgRNA 靶向点）
    length: int                   # 基因间区长度
    score: float
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    #: 落点相对区间中心的位置（0 = 正中）
    off_center: int = 0
    #: 两侧最近的基因
    left_gene: str = ""
    right_gene: str = ""
    gc: float = 0.0
    #: 是否靠近可移动元件
    near_mobile: bool = False

    def as_dict(self) -> dict:
        return {
            "start": self.start, "end": self.end, "insert_pos": self.insert_pos,
            "length": self.length, "score": round(self.score, 1),
            "off_center": self.off_center, "gc": round(self.gc, 3),
            "left_gene": self.left_gene, "right_gene": self.right_gene,
            "near_mobile": self.near_mobile,
            "reasons": self.reasons, "warnings": self.warnings,
        }


# ---------------------------------------------------------------------------
# 评分权重（显式、取整、可审阅）
# ---------------------------------------------------------------------------

W_CENTER = 40.0            # 越靠基因间区中心越好（避免误伤邻近元件）
W_GAP_LENGTH = 30.0        # 区间越长越好（但超过阈值后不再加分）
GAP_LENGTH_SATURATION = 400
W_GC_MATCH = 20.0          # GC 与基因组均值接近
W_MOBILE_PENALTY = -50.0   # 靠近可移动元件
W_LEFT_PROXIMITY = -4.0    # 每 bp 离左侧基因末端太近
W_RIGHT_PROXIMITY = -4.0
SAFE_MARGIN = 40           # 建议插入点距基因边界至少留这么多 bp

#: 可移动元件 / 易重排特征的关键词（按注释 feature key 与 product/note 匹配）
_MOBILE_KEYS = {"mobile_element", "repeat_region"}
_MOBILE_WORDS = (
    "transposase", "transposon", "insertion sequence", "is element",
    "phage", "prophage", "recombinase", "integrase", "resolvase",
    "转座酶", "转座子", "噬菌体", "重组酶", "整合酶",
)

#: rRNA / tRNA 邻近区通常不适合（转录活跃且结构保守）
_STRUCTURAL_KEYS = {"rRNA", "tRNA", "tmRNA", "ncRNA"}


def _gc(seq: str) -> float:
    if not seq:
        return 0.0
    g = seq.count("G") + seq.count("C")
    return g / len(seq)


def _is_mobile(feat: Any) -> bool:
    if getattr(feat, "key", "") in _MOBILE_KEYS:
        return True
    text = " ".join(str(v) for v in getattr(feat, "qualifiers", {}).values()).lower()
    text += " " + str(getattr(feat, "key", "")).lower()
    return any(w in text for w in _MOBILE_WORDS)


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def find_neutral_sites(seq: str, features: list[Any], *,
                       min_gap: int = 120,
                       max_candidates: int = 12,
                       avoid: tuple[int, int] | None = None,
                       require_annotations: bool = True) -> list[NeutralSite]:
    """
    在基因组中查找候选中性位点。

    参数:
      seq                 基因组序列
      features            GenBank 特征表（**必需**：没有注释就无法判断"哪里不是基因"）
      min_gap             最小基因间区长度（bp）
      max_candidates      最多返回多少个
      avoid               额外避开的区间 (start, end)，如用户已知的敏感区
      require_annotations True 时若无注释直接报错（推荐）

    返回按分数降序的候选列表。

    ⚠️ 「中性位点」的判定**完全依赖注释的完整性**：若注释漏掉了某个基因，
    工具会把它当成基因间区。这不是可以靠算法弥补的 —— 报告里会明确标注
    这一点，并建议用户对照已知操纵子/必需基因清单人工复核。
    """
    if not features:
        if require_annotations:
            raise ValueError(
                "无法查找中性位点：**基因组文件没有注释**。\n"
                "「中性位点」的定义就是「不在任何基因内部」，这需要注释才能判断。\n"
                "纯 FASTA 只有序列，无法区分基因与基因间区。\n"
                "请改用 GenBank 格式（.gb/.gbk/.gbff）。")
        raise ValueError("需要注释特征表")

    n = len(seq)
    genome_gc = _gc(seq)

    # ---- 1) 收集所有"被占据"的区间（含 rRNA/tRNA，它们也是功能区）----
    occupied: list[tuple[int, int]] = []
    for f in features:
        if getattr(f, "key", "") in ("source",):
            continue
        # 多段特征（join）全部算占用
        parts = getattr(f, "parts", None) or [(f.start, f.end)]
        for a, b in parts:
            occupied.append((max(0, a), min(n, b)))
    if avoid:
        occupied.append((max(0, avoid[0]), min(n, avoid[1])))
    occupied.sort()

    # ---- 2) 扫出基因间区 ----
    gaps: list[tuple[int, int]] = []
    cur = 0
    for a, b in occupied:
        if a > cur:
            gaps.append((cur, a))
        cur = max(cur, b)
    if cur < n:
        gaps.append((cur, n))

    # ---- 3) 为每个基因间区找左右邻居基因名 + 评分 ----
    feats_by_pos = sorted(
        [(getattr(f, "start", 0), getattr(f, "end", 0), f) for f in features
         if getattr(f, "key", "") not in ("source",)],
        key=lambda t: t[0],
    )

    def neighbor(pos: int, side: str) -> tuple[str, int]:
        """返回 (基因名, 距离)。"""
        best = ("", 10 ** 9)
        for a, b, f in feats_by_pos:
            if side == "left" and b <= pos:
                d = pos - b
            elif side == "right" and a >= pos:
                d = a - pos
            else:
                continue
            if d < best[1]:
                nm = getattr(f, "name", lambda: "")()
                if not nm:
                    nm = getattr(f, "qualifiers", {}).get("locus_tag", "")
                best = (nm or "?", d)
        return best

    def mobile_within(a: int, b: int, radius: int = 2000) -> bool:
        for f in features:
            if not _is_mobile(f):
                continue
            fs, fe = getattr(f, "start", 0), getattr(f, "end", 0)
            if fs < b + radius and fe > a - radius:
                return True
        return False

    out: list[NeutralSite] = []
    for a, b in gaps:
        gap_len = b - a
        if gap_len < min_gap:
            continue

        center = (a + b) // 2
        off = 0
        score = 0.0
        reasons: list[str] = []
        warns: list[str] = []

        # (1) 中心性
        max_off = max(1, gap_len // 2)
        off_frac = abs(off) / max_off
        center_score = W_CENTER * (1.0 - off_frac)
        score += center_score

        # (2) 区间长度（饱和）
        eff = min(gap_len, GAP_LENGTH_SATURATION)
        score += W_GAP_LENGTH * (eff / GAP_LENGTH_SATURATION)
        reasons.append(f"基因间区 {gap_len} bp")

        # (3) GC 与基因组均值接近
        g = _gc(seq[a:b])
        gc_delta = abs(g - genome_gc)
        score += W_GC_MATCH * max(0.0, 1.0 - gc_delta / 0.15)
        if gc_delta > 0.10:
            warns.append(f"GC {g:.1%} 与基因组均值 {genome_gc:.1%} 相差 {gc_delta:.1%}")

        # (4) 可移动元件
        near_mob = mobile_within(a, b)
        if near_mob:
            score += W_MOBILE_PENALTY
            warns.append("邻近可移动元件（转座酶/前噬菌体/重组酶）—— 易发生重排")

        # (5) 距两侧基因过近
        lg, ld = neighbor(a, "left")
        rg, rd = neighbor(b, "right")
        # 插入点距边界的安全余量
        inset = abs(center - a)
        if inset < SAFE_MARGIN:
            score += W_LEFT_PROXIMITY * (SAFE_MARGIN - inset)
            warns.append(f"距左侧基因仅 {inset} bp（建议 ≥{SAFE_MARGIN} bp）")
        if abs(b - center) < SAFE_MARGIN:
            score += W_RIGHT_PROXIMITY * (SAFE_MARGIN - abs(b - center))
            warns.append(f"距右侧基因仅 {abs(b - center)} bp（建议 ≥{SAFE_MARGIN} bp）")

        # (6) 邻近结构性 RNA → 提示
        for f in features:
            if getattr(f, "key", "") in _STRUCTURAL_KEYS:
                fs, fe = getattr(f, "start", 0), getattr(f, "end", 0)
                if fs < b + 300 and fe > a - 300:
                    warns.append(f"邻近 {f.key}（{getattr(f, 'name', lambda: '')()}）"
                                 "—— 转录活跃且结构保守")
                    score -= 15.0
                    break

        if inset >= SAFE_MARGIN and abs(b - center) >= SAFE_MARGIN:
            reasons.append("两侧均有安全余量")
        if not near_mob:
            reasons.append("不邻近可移动元件")

        out.append(NeutralSite(
            start=a, end=b, insert_pos=center, length=gap_len,
            score=score, reasons=reasons, warnings=warns,
            off_center=off, left_gene=lg, right_gene=rg, gc=g,
            near_mobile=near_mob,
        ))

    out.sort(key=lambda s: -s.score)
    return out[:max_candidates]


# ---------------------------------------------------------------------------
# 靶向 spacer 设计（要求插入位点落在中性位点区域内）
# ---------------------------------------------------------------------------

def sgrnas_for_neutral_site(seq: str, site: NeutralSite, *,
                            window: int = 60, top: int = 5) -> list[dict]:
    """
    为某个中性位点找 STIB sgRNA，使**预测插入位点落在该区间内**。

    做法：在插入点附近 `window` bp 内扫描 GTN PAM，
    对每个候选计算预测插入位点（PAM 起点 + 偏移），只保留落在
    [site.start + margin, site.end - margin] 内的。

    返回 [{spacer, pam, pam_start, insert_site, distance_to_center}]。
    """
    from .sequence import scan_sgrnas, PAM_LEN, INSERT_OFFSET_FROM_PAM3

    lo = max(0, site.start)
    hi = min(len(seq), site.end)
    margin = 20
    ok_lo, ok_hi = lo + margin, hi - margin

    cands = scan_sgrnas(seq, max(0, site.insert_pos - window),
                        min(len(seq), site.insert_pos + window), flank=0)
    out: list[dict] = []
    for c in cands:
        # 预测插入位点 = PAM 起点 + PAM 长度 + 名义偏移
        pred = c.pam_start + PAM_LEN + INSERT_OFFSET_FROM_PAM3
        if not (ok_lo <= pred <= ok_hi):
            continue
        out.append({
            "spacer": c.spacer, "pam": c.pam, "pam_start": c.pam_start,
            "strand": c.strand, "insert_site": pred,
            "distance_to_center": abs(pred - site.insert_pos),
            "distance_to_edge": min(pred - lo, hi - pred),
        })
    out.sort(key=lambda d: d["distance_to_center"])
    return out[:top]
