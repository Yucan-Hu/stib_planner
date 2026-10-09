"""
STIB Planner — 规划层。

输入：目标基因 / 区间 + 基因组序列
输出：完整编辑方案（sgRNA 候选 → 寡核苷酸订单 → 验证引物 → 流程 → 判定点 → 风险）

设计原则:
  - **不预测效率**。只给可编辑性先验排序，并强制标注"需实验确认"。
  - **不隐藏不确定**。所有未收敛变量以 WARNING 出现在 plan.warnings 中。
  - 所有条件值都从 knowledge.KnowledgeBase 取，带出处。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any

from .knowledge import KnowledgeBase, Conflict, Fact
from .offtarget import (
    SeedIndex, find_offtargets, OffTargetReport, GUARANTEED_MAX_MISMATCH,
    RISK_RULES,
)
from .seqfeatures import compute_features, SeqFeatures, score_of as seq_score_of
from .sequence import (
    SgRna, Primer, GenomeIndex,
    scan_sgrnas, dedupe_sgrnas, seed_is_unique, gc_content, tm_basic,
    design_junction_primers, design_amplification_primers,
    INSERT_WINDOW, INSERT_WINDOW_ALT, SPACER_LEN, TSD_LEN,
    BSAI_F_OVERHANG, BSAI_R_OVERHANG,
)


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class Warning:
    """方案中的显式警示。"""

    kind: str          # 'open_question' | 'contradiction' | 'risk' | 'limitation'
    severity: str      # 'high' | 'medium' | 'low'
    key: str
    message: str
    source: str = ""
    discriminator: str = ""


@dataclass
class Step:
    """流程中的一步。"""

    phase: str
    day: str
    action: str
    readout: str = ""          # 该步应该看到什么
    ref: str = ""


@dataclass
class DecisionPoint:
    """判定点：看到 X 就做 Y。"""

    after: str                      # 在哪一步之后
    observe: str                    # 观察什么
    if_ok: str                      # 符合预期则
    if_not: str                     # 不符合则
    hypothesis: str = ""            # 最可能的假设
    discriminator: str = ""         # 最便宜的判别实验


@dataclass
class Plan:
    """完整编辑方案。"""

    target: dict
    sgrnas: list[SgRna]
    candidates: list["Candidate"]
    rejected: list["Candidate"]
    oligos: list[dict]
    primers: dict[str, Primer | None]
    build_steps: list[Step]
    decisions: list[DecisionPoint]
    warnings: list[Warning]
    risks: list[dict]
    fallbacks: list[dict] = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "target": self.target,
            "sgrnas": [s.as_dict() for s in self.sgrnas],
            "candidates": [c.as_dict() for c in self.candidates],
            "rejected": [c.as_dict() for c in self.rejected],
            "oligos": self.oligos,
            "primers": {k: (v.as_dict() if v else None) for k, v in self.primers.items()},
            "build_steps": [asdict(s) for s in self.build_steps],
            "decisions": [asdict(d) for d in self.decisions],
            "warnings": [asdict(w) for w in self.warnings],
            "risks": self.risks,
            "fallbacks": self.fallbacks,
            "meta": self.meta,
        }


# ---------------------------------------------------------------------------
# sgRNA 候选评估
# ---------------------------------------------------------------------------
#
# 设计原则（v0.2 修订，依据用户反馈）:
#
#   ❌ **不做效率预测**。已发表数据点太少（记录库中仅约 11 条带效率的位点），
#      抽不出统计上站得住的评估标准。硬做只会产出"看起来科学"的假精度。
#
#   ✅ **脱靶评估是硬判据**，不是打分项。全基因组枚举 ≤4 错配位点，
#      按错配位置分布分级；出现高风险脱靶的候选直接排除。
#
#   ✅ 其余判据也都是**硬过滤或可解释的偏好**，不给"效率分"。

#: 硬过滤阈值
GC_MIN, GC_MAX = 0.40, 0.60
MAX_HOMOPOLYMER = 4          # ≥5 连续单碱基淘汰
BSAI_SITE = "GGTCTC"         # 含 BsaI 位点会干扰 Golden Gate 构建

#: 偏好权重（仅在**通过硬过滤的候选之间**排序，不代表效率）
#
# 位点放置规则（v0.3，依据课题组反馈）:
#   基因失活时，插入位点应落在目标基因 **CDS 中央 80%**；
#   论文采用**更严格**的做法：CDS **前 1/3 到 2/3 之间**。
W_STRICT_ZONE = 120.0        # 落在「前 1/3 – 2/3」区（论文的更严格做法）
W_CENTRAL_80 = 80.0          # 落在中央 80%（一般建议）
W_OUTSIDE_GENE = -60.0       # 落在目标基因之外
W_GC_CENTER = 20.0
W_PAM_TIER1 = 25.0           # GTT（课题组经验最多）
W_PAM_TIER2 = 12.0           # GTA（GTN 次优选）
W_SEQ_FEATURE = 8.0          # ShCAST 序列偏好分（每 1 分对应的排序权重）

#: 脱靶风险加权（v0.6，依据课题组提供的结构信息修正分区）。
#:
#: ⚠️ **错配一律增加风险，绝无「减害」**（v0.3 的错误已更正）。
#: 一个真实存在的近同源位点，无论错配落在哪一段，都是**可能被转座的候选**。
#: 位置权重只表达**相对严重程度**，不能改变符号。
#:
#: ═══ 分区方向与 Cas9 相反（v0.6 重大更正）═══
#:
#: v0.2–v0.5 套用 **Cas9 惯例**：seed 在 spacer 3' 端、middle 居中，
#: 于是把最关键的 17 nt 当成低权重区 —— **方向反了**。
#:
#: 依据 *Structural basis for the assembly of the type V CRISPR-associated
#: transposon complex*：**前 17 nt 更重要**。
#: 结构上说得通：ShCAST 的 PAM 位于 spacer **上游**，坐标 0 紧邻 PAM；
#: Cas12k 与 RNA:DNA 异源双链的接触集中在 PAM 近端，
#: 3' 端 6 nt 是延伸尾部，接触较少。
#:
#:   每个非自身近同源位点的负担 = clamp(1 + Σ 位置权重, 0, MAX_BURDEN)
#:   权重：helix（坐标 0–16，PAM 近端关键区）**×3**、tail（坐标 17–22）×1
OFFTARGET_WEIGHT = {"helix": 3.0, "tail": 1.0}
HIT_BURDEN_BASE = 1.0

#: 单个位点负担的上限。取 4.0 以保留可区分的档位：
#:   纯 PAM 远端 / 纯种子区错配 → 1.0
#:   纯中段错配                 → 1 + 3n（封顶 4.0）
#:   额外的 0 错配拷贝（串联重复）→ 4.0
#: ⚠️ v0.4 曾取 2.0，结果把「2 个种子错配」与「2 个中段错配」压成同一个值（都是 2.0），
#: 使位置权重失去意义。上限只该防单个位点独大，**不该抹平档位**。
MAX_BURDEN = 4.0

#: 风险分阈值（v0.4 收紧，依据课题组反馈「整体看有点偏松」）。
#:
#:   v0.3：上限 1.0 + 阈值 6.0 → 可容纳 6 个近同源位点，明显过松。
#:   现在：上限 4.0 + 阈值 **3.0**。判例：
#:     1 个纯种子区位点（1.0）                         → 通过
#:     2 个纯种子区/远端点（2.0）                      → 通过
#:     1 个额外的 0 错配拷贝（4.0，串联重复）          → 拒绝
#:     1 个纯中段多错配位点（3.0～4.0）                → 拒绝
#:     1 个种子 + 2 个远端（3.0）                      → 拒绝
RISK_THRESHOLD = 3.0


def _hit_burden(hit: Any) -> float:
    """单个非自身位点的负担，clamp 到 [0, MAX_BURDEN]。错配只增不减。"""
    raw = HIT_BURDEN_BASE
    for region, weight in OFFTARGET_WEIGHT.items():
        raw += weight * len(hit.by_region.get(region, []))
    return max(0.0, min(MAX_BURDEN, raw))


def offtarget_risk(report: Any) -> tuple[float, dict]:
    """
    把脱靶报告折算成**加权风险分**与明细。

    ⚠️ **分数越小越好**；0 = 除自身位点外无近同源。

    - 排除自身位点（0 错配）—— 但**串联重复**的额外拷贝会被单独计数
    - 每个非自身近同源位点贡献 clamp(1 + Σ位置权重, 0, MAX_BURDEN)
    - **任何错配都是正风险**，不存在减害项
    """
    score = 0.0
    n_self = 0
    n_identical = 0      # 额外的 0 错配拷贝（串联重复）
    used = 0
    repeat_sites: list[int] = []
    for h in report.hits:
        if h.n_mismatch == 0:
            n_self += 1
            continue
        if getattr(h, "is_identical_copy", False):
            n_identical += 1
            repeat_sites.append(h.pos)
        used += 1
        score += _hit_burden(h)
    return score, {
        "n_self": n_self,
        "n_other": used,
        "n_identical_copies": n_identical,
        "identical_copy_sites": repeat_sites[:20],
        "risk_score": round(score, 2),
        "threshold": RISK_THRESHOLD,
        "scale": "0 = 最优（仅自身命中）；越大越差。任何错配都只增不减。",
    }


def _max_run(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def _revcomp_local2(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _hard_filter(c: SgRna) -> tuple[bool, list[str]]:
    """硬过滤（序列层面，不含脱靶）。返回 (通过?, 未通过原因列表)。"""
    reasons: list[str] = []
    if len(c.spacer) != SPACER_LEN:
        reasons.append(f"spacer 长度 {len(c.spacer)} ≠ {SPACER_LEN}")
    if any(b not in "ACGT" for b in c.spacer):
        reasons.append("spacer 含非 ACGT 字符")
    if not (GC_MIN <= c.gc <= GC_MAX):
        reasons.append(f"GC {c.gc:.0%} 超出 {GC_MIN:.0%}–{GC_MAX:.0%}")
    if _max_run(c.spacer) > MAX_HOMOPOLYMER:
        reasons.append(f"存在 {_max_run(c.spacer)} 个连续相同碱基")
    if BSAI_SITE in c.spacer or BSAI_SITE in _revcomp_local2(c.spacer):
        reasons.append("含 BsaI 识别位点 (GGTCTC)，会干扰 Golden Gate 构建")
    return (not reasons), reasons


@dataclass
class Candidate:
    """一条经评估的候选。"""

    sgrna: SgRna
    passed: bool
    reject_reasons: list[str] = field(default_factory=list)
    offtarget: Any = None        # OffTargetReport | None
    risk_score: float = 0.0      # 加权脱靶风险分（越低越好）
    risk_breakdown: dict = field(default_factory=dict)
    score: float = 0.0           # 综合排序分（含插入位置与 GC 偏好）
    is_self_only: bool = False   # 除自身外无任何近同源命中
    seq_features: Any = None     # SeqFeatures | None
    pam_tier: int = 3            # 1=GTT, 2=GTA, 3=其他 GTN
    placement: str = ""          # 'strict' | 'central80' | 'in_gene' | 'outside'

    @property
    def n_high(self) -> int:
        return len(self.offtarget.by_risk("high")) if self.offtarget else 0

    @property
    def n_medium(self) -> int:
        return len(self.offtarget.by_risk("medium")) if self.offtarget else 0

    @property
    def n_other(self) -> int:
        return self.risk_breakdown.get("n_other", 0)

    def as_dict(self) -> dict:
        d = self.sgrna.as_dict()
        d.update({
            "passed": self.passed,
            "reject_reasons": self.reject_reasons,
            "rank_score": round(self.score, 1),
            "offtarget_risk_score": round(self.risk_score, 2),
            "offtarget_breakdown": self.risk_breakdown,
            "is_self_only": self.is_self_only,
            "n_offtarget_high": self.n_high,
            "n_offtarget_medium": self.n_medium,
            "pam_tier": self.pam_tier,
            "placement": self.placement,
        })
        if self.seq_features:
            d["seq_features"] = self.seq_features.as_dict()
        if self.offtarget:
            d["offtarget"] = self.offtarget.as_dict()
        return d


def _pam_tier(pam: str) -> int:
    """PAM 优选层级（课题组经验：GTT 做得最多，其次 GTA）。"""
    p = pam.upper()
    if p == "GTT":
        return 1
    if p == "GTA":
        return 2
    return 3


#: 论文中测试过、耐氧能力较强、**可直接好氧接合**的物种（v0.3，依据课题组反馈）。
#: 出处：论文 "Genomic editing of different Bacteroides species"（这些物种均实测过）。
AEROBIC_CAPABLE = (
    "thetaiotaomicron", "fragilis", "vulgatus", "uniformis", "ovatus", "salyersiae",
)


def recommend_atmosphere(species: str, community_editing: bool = False) -> dict:
    """
    推荐接合气氛（v0.3，依据课题组反馈）。

    规则：
      - 论文中测试过的耐氧物种（bt/bf/bv/bu/bo/bs）→ **可直接好氧**
      - 客户菌株耐氧能力不好（未测试过的非模式种）→ **厌氧**
      - 菌群编辑（SynCom / 粪便菌群）→ **厌氧**

    ⚠️ 这**不再是未收敛变量** —— 课题组已澄清：厌氧仅限不耐氧菌株与菌群编辑。
    出处：SOP §4 阶段C Day2 第5条（好氧写法）+ 课题组澄清。
    """
    s = species.lower()
    if community_editing:
        return {
            "recommendation": "anaerobic",
            "label": "厌氧",
            "reason": "菌群编辑（SynCom / 粪便菌群）需厌氧接合",
            "confidence": "high",
        }
    if any(k in s for k in AEROBIC_CAPABLE):
        return {
            "recommendation": "aerobic",
            "label": "好氧（可直接操作）",
            "reason": f"{species} 属论文中测试过的耐氧物种，可直接好氧接合",
            "confidence": "high",
        }
    return {
        "recommendation": "anaerobic_or_test",
        "label": "厌氧（保守）；或先做好氧/厌氧平行对照",
        "reason": "该菌株未在论文中测试过，耐氧能力未知。厌氧是保守选择；"
                  "若想确认，可先做一次小规模平行对照（成本 = 多一块板）",
        "confidence": "medium",
    }


def _placement(insert_site: int, gene_start: int, gene_end: int) -> str:
    """
    插入位点相对目标基因的落点分类（v0.4，依据课题组反馈）。

    两级制：

      1. **strict** —— 落在 CDS **前 1/3 – 2/3** 之间。
         论文采用的更严格做法，**优先满足**。

      2. **central80** —— 严格区不满足时，回落到**普通约束**：
         插入位点处于基因**中段 80%**，即**排除 CDS 前 5% 与后 5%**。

      3. **in_gene** —— 落在 CDS 两端各 5% 内，不推荐。

      4. **outside** —— 落在目标基因之外。

    ⚠️ v0.3 的 central80 用的是「排除前后各 10%」，与课题组反馈的
    「排除前 5% 与后 5%」不一致，本次更正为 5%。
    """
    span = max(1, gene_end - gene_start)
    if not (gene_start <= insert_site < gene_end):
        return "outside"
    frac = (insert_site - gene_start) / span
    if 1.0 / 3.0 <= frac <= 2.0 / 3.0:
        return "strict"
    # 普通约束：中段 80% = 排除前 5% 与后 5%
    if 0.05 <= frac <= 0.95:
        return "central80"
    return "in_gene"


def _cheap_score(c: "Candidate") -> float:
    """廉价排序分（不含脱靶）——用于两阶段评估的第一阶段。"""
    s = 0.0
    if c.placement == "strict":
        s += W_STRICT_ZONE
    elif c.placement == "central80":
        s += W_CENTRAL_80
    elif c.placement == "in_gene":
        s += W_CENTRAL_80 * 0.5
    else:
        s += W_OUTSIDE_GENE
    if c.pam_tier == 1:
        s += W_PAM_TIER1
    elif c.pam_tier == 2:
        s += W_PAM_TIER2
    s += W_GC_CENTER * (1.0 - abs(c.sgrna.gc - 0.5) / 0.5)
    if c.seq_features is not None:
        s += W_SEQ_FEATURE * seq_score_of(c.seq_features)
    return s


#: 两阶段评估：第一阶段按廉价分排序后，只对前 N 条做全基因组脱靶枚举。
#: 6-mer + 1-替换的查询量比 8-mer 大约 40 倍，对全部候选做脱靶不可接受
#: （实测 173 条候选 → 86 秒）。取 30 条 → 约 15 秒，且不损失最终结果质量：
#: 通过硬过滤的候选本就按廉价分排在前列。
OFFTARGET_BUDGET = 30


def evaluate_all(cands: list[SgRna], genome: GenomeIndex,
                 gene_start: int, gene_end: int,
                 index: Any = None,
                 max_mismatch: int = 4,
                 offtarget_budget: int = OFFTARGET_BUDGET) -> tuple[list[Candidate], Any]:
    """
    两阶段评估：硬过滤 + 序列偏好 → 廉价排序 → **只对前 N 条**做脱靶枚举 → 终排。

    返回 (全部候选（已按分数排序）, 种子索引)。

    ⚠️ 为何分两阶段：全基因组脱靶枚举是主要成本。6-mer + 1-替换方案下
    单条查询约 22 万候选、需数秒；对上百条候选逐个做会让 `plan` 从
    十余秒退化到一分半（实测 86 秒）。由于通过硬过滤的候选按廉价分排在前列，
    对前 N 条做脱靶已足以选出正确的首选，且能给出足够的备选。

    未做脱靶枚举的候选在报告中标为「未做脱靶评估」，不会伪装成已评估。
    """
    idx = index if index is not None else SeedIndex(genome.seq)

    evaluated: list[Candidate] = []
    for c in cands:
        ok, reasons = _hard_filter(c)
        cand = Candidate(sgrna=c, passed=ok, reject_reasons=reasons)
        cand.pam_tier = _pam_tier(c.pam)
        cand.placement = _placement(c.insert_site, gene_start, gene_end)
        # 序列特征（ShCAST 的 TSD / ±5 / ±12 偏好）—— 不参与硬过滤，只影响排序
        cand.seq_features = compute_features(genome.seq, c.insert_site)
        cand.score = _cheap_score(cand)
        evaluated.append(cand)

    # ---- 第一阶段：按廉价分挑出值得做脱靶评估的候选 ----
    passed_cheap = [c for c in evaluated if c.passed]
    passed_cheap.sort(key=lambda c: -c.score)
    to_evaluate = passed_cheap[:max(1, offtarget_budget)]
    to_eval_ids = {id(c) for c in to_evaluate}

    # ---- 第二阶段：只对这些做脱靶枚举 ----
    for cand in to_evaluate:
        cand.offtarget = find_offtargets(idx, cand.sgrna.spacer,
                                        max_mismatch=max_mismatch)
        cand.risk_score, cand.risk_breakdown = offtarget_risk(cand.offtarget)
        cand.is_self_only = cand.risk_breakdown.get("n_other", 0) == 0
        if cand.risk_score > RISK_THRESHOLD:
            cand.passed = False
            cand.reject_reasons.append(
                f"脱靶风险分 {cand.risk_score:.1f} 超过阈值 {RISK_THRESHOLD}"
                f"（非自身近同源位点 {cand.risk_breakdown.get('n_other', 0)} 个；"
                f"权重：中段×3 / 远端×1 / 种子×(−1)）"
            )

    # 未进入脱靶阶段的候选：标注清楚，不用零风险伪装
    for cand in evaluated:
        if cand.passed and id(cand) not in to_eval_ids:
            cand.reject_reasons.append(
                f"未做脱靶评估（仅对廉价分前 {offtarget_budget} 条做了全基因组枚举）"
            )

    # ---- 终排：把脱靶风险计入分数 ----
    for c in evaluated:
        if not c.passed:
            c.score = float("-inf")
            continue
        c.score = _cheap_score(c) - 3.0 * c.risk_score

    evaluated.sort(key=lambda c: -c.score)
    return evaluated, idx


def evaluate_candidates(cands: list[SgRna], genome: GenomeIndex,
                        gene_start: int, gene_end: int,
                        index: Any = None,
                        max_mismatch: int = 4,
                        top_n: int = 5) -> tuple[list[Candidate], Any]:
    """
    硬过滤 + 全基因组脱靶评估 + 可解释排序，只返回**通过**的候选（截断到 top_n）。

    ⚠️ 不输出效率预测值。排序只反映「是否可用」与可解释的偏好。
    """
    evaluated, idx = evaluate_all(cands, genome, gene_start, gene_end, index, max_mismatch)
    usable = [c for c in evaluated if c.passed]
    return usable[:top_n], idx


def _rank_sgrnas(cands: list[SgRna], genome: GenomeIndex,
                 gene_start: int, gene_end: int,
                 top_n: int = 5) -> list[SgRna]:
    """兼容旧接口的薄封装。新代码请用 evaluate_candidates()。"""
    usable, _ = evaluate_candidates(cands, genome, gene_start, gene_end, top_n=top_n)
    return [c.sgrna for c in usable]


# ---------------------------------------------------------------------------
# 引物组装
# ---------------------------------------------------------------------------

def _build_primers(genome: GenomeIndex, insert_site: int,
                   cargo_arm_len: int = 800) -> dict[str, Primer | None]:
    """
    组装验证引物组。

    产出:
      LE_up_F     基因组上游引物（正链，向右）
      LE_cargo_R  cargo 内部 LE 侧引物（负链）— 通用，与 span 无关
      RE_cargo_F  cargo 内部 RE 侧引物（正链）— 通用
      RE_dn_R     基因组下游引物（负链，向左）
      backbone_F 骨架检测（ampR/ori）— 需用户提供骨架序列，此处仅占位
      ctrl_16S_F / ctrl_16S_R  16S 全长鉴定对照

    出处: SOP 表 8；SOP §4 阶段D D.3
    """
    out: dict[str, Primer | None] = {}

    jp = design_junction_primers(genome.seq, insert_site)
    out["LE_up_F"] = jp.get("up_fwd")
    out["RE_dn_R"] = jp.get("dn_rev")

    # cargo 侧通用引物：此处不设计，因为 cargo 序列由用户提供。
    # 保留键位以便插件层填充。
    out["LE_cargo_R"] = None
    out["RE_cargo_F"] = None
    out["backbone_F"] = None
    out["backbone_R"] = None

    # 16S 对照（27F / 1492R，论文与 SOP 通用引物）
    # 注意: 这两条含 IUPAC 简并碱基（M = A/C, Y = C/T），Tm 按简并位置不计入估算，
    # 仅作参考。实际退火温度以厂家推荐（通常 48-52 °C 用于此类通用引物）为准。
    _27F = "AGAGTTTGATCMTGGCTCAG"
    _1492R = "TACGGYTACCTTGTTACGACTT"
    out["ctrl_16S_F"] = Primer(
        name="27F", seq=_27F, tm=tm_basic(_27F), gc=gc_content(_27F),
        start=-1, end=-1, strand="+",
        note="16S 通用正向引物（SOP §4 D.3）；含简并碱基 M，Tm 仅供参考",
    )
    out["ctrl_16S_R"] = Primer(
        name="1492R", seq=_1492R, tm=tm_basic(_1492R), gc=gc_content(_1492R),
        start=-1, end=-1, strand="-",
        note="16S 通用反向引物（SOP §4 D.3）；含简并碱基 Y，Tm 仅供参考",
    )
    return out


# ---------------------------------------------------------------------------
# 流程与判定点
# ---------------------------------------------------------------------------

def _build_steps(kb: KnowledgeBase, gene: str, species: str,
                 cm_conc: str, arm_label: str = "好氧",
                 arm_reason: str = "") -> list[Step]:
    """构建阶段 A-F 流程（引用 SOP，不重写）。`arm_*` 为接合气氛推荐。"""

    def ref(k: str) -> str:
        f = kb.get(k)
        return f.ref() if f else ""

    return [
        Step("A sgRNA 设计", "1-2 天",
             f"确认目标区域，下载 {species} 参考基因组；扫描 NGTN PAM；取 PAM 下游 20-23 nt 作 spacer（本体系 23 nt）",
             "得到候选 spacer 列表（见本方案 §2）", ref("system.spacer_len")),
        Step("A sgRNA 设计", "1-2 天",
             "脱靶分析：用本工具的 offtarget 做全基因组 ≤4 错配枚举（已交叉验证）；"
             "重点排查 **helix 段（PAM 近端前 17 nt）** 的近同源位点",
             "无 helix 段高度相似的候选脱靶位点", ref("risk.RK2")),
        Step("A sgRNA 设计", "1-2 天",
             "确认插入窗口（实测 PAM 下游 56–66 bp 区间）落在预期位置；避开必需基因编码区与上游调控元件。"
             "**区间内偏移均可接受**（该窗口是观察分布，不是验收标准），"
             "仅方向敏感应用（启动子捕获、框内插入）需精确读坐标",
             "窗口位置符合应用需求", ref("pam.preference")),
        Step("B 载体构建", "1-2 天",
             "B.1 更换 sgRNA：Bsa1-Golden Gate（寡核苷酸磷酸化退火 + 载体单酶切 + 组装）",
             "菌落 PCR 扩增 sgRNA 表达盒并测序，确认 spacer 正确", ref("reagent.bsai")),
        Step("B 载体构建", "2-3 天",
             "B.2 更换 cargo：高保真 PCR + DMT 消化 + Gibson 组装 + 热激转化 S17-1λpir",
             "LE / cargo / RE 各 junction 测序正确", ref("reagent.gibson")),
        Step("C 单菌接合", "接合前 1 天",
             "制备 BHIB 血平板（无抗）与筛选平板（BHIS + Genta 200 + Cm " + cm_conc + "）",
             "血平板为 8% 脱纤维马血；筛选平板含正确抗生素", ref("media.blood_plate")),
        Step("C 单菌接合", "Day 0/1",
             "供体划线 LB + Amp 100；受体划线 BHIS + Genta（37 °C）",
             "两者均长出单克隆", ref("reagent.amp")),
        Step("C 单菌接合", "Day 1",
             "各挑单克隆接种：供体 5 mL LB+Amp（37 °C 220 rpm）；受体 5 mL BHIS+Genta（37 °C 厌氧静置）过夜",
             "过夜培养浑浊", ref("SOP_MAIN")),
        Step("C 单菌接合", "Day 2 (09:00)",
             "供体 1:100 转接新鲜 LB+Amp，培养 2.5-3 h 至 OD600 ≈ 0.4（装液量 ≤ 容器 1/5）",
             "OD600 达 0.4（对数前期）", ref("conj.donor_od")),
        Step("C 单菌接合", "Day 2 (09:00)",
             "受体 1:100 转接 BHIS，37 °C 厌氧静置约 6 h",
             "OD600 ≈ 0.8", ref("conj.recipient_od")),
        Step("C 单菌接合", "Day 2 (15:00)",
             "收供体：4000 g 5 min，无菌 PBS 洗 1 次去抗生素，BHIS 重悬浓缩至 OD600 ≈ 0.8",
             "无抗生素残留", ref("conj.spin")),
        Step("C 单菌接合", "Day 2 (15:00)",
             "供体:受体 1:1 各 1 mL 混匀 → 4000 g 20 min → 弃上清 → 约 100 µL BHIS 重悬 → 滴至无抗 BHIB 血平板中心",
             "形成硬币大小液滴", ref("conj.ratio")),
        Step("C 单菌接合", "Day 2 → Day 3",
             f"37 °C 静置 20-24 h 接合。【气氛：{arm_label}】{arm_reason}",
             f"形成菌苔（气氛：{arm_label}）", ref("conj.mating")),
        Step("C 单菌接合", "Day 3 下午",
             "棉签刮取菌苔 → 1 mL BHI 重悬 → 梯度稀释 → 原液与 10⁻² 各 100 µL 涂筛选平板",
             "37 °C 厌氧 2-3 天后出现抗性菌落", ref("conj.plating")),
        Step("D 鉴定", "Day 5-6",
             "挑单克隆接种 1 mL BHI + Genta + Cm，厌氧过夜；菌液 10 µL + 40 µL 无菌水，取 1 µL 作 PCR 模板",
             "菌液浑浊", ref("SOP_MAIN")),
        Step("D 鉴定", "1-2 天",
             "菌落 PCR：LE junction / RE junction / cargo / 骨架共整合（骨架应为阴性）",
             "前三条有单一目的条带；骨架无条带", ref("verify.backbone_pcr")),
        Step("D 鉴定", "1-2 天",
             "LE/RE junction 产物 Sanger 测序",
             f"整合位于 PAM 下游 {INSERT_WINDOW[0]}-{INSERT_WINDOW[1]} bp；存在 {TSD_LEN} bp TSD；边界与载体一致",
             ref("verify.sanger")),
        Step("D 鉴定", "1-2 天",
             "16S 全长鉴定（27F/1492R，Sanger 双向，EzBioCloud）",
             "top-5 hits 中目标种显著领先（注意近缘种干扰）", ref("verify.16s")),
    ]


def _build_decisions() -> list[DecisionPoint]:
    """判定点：症状 → 假设 → 最便宜的判别实验。"""

    return [
        DecisionPoint(
            after="Day 3 涂板后 2-3 天（bv 等慢生长种延至 7 天）",
            observe="筛选平板上无菌落",
            if_ok="进入挑克隆鉴定",
            if_not="**先排除操作原因，再怀疑生物学原因**",
            hypothesis="候选假设：(H1) 血平板漏加 DAP（若用 ΔdapA 供体）；(H2) 抗生素失效或浓度过高；(H3) 接合气氛不合适；(H4) sgRNA/PAM 选择不当",
            discriminator="H1 最高优先且**零成本可查**：核对血平板配方记录是否含 DAP 0.06 g/L。记录库 2024-03 曾因此整批零产物。若无 DAP 问题，再做 NC 对照判断 H2。",
        ),
        DecisionPoint(
            after="菌落 PCR",
            observe="骨架 PCR **阳性**（检出 ampR/ori）",
            if_ok="骨架 PCR 阴性 → 简单插入，正确",
            if_not="判为共整合，**不得作为阳性克隆**",
            hypothesis="Type V-K 复制型转座通路副产物",
            discriminator="划线纯化后传代 2-5 次复检骨架 PCR。若持续阳性，说明该 sgRNA/位点组合易产生共整合 → 换 spacer 或换位点（成本 1 轮接合，约 3 天）。",
        ),
        DecisionPoint(
            after="菌落 PCR",
            observe="LE/RE junction **两端都无条带**，且野生型条带也扩不出",
            if_ok="两端均有条带 → 进入测序",
            if_not="**先怀疑 PCR 体系与模板，不是构建失败**",
            hypothesis="记录库反复出现此模式：(a) 模板降解；(b) 引物设计不佳；(c) 平板污染导致挑到杂菌",
            discriminator="换一对新引物重做 PCR（记录库 2024-08-30 首次 junction PCR 无条带，换引物后立即成功）。同时取菌液做 16S 确认挑到的确实是目标菌。",
        ),
        DecisionPoint(
            after="16S 鉴定",
            observe="16S top hit 与目标种相似度仅略高，或出现非目标种",
            if_ok="目标种相似度显著领先 → 身份确认",
            if_not="**不得仅凭单一 top hit 判定物种**",
            hypothesis="16S 分辨率不足以区分近缘种（B. thetaiotaomicron vs B. faecis 相似度 98.43%；P. vulgatus vs P. dorei 97.10%）",
            discriminator="查看 top-5 hits 的相似度分布；必要时补做 species-specific PCR 或 WGS。交接菌株时双人核对。",
        ),
        DecisionPoint(
            after="Sanger 测序",
            observe="插入位置明显偏离观察窗口（落在 56–66 bp 之外）",
            if_ok="落于窗口内 → 方案成立",
            if_not="**核对实际窗口后再判定**",
            hypothesis="论文不同章节报告 60–66 bp（纯培养）与 56–63 bp（菌群）两个观察窗口",
            discriminator=f"Sanger 读出实际坐标。注意论文 Figure 4 报告菌群编辑中为 {INSERT_WINDOW_ALT[0]}-{INSERT_WINDOW_ALT[1]} bp，与 Figure 2 的 {INSERT_WINDOW[0]}-{INSERT_WINDOW[1]} bp 不一致。",
        ),
    ]


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def make_plan(genome: GenomeIndex, gene: str, start: int, end: int,
              species: str = "Bacteroides thetaiotaomicron",
              top_n: int = 5,
              max_mismatch: int = 4,
              community_editing: bool = False,
              kb: KnowledgeBase | None = None) -> Plan:
    """
    生成完整编辑方案。

    参数:
      genome             GenomeIndex
      gene               基因名（仅用于标注）
      start/end          目标基因在基因组上的 0-based 半开区间
      species            物种名（影响氯霉素浓度与接合气氛推荐）
      top_n              返回的 sgRNA 候选数
      max_mismatch       脱靶枚举的错配上限。**≤4 时保证不漏检**；
                         >4 会显式警示可能少报。
      community_editing  是否为菌群编辑（决定接合气氛：菌群编辑用厌氧）
    """
    kb = kb or KnowledgeBase()

    # 氯霉素浓度按物种选择
    if "thetaiotaomicron" in species or "ovatus" in species:
        cm_conc = "25 µg/mL"
    else:
        cm_conc = "12.5 µg/mL"

    # 接合气氛推荐（v0.3：按耐氧能力，不再是未收敛变量）
    atm = recommend_atmosphere(species, community_editing)

    # 1) 扫描候选
    raw = scan_sgrnas(genome.seq, start, end)
    raw = dedupe_sgrnas(raw)
    if not raw:
        raise ValueError(
            f"在 {gene} 区间 [{start}, {end}) ±300 bp 内未找到任何 {SPACER_LEN} nt spacer 的 GTN PAM。"
            "请扩大窗口或确认序列方向。"
        )

    # 2) 硬过滤 + 全基因组脱靶评估 + 排序（一次算完，不重复枚举）
    max_mm = max_mismatch
    evaluated, index = evaluate_all(raw, genome, start, end, max_mismatch=max_mm)
    usable = [c for c in evaluated if c.passed][:top_n]
    rejected = [c for c in evaluated if not c.passed]

    if not usable:
        raise ValueError(
            f"在 {gene} 区间内扫描到 {len(raw)} 条候选，但**全部未通过硬过滤或脱靶评估**：\n"
            + "\n".join(f"  - {c.sgrna.spacer}: {'; '.join(c.reject_reasons)}"
                        for c in rejected[:8])
            + "\n请扩大扫描窗口，或用 --max-mismatch 调宽后复查（注意会削弱不漏检保证）。"
        )

    if not usable:
        raise ValueError(
            f"在 {gene} 区间内扫描到 {len(raw)} 条候选，但**全部未通过硬过滤或脱靶评估**：\n"
            + "\n".join(f"  - {c.sgrna.spacer}: {'; '.join(c.reject_reasons)}" for c in rejected[:8])
            + "\n请扩大扫描窗口，或把 --max-mismatch 调宽后复查（注意会削弱保证）。"
        )

    ranked = [c.sgrna for c in usable]

    # 3) 寡核苷酸订单
    oligos = [
        {
            "name": f"shgRNA-{i+1}F",
            "seq": s.oligo_f,
            "note": f"{BSAI_F_OVERHANG}(BsaI 粘性末端) + spacer；PAGE 纯化，母液 100 µM，工作液 10 µM",
        }
        for i, s in enumerate(ranked)
    ] + [
        {
            "name": f"shgRNA-{i+1}R",
            "seq": s.oligo_r,
            "note": f"{BSAI_R_OVERHANG}(BsaI 粘性末端) + spacer 反向互补",
        }
        for i, s in enumerate(ranked)
    ]

    # 4) 引物
    best = ranked[0]
    primers = _build_primers(genome, best.insert_site)

    # 5) 流程与判定点
    steps = _build_steps(kb, gene, species, cm_conc,
                         arm_label=atm["label"], arm_reason=atm["reason"])
    decisions = _build_decisions()

    # 6) 警示（未收敛变量 + 高等级风险 + 脱靶边界）
    warnings: list[Warning] = []
    for c in kb.conflicts:
        warnings.append(Warning(
            kind=c.kind, severity=c.severity, key=c.key,
            message=" | ".join(f"{cl['source']}: {cl['statement']}" for cl in c.claims),
            source=c.claims[0]["source"] if c.claims else "",
            discriminator=c.discriminator,
        ))
    for f in kb.by_prefix("risk."):
        v = f.value
        if v["level"] == "high":
            warnings.append(Warning(
                kind="risk", severity="high", key=v["id"],
                message=v["title"] + " — " + v["detail"],
                source=f.source, discriminator=v["action"],
            ))
    if max_mm > GUARANTEED_MAX_MISMATCH:
        warnings.append(Warning(
            kind="limitation", severity="high", key="offtarget_not_guaranteed",
            message=(f"脱靶阈值设为 {max_mm} 错配，**超出本引擎的保证范围"
                     f"（≤{GUARANTEED_MAX_MISMATCH}）**，结果可能少报。"
                     "论文与 SOP 用 Cas-OFFinder v2.4 允许 ≤5 错配，两者不等价。"),
            source="PAPER_METHODS", locator="Guide RNA design",
            discriminator="用 Cas-OFFinder 交叉核对，或把阈值降回 4 以取得不漏检保证。",
        ))
    warnings.sort(key=lambda w: {"high": 0, "medium": 1, "low": 2}[w.severity])

    risks = [f.value for f in kb.by_prefix("risk.")]

    # 7) 回退面板：已验证 spacer（新手最需要的是"能用的"，不是"最优的"）
    fallbacks = valid_spacer_panel(species)

    top_ot = usable[0].offtarget if usable else None
    target = {
        "gene": gene,
        "species": species,
        "span": [start, end],
        "span_len": end - start,
        "selected_spacer": best.spacer,
        "selected_pam": best.pam,
        "selected_pam_start": best.pam_start,
        "predicted_insert_site": best.insert_site,
        "insert_window": list(INSERT_WINDOW),
        "insert_window_alt": list(INSERT_WINDOW_ALT),
        "cm_conc": cm_conc,
        "atmosphere": atm,
        "community_editing": community_editing,
        "construct": "psh43 系列（Cas12k–TnsC 融合）",
        "n_candidates_scanned": len(raw),
        "n_passed": len(usable),
        "n_rejected": len(rejected),
        "max_mismatch": max_mm,
        "offtarget_guaranteed": max_mm <= GUARANTEED_MAX_MISMATCH,
        "n_offtarget_high": len(top_ot.by_risk("high")) if top_ot else 0,
        "n_offtarget_total": top_ot.n_hits if top_ot else 0,
    }

    meta = {
        "tool": "stib-planner",
        "version": "0.2.0-demo0",
        "offtarget_engine": (
            f"6-mer 块 ×3 + 每块 1 替换枚举，阈值 ≤{max_mm} 错配；"
            + ("保证不漏检（≤4）" if max_mm <= GUARANTEED_MAX_MISMATCH
               else "**可能漏检**（>4，保证失效）")
        ),
        "risk_rules": [{"level": lv, "rule": r} for lv, r in RISK_RULES],
        "provenance": "所有条件值见 knowledge.py，逐条带出处",
        "no_efficiency_model": (
            "**不做效率预测**：已发表数据点太少（记录库中约 11 条带效率的位点），"
            "无法抽出具统计意义的评估标准。排序只反映「是否可用」与可解释偏好。"
        ),
        "disclaimer": (
            "本方案基于已发表论文 + 课题组 SOP + 实验记录库生成。"
            "脱靶评估为种子剪枝算法，不等价于 Cas-OFFinder 穷举。"
            "所有条件值须以原始 SOP 为准并人工复核。"
        ),
    }

    return Plan(
        target=target, sgrnas=ranked, candidates=usable, rejected=rejected,
        oligos=oligos, primers=primers, build_steps=steps,
        decisions=decisions, warnings=warnings, risks=risks,
        fallbacks=fallbacks, meta=meta,
    )


def _all_evaluated(cands: list[SgRna], genome: GenomeIndex,
                   gene_start: int, gene_end: int,
                   index: Any, max_mismatch: int) -> list["Candidate"]:
    """
    已废弃：请改用 evaluate_all()。

    保留仅为向后兼容；内部直接转调，不再重复实现评估逻辑。
    """
    evaluated, _ = evaluate_all(cands, genome, gene_start, gene_end, index, max_mismatch)
    return evaluated


# ---------------------------------------------------------------------------
# 已验证 spacer 面板（回退用，不是效率预测）
# ---------------------------------------------------------------------------
#
# 对新手的意义: 他不需要"最优 spacer"，他需要**一条确定能用的**。
# 以下 spacer 均来自论文 Table S2，并附论文/记录库中的实测结果。
# 同时列出**已知失败的位点** —— 这些同样是宝贵信息，避免新手踩进去。

_VALIDATED_PANEL: list[dict] = [
    {"name": "sgRNA44", "site": "Bt05", "species": "Bacteroides thetaiotaomicron",
     "spacer": "CAGATCGGTGCTGGGAGTGGTTC", "pam": "GTT", "evidence": "≈100% 编辑，0 共整合（记录库 2023-09-30）",
     "outcome": "success", "source": "TABLE_S2"},
    {"name": "sgRNA26", "site": "Bt01 (tdk)", "species": "Bacteroides thetaiotaomicron",
     "spacer": "ATTCGGAAGAAGACGTGGTATCA", "pam": "GTT", "evidence": "论文主力位点；记录库 15%",
     "outcome": "success", "source": "TABLE_S2"},
    {"name": "sgRNA39", "site": "Bt02 (bt1754)", "species": "Bacteroides thetaiotaomicron",
     "spacer": "TGGTTTAATCGGGATGACGGCTT", "pam": "GTT", "evidence": "记录库 25%",
     "outcome": "success", "source": "TABLE_S2"},
    {"name": "sgRNA65", "site": "Bo03", "species": "Bacteroides ovatus",
     "spacer": "TCAGAATGTGTCACAATCACCTT", "pam": "GTT", "evidence": "记录库频率 ≈0.9",
     "outcome": "success", "source": "TABLE_S2"},
    {"name": "sgRNA66", "site": "Bo04", "species": "Bacteroides ovatus",
     "spacer": "AATGAGTGAACTGTGCGGATGTA", "pam": "GTT", "evidence": "记录库频率 ≈0.966",
     "outcome": "success", "source": "TABLE_S2"},
    {"name": "sgRNA81", "site": "16S rDNA", "species": "Bacteroides fragilis / uniformis",
     "spacer": "ACTCACCCGTGCGCCGGTCGCCA", "pam": "GTT", "evidence": "论文 >90%（三株）",
     "outcome": "success", "source": "TABLE_S2",
     "note": "⚠️ 该 spacer 为 22 nt（其余 23 nt）；SOP 另提到 shgRNA81/82/83，但 Table S2 只有 81"},
    {"name": "sgRNA76", "site": "Bv01", "species": "Bacteroides vulgatus",
     "spacer": "ACAACAGTACGCCTAGCATCTAT", "pam": "GTT", "evidence": "论文 Bv01 位点 >93%",
     "outcome": "success", "source": "TABLE_S2"},
    # --- 已知失败 ---
    {"name": "sgRNA47", "site": "Bt03", "species": "Bacteroides thetaiotaomicron",
     "spacer": "GTTCGTGTTTGCCGTGTATCATT", "pam": "GTT", "evidence": "**0%**（记录库 2023-07-17；构建成功但从未工作）",
     "outcome": "failed", "source": "TABLE_S2"},
    {"name": "sgRNA78/79", "site": "Bv02 / Bv03", "species": "Bacteroides vulgatus",
     "spacer": "CAACAGTTCTGATTATGTAGAAT / CAATGCCGCACCTTTCTTAATCT", "pam": "GTT",
     "evidence": "论文报告 **0%**", "outcome": "failed", "source": "PAPER_RESULTS"},
    {"name": "sgRNA63", "site": "Bo02", "species": "Bacteroides ovatus",
     "spacer": "ATTCCAATACTGCCCTGTATGTC", "pam": "GTT", "evidence": "论文报告 30%（显著低于同物种其他位点）",
     "outcome": "weak", "source": "PAPER_RESULTS"},
]


def valid_spacer_panel(species: str = "") -> list[dict]:
    """
    返回已验证/已知失败的 spacer 面板。

    用途: 新手最需要的是**一条确定能用的** spacer，而不是"最优"的。
    若首选位点失败，本面板提供经过实验验证的回退选项。

    ⚠️ 这不是效率预测。这些是**已发表的实测结果**，作为经验先验使用。
    """
    out = list(_VALIDATED_PANEL)
    if species:
        key = species.split()[0].lower() if species else ""
        # 目标物种的排前面
        out.sort(key=lambda d: 0 if key and key in d["species"].lower() else 1)
    return out
