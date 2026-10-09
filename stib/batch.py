"""
STIB Planner — 批量设计模式。

场景（用户反馈第 7 条）：用户可能不只是编辑单个基因，而是想**在全基因组范围内
筛选更多靶向基因**（例如做全基因组筛选、或为某个表型找靶点）。

批量模式与单基因模式的差别:
  - 单基因：目标区间小，候选少，可以做**全基因组脱靶枚举**（贵但值得）
  - 批量  ：目标区间大（整个基因组），候选成千上万 → **先做廉价过滤与打分**，
            再对 Top-N 做脱靶枚举

分两步走的原因：全基因组脱靶枚举是 6-mer 块查表（含替换变体）+ 逐位校验，虽然单条几秒，
但数千条相乘就不可接受。所以批量模式**默认不做脱靶枚举**，只做序列层面筛选，
并在输出中明确标注「未做脱靶评估，进入名单后需单独评估」。
"""

from __future__ import annotations

from typing import Any

from .knowledge import KnowledgeBase
from .planner import (
    Candidate, _pam_tier, _hard_filter, offtarget_risk, RISK_THRESHOLD,
    W_PAM_TIER1, W_PAM_TIER2, W_GC_CENTER, W_SEQ_FEATURE,
)
from .seqfeatures import compute_features, score_of as seq_score_of
from .sequence import GenomeIndex, scan_sgrnas, dedupe_sgrnas
from .offtarget import SeedIndex, find_offtargets, GUARANTEED_MAX_MISMATCH


def scan_genome(genome: GenomeIndex,
                pane_size: int = 200_000,
                min_spacing: int = 0,
                max_mismatch: int = 4,
                do_offtarget: bool = False,
                offtarget_top: int = 50,
                progress: bool = False) -> list[Candidate]:
    """
    全基因组扫描候选 sgRNA。

    参数:
      pane_size        分块扫描的块大小（避免一次性构造巨量候选）
      min_spacing      同一位点附近只保留一条（0 = 不过滤）
      do_offtarget     是否为最终 Top-N 做脱靶枚举
      offtarget_top    对前多少条做脱靶枚举
      progress         是否打印进度

    返回按分数排序的 Candidate 列表（未做脱靶的条目 risk_score 为 0）。
    """
    from .planner import Candidate as C

    n = len(genome)
    pool: list[Candidate] = []
    seen_sites: list[int] = []

    for start in range(0, n, pane_size):
        end = min(n, start + pane_size)
        raw = dedupe_sgrnas(scan_sgrnas(genome.seq, start, end, flank=0))
        for sg in raw:
            ok, reasons = _hard_filter(sg)
            if not ok:
                continue
            # 位点间距去重
            if min_spacing > 0 and any(abs(sg.insert_site - s) < min_spacing
                                       for s in seen_sites):
                continue
            seen_sites.append(sg.insert_site)

            cand = C(sgrna=sg, passed=True,
                     pam_tier=_pam_tier(sg.pam),
                     placement="outside",      # 批量模式无目标基因，落点不适用
                     seq_features=compute_features(genome.seq, sg.insert_site))
            s = 0.0
            if cand.pam_tier == 1:
                s += W_PAM_TIER1
            elif cand.pam_tier == 2:
                s += W_PAM_TIER2
            s += W_GC_CENTER * (1.0 - abs(sg.gc - 0.5) / 0.5)
            if cand.seq_features is not None:
                s += W_SEQ_FEATURE * seq_score_of(cand.seq_features)
            cand.score = s
            pool.append(cand)
        if progress:
            print(f"  扫描 {end:,}/{n:,} … 累计 {len(pool)} 条")

    pool.sort(key=lambda c: -c.score)

    # 对 Top-N 做脱靶枚举
    if do_offtarget and pool:
        idx = SeedIndex(genome.seq)
        for cand in pool[:offtarget_top]:
            rep = find_offtargets(idx, cand.sgrna.spacer, max_mismatch=max_mismatch)
            cand.offtarget = rep
            cand.risk_score, cand.risk_breakdown = offtarget_risk(rep)
            cand.is_self_only = cand.risk_breakdown.get("n_other", 0) == 0
            if cand.risk_score > RISK_THRESHOLD:
                cand.passed = False
                cand.reject_reasons.append(
                    f"脱靶风险分 {cand.risk_score:.1f} 超过阈值 {RISK_THRESHOLD}"
                )
            cand.score -= 3.0 * cand.risk_score
        pool.sort(key=lambda c: -c.score)

    return pool


def batch_plan(genome: GenomeIndex, gene: str, start: int, end: int,
               species: str = "Bacteroides thetaiotaomicron",
               top_n: int = 24, cm_conc: str = "25 µg/mL") -> Any:
    """
    在一个**较宽的区间**内批量出候选，并包装成可视化层需要的 Plan 形状。

    用途：为某个基因簇 / 基因区域一次性给出多条备选 sgRNA，
    或为筛选实验准备候选池。
    """
    from .planner import (Plan, evaluate_all, _build_steps, _build_decisions,
                          Warning, recommend_atmosphere)
    from .offtarget import RISK_RULES
    from .sequence import INSERT_WINDOW, INSERT_WINDOW_ALT

    kb = KnowledgeBase()
    raw = dedupe_sgrnas(scan_sgrnas(genome.seq, start, end, flank=0))
    if not raw:
        raise ValueError(f"在 [{start}, {end}) 内未找到候选 sgRNA")
    evaluated, _idx = evaluate_all(raw, genome, start, end)
    usable = [c for c in evaluated if c.passed][:top_n]
    rejected = [c for c in evaluated if not c.passed]
    if not usable:
        raise ValueError("批量模式：全部候选未通过筛选，请扩大区间或放宽阈值")

    atm = recommend_atmosphere(species)
    best = usable[0].sgrna
    steps = _build_steps(kb, gene, species, cm_conc,
                         arm_label=atm["label"], arm_reason=atm["reason"])
    warnings = [
        Warning(kind=c.kind, severity=c.severity, key=c.key,
                message=" | ".join(f"{cl['source']}: {cl['statement']}"
                                   for cl in c.claims),
                discriminator=c.discriminator)
        for c in kb.conflicts
    ]

    target = {
        "gene": gene, "species": species,
        "span": [start, end], "span_len": end - start,
        "selected_spacer": best.spacer, "selected_pam": best.pam,
        "selected_pam_start": best.pam_start,
        "predicted_insert_site": best.insert_site,
        "insert_window": list(INSERT_WINDOW),
        "insert_window_alt": list(INSERT_WINDOW_ALT),
        "cm_conc": cm_conc, "atmosphere": atm,
        "construct": "psh43 系列（Cas12k–TnsC 融合）",
        "n_candidates_scanned": len(raw),
        "n_passed": len(usable), "n_rejected": len(rejected),
        "max_mismatch": 4, "offtarget_guaranteed": True,
    }
    meta = {
        "tool": "stib-planner", "version": "0.3.0-demo0",
        "risk_rules": [{"level": lv, "rule": r} for lv, r in RISK_RULES],
        "genome_len": len(genome), "mode": "batch",
        "no_efficiency_model": (
            "**不做效率预测**：已发表数据点太少，无法抽出具统计意义的评估标准。"
            "排序只反映「是否可用」与可解释偏好。"
        ),
        "disclaimer": (
            "批量模式：在较宽区间内一次给出多条候选。"
            "脱靶评估为 6-mer 块 + 每块 1 替换枚举（≤4 不漏检，已交叉验证），"
            "不等价于 Cas-OFFinder 穷举。"
            "所有条件值须以原始 SOP 为准并人工复核。"
        ),
    }
    return Plan(
        target=target, sgrnas=[c.sgrna for c in usable],
        candidates=usable, rejected=rejected,
        oligos=[], primers={}, build_steps=steps,
        decisions=_build_decisions(), warnings=warnings,
        risks=[f.value for f in kb.by_prefix("risk.")],
        fallbacks=[], meta=meta,
    )
