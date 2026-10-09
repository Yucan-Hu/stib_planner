#!/usr/bin/env python3
"""
STIB Planner CLI。

用法:
  python -m stib.cli selftest
  python -m stib.cli check
  python -m stib.cli facts [--prefix risk.]
  python -m stib.cli fetch NC_004663.1 [--out data/]
  python -m stib.cli scan --genome data/NC_004663.1.fasta --motif <SEQ>
  python -m stib.cli plan --genome data/NC_004663.1.fasta \
      --gene tdk --motif ATTCGGAAGAAGACGTGGTATCA \
      --species "Bacteroides thetaiotaomicron" \
      --out plan_tdk.md [--json plan_tdk.json]

设计说明: 本工具**不预测效率**、**不做全基因组脱靶搜索**，见 README §7。
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.parse
import urllib.request

from .knowledge import KnowledgeBase, SOURCES, validate
from .offtarget import SeedIndex, find_offtargets, GUARANTEED_MAX_MISMATCH
from .planner import make_plan, valid_spacer_panel, offtarget_risk, RISK_THRESHOLD
from .render import render_markdown, render_json
from .visual import render_html
from .protocol import render_protocol
from .genome import (load_genome_any, resolve_region, RegionError,
                     parse_genbank, detect_format)
from .batch import batch_plan
from .sequence import (
    GenomeIndex, read_fasta, write_fasta, scan_sgrnas, dedupe_sgrnas,
    selftest as seq_selftest, revcomp, SPACER_LEN, INSERT_WINDOW,
)

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


# ---------------------------------------------------------------------------
# 子命令
# ---------------------------------------------------------------------------

def cmd_selftest(args: argparse.Namespace) -> int:
    """自检：论文已发表序列基准 + 知识库校验。"""
    print("=" * 68)
    print("STIB Planner — selftest")
    print("=" * 68)

    print("\n[1] 序列层（基于 Table S2 已发表序列）\n")
    rows = seq_selftest()
    w = max(len(r[0]) for r in rows)
    n_ok = 0
    for name, ok, note in rows:
        flag = "PASS" if ok else "FAIL"
        n_ok += 1 if ok else 0
        print(f"  [{flag}] {name:<{w}}  {note}")
    print(f"\n  {n_ok}/{len(rows)} passed")

    print("\n[2] 知识库层（来源锚定 / 冲突显式）\n")
    kb = KnowledgeBase()
    probs = validate(kb)
    print(f"  facts      : {len(kb.facts)}")
    print(f"  conflicts  : {len(kb.conflicts)}"
          f"  (open_question={len(kb.conflicts_of('open_question'))},"
          f" contradiction={len(kb.conflicts_of('contradiction'))})")
    print(f"  sources    : {len(SOURCES)}")
    for p in probs:
        print(f"  [FAIL] {p}")
    print(f"\n  {'PASS' if not probs else 'FAIL'}: {len(probs)} problems")

    print("\n[3] 脱靶引擎（合成底物，错配数与位置均为已知设定）\n")
    off_rows = _offtarget_selftest()
    w2 = max(len(r[0]) for r in off_rows)
    n_off = 0
    for name, ok, note in off_rows:
        flag = "PASS" if ok else "FAIL"
        n_off += 1 if ok else 0
        print(f"  [{flag}] {name:<{w2}}  {note}")
    print(f"\n  {n_off}/{len(off_rows)} passed")

    total_ok = (n_ok == len(rows)) and not probs and (n_off == len(off_rows))
    print("\n" + "=" * 68)
    print("RESULT:", "ALL PASS" if total_ok else "FAILURES PRESENT")
    print("=" * 68)
    return 0 if total_ok else 1


def _offtarget_selftest() -> list[tuple[str, bool, str]]:
    """
    脱靶引擎自检。

    构造一个已知底物：埋入一条 23 nt spacer，再埋入几个**我手工设定错配数与
    错配位置**的变体。因此期望命中数是可以手算核对的。

    分区（23 nt）: distal = 0-7, middle = 8-14, seed = 15-22

    ⚠️ 底物设计的两条纪律（v0.3 踩过坑）：
      1. **每个植入片段必须唯一**。早期版本埋了多份完全相同的拷贝，
         结果「0 错配命中数」远超预期 —— 那不是引擎错，是底物里真有多个位点。
      2. **填充序列不得与植入片段的前缀撞车**。早期填充以 `TTGACC…` 开头，
         与用 `GAC` 打头的变体形成额外 1 错配位点。
         现在的填充由确定性伪随机生成，且每个植入片段两侧用独立的分隔串。
    """
    rows: list[tuple[str, bool, str]] = []

    SP = "ATTCGGAAGAAGACGTGGTATCA"          # 23 nt（论文 sgRNA26）
    FLIP = {"A": "C", "C": "G", "G": "T", "T": "A"}

    def mutate(seq: str, idxs) -> str:
        s = list(seq)
        for i in idxs:
            s[i] = FLIP[s[i]]
        return "".join(s)

    # 确定性伪随机填充（避开植入片段的前缀，保证位点唯一）
    def filler(n: int, seed: int) -> str:
        out: list[str] = []
        i = seed
        while len(out) < n:
            i = (i * 1103515245 + 12345) & 0x7FFFFFFF
            out.append("ACGT"[(i >> 16) & 3])
        return "".join(out)

    # 变体：按 v0.6 的 helix/tail 分区选取位置
    #   helix = 坐标 0–16（PAM 近端关键区）；tail = 坐标 17–22（外周）
    v1 = mutate(SP, [11])                   # helix 内 1 错配
    v2 = mutate(SP, [10, 11, 12])           # helix 内 3 错配
    v3 = mutate(SP, [2])                    # helix 内 1 错配（另一位点）
    v4 = mutate(SP, [20])                   # tail 内 1 错配

    seqs = [SP, v1, v2, v3, v4]
    # 唯一性自检：任何两段不得相同
    assert len(set(seqs)) == len(seqs), "自检底物含重复片段"

    genome = "".join(filler(60, 100 + 7 * k) + s for k, s in enumerate(seqs)) + filler(60, 999)
    idx = SeedIndex(genome)

    rows.append(("索引非空", idx.n_seeds > 0,
                 f"seeds={idx.n_seeds}, positions={idx.n_positions}"))

    rep = find_offtargets(idx, SP, max_mismatch=1)
    zero = [h for h in rep.hits if h.n_mismatch == 0]
    one = [h for h in rep.hits if h.n_mismatch == 1]
    rows.append(("完整匹配唯一", len(zero) == 1,
                 f"0 错配命中 {len(zero)} 个（期望 1，底物含 {len(set(seqs))} 段唯一序列）"))
    rows.append(("1 错配变体全部找到", len(one) == 3,
                 f"1 错配命中 {len(one)} 个（期望 3：helix×2 + tail×1）"))

    def region_of(hit) -> str:
        for k in ("helix", "tail"):
            if hit.by_region.get(k):
                return k
        return "none"

    regions = sorted(region_of(h) for h in one)
    rows.append(("分区判定正确（helix/tail）",
                 regions == ["helix", "helix", "tail"],
                 f"分区 = {regions}（期望 helix×2 + tail×1，"
                 f"helix = PAM 近端前 17 nt）"))

    hx_hits = [h for h in one if region_of(h) == "helix"]
    rows.append(("helix 1 错配 → high",
                 bool(hx_hits) and all(h.risk == "high" for h in hx_hits),
                 f"helix {len(hx_hits)} 个，风险 {[h.risk for h in hx_hits]}"))

    # 负链：单独造一个只含反向互补的底物，避免被正链命中干扰
    idx_rc = SeedIndex(filler(60, 555) + revcomp(SP) + filler(60, 777))
    rep_rc = find_offtargets(idx_rc, SP, max_mismatch=0)
    rows.append(("负链命中被找到", any(h.strand == "-" for h in rep_rc.hits),
                 f"命中 {rep_rc.n_hits} 个（0 错配）"))

    rep3 = find_offtargets(idx, SP, max_mismatch=3)
    v2_hits = [h for h in rep3.hits
               if h.n_mismatch == 3 and len(h.by_region.get("helix", [])) == 3]
    rows.append(("3 错配（全 helix）被找到且判 high",
                 bool(v2_hits) and all(h.risk == "high" for h in v2_hits),
                 f"命中 {len(v2_hits)} 个，风险 {[h.risk for h in v2_hits]}"))

    rep_lo = find_offtargets(idx, SP, max_mismatch=0)
    rows.append(("阈值 0 → 仅完整匹配", rep_lo.n_hits == 1,
                 f"n_hits={rep_lo.n_hits}（期望 1）"))

    rows.append(("≤4 标注为保证", rep3.guaranteed is True, f"max_mm=3 → {rep3.guaranteed}"))
    rep5 = find_offtargets(idx, SP, max_mismatch=5)
    rows.append((">4 标注为不保证", rep5.guaranteed is False, f"max_mm=5 → {rep5.guaranteed}"))

    try:
        _ = find_offtargets(idx, "A" * 23, max_mismatch=4)
        rows.append(("低复杂度 spacer 不崩溃", True, "A" * 23))
    except Exception as exc:  # noqa: BLE001
        rows.append(("低复杂度 spacer 不崩溃", False, repr(exc)))

    with_pam = [h for h in rep.hits if h.pam_nearby]
    rows.append(("记录邻近 PAM 序列", len(with_pam) > 0,
                 f"{len(with_pam)}/{rep.n_hits} 个命中带邻近序列"))

    return rows


def cmd_check(args: argparse.Namespace) -> int:
    """只校验知识库。"""
    kb = KnowledgeBase()
    probs = validate(kb)
    for p in probs:
        print(p)
    print(f"\n{len(kb.facts)} facts, {len(kb.conflicts)} conflicts, {len(probs)} problems")
    return 1 if probs else 0


def cmd_facts(args: argparse.Namespace) -> int:
    """列出知识条目及其出处。"""
    kb = KnowledgeBase()
    facts = kb.by_prefix(args.prefix) if args.prefix else kb.facts
    print(f"{len(facts)} facts\n")
    for f in facts:
        val = f.value
        if isinstance(val, dict):
            val = json.dumps(val, ensure_ascii=False)
        elif isinstance(val, str) and len(val) > 90:
            val = val[:87] + "..."
        print(f"  {f.key}")
        print(f"      value : {val}")
        print(f"      prov  : {f.ref()}")
        print(f"      level : {f.evidence}")
        if f.note:
            print(f"      note  : {f.note}")
        print()
    return 0


def cmd_conflicts(args: argparse.Namespace) -> int:
    """列出所有矛盾与未收敛变量（Demo-0 的核心产出之一）。"""
    kb = KnowledgeBase()
    for c in kb.conflicts:
        icon = {"high": "!!!", "medium": " ! ", "low": "   "}[c.severity]
        print(f"[{icon}] {c.key}  ({c.kind}, severity={c.severity})")
        for cl in c.claims:
            print(f"        - [{cl['source']} / {cl['locator']}] {cl['statement']}")
        if c.guidance:
            print(f"        => 建议: {c.guidance}")
        if c.discriminator:
            print(f"        => 判别实验: {c.discriminator}")
        print()
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    """从 NCBI 取回基因组并缓存为 FASTA。"""
    acc = args.accession
    outdir = args.out or "data"
    os.makedirs(outdir, exist_ok=True)
    dest = os.path.join(outdir, f"{acc}.fasta")

    if os.path.exists(dest) and not args.force:
        seqs = read_fasta(dest)
        n = sum(len(v) for v in seqs.values())
        print(f"已存在，跳过下载: {dest}  ({len(seqs)} 条序列, {n:,} bp)")
        return 0

    params = urllib.parse.urlencode({
        "db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text",
    })
    url = f"{NCBI_EFETCH}?{params}"
    print(f"下载 {acc} ...\n  {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "stib-planner/0.1 (research)"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    if not text.lstrip().startswith(">"):
        print("ERROR: 返回内容不是 FASTA。前 300 字符：")
        print(text[:300])
        return 2

    tmp = dest + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    seqs = read_fasta(tmp)
    os.replace(tmp, dest)
    n = sum(len(v) for v in seqs.values())
    print(f"完成: {dest}  ({len(seqs)} 条序列, {n:,} bp)")
    return 0


def _load_genome(path: str) -> GenomeIndex:
    """只加载序列（不关心注释）。见 _load_genome_ex 取得注释。"""
    g, _feats, _meta = _load_genome_ex(path)
    return g


def _load_genome_ex(path: str):
    """
    加载基因组，自动识别 FASTA / GenBank。

    返回 (GenomeIndex, features, meta)。

    **GenBank 才有 features** —— 这是按基因名定位的前提（基因名不在碱基序列里）。
    FASTA 时 features 为空列表，meta["has_annotations"] 为 False。
    """
    sid, seq, feats, meta = load_genome_any(path)
    if not seq:
        raise SystemExit(f"ERROR: {path} 中未读到序列")
    fmt = meta.get("format", "?")
    extra = ""
    if fmt == "genbank":
        extra = f"，{meta.get('n_cds', 0)} 个编码特征"
    print(f"基因组: {sid}  {len(seq):,} bp  [{fmt}{extra}]")
    return GenomeIndex(seq, sid), feats, meta


def cmd_scan(args: argparse.Namespace) -> int:
    """扫描一个区域的候选 sgRNA（不生成完整方案）。"""
    g, feats, _meta = _load_genome_ex(args.genome)
    span = _resolve_span(g, args, feats)
    _s, _e, solver = _resolve_span_ex(g, args, feats)
    print(f"定位方式    : {solver.get('note', '?')}")
    cands = dedupe_sgrnas(scan_sgrnas(g.seq, span[0], span[1]))
    print(f"扫描区间    : {span[0]:,} – {span[1]:,}")
    print(f"候选数      : {len(cands)}  (GTN PAM, spacer {SPACER_LEN} nt)\n")
    print(f"{'#':>3}  {'spacer':<26} {'PAM':<4} {'链':<2} {'PAM位置':>12}  {'预测插入':>12}  GC")
    for i, c in enumerate(cands[: args.top], 1):
        print(f"{i:>3}  {c.spacer:<26} {c.pam:<4} {c.strand:<2} "
              f"{c.pam_start:>12,}  {c.insert_site:>12,}  {c.gc:.0%}")
    if len(cands) > args.top:
        print(f"\n... 另有 {len(cands) - args.top} 条未显示（用 --top 调整）")
    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    """生成完整编辑方案。"""
    g, feats, _meta = _load_genome_ex(args.genome)
    span = _resolve_span(g, args, feats)
    _s, _e, solver = _resolve_span_ex(g, args, feats)
    print(f"目标区间    : {_s:,} – {_e:,}（{_e - _s:,} bp）")
    print(f"定位方式    : {solver.get('note', '?')}")
    kb = KnowledgeBase()

    plan = make_plan(
        genome=g, gene=args.gene, start=span[0], end=span[1],
        species=args.species, top_n=args.top,
        max_mismatch=args.max_mismatch, community_editing=args.community, kb=kb,
    )
    plan.meta["genome_len"] = len(g)     # 供可视化层画全基因组视图

    # --detailed：输出**可照着做**的定制实验协议（含体系体积、程序、时间轴）
    md = render_protocol(plan, kb, species=args.species) if args.detailed \
        else render_markdown(plan)

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(md)
        print(f"方案已写出: {args.out}")
    else:
        print(md)

    if args.json:
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_json(plan))
        print(f"JSON 已写出: {args.json}")

    if args.html:
        os.makedirs(os.path.dirname(os.path.abspath(args.html)), exist_ok=True)
        with open(args.html, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_html(plan))
        print(f"HTML 可视化已写出: {args.html}（双击即可打开，无外部依赖）")

    # 终端摘要
    t = plan.target
    print("\n--- 摘要 ---")
    print(f"目标基因        : {t['gene']}")
    print(f"推荐 spacer     : {t['selected_spacer']}  (PAM {t['selected_pam']})")
    print(f"预测插入位点    : {t['predicted_insert_site']:,}")
    print(f"氯霉素浓度      : {t['cm_conc']}")
    print(f"候选通过/淘汰   : {t['n_passed']} / {t['n_rejected']}（共扫描 {t['n_candidates_scanned']}）")
    ot_note = "保证不漏检" if t["offtarget_guaranteed"] else "**可能漏检**"
    print(f"脱靶阈值        : ≤{t['max_mismatch']} 错配（{ot_note}）")
    print(f"首选落点        : "
          + ({"strict": "★严格区(CDS 前1/3-2/3)", "central80": "中段80%(排除前后5%)",
              "in_gene": "基因端部5%", "outside": "基因外"}.get(plan.candidates[0].placement, "-")
             if plan.candidates else "n/a"))
    print(f"首选脱靶        : 风险分 "
          + (f"{plan.candidates[0].risk_score:.1f}" if plan.candidates else "n/a")
          + f"（阈值 {RISK_THRESHOLD}）｜ 非自身近同源 "
          + (f"{plan.candidates[0].n_other} 处" if plan.candidates else "n/a"))
    print(f"高等级警示      : {sum(1 for w in plan.warnings if w.severity == 'high')} 条")
    print(f"风险条目        : {len(plan.risks)} 条")
    print(f"判定点          : {len(plan.decisions)} 个")
    print(f"回退面板        : {len(plan.fallbacks)} 条已验证 spacer")
    return 0


def cmd_batch(args: argparse.Namespace) -> int:
    """
    批量设计模式：在一个较宽区间内一次给出多条候选。

    用途：为某个基因簇/区域准备候选池，或做小规模筛选。
    """
    g, feats, _meta = _load_genome_ex(args.genome)
    span = _resolve_span(g, args, feats)
    _s, _e, solver = _resolve_span_ex(g, args, feats)
    print(f"定位方式    : {solver.get('note', '?')}")
    plan = batch_plan(g, gene=args.gene, start=span[0], end=span[1],
                      species=args.species, top_n=args.top)
    plan.meta["genome_len"] = len(g)

    t = plan.target
    print(f"批量模式：{args.gene}")
    print(f"区间            : {span[0]:,} – {span[1]:,}（{span[1]-span[0]:,} bp）")
    print(f"候选 扫描/通过/淘汰: {t['n_candidates_scanned']} / "
          f"{t['n_passed']} / {t['n_rejected']}")
    print()
    print(f"{'#':>3}  {'spacer':<26}{'PAM':<5}{'落点':<10}{'GC':<5}"
          f"{'GCWGC':<7}{'±5A/T':<7}{'风险分':<7}{'总分'}")
    for i, c in enumerate(plan.candidates, 1):
        sf = c.seq_features.as_dict() if c.seq_features else {}
        place = {"strict": "★严格区", "central80": "中央80%",
                 "in_gene": "基因端部", "outside": "基因外"}.get(c.placement, c.placement)
        print(f"{i:>3}  {c.sgrna.spacer:<26}{c.sgrna.pam:<5}{place:<10}"
              f"{c.sgrna.gc:>4.0%} {('✓' if sf.get('tsd_preferred_GCWGC') else '—'):<7}"
              f"{('✓' if sf.get('at_pm5_ok') else '—'):<7}"
              f"{c.risk_score:<7.2f}{c.score:.0f}")

    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_markdown(plan))
        print(f"\nMarkdown 已写出: {args.out}")
    if args.json:
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_json(plan))
        print(f"JSON 已写出: {args.json}")
    if args.html:
        with open(args.html, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_html(plan))
        print(f"HTML 可视化已写出: {args.html}")
    return 0


def cmd_offtarget(args: argparse.Namespace) -> int:
    """单独评估一条或多条 spacer 的全基因组脱靶。"""
    g, _feats, _meta = _load_genome_ex(args.genome)
    print("建索引中 ...", end="", flush=True)
    idx = SeedIndex(g.seq)
    print(f" 完成（{idx.n_seeds:,} 个种子 / {idx.n_positions:,} 个位置）\n")

    spacers = [s.strip().upper() for s in (args.spacer or []) if s.strip()]
    if args.panel:
        kb = KnowledgeBase()
        spacers += [d["spacer"].split(" / ")[0] for d in valid_spacer_panel()]
    if not spacers:
        print("ERROR: 需提供 --spacer 或 --panel")
        return 2

    rc = {"high": "🔴", "medium": "🟠", "low": "🟡"}
    for sp in spacers:
        rep = find_offtargets(idx, sp, max_mismatch=args.max_mismatch,
                              max_report=args.top)
        # ⚠️ 候选**必然精确命中它自己**（0 错配）。这是正常的，不是脱靶。
        # 风险判定与显示都必须排除自身位点。
        self_hits = [h for h in rep.hits if h.n_mismatch == 0]
        other = [h for h in rep.hits if h.n_mismatch > 0]
        n_copy = len(self_hits)
        n_hi = sum(1 for h in other if h.risk == "high")
        n_mid = sum(1 for h in other if h.risk == "medium")
        n_lo = sum(1 for h in other if h.risk == "low")

        risk_score, _bd = offtarget_risk(rep)
        if not other:
            verdict = "✅ 干净（除自身位点外无近同源）"
        elif n_hi or risk_score > RISK_THRESHOLD:
            verdict = f"🔴 不建议（风险分 {risk_score:.1f} > {RISK_THRESHOLD}）"
        elif n_mid:
            verdict = f"🟠 可用但需注意（风险分 {risk_score:.1f}）"
        else:
            verdict = f"🟡 可接受（风险分 {risk_score:.1f}）"

        print(f"{sp}")
        copy_note = "（仅 1 处）" if n_copy <= 1 else f"（**{n_copy} 处，串联重复**）"
        print(f"   精确命中 {n_copy} 处{copy_note}｜ "
              f"非自身近同源 {len(other)} 处 "
              f"(高 {rc['high']}{n_hi} / 中 {rc['medium']}{n_mid} / 低 {rc['low']}{n_lo})")
        print(f"   风险分 {risk_score:.1f}（阈值 {RISK_THRESHOLD}）｜ "
              f"错配阈值 ≤{rep.max_mismatch}｜ 候选扫描 {rep.n_candidates_scanned:,}")
        print(f"   判定: {verdict}")
        for h in other[: args.top]:
            dist = ", ".join(f"{k}:{len(v)}" for k, v in h.by_region.items() if v)
            print(f"       {h.pos:>12,} {h.strand} 错配{h.n_mismatch}  [{dist}]  {h.risk}")
        for w in rep.warnings:
            print(f"   [WARN] {w}")
        print()
    return 0



def cmd_genes(args: argparse.Namespace) -> int:
    """列出/搜索基因组注释里的基因（需要 GenBank）。"""
    g, feats, meta = _load_genome_ex(args.genome)
    if not feats:
        raise SystemExit(
            "ERROR: 当前基因组文件没有注释，无法列出基因。\n"
            "纯 FASTA 只有碱基序列。请改用 GenBank 格式（.gb/.gbk/.gbff）。")
    cds = [f for f in feats if f.key in ("CDS", "tRNA", "rRNA", "ncRNA")]
    if args.search:
        hits = [f for f in cds
                if any(args.search.lower() in f.qualifiers.get(k, "").lower()
                       for k in ("gene", "locus_tag", "old_locus_tag", "product"))]
    else:
        hits = cds
    # 基因组基本信息已由 _load_genome_ex 打印，此处只报编码特征数
    print(f"编码特征    : {len(cds)} 个")
    if args.search:
        print(f"搜索        : \"{args.search}\" → {len(hits)} 条")
    print()
    print(f"{'基因':<14}{'locus_tag':<14}{'位置':>20}  {'链':<3} product")
    for f in hits[: args.top]:
        gene = f.qualifiers.get("gene", "")
        tag = f.qualifiers.get("locus_tag", "")
        prod = f.qualifiers.get("product", "")[:44]
        loc = f"{f.start + 1:,}–{f.end:,}"
        print(f"{gene:<14}{tag:<14}{loc:>20}  {f.strand:<3} {prod}")
    if len(hits) > args.top:
        print(f"... 另有 {len(hits) - args.top} 条（用 --top 调整）")
    return 0


def cmd_neutral(args: argparse.Namespace) -> int:
    """
    查找中性位点（用于插入大片段 cargo），并为每个位点给出 sgRNA。

    ⚠️ 「中性位点」的定义是「不在任何基因内部」—— **必须有注释才能判断**。
    """
    from .neutsche import find_neutral_sites, sgrnas_for_neutral_site, KNOWN_NEUTRAL_SITES

    g, feats, meta = _load_genome_ex(args.genome)
    if not feats:
        raise SystemExit(
            "ERROR: 无法查找中性位点：当前基因组文件**没有注释**。\n"
            "「中性位点」= 不在任何基因内部，这需要注释才能判断；\n"
            "纯 FASTA 无法区分基因与基因间区。请改用 GenBank 格式。")

    print()
    print("⚠️ 注意：论文中**已验证**的中性位点（下方列出）价值高于本工具从序列特征"
          "预测的候选。")
    for k in KNOWN_NEUTRAL_SITES:
        print(f"     {k['label']:<18} {k['pos']:>12,}  {k['strand']}  {k['note']}")
    print()

    try:
        sites = find_neutral_sites(g.seq, feats, min_gap=args.min_gap,
                                   max_candidates=args.top)
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}") from None

    print(f"中性位点候选: {len(sites)} 个（基因间区 ≥{args.min_gap} bp）")
    print()
    for i, site in enumerate(sites, 1):
        print(f"#{i}  [{site.start:,}–{site.end:,}]  {site.length} bp  "
              f"评分 {site.score:.1f}")
        print(f"     左邻 {site.left_gene or '?'}  右邻 {site.right_gene or '?'}  "
              f"GC {site.gc:.0%}")
        for r in site.reasons:
            print(f"     + {r}")
        for w in site.warnings:
            print(f"     ! {w}")
        sgs = sgrnas_for_neutral_site(g.seq, site, window=args.window, top=args.sg)
        if sgs:
            print(f"     可用 sgRNA（预测插入位点落在该区间内）:")
            for d in sgs:
                print(f"       {d['spacer']}  {d['pam']}@{d['pam_start']:,}  "
                      f"insert={d['insert_site']:,}  距中心 {d['distance_to_center']} bp")
        else:
            print("     ! 该区间内未找到合适的 sgRNA（可能需要放宽 window）")
        print()

    if args.out:
        lines = ["# 中性位点候选\n",
                 f"基因组: {g.id}（{len(g):,} bp）\n\n",
                 "## 论文已验证的位点（优先考虑）\n\n"]
        for k in KNOWN_NEUTRAL_SITES:
            lines.append(f"- **{k['label']}** @ {k['pos']:,} ({k['strand']})：{k['note']}\n")
        lines.append("\n## 序列特征预测的候选\n\n")
        lines.append("| # | 区间 | 长度 | 评分 | 左邻 | 右邻 | GC | 告警 |\n")
        lines.append("|---|---|---|---|---|---|---|---|\n")
        for i, site in enumerate(sites, 1):
            warns = "；".join(site.warnings) or "—"
            lines.append(f"| {i} | {site.start:,}–{site.end:,} | {site.length} | "
                         f"{site.score:.1f} | {site.left_gene or '?'} | "
                         f"{site.right_gene or '?'} | {site.gc:.0%} | {warns} |\n")
        lines.append("\n> ⚠️ 这些是**序列特征预测**的候选，不是已验证位点。\n")
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("".join(lines))
        print(f"已写出: {args.out}")
    return 0


def cmd_cargo(args: argparse.Namespace) -> int:
    """
    为 cargo 设计 Gibson 组装引物，并输出**保留骨架全部注释**的可视化图谱。

    两种策略（见 stib/cargo_assembly.py）：
      keep_cmr（默认）   保留 PcfxA-CmR，新 cargo 插在 LE 之后、筛选盒之前
      replace-selection  用 cargo 替换整段 LE–RE 内容（失去氯霉素筛选）
    """
    from .cargo import read_cargo_fasta, check_cargo  # noqa: F401
    from .cargo_assembly import (STRATEGY_KEEP_CMR, STRATEGY_REPLACE_SELECTION,
                                 build_construct_map, cargo_report)
    from .construct import load_spec, features_from_snapgene
    from .buildout import read_snapgene, roundtrip_check, write_snapgene_file, \
        write_genbank, write_fasta

    spec = load_spec(args.spec)

    if args.list_sites:
        print(f"载体规范: {args.spec}")
        print(f"载体    : {spec.meta.get('name')}  "
              f"{spec.meta.get('backbone_length', 0):,} bp")
        print()
        print("可选的 cargo 插入位点：")
        for s_ in spec.cargo_slot.sites:
            bl = spec.blockers_in(*s_.region)
            gaps = spec.free_gaps(*s_.region)
            mx = max((b - a for a, b in gaps), default=0)
            print(f"\n  【{s_.id}】 {s_.label}")
            print(f"      区间        : {s_.region[0]:,}–{s_.region[1]:,}"
                  f"（{s_.region[1] - s_.region[0]:,} bp）")
            print(f"      cargo 方向  : {'负链（与 CmR 同向）' if s_.cargo_strand < 0 else '正链'}")
            if getattr(s_, "replace_all", False):
                print(f"      语义        : **整段替换** —— cargo 长度"
                      "不受区间长度限制；区间内元件会被一并替换")
                if bl:
                    print(f"      将被替换掉的: " + "；".join(
                        f"{b['name']}（{b['end'] - b['start']} bp）" for b in bl))
            elif bl:
                print(f"      区间内需保留: " + "；".join(
                    f"{b['name']}（{b['end'] - b['start']} bp）" for b in bl))
                print(f"      可用空档    : " + " │ ".join(
                    f"{a:,}–{b:,}={b - a} bp" for a, b in gaps))
            else:
                print(f"      可用空档    : 全区间（无内置元件）")
            print(f"      说明        : {' '.join(s_.note.split())}")
        return 0

    print(f"载体规范    : {args.spec}")
    print(f"载体        : {spec.meta.get('name')}  "
          f"{spec.meta.get('backbone_length', 0):,} bp  "
          f"{spec.meta.get('topology')}")
    print(f"cargo 槽位  : {spec.cargo_slot.start:,}–{spec.cargo_slot.end:,}"
          f"（{spec.cargo_slot.length:,} bp）")
    print(f"转座子臂    : LE {spec.le.core_start:,}–{spec.le.core_end:,} ； "
          f"RE {spec.re.core_start:,}–{spec.re.core_end:,}")
    cm_region = spec.cmr_region()
    if cm_region:
        print(f"筛选盒      : CmR+PcfxA {cm_region[0]:,}–{cm_region[1]:,}")

    # ---- 读载体（序列 + 全部注释）----
    bb = read_snapgene(open(args.backbone, "rb").read())
    if not bb["seq"]:
        raise SystemExit("ERROR: 载体文件里没有序列")
    bbseq = bb["seq"].upper()
    feats = features_from_snapgene(bb)
    print(f"载体注释    : {len(feats)} 个特征（将全部保留到输出图谱）")

    # ---- 读 cargo ----
    cname, cargo = read_cargo_fasta(args.cargo)
    print(f"cargo       : {args.cargo}  名={cname}  {len(cargo):,} bp")

    chk = check_cargo(cargo, name=cname)
    print()
    print("=== cargo 校验 ===")
    print(f"  ok={chk.ok}  {chk.stats}")
    for e in chk.errors:
        print("  ERROR:", e)
    for w in chk.warnings:
        print("  WARN :", w)
    if not chk.ok:
        raise SystemExit("cargo 未通过硬性校验，停止。")

    # ---- 策略 ----
    if args.replace_selection:
        strategy = STRATEGY_REPLACE_SELECTION
    else:
        strategy = STRATEGY_KEEP_CMR
    rr = None
    if args.replace_region:
        rr = tuple(int(x) for x in args.replace_region)

    cm, plan, warns = build_construct_map(
        spec, bbseq, feats, cargo=cargo, cargo_name=cname,
        strategy=strategy, overlap_len=args.overlap_len, flank=args.flank,
        replace_region=rr, site_id=args.site, cargo_strand=args.cargo_strand)

    print()
    print("=== 落点与筛选标记 ===")
    for w in warns:
        print(f"  {w}")

    print()
    print("=== 组装引物 ===")
    for pr in plan.primers:
        print(f"  {pr.name}")
        print(f"    序列 5'->3' : {pr.seq}")
        print(f"    overlap(5') : {pr.overlap}")
        print(f"    退火段(3')  : {pr.anneal}")
        print(f"    Tm 估算     : {pr.tm_full:.1f} °C（全长，估算值）")

    print()
    print("=== 图谱（保留骨架全部注释）===")
    print(f"  总长 {len(cm.seq):,} bp；特征 {len(cm.features)} 个；"
          f"引物 {len(cm.primers)} 条")
    has_cmr = any(f.name == "CmR" for f in cm.features)
    print(f"  CmR 仍在图谱中: {'是 ✓（可继续用氯霉素筛选）' if has_cmr else '否 ✗'}")

    print()
    print("=== 图谱往返校验 ===")
    rt_ok, diffs = roundtrip_check(cm)
    print(f"  一致: {rt_ok}")
    for d in diffs:
        print("  DIFF:", d)

    # ---- 输出 ----
    os.makedirs(args.out, exist_ok=True)
    base = os.path.join(args.out, cm.name)
    write_snapgene_file(base + ".dna", cm)
    with open(base + ".gb", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(write_genbank(cm))
    with open(base + ".fa", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(write_fasta(cm))
    with open(base + ".md", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(cargo_report(spec, cm, plan, warns, cm.comments))

    print()
    print("=== 产出 ===")
    for ext in (".dna", ".gb", ".fa", ".md"):
        print(f"  {os.path.abspath(base + ext)}")
    print()
    print("  → .dna 用 SnapGene 打开查看 cargo / LE / RE / 引物位置")
    return 0


def _resolve_span(g: GenomeIndex, args: argparse.Namespace,
                  features: list | None = None) -> tuple[int, int]:
    """解析目标区间（兼容旧调用）。见 _resolve_span_ex。"""
    s, e, _info = _resolve_span_ex(g, args, features)
    return s, e


def _resolve_span_ex(g: GenomeIndex, args: argparse.Namespace,
                     features: list | None = None) -> tuple[int, int, dict]:
    """
    解析目标区间，返回 (start, end, resolver_info)。

    四种模式，优先级：坐标 > 序列片段 > 基因名

      --start/--end   直接坐标
      --motif SEQ     在基因组中定位该序列（FASTA 也可用）
      --gene NAME     按基因名定位（**需要 GenBank 注释**）

    resolver_info 记录本次用的是什么定位方式 —— 报告里如实标注，
    因为使用者需要知道这个区间是查到的还是他给的。
    """
    try:
        s, e, info = resolve_region(
            g.seq, features or [],
            gene=getattr(args, "gene", None),
            start=getattr(args, "start", None),
            end=getattr(args, "end", None),
            motif=getattr(args, "motif", None),
            pad=int(getattr(args, "pad", 500) or 500),
        )
    except RegionError as exc:
        raise SystemExit(f"ERROR: {exc}") from None
    return s, e, info


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stib",
        description="STIB Planner — 拟杆菌 STIB 基因编辑方案生成器（Demo-0）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("selftest", help="自检（序列层 + 知识库层）")
    sp.set_defaults(func=cmd_selftest)

    sp = sub.add_parser("check", help="校验知识库")
    sp.set_defaults(func=cmd_check)

    sp = sub.add_parser("facts", help="列出知识条目及出处")
    sp.add_argument("--prefix", help="按 key 前缀过滤，如 risk. 或 reagent.")
    sp.set_defaults(func=cmd_facts)

    sp = sub.add_parser("conflicts", help="列出矛盾与未收敛变量")
    sp.set_defaults(func=cmd_conflicts)

    sp = sub.add_parser("fetch", help="从 NCBI 取回基因组")
    sp.add_argument("accession", help="如 NC_004663.1")
    sp.add_argument("--out", help="输出目录（默认 data/）")
    sp.add_argument("--force", action="store_true", help="覆盖已存在文件")
    sp.set_defaults(func=cmd_fetch)

    sp = sub.add_parser("scan", help="扫描候选 sgRNA")
    sp.add_argument("--genome", required=True, help="基因组文件（FASTA 或 GenBank）")
    sp.add_argument("--start", type=int)
    sp.add_argument("--end", type=int)
    sp.add_argument("--gene", help="按基因名定位（需 GenBank 格式基因组）")
    sp.add_argument("--motif", help="定位该序列并以其为中心取区间")
    sp.add_argument("--pad", type=int, default=500, help="motif 模式下的两侧扩展量")
    sp.add_argument("--top", type=int, default=20, help="显示条数")
    sp.set_defaults(func=cmd_scan)

    sp = sub.add_parser("plan", help="生成完整编辑方案")
    sp.add_argument("--genome", required=True, help="基因组文件（FASTA 或 GenBank）")
    sp.add_argument("--gene", required=True,
                    help="目标基因名。GenBank 基因组时**同时用于定位**；"
                         "FASTA 时仅作标注，需另给 --motif 或 --start/--end")
    sp.add_argument("--species", default="Bacteroides thetaiotaomicron")
    sp.add_argument("--start", type=int)
    sp.add_argument("--end", type=int)
    sp.add_argument("--motif", help="定位该序列并以其为中心取区间（FASTA 也能用）")
    sp.add_argument("--pad", type=int, default=500)
    sp.add_argument("--top", type=int, default=5, help="返回候选数")
    sp.add_argument("--max-mismatch", type=int, default=4,
                    help=f"脱靶枚举错配上限。≤{GUARANTEED_MAX_MISMATCH} 保证不漏检；"
                         ">4 可能少报（会显式警示）")
    sp.add_argument("--out", help="输出 Markdown 路径")
    sp.add_argument("--json", help="同时输出 JSON 路径")
    sp.add_argument("--html", help="同时输出 HTML 可视化（双击即可打开）")
    sp.add_argument("--community", action="store_true",
                    help="菌群编辑模式（接合气氛转为厌氧推荐）")
    sp.add_argument("--detailed", action="store_true",
                    help="输出**完整定制实验协议**（含 Golden Gate/Gibson 体系体积、"
                         "程序、时间轴、核对清单），而非决策摘要")
    sp.set_defaults(func=cmd_plan)

    sp = sub.add_parser("batch", help="批量设计：一个较宽区间内给出多条候选")
    sp.add_argument("--genome", required=True, help="基因组文件（FASTA 或 GenBank）")
    sp.add_argument("--gene", required=True,
                    help="区域名。GenBank 基因组时同时用于定位；FASTA 时仅作标注")
    sp.add_argument("--species", default="Bacteroides thetaiotaomicron")
    sp.add_argument("--start", type=int)
    sp.add_argument("--end", type=int)
    sp.add_argument("--motif", help="定位该序列并以其为中心取区间")
    sp.add_argument("--pad", type=int, default=2000,
                    help="motif 模式下的两侧扩展（批量模式默认更大）")
    sp.add_argument("--top", type=int, default=24, help="返回候选数")
    sp.add_argument("--out", help="输出 Markdown 路径")
    sp.add_argument("--json", help="输出 JSON 路径")
    sp.add_argument("--html", help="输出 HTML 可视化")
    sp.set_defaults(func=cmd_batch)

    sp = sub.add_parser("genes", help="列出/搜索基因（需 GenBank 注释）")
    sp.add_argument("--genome", required=True, help="基因组文件（GenBank）")
    sp.add_argument("--search", help="按关键词搜索基因名/locus_tag/product")
    sp.add_argument("--top", type=int, default=40, help="显示条数")
    sp.set_defaults(func=cmd_genes)

    sp = sub.add_parser("neutral", help="查找中性位点（插入大片段用，需 GenBank 注释）")
    sp.add_argument("--genome", required=True, help="基因组文件（GenBank）")
    sp.add_argument("--min-gap", type=int, default=120, help="最小基因间区长度")
    sp.add_argument("--window", type=int, default=120,
                    help="在插入点两侧多大范围内找 sgRNA")
    sp.add_argument("--sg", type=int, default=3, help="每个位点显示几条 sgRNA")
    sp.add_argument("--top", type=int, default=8, help="候选位点数")
    sp.add_argument("--out", help="输出 Markdown 路径")
    sp.set_defaults(func=cmd_neutral)

    sp = sub.add_parser("cargo", help="为 cargo 设计 Gibson 组装引物 + 输出可视化图谱")
    sp.add_argument("--cargo", required=True, help="cargo 的 FASTA 文件")
    sp.add_argument("--backbone", required=True, help="载体文件（.dna / .gb）")
    _here = os.path.dirname(os.path.abspath(__file__))
    sp.add_argument("--spec",
                    default=os.path.join(_here, "data", "construct_psh43.yml"),
                    help="构建规范 YAML（默认 stib/data/construct_psh43.yml）")
    sp.add_argument("--overlap-len", type=int, default=22, help="Gibson overlap 长度")
    sp.add_argument("--flank", type=int, default=30,
                    help="从槽位两侧取多少 bp 作为 overlap 依据")
    sp.add_argument("--site", help="插入位点 id（见规范 cargo_slot.sites）；默认 after_le")
    sp.add_argument("--cargo-strand", type=int, choices=[1, -1],
                    help="cargo 所在链；-1 = 与 CmR 同向（cargo 需表达时建议）")
    sp.add_argument("--list-sites", action="store_true", help="只列出可选插入位点")
    sp.add_argument("--replace-selection", action="store_true",
                    help="连 PcfxA-CmR 筛选盒一起替换（会失去氯霉素筛选）")
    sp.add_argument("--replace-region", nargs=2, type=int, metavar=("START", "END"),
                    help="手动指定被替换的载体区间（0-based 半开）；高级用法")
    sp.add_argument("--out", default="out", help="输出目录")
    sp.set_defaults(func=cmd_cargo)

    sp = sub.add_parser("offtarget", help="单独评估 spacer 的全基因组脱靶")
    sp.add_argument("--genome", required=True, help="基因组 FASTA")
    sp.add_argument("--spacer", action="append",
                    help="待评估的 23 nt spacer（可重复指定）")
    sp.add_argument("--panel", action="store_true",
                    help="评估内置的已验证 spacer 面板")
    sp.add_argument("--max-mismatch", type=int, default=4)
    sp.add_argument("--top", type=int, default=8, help="每条最多显示几个命中")
    sp.set_defaults(func=cmd_offtarget)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:
        return 0
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
