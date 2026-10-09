"""
STIB Planner — 可视化层。

产出**自包含的单文件 HTML 报告**（内联 SVG + CSS，无外部依赖、无 JS 依赖、可离线打开）。

设计理由:
  用户在命令行看到的是一堆坐标数字，难以建立空间直觉。一张图能立刻回答
  「这个位点落在基因哪里」「候选之间相隔多远」「脱靶位点分布在哪儿」。

⚠️ 刻意不使用任何 CDN / 外部字体 / 前端框架 —— 目标是**双击就能看**，
   在离线的实验室电脑上也能打开。
"""

from __future__ import annotations

import html
from typing import Any

# ---------------------------------------------------------------------------
# 绘制原语
# ---------------------------------------------------------------------------

#: 配色（与 Markdown 报告的风险图标一致）
C = {
    "bg": "#ffffff",
    "ink": "#1f2328",
    "muted": "#6b7280",
    "line": "#d0d7de",
    "gene": "#3b82f6",
    "gene_light": "#dbeafe",
    "strict": "#16a34a",
    "central": "#f59e0b",
    "outside": "#9ca3af",
    "high": "#dc2626",
    "medium": "#f59e0b",
    "low": "#eab308",
    "ok": "#16a34a",
    "pam1": "#16a34a",
    "pam2": "#f59e0b",
    "pam3": "#9ca3af",
}


def _esc(s: Any) -> str:
    return html.escape(str(s))


def _scale(v: float, lo: float, hi: float, width: float) -> float:
    if hi <= lo:
        return 0.0
    return (v - lo) / (hi - lo) * width


# ---------------------------------------------------------------------------
# 基因/位点图
# ---------------------------------------------------------------------------

def locus_svg(candidates: list[Any], gene_start: int, gene_end: int,
              genome_len: int, width: int = 900, pad: int = 60,
              window: int | None = None) -> str:
    """
    画目标基因与候选插入位点的相对位置。

    window: 若给定，则只画基因区间 ±window 的范围（放大视图）；否则画全基因组。
    """
    if window is not None:
        lo = max(0, gene_start - window)
        hi = min(genome_len, gene_end + window)
    else:
        lo, hi = 0, genome_len

    track_w = width - 2 * pad
    row_h = 26

    parts: list[str] = []
    h = 90 + row_h * max(1, len(candidates)) + 40
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{h}" '
        f'viewBox="0 0 {width} {h}" font-family="system-ui,-apple-system,Segoe UI,'
        f'sans-serif" font-size="12">'
    )
    parts.append(f'<rect width="{width}" height="{h}" fill="{C["bg"]}"/>')

    y0 = 50
    # 基因组轴线
    parts.append(f'<line x1="{pad}" y1="{y0}" x2="{pad + track_w}" y2="{y0}" '
                 f'stroke="{C["line"]}" stroke-width="2"/>')

    # 基因体
    gx1 = pad + _scale(gene_start, lo, hi, track_w)
    gx2 = pad + _scale(gene_end, lo, hi, track_w)
    parts.append(f'<rect x="{gx1:.1f}" y="{y0 - 9}" width="{max(2, gx2 - gx1):.1f}" '
                 f'height="18" rx="3" fill="{C["gene_light"]}" stroke="{C["gene"]}" '
                 f'stroke-width="1.5"/>')

    # 严格区（前 1/3 – 2/3）
    span = max(1, gene_end - gene_start)
    s1 = gene_start + span / 3.0
    s2 = gene_start + span * 2.0 / 3.0
    sx1 = pad + _scale(s1, lo, hi, track_w)
    sx2 = pad + _scale(s2, lo, hi, track_w)
    parts.append(f'<rect x="{sx1:.1f}" y="{y0 - 9}" width="{max(1, sx2 - sx1):.1f}" '
                 f'height="18" fill="{C["strict"]}" opacity="0.22"/>')

    # 标签
    parts.append(f'<text x="{pad}" y="26" fill="{C["ink"]}" font-size="14" '
                 f'font-weight="600">目标基因与候选插入位点</text>')
    parts.append(f'<text x="{pad}" y="42" fill="{C["muted"]}">'
                 f'区间 {lo:,} – {hi:,}'
                 f'{"（放大视图）" if window is not None else "（全基因组）"}'
                 f'　绿底 = 前 1/3–2/3 严格区</text>')

    # 候选
    for i, c in enumerate(candidates):
        yy = y0 + 34 + i * row_h
        pos = c.sgrna.insert_site
        x = pad + _scale(pos, lo, hi, track_w)
        col = {"strict": C["strict"], "central80": C["central"],
               "in_gene": C["outside"], "outside": C["outside"]}.get(c.placement, C["outside"])
        # 从轴线到候选行画引导线
        parts.append(f'<line x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{yy - 5}" '
                     f'stroke="{col}" stroke-width="1" stroke-dasharray="2,2" opacity="0.6"/>')
        parts.append(f'<circle cx="{x:.1f}" cy="{yy - 5}" r="5" fill="{col}"/>')
        # spacer 文本（截断显示）
        sp = c.sgrna.spacer
        pam_s = {1: "①", 2: "②"}.get(c.pam_tier, "③")
        label = f"{i+1}. {sp[:10]}…{sp[-4:]}  {pam_s}{c.sgrna.pam}  @{pos:,}"
        tx = x + 10 if x < pad + track_w * 0.6 else x - 10
        anchor = "start" if x < pad + track_w * 0.6 else "end"
        parts.append(f'<text x="{tx:.1f}" y="{yy}" fill="{C["ink"]}" '
                     f'text-anchor="{anchor}">{_esc(label)}</text>')

    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# 候选评分条形图
# ---------------------------------------------------------------------------

def scores_svg(candidates: list[Any], width: int = 900, pad: int = 60) -> str:
    """画候选的综合排序分（横向条形图），并标注构成要素。"""
    if not candidates:
        return ""
    row_h = 30
    h = 70 + row_h * len(candidates) + 30
    track_w = width - 2 * pad - 180

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{h}" '
        f'viewBox="0 0 {width} {h}" font-family="system-ui,-apple-system,Segoe UI,'
        f'sans-serif" font-size="12">',
        f'<rect width="{width}" height="{h}" fill="{C["bg"]}"/>',
        f'<text x="{pad}" y="26" fill="{C["ink"]}" font-size="14" font-weight="600">'
        f'候选综合排序分（相对值，不代表效率）</text>',
        f'<text x="{pad}" y="44" fill="{C["muted"]}">'
        f'构成：落点位置 + PAM 优选 + GC + ShCAST 序列偏好 − 脱靶风险</text>',
    ]

    mx = max((c.score for c in candidates), default=1.0) or 1.0
    mn = min((c.score for c in candidates), default=0.0)
    rng = max(1e-6, mx - mn)

    for i, c in enumerate(candidates):
        y = 62 + i * row_h
        bl = pad + _scale(c.score, mn, mx, track_w)
        col = C["pam1"] if c.pam_tier == 1 else (C["pam2"] if c.pam_tier == 2 else C["pam3"])
        parts.append(f'<rect x="{pad}" y="{y}" width="{max(2, bl):.1f}" height="16" '
                     f'rx="3" fill="{col}" opacity="0.85"/>')
        sp = c.sgrna.spacer
        parts.append(f'<text x="{pad}" y="{y + 12}" fill="#ffffff" font-size="11" '
                     f'font-weight="600">　{_esc(sp[:14])}…</text>')
        parts.append(f'<text x="{pad + track_w + 12}" y="{y + 12}" fill="{C["ink"]}">'
                     f'{c.score:.0f}　GC {c.sgrna.gc:.0%}　risk {c.risk_score:.1f}</text>')

    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# 完整 HTML 报告
# ---------------------------------------------------------------------------

def render_html(plan: Any, batch: list[Any] | None = None) -> str:
    """
    产出自包含 HTML 报告。

    plan  : planner.Plan
    batch : 可选的批量候选列表（planner.Candidate），用于附加一节批量视图
    """
    t = plan.target
    gene_start, gene_end = t["span"]
    genome_len = plan.meta.get("genome_len") or max(gene_end * 4, gene_end + 10000)

    parts: list[str] = []
    parts.append("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'>")
    parts.append(f"<title>STIB 编辑方案 — {_esc(t['gene'])}</title>")
    parts.append(
        "<style>"
        "body{margin:0;padding:32px;background:#f6f8fa;color:#1f2328;"
        "font-family:system-ui,-apple-system,'Segoe UI',sans-serif;line-height:1.6}"
        ".wrap{max-width:1000px;margin:0 auto}"
        "h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 10px;"
        "border-bottom:1px solid #d0d7de;padding-bottom:6px}"
        ".card{background:#fff;border:1px solid #d0d7de;border-radius:8px;padding:18px;"
        "margin-bottom:16px;overflow-x:auto}"
        "table{border-collapse:collapse;width:100%;font-size:13px}"
        "th,td{border:1px solid #d0d7de;padding:6px 8px;text-align:left}"
        "th{background:#f6f8fa;font-weight:600}"
        "code{background:#f6f8fa;padding:1px 5px;border-radius:4px;font-size:12px}"
        ".warn{border-left:4px solid #dc2626;background:#fef2f2;padding:12px 14px;"
        "border-radius:6px;margin-bottom:10px}"
        ".warn.med{border-left-color:#f59e0b;background:#fffbeb}"
        ".warn.low{border-left-color:#eab308;background:#fefce8}"
        ".muted{color:#6b7280;font-size:12px}"
        ".pill{display:inline-block;padding:1px 8px;border-radius:10px;font-size:11px;"
        "border:1px solid #d0d7de;background:#f6f8fa}"
        "</style></head><body><div class='wrap'>")

    parts.append(f"<h1>STIB 编辑方案 — {_esc(t['gene'])}</h1>")
    parts.append(f"<p class='muted'>{_esc(t['species'])}　｜　"
                 f"工具版本 {_esc(plan.meta.get('version',''))}</p>")

    # 摘要
    atm = t.get("atmosphere") or {}
    parts.append("<div class='card'><table>")
    rows = [
        ("推荐 spacer", f"<code>{_esc(t['selected_spacer'])}</code>"),
        ("PAM", f"<code>{_esc(t['selected_pam'])}</code> @ {t['selected_pam_start']:,}"),
        ("预测插入位点", f"<b>{t['predicted_insert_site']:,}</b>"),
        ("接合气氛", _esc(atm.get("label", ""))),
        ("构建", _esc(t.get("construct", "psh43 系列"))),
        ("氯霉素", _esc(t["cm_conc"])),
        ("候选 扫描/通过/淘汰",
         f"{t['n_candidates_scanned']} / {t['n_passed']} / {t['n_rejected']}"),
    ]
    for k, v in rows:
        parts.append(f"<tr><th style='width:180px'>{_esc(k)}</th><td>{v}</td></tr>")
    parts.append("</table></div>")

    # 图：位点分布（放大 + 全基因组）
    parts.append("<h2>① 候选在目标基因上的位置（放大）</h2>")
    parts.append("<div class='card'>")
    parts.append(locus_svg(plan.candidates, gene_start, gene_end, genome_len,
                           window=max(500, (gene_end - gene_start))))
    parts.append("</div>")

    parts.append("<h2>② 候选综合排序分</h2>")
    parts.append("<div class='card'>")
    parts.append(scores_svg(plan.candidates))
    parts.append("</div>")

    # 候选表
    parts.append("<h2>③ 候选明细</h2>")
    parts.append("<div class='card'><table>")
    parts.append("<tr><th>#</th><th>spacer</th><th>PAM</th><th>优选</th><th>落点</th>"
                 "<th>GC</th><th>序列偏好</th><th>风险分</th><th>非自身近同源</th></tr>")
    for i, c in enumerate(plan.candidates, 1):
        sf = c.seq_features.as_dict() if c.seq_features else {}
        bits = []
        if sf.get("tsd_preferred_GCWGC"):
            bits.append("GCWGC")
        if sf.get("at_pm5_ok"):
            bits.append("±5 A/T")
        pam_s = {1: "① GTT", 2: "② GTA"}.get(c.pam_tier, "③ 其他")
        place_s = {"strict": "★严格区", "central80": "中央80%",
                   "in_gene": "基因端部", "outside": "基因外"}.get(c.placement, c.placement)
        risk_s = "✅ 0" if c.is_self_only else f"{c.risk_score:.2f}"
        other = "✅ 0" if c.n_other == 0 else str(c.n_other)
        parts.append(
            f"<tr><td>{i}{' ★' if i == 1 else ''}</td>"
            f"<td><code>{_esc(c.sgrna.spacer)}</code></td>"
            f"<td><code>{_esc(c.sgrna.pam)}</code></td><td>{pam_s}</td>"
            f"<td>{place_s}</td><td>{c.sgrna.gc:.0%}</td>"
            f"<td>{'＋'.join(bits) if bits else '—'}</td>"
            f"<td>{risk_s}</td><td>{other}</td></tr>")
    parts.append("</table></div>")

    # 高等级警示
    highs = [w for w in plan.warnings if w.severity == "high"]
    if highs:
        parts.append("<h2>④ 高等级警示</h2><div class='card'>")
        for w in highs:
            parts.append(f"<div class='warn'><b>{_esc(w.key)}</b><br>"
                         + "<br>".join(_esc(s) for s in w.message.split(" | "))
                         + (f"<br><span class='muted'>建议：{_esc(w.discriminator)}</span>"
                            if w.discriminator else "")
                         + "</div>")
        parts.append("</div>")

    # 寡核苷酸
    parts.append("<h2>⑤ 寡核苷酸订单</h2><div class='card'><table>")
    parts.append("<tr><th>名称</th><th>序列 (5'→3')</th><th>说明</th></tr>")
    for o in plan.oligos:
        parts.append(f"<tr><td>{_esc(o['name'])}</td><td><code>{_esc(o['seq'])}</code></td>"
                     f"<td class='muted'>{_esc(o['note'])}</td></tr>")
    parts.append("</table></div>")

    # 脱靶明细（首选）
    if plan.candidates and plan.candidates[0].offtarget:
        ot = plan.candidates[0].offtarget
        other = [h for h in ot.hits if h.n_mismatch > 0]
        parts.append("<h2>⑥ 首选候选的脱靶情况</h2><div class='card'>")
        if not other:
            parts.append("<p>✅ <b>除自身位点外无任何近同源</b>（这是最好的情况）。</p>"
                         "<p class='muted'>任何 spacer 都必然精确命中它自己，"
                         "这一条不计入风险。6.26 Mb 基因组上 ≤4 错配的随机期望命中数约 0.003。</p>")
        else:
            parts.append("<table><tr><th>位置</th><th>链</th><th>错配数</th>"
                         "<th>错配位置分布</th><th>风险</th></tr>")
            for h in other[:20]:
                dist = ", ".join(f"{k}:{len(v)}" for k, v in h.by_region.items() if v)
                parts.append(f"<tr><td>{h.pos:,}</td><td>{h.strand}</td>"
                             f"<td>{h.n_mismatch}</td><td>{_esc(dist)}</td>"
                             f"<td>{_esc(h.risk)}</td></tr>")
            parts.append("</table>")
        parts.append("</div>")

    # 批量视图
    if batch:
        parts.append("<h2>⑦ 批量候选总览（全基因组）</h2><div class='card'>")
        parts.append(locus_svg(batch, gene_start, gene_end, genome_len, window=None))
        parts.append("<table><tr><th>#</th><th>spacer</th><th>PAM</th>"
                     "<th>插入位点</th><th>GC</th><th>GCWGC</th><th>±5 A/T</th></tr>")
        for i, b in enumerate(batch[:60], 1):
            sf = b.seq_features.as_dict() if b.seq_features else {}
            parts.append(f"<tr><td>{i}</td><td><code>{_esc(b.sgrna.spacer)}</code></td>"
                         f"<td><code>{_esc(b.sgrna.pam)}</code></td>"
                         f"<td>{b.sgrna.insert_site:,}</td><td>{b.sgrna.gc:.0%}</td>"
                         f"<td>{'✓' if sf.get('tsd_preferred_GCWGC') else '—'}</td>"
                         f"<td>{'✓' if sf.get('at_pm5_ok') else '—'}</td></tr>")
        parts.append("</table></div>")

    # 能力边界
    parts.append("<h2>⑧ 能力边界</h2><div class='card'><ul>")
    for line in [
        "<b>不预测编辑效率</b> —— 已发表数据点太少（约 11 条），抽不出具统计意义的模型。",
        "脱靶评估为 <b>6-mer 块 ×3 + 每块 1 替换枚举</b>：≤4 错配保证不漏检（已用暴力穷举交叉验证）；<b>不等价于 Cas-OFFinder 穷举</b>，>4 时保证失效。",
        "序列偏好（GCWGC / ±5 A/T）来自 ShCAST 机制研究，"
        "<b>AT 富集区相对整合位点的精确偏移尚未核实</b>，故只做定性偏好、不参与硬过滤。",
        "未设计 cargo 内部序列与骨架引物（需用户提供序列）。",
        "<b>未做湿实验</b> —— 本报告的正确性上限是「忠实于已发表方法与记录库」。",
    ]:
        parts.append(f"<li>{line}</li>")
    parts.append("</ul></div>")

    parts.append(f"<p class='muted'>生成于 {_esc(plan.meta.get('tool',''))} "
                 f"{_esc(plan.meta.get('version',''))}</p>")
    parts.append("</div></body></html>")
    return "".join(parts)
