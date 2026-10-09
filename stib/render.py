"""
STIB Planner — 渲染层。

把 Plan 渲染成人类可读的 Markdown 方案（7 节），
或机器可读的 JSON（供 DSH 插件的 tool 消费）。
"""

from __future__ import annotations

import json
from datetime import datetime

from .planner import Plan, Warning, RISK_THRESHOLD, MAX_BURDEN
from .seqfeatures import INSERT_OFFSET_DEFAULT


SEV_ICON = {"high": "🔴", "medium": "🟠", "low": "🟡"}
KIND_LABEL = {
    "open_question": "未收敛变量",
    "contradiction": "材料矛盾",
    "risk": "风险",
    "limitation": "能力边界",
}


def render_markdown(plan: Plan, title: str | None = None) -> str:
    t = plan.target
    L: list[str] = []

    L.append(f"# STIB 编辑方案 — {t['gene']}")
    L.append("")
    L.append(f"**目标物种**：{t['species']}　　**生成时间**：{datetime.now():%Y-%m-%d %H:%M}")
    L.append("")
    L.append(f"> {plan.meta['disclaimer']}")
    L.append("")

    # ---- 0. 警示前置 ----
    high = [w for w in plan.warnings if w.severity == "high"]
    if high:
        L.append("## ⚠️ 开始前必读（高等级警示）")
        L.append("")
        for w in high:
            L.append(f"### {SEV_ICON[w.severity]} [{KIND_LABEL.get(w.kind, w.kind)}] {w.key}")
            L.append("")
            for seg in w.message.split(" | "):
                L.append(f"- {seg}")
            L.append("")
            if w.discriminator:
                L.append(f"**建议动作 / 判别实验**：{w.discriminator}")
                L.append("")

    # ---- 1. 目标摘要 ----
    L.append("## 1. 目标摘要")
    L.append("")
    L.append("| 项 | 值 |")
    L.append("|---|---|")
    L.append(f"| 目标基因 | `{t['gene']}` |")
    L.append(f"| 基因组区间 | {t['span'][0]:,} – {t['span'][1]:,}（{t['span_len']:,} bp） |")
    L.append(f"| 选定 spacer | `{t['selected_spacer']}` |")
    L.append(f"| 对应 PAM | `{t['selected_pam']}` @ {t['selected_pam_start']:,} |")
    L.append(f"| **预测插入位点** | **{t['predicted_insert_site']:,}**"
             f"（按名义偏移 PAM 3' 端下游 {INSERT_OFFSET_DEFAULT} bp 计算；"
             f"实测落在 **{t['insert_window'][0]}–{t['insert_window'][1]} bp** 区间内均可接受） |")
    L.append(f"| 氯霉素筛选浓度 | {t['cm_conc']} |")
    atm = t.get("atmosphere") or {}
    if atm:
        conf = {"high": "✅", "medium": "🟠", "low": "🟡"}.get(atm.get("confidence", ""), "")
        L.append(f"| **接合气氛** | **{atm.get('label', '')}** {conf} — {atm.get('reason', '')} |")
    L.append(f"| 构建 | {t.get('construct', 'psh43 系列')} |")
    L.append(f"| 候选 扫描/通过/淘汰 | {t['n_candidates_scanned']} / "
             f"{t['n_passed']} / {t['n_rejected']} |")
    L.append(f"| 脱靶阈值 | ≤{t['max_mismatch']} 错配"
             f"（{'保证不漏检' if t['offtarget_guaranteed'] else '**可能漏检**'}） |")
    L.append("")
    w0, w1 = t["insert_window"]
    a0, a1 = t["insert_window_alt"]
    L.append(f"> 插入窗口论文两处写法略有差异（Figure 2 为 {w0}–{w1} bp，"
             f"Figure 4 菌群为 {a0}–{a1} bp）。**一定范围内的偏移可接受**；"
             "仅方向敏感应用（启动子捕获、框内插入）需精确确认。")
    L.append("")
    L.append("**排序依据是编辑位点落在基因内 + GC 趋于居中 + 脱靶风险，"
             "不是效率预测。** 效率须实验确认，本工具不做效率预测。")
    L.append("")

    # ---- 2. sgRNA 候选（含脱靶评估）----
    L.append("## 2. sgRNA 候选与脱靶评估")
    L.append("")
    L.append("### 2.1 通过筛选的候选")
    L.append("")
    L.append("| # | spacer (5'→3') | PAM | 优选 | 落点 | GC | 序列偏好 | 风险分 | 非自身近同源 |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for i, c in enumerate(plan.candidates, 1):
        s = c.sgrna
        star = " ★" if i == 1 else ""
        risk = c.risk_score
        if c.is_self_only:
            risk_s = "✅ 0"
        elif risk <= 1.0:
            risk_s = f"🟡 {risk:.2f}"
        elif risk <= 3.0:
            risk_s = f"🟠 {risk:.2f}"
        else:
            risk_s = f"🔴 {risk:.2f}"
        other_s = "✅ 0（仅自身）" if c.n_other == 0 else f"{c.n_other} 个"
        pam_s = {1: "① GTT", 2: "② GTA"}.get(c.pam_tier, "③ 其他")
        place_s = {"strict": "★严格区", "central80": "中段80%",
                   "in_gene": "基因端部5%", "outside": "基因外"}.get(c.placement, c.placement)
        sf = c.seq_features.as_dict() if c.seq_features else {}
        bits = []
        if sf.get("tsd_preferred_GCWGC"):
            bits.append("GCWGC")
        if sf.get("at_pm5_ok"):
            bits.append("±5 A/T")
        seq_s = "＋".join(bits) if bits else "—"
        L.append(f"| {i}{star} | `{s.spacer}` | `{s.pam}` | {pam_s} | {place_s} | "
                 f"{s.gc:.0%} | {seq_s} | {risk_s} | {other_s} |")
    L.append("")
    L.append("★ = 推荐首选。**落点两级制**：「★严格区」= 插入位点落在 CDS **前 1/3 – 2/3**"
             "（论文做法，优先满足）；严格区不满足时回落到「中段80%」"
             "（**排除 CDS 前 5% 与后 5%**）的普通约束。")
    L.append("")
    L.append("**PAM 优选**：① GTT（课题组经验最多）→ ② GTA → ③ 其他 GTN。")
    L.append("")
    L.append("**序列偏好**：ShCAST 的整合位点不只由 sgRNA 决定，DNA 序列本身也有偏好 —— "
             "TSD（5 bp）内含 **GCWGC**（TnsB 偏好的插入点）、TSD 中心 **±5 bp** 偏好 **(A/T)**。"
             "这两项参与排序，但不作硬性要求。")
    L.append("")
    L.append("**风险分为加权值**：每个非自身近同源位点贡献 "
             "clamp(1 + Σ位置权重, 0, " + str(MAX_BURDEN) + ")；"
             "**helix ×3**（PAM 近端前 17 nt 关键区）、tail ×1（3' 端外周 6 nt）。"
             "**任何错配都只增不减** —— 近同源位点无论错配在哪一段都是可能被转座的候选，"
             "不存在「减害」。\n"
             f"**分数越小越好**（0 = 最优），超过阈值 {RISK_THRESHOLD} 的候选已被淘汰。")
    L.append("")
    L.append("> **注意**：候选在基因组中**必然精确命中它自己**，这一条不计入风险。"
             "因此「非自身近同源 0 个」是**最好的情况**，表示除自身位点外无近同源。")
    L.append("")

    # 脱靶引擎的能力边界
    t = plan.target
    L.append("### 2.2 脱靶评估说明（能力边界，务必阅读）")
    L.append("")
    L.append(f"- **引擎**：6-mer 块 ×3（互不重叠）+ 每块 1 处替换枚举，全基因组正负链枚举，阈值 **≤{t['max_mismatch']} 错配**")
    if t["offtarget_guaranteed"]:
        L.append("- **保证**：≤4 错配**不漏检**（3 块 × 每块 1 替换 → 需 6 个错配才能弄脏全部块，故 ≤5 必有命中块；本工具保守声明 ≤4）。**已用暴力穷举交叉验证**")
    else:
        L.append(f"- 🔴 **警告**：阈值 {t['max_mismatch']} > 4，**鸽巢保证失效，可能漏检**。"
                 "如需完整 ≤5 错配枚举，请用 Cas-OFFinder v2.4 交叉核对。")
    L.append("- **不等价声明**：本引擎是种子剪枝，与 Cas-OFFinder 的穷举**不等价**，"
             "可能少报。论文与 SOP 用的是 Cas-OFFinder ≤5 错配。")
    L.append("- **不建模**：脱靶位点是否真会发生转座，还取决于局部染色质状态、"
             "靶免疫（target immunity）等本引擎不涉及的因素。")
    L.append("")
    L.append("### 2.3 风险分级规则（按错配位置，CAST 特有）")
    L.append("")
    L.append("按 **helix（PAM 近端前 17 nt）/ tail（3' 端外周 6 nt）** 分区。\n"
             "依据结构研究：ShCAST 的 PAM 在 spacer **上游**，Cas12k 与 RNA:DNA "
             "异源双链的接触集中在 **前 17 nt**，3' 端 6 nt 为延伸尾部。\n"
             "因此同样错配数下，**错配落在前 17 nt 的位点更危险**。\n"
             "（v0.6 更正 —— 此前误套 Cas9 的 3' 端种子区惯例，方向相反）")
    L.append("")
    L.append("| 等级 | 规则 |\n")
    L.append("|---|---|\n")
    # fallback：plan.meta["risk_rules"] 总会存在（render_json 也依赖它），
    # 这里是防御性默认值 —— 必须与 offtarget.RISK_RULES 保持一致，
    # 否则万一 meta 缺失会打印出过时的规则。
    from .offtarget import RISK_RULES as _RR
    for rule in plan.meta.get("risk_rules", []) or [
        {"level": lv, "rule": r} for lv, r in _RR
    ]:
        L.append(f"| {rule['level']} | {rule['rule']} |\n")
    L.append("\n")
    L.append("**出现高风险脱靶的候选已被自动排除**，不进入上表。\n")
    L.append("\n")

    # 被淘汰的候选（透明化）
    if plan.rejected:
        L.append("### 2.4 被淘汰的候选（及原因）")
        L.append("")
        L.append("| spacer | 淘汰原因 |")
        L.append("|---|---|")
        for c in plan.rejected[:15]:
            L.append(f"| `{c.sgrna.spacer}` | {'; '.join(c.reject_reasons)} |")
        if len(plan.rejected) > 15:
            L.append(f"| *…另有 {len(plan.rejected) - 15} 条* | |")
        L.append("")

    # 首选候选的脱靶明细
    if plan.candidates and plan.candidates[0].offtarget:
        ot = plan.candidates[0].offtarget
        self_hits = [h for h in ot.hits if h.n_mismatch == 0]
        other = [h for h in ot.hits if h.n_mismatch > 0]
        L.append("### 2.5 首选候选的脱靶明细")
        L.append("")
        n_copy = sum(1 for h in self_hits if getattr(h, "is_identical_copy", False))
        L.append(f"- 全基因组命中总数（≤{ot.max_mismatch} 错配）：**{ot.n_hits}**"
                 f"（**精确命中 {len(self_hits)} 处**；**非自身近同源 {len(other)} 处**）")
        if len(self_hits) > 1:
            L.append(f"- ⚠️ **该位点存在 {len(self_hits)} 份完全相同的拷贝**"
                     "（0 错配）—— 这是**串联重复区**的特征（如 rRNA 多拷贝、"
                     "重复基因家族）。这些拷贝不计入脱靶风险，"
                     "但**插入可能落在任一份拷贝上**，鉴定时必须覆盖全部拷贝位点。")
            L.append("")
            L.append("| 完全相同拷贝位置 | 链 |")
            L.append("|---|---|")
            for h in self_hits[:12]:
                L.append(f"| {h.pos:,} | {h.strand} |")
            if len(self_hits) > 12:
                L.append(f"| *…另有 {len(self_hits) - 12} 处* | |")
            L.append("")
        if other:
            L.append("")
            L.append("| 位置 | 链 | 错配数 | 错配位置分布 | 风险 |")
            L.append("|---|---|---|---|---|")
            for h in other[:10]:
                dist = ", ".join(f"{k}:{len(v)}" for k, v in h.by_region.items() if v)
                L.append(f"| {h.pos:,} | {h.strand} | {h.n_mismatch} | {dist} | {h.risk} |")
        else:
            L.append("")
            L.append("✅ **除自身位点外无任何近同源** —— 这是最好的情况。")
            L.append("")
            L.append("> 任何 spacer 都**必然精确命中它自己**（0 错配），这一条不计入风险。"
                     "在 6.26 Mb 基因组上，≤4 错配的随机期望命中数约 0.003，"
                     "因此「仅自身命中」是最常见的结果，属正常。")
        L.append("")
        for w in ot.warnings:
            L.append(f"> ⚠️ {w}")
        L.append("")


    # ---- 3. 寡核苷酸订单 ----
    L.append("## 3. 寡核苷酸订单（可直接下单）")
    L.append("")
    L.append("| 名称 | 序列 (5'→3') | 说明 |")
    L.append("|---|---|---|")
    for o in plan.oligos:
        L.append(f"| {o['name']} | `{o['seq']}` | {o['note']} |")
    L.append("")
    L.append("**合成要求**：PAGE 纯化，母液 100 µM，工作液 10 µM。")
    L.append("")

    # ---- 4. 验证引物 ----
    L.append("## 4. 验证引物")
    L.append("")
    p = plan.primers
    L.append("| 名称 | 序列 (5'→3') | 长度 | Tm≈ | GC | 位置 | 用途 |")
    L.append("|---|---|---|---|---|---|---|")
    order = ["LE_up_F", "LE_cargo_R", "RE_cargo_F", "RE_dn_R", "backbone_F", "backbone_R",
             "ctrl_16S_F", "ctrl_16S_R"]
    for k in order:
        pr = p.get(k)
        if pr is None:
            L.append(f"| {k} | *（需提供 cargo / 骨架序列后自动设计）* | — | — | — | — | — |")
        else:
            pos = "基因组" if pr.start >= 0 else "通用"
            L.append(f"| {pr.name} | `{pr.seq}` | {pr.length} | {pr.tm:.1f} | {pr.gc:.0%} | {pos} | {pr.note or k} |")
    L.append("")
    L.append("**PCR 组合方案**（出处：SOP 表 8）")
    L.append("")
    L.append("| 检测 | 引物组合 | 预期结果 |")
    L.append("|---|---|---|")
    L.append("| LE junction | `LE_up_F` + `LE_cargo_R` | 单一目的条带 |")
    L.append("| RE junction | `RE_cargo_F` + `RE_dn_R` | 单一目的条带 |")
    L.append("| cargo | cargo 内部引物 | 单一目的条带 |")
    L.append("| **骨架共整合** | `backbone_F` + `backbone_R` | **无条带（阴性）** |")
    L.append("| 物种鉴定 | `ctrl_16S_F` + `ctrl_16S_R` | ~1.5 kb 条带 |")
    L.append("")
    L.append("> Tm 为**估算值**（≤13 nt 用 Wallace 规则，更长用 GC 校正公式），"
             "不可直接作为退火温度。首次使用请做梯度 PCR 确认。")
    L.append("")
    L.append("> `LE_cargo_R` / `RE_cargo_F` / 骨架引物需要 cargo 与骨架序列才能设计 —— "
             "请在提供 cargo 序列后调用 `plan` 重跑，或手工按 SOP 表 5 的 GA725/GA778 体系设计。")
    L.append("")

    # ---- 5. 构建与接合流程 ----
    L.append("## 5. 构建与接合流程")
    L.append("")
    L.append("> 流程细节以课题组既有 SOP 为准（`039_STIB 使用protocol.md`）。下表为**检查清单**，"
             "每步附「应该看到什么」。")
    L.append("")
    cur = None
    for s in plan.build_steps:
        if s.phase != cur:
            cur = s.phase
            L.append(f"### 阶段 {cur}")
            L.append("")
            L.append("| 时间 | 操作 | 应该看到 | 出处 |")
            L.append("|---|---|---|---|")
        L.append(f"| {s.day} | {s.action} | {s.readout} | {s.ref} |")
    L.append("")
    L.append("**时间预估**：单菌编辑自收到寡核苷酸起约 **1 周**；菌群编辑约 **2–3 周**。")
    L.append("")

    # ---- 6. 判定点与停止规则 ----
    L.append("## 6. 判定点与停止规则")
    L.append("")
    L.append("> 这一节是 SOP 没有的部分：SOP 给的是「可能原因」并列列表，"
             "这里给的是**症状 → 假设 → 最便宜的判别实验**。")
    L.append("")
    for i, d in enumerate(plan.decisions, 1):
        L.append(f"### 6.{i} 观察到：{d.observe}")
        L.append("")
        L.append(f"- **时点**：{d.after}")
        L.append(f"- **符合预期则**：{d.if_ok}")
        L.append(f"- **否则**：{d.if_not}")
        if d.hypothesis:
            L.append(f"- **假设**：{d.hypothesis}")
        if d.discriminator:
            L.append(f"- **判别实验**：{d.discriminator}")
        L.append("")

    # ---- 7. 风险清单 ----
    L.append("## 7. 风险清单")
    L.append("")
    for r in sorted(plan.risks, key=lambda x: {"high": 0, "medium": 1, "low": 2}[x["level"]]):
        L.append(f"### {SEV_ICON[r['level']]} {r['id']} · {r['title']}")
        L.append("")
        L.append(f"- **机制/表现**：{r['detail']}")
        L.append(f"- **如何发现**：{r['detect']}")
        L.append(f"- **如何处置**：{r['action']}")
        if r.get("note"):
            L.append(f"- **备注**：{r['note']}")
        L.append("")

    # ---- 8. 回退面板（新手最需要的：一条确定能用的）----
    if plan.fallbacks:
        L.append("## 8. 回退方案：已验证的 spacer 面板")
        L.append("")
        L.append("> **如果首选位点失败，不要从零开始设计。** 下表是论文与记录库中"
                 "**已实测**的 spacer。新手需要的不是「最优」的，而是**一条确定能用的**。")
        L.append("")
        L.append("| 名称 | 位点 | 适用物种 | spacer | 实测结果 | 状态 |")
        L.append("|---|---|---|---|---|---|")
        for d in plan.fallbacks:
            icon = {"success": "✅", "weak": "🟠", "failed": "🔴"}.get(d["outcome"], "•")
            L.append(f"| {d['name']} | {d['site']} | {d['species']} | `{d['spacer']}` | "
                     f"{d['evidence']} | {icon} |")
        L.append("")
        L.append("**图例**：✅ 已验证可用 / 🟠 可用但偏弱 / 🔴 **已知失败，避开**")
        L.append("")
        L.append("> ⚠️ 这是**经验先验，不是效率预测**。样本量小，仅供回退参考。")
        L.append("> 🔴 的条目尤其重要 —— 它们能让你避免重复别人已经踩过的坑。")
        L.append("")

    # ---- 9. 剩余警示 ----
    rest = [w for w in plan.warnings if w.severity != "high"]
    if rest:
        L.append("## 9. 其他警示")
        L.append("")
        for w in rest:
            L.append(f"- {SEV_ICON[w.severity]} **[{KIND_LABEL.get(w.kind, w.kind)}] {w.key}**："
                     + " | ".join(w.message.split(" | ")))
        L.append("")

    # ---- 附：本方案未做的事 ----
    L.append("---")
    L.append("")
    L.append("## 附：本方案**未**做的事（能力边界）")
    L.append("")
    L.append("1. **未预测编辑效率** —— 已发表数据点太少（约 11 条带效率的位点），"
             "抽不出具统计意义的评估标准。硬做只会产出假精度。")
    L.append("2. **脱靶评估不等价于 Cas-OFFinder 穷举** —— 本引擎是 6-mer 块 + 替换枚举，"
             "≤4 错配不漏检（已交叉验证），但 >4 时保证失效（引擎会警示）。"
             "论文用 Cas-OFFinder ≤5 错配。")
    L.append("3. **未设计 cargo 内部序列** —— 需你提供目的片段。")
    L.append("4. **未设计骨架 / LE / RE 引物** —— 需你提供骨架与 cargo 序列。")
    L.append("5. **未验证引物二聚体/发夹结构** —— Tm 为估算值，须梯度 PCR 确认。")
    L.append("6. **未做湿实验** —— 本方案的正确性上限是「忠实于已发表方法与记录库」，"
             "不等于在你的实验室、你的菌株上会成功。")
    L.append("")
    L.append(f"*工具版本 {plan.meta['version']}*")
    L.append("")

    return "\n".join(L)


def render_json(plan: Plan) -> str:
    return json.dumps(plan.as_dict(), ensure_ascii=False, indent=2)
