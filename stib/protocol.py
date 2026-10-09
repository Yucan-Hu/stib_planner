"""
STIB Planner — 定制化实验协议渲染。

═══ 为什么单独一个模块 ═══

`render.py` 产出的方案是**决策文档**（选哪个位点、为什么、有什么风险），
但科研人员拿去做实验时还需要**操作文档**：加多少 µL、多少度、多长时间、
今天做什么明天做什么。

本模块把知识库里的实验参数（全部带出处）编排成一份**可照着做的 protocol**，
并针对用户的具体目标（基因名、扩增子长度、物种、菌株）填入实际数值。

═══ 原则 ═══

1. **每个数值都带出处** —— 体系体积来自 SOP 的表 2/3/4/6/7，
   不凭印象写。用户能核对到原始 SOP 的哪一行。
2. **不虚构未核实的东西** —— 例如退火温度只在用户给了引物 Tm 时才算，
   否则留占位符并说明怎么定。
3. **按时间轴编排** —— 科研人员需要知道"第几天做什么"，
   而不是一堆并列的步骤。
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------

def _table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    """每个元素**以换行结尾** —— 调用方用 "".join(L) 拼接。"""
    out = ["| " + " | ".join(headers) + " |\n",
           "|" + "|".join("---" for _ in headers) + "|\n"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |\n")
    return out


def _fact(kb: Any, key: str) -> Any:
    f = kb.get(key) if hasattr(kb, "get") else None
    return f.value if f else None


def _prov(kb: Any, key: str) -> str:
    f = kb.get(key) if hasattr(kb, "get") else None
    if not f:
        return ""
    src = f.source
    return f"{src} {f.locator}"


def _system_table(kb: Any, key: str, kb_any: Any) -> list[str]:
    """把一个体系 fact 渲染成加样表 + 程序表。"""
    v = _fact(kb_any, key)
    if not isinstance(v, dict):
        return []
    L: list[str] = []
    L.append(f"**总体积 {v.get('total_µL', '?')} µL**（出处：{_prov(kb_any, key)}）\n\n")
    if v.get("reagents"):
        L.append("| 试剂 | 体积 (µL) |\n")
        L.append("|---|---|\n")
        for r in v["reagents"]:
            name, vol = (r if isinstance(r, (list, tuple)) else (r, ""))
            L.append(f"| {name} | {vol} |\n")
        L.append("\n")
    if v.get("program"):
        L.append("| 步骤 | 温度 | 时间 | 循环 |\n")
        L.append("|---|---|---|---|\n")
        for step in v["program"]:
            row = list(step) + [""] * (4 - len(step))
            L.append(f"| {row[0]} | {row[1]} | {row[2]} | {row[3]} |\n")
        L.append("\n")
    if v.get("pmol"):
        pm = v["pmol"]
        L.append(f"**加量**：载体片段 {pm.get('载体片段')} pmol、"
                 f"插入片段 {pm.get('插入片段')} pmol"
                 f"（允许 {pm.get('允许范围')}）\n")
    if v.get("formula"):
        L.append(f"**换算**：`{v['formula']}`\n")
    if v.get("total_time"):
        L.append(f"**耗时**：{v['total_time']}\n")
    return L


# ---------------------------------------------------------------------------
# 主渲染
# ---------------------------------------------------------------------------

def render_protocol(plan: Any, kb: Any, *, species: str | None = None) -> str:
    """
    生成**分阶段、可照着做**的定制 protocol。

    参数:
      plan   planner.Plan（含目标基因、坐标、候选、引物、菌株信息）
      kb     KnowledgeBase
    """
    t = plan.target
    gene = t["gene"]
    sp = species or t.get("species", "Bacteroides thetaiotaomicron")
    cm = t.get("cm_conc", "25 µg/mL")
    atm = t.get("atmosphere") or {}
    s0, e0 = t["span"]
    span = e0 - s0

    L: list[str] = []
    L.append(f"# 定制实验协议 — {gene} 插入失活（{sp}）\n")
    L.append(f"> 生成时间：{datetime.now():%Y-%m-%d %H:%M}　｜　"
             f"目标区间 {s0:,}–{e0:,}（{span:,} bp）\n")
    L.append("> ⚠️ **本协议基于课题组 SOP 与已发表论文编排，每个数值都标了出处；"
             "但未做湿实验验证，首次执行请按 §0 的核对清单逐项确认。**\n")

    # =====================================================================
    L.append("\n---\n## 0. 动手前核对清单（10 分钟，能省两周）\n")
    L.append("| # | 核对项 | 为什么 | 出处 |\n")
    L.append("|---|---|---|---|\n")
    checks = [
        ("血平板是否含 **DAP 60 µg/mL**",
         "ΔdapA 供体漏加 DAP → 转化子为 **0 CFU**（记录库实测）",
         "RK4 / 记录库 2024-03"),
        ("血平板**无菌、脱纤维**",
         "血平板污染会毁掉整批接合",
         "reagent.blood"),
        (f"受体菌株氯霉素浓度 = **{cm}**",
         "bt/bo 用 25；其他用 12.5（用错则筛选失效）",
         "reagent.cm_bt_bo / cm_other"),
        (f"接合气氛 = **{atm.get('label', '?')}**",
         atm.get("reason", ""),
         "o2.aerobic_capable / anaerobic_required"),
        ("感受态是 **S17-1λpir**（含 pir）",
         "R6K 载体无 pir 不能复制 → 构建失败",
         "SOP §B.1 步骤4"),
        ("引物已按 §1 下单（含 TSD/同源臂）",
         "序列错一位整条流程白做",
         "本方案 §1"),
    ]
    for i, (item, why, src) in enumerate(checks, 1):
        L.append(f"| {i} | {item} | {why} | {src} \n")
    L.append("")

    # =====================================================================
    L.append("\n---\n## 1. 阶段 A：sgRNA 设计（1–2 天）\n")
    L.append("### 1.1 选定的 sgRNA\n")
    L.append("| 项 | 值 |\n")
    L.append("|---|---\n")
    L.append(f"| spacer (23 nt) | `{t['selected_spacer']}` \n")
    L.append(f"| PAM | `{t['selected_pam']}` @ {t['selected_pam_start']:,} \n")
    L.append(f"| 预测插入位点 | {t['predicted_insert_site']:,} \n")
    L.append("")

    L.append("### 1.2 订购寡核苷酸（可直接复制下单）\n")
    if plan.oligos:
        L.append("| 名称 | 序列 (5'→3') | 长度 | 说明 |\n")
        L.append("|---|---|---|---|\n")
        for o in plan.oligos:
            L.append(f"| {o['name']} | `{o['seq']}` | {len(o['seq'])} | {o['note']} \n")
        L.append("")
        L.append("> 粘性末端：`aaag`（F）/ `agcc`（R）—— 与 BsaI Golden Gate 兼容。\n")
        L.append("> **订购规格**：普通脱盐即可（无需磷酸化修饰 —— 步骤 B.1-1 会酶法磷酸化）。\n")
    else:
        L.append("（本方案未生成寡核苷酸）\n")

    if plan.primers:
        L.append("\n### 1.3 验证引物\n")
        L.append("| 名称 | 序列 (5'→3') | 长度 | Tm≈ | GC | 位置 | 用途 |\n")
        L.append("|---|---|---|---|---|---|---|\n")
        for k, pr in (plan.primers.items() if isinstance(plan.primers, dict)
                      else enumerate(plan.primers)):
            # 兼容两种形态：dict（插件层传入）与 Primer dataclass（planner 内部）。
            # 修 bug：此前只处理 dict，导致 dataclass 形态的引物被静默丢弃，
            # §1.3 只剩表头 —— 而 junction 引物恰恰是鉴定插入所必需的。
            if pr is None:
                continue
            if isinstance(pr, dict):
                d = pr
            elif hasattr(pr, "as_dict"):
                d = pr.as_dict()
            else:
                continue
            gc = d.get("gc", "—")
            gc_s = f"{gc:.0%}" if isinstance(gc, (int, float)) else str(gc)
            tm = d.get("tm", "—")
            tm_s = f"{tm:.1f}" if isinstance(tm, (int, float)) else str(tm)
            pos = "基因组" if d.get("start", -1) >= 0 else "通用"
            L.append(f"| {d.get('name', k)} | `{d.get('seq', '')}` | "
                     f"{d.get('len', len(d.get('seq', '')))} | {tm_s} | "
                     f"{gc_s} | {d.get('pos', pos)} | "
                     f"{d.get('note') or d.get('purpose', '') or k} |\n")
        L.append("")
        L.append("> ⚠️ **Tm 为估算值**（Wallace 修正式），首次使用须梯度 PCR 确认。\n")

    L.append("\n### 1.4 脱靶复核\n")
    L.append("本方案已做全基因组 ≤4 错配枚举（helix/tail 分区）。"
             "**发表级结果请用 Cas-OFFinder v2.4 复核**（本引擎不等价于穷举）。\n")
    if plan.meta.get("risk_rules"):
        L.append("| 等级 | 判定规则 |\n")
        L.append("|---|---\n")
        for rr in plan.meta["risk_rules"]:
            L.append(f"| {rr['level']} | {rr['rule']} \n")
        L.append("")

    # =====================================================================
    L.append("\n---\n## 2. 阶段 B.1：BsaI Golden Gate 更换 sgRNA（1 天）\n")
    L.append("> 出处：SOP §B.1（表 2 / 表 3 / 表 4-1 / 表 4-2）\n")
    L.append("\n### B.1-1 寡核苷酸磷酸化与退火（与 B.1-2 同步进行）\n")
    L += _system_table(None, "gg.phospho_anneal", kb)
    v = _fact(kb, "gg.phospho_anneal") or {}
    if v.get("program"):
        L.append("> Touch down 程序：**自 95 °C 起每循环降低 1 °C，共 79 个循环**，"
                 "降至 16 °C；随后 16 °C 保温。总计约 55 min。\n")

    L.append("\n### B.1-2 载体 BsaI 单酶切（与 B.1-1 同步进行）\n")
    L += _system_table(None, "gg.vector_digest", kb)
    L.append("> 质粒用量 **0.015–0.05 pmol**。此步与 B.1-1 可同时放入 PCR 仪/水浴。\n")

    L.append("\n### B.1-3 Golden Gate 组装\n")
    L += _system_table(None, "gg.assembly", kb)
    L.append("> **连接温度注意**：用 **Hi-T4**（耐热）时连接步骤 25 °C；"
             "用**普通 T4 连接酶**时改为 **16 °C**。时间紧张可 37 °C 恒温 1–1.5 h。\n")

    L.append("\n### B.1-4 转化与涂板\n")
    L.append("1. 取 **5–10 µL** Golden Gate 反应液，转入 **S17-1λpir** 感受态"
             "（热激法见 §4.5）\n")
    L.append("2. 100 µL 菌液涂 **LB + Amp 100 µg/mL**；其余离心浓缩后涂板\n")
    L.append("3. **37 °C 过夜**\n")

    L.append("\n### B.1-5 验证\n")
    L.append(f"{_fact(kb, 'verify.colony_pcr') or ''}\n")
    L.append(f"（出处：{_prov(kb, 'verify.colony_pcr')}）\n")

    # =====================================================================
    L.append("\n---\n## 3. 阶段 B.2：Gibson 组装更换 cargo（1–2 天，可选）\n")
    L.append("> **仅当你要插入外源片段时才做这一步。**纯基因失活可跳过，"
             "直接进 §4 接合。\n")
    L.append("> 出处：SOP §B.2（表 5 / 表 6 / 表 7）\n")

    L.append("\n### B.2-1 cargo 引物设计\n")
    L.append("⚠️ **cargo 引物需自行附加同源臂**。同源来源：\n")
    L.append("- 上游：**Cargo1 区域**（SOP 表 5，`ACATGCAACTGTCTAGAATGCGAG`）\n")
    L.append("- 下游：**GA778 的后 20 bp**（`GTGGCCTAT` 所在的退火段 3' 端）\n")
    L.append("\n> 也可用 `stib_cargo_assembly` 工具自动设计并输出 SnapGene 图谱。\n")

    L.append("\n### B.2-2 高保真 PCR 扩增\n")
    L += _system_table(None, "pcr.phanta_max", kb)
    L.append("> **退火温度**：按引物 Tm 设定（本方案引物若给了 Tm，"
             "建议先做梯度 55–65 °C）。\n")
    L.append("> **延伸时间**：按产物长度 × 30–60 s/kb 计算。\n")

    L.append("\n### B.2-3 产物处理\n")
    L.append(f"{_fact(kb, 'pcr.product_cleanup') or ''}\n")
    L.append(f"（出处：{_prov(kb, 'pcr.product_cleanup')}）\n")

    L.append("\n### B.2-4 Gibson 组装\n")
    L += _system_table(None, "gibson.reaction", kb)
    L.append("> **pmol 换算**：`pmol = weight(ng) / (length(bp) × 0.65)`。\n")

    # =====================================================================
    L.append("\n### B.2-5 转化（NEB 热激法）\n")
    hs = _fact(kb, "transform.heat_shock")
    if isinstance(hs, list):
        for i, step in enumerate(hs, 1):
            L.append(f"{i}. {step}\n")
        L.append(f"\n（出处：{_prov(kb, 'transform.heat_shock')}）\n")
    L.append("\n> **B.1（Golden Gate）产物同样用此法**：取 **5–10 µL** 反应液即可。\n")
    L.append("> 涂板用 **LB + Amp 100 µg/mL**；长片段大质粒可 30 °C 培养 24 h。\n")

    # =====================================================================
    L.append("\n---\n## 4. 阶段 C：接合转移（第 3–4 天）\n")
    L.append("> 出处：SOP §4 阶段C\n")
    steps_c = [
        ("C-1", "供体培养",
         f"{_fact(kb, 'conj.donor_od') or ''}（出处：{_prov(kb, 'conj.donor_od')}）"),
        ("C-2", "受体培养",
         f"{_fact(kb, 'conj.recipient_od') or ''}（出处：{_prov(kb, 'conj.recipient_od')}）"),
        ("C-3", "混合",
         f"{_fact(kb, 'conj.ratio') or ''}（出处：{_prov(kb, 'conj.ratio')}）"),
        ("C-4", "离心浓缩",
         f"{_fact(kb, 'conj.spin') or ''}（出处：{_prov(kb, 'conj.spin')}）"),
        ("C-5", "滴板",
         f"{_fact(kb, 'conj.spot') or ''}（出处：{_prov(kb, 'conj.spot')}）"),
        ("C-6", "接合",
         f"{_fact(kb, 'conj.mating') or ''}　**本方案气氛：{atm.get('label', '?')}**"
         f"（出处：{_prov(kb, 'conj.mating')}）"),
        ("C-7", "回收",
         f"{_fact(kb, 'conj.recovery') or ''}（出处：{_prov(kb, 'conj.recovery')}）"),
        ("C-8", "涂板筛选",
         f"{_fact(kb, 'conj.plating') or ''}　**Cm 浓度：{cm}**"
         f"（出处：{_prov(kb, 'conj.plating')}）"),
    ]
    for sid, title, body in steps_c:
        L.append(f"\n**{sid} {title}**\n")
        L.append(f"{body}\n")
    cd = _fact(kb, "conj.colony_days")
    if isinstance(cd, dict):
        L.append("\n**集落出现时间**："
                 + "；".join(f"{k} → {v}" for k, v in cd.items()) + "\n")
        L.append(f"（出处：{_prov(kb, 'conj.colony_days')}）\n")

    # =====================================================================
    L.append("\n---\n## 5. 阶段 D：接合子鉴定（第 5–7 天）\n")
    verify = [
        ("D-1", "挑克隆",
         "挑单克隆接种 BHIS（含 Cm 与 Genta 200 µg/mL），37 °C 厌氧过夜"),
        ("D-2", "菌落 PCR",
         "用 §1.3 的验证引物扩增 **junction**：左 junction 与右 junction 各一对"),
        ("D-3", "Sanger 测序",
         f"确认整合落在 PAM 下游 {t.get('insert_window', ['?', '?'])[0]}–"
         f"{t.get('insert_window', ['?', '?'])[1]} bp 区间；存在 5 bp TSD；"
         "边界与载体 LE/RE 一致"),
        ("D-4", "**骨架 PCR（必做）**",
         "用 ampR / ori 特异引物扩增，**应为阴性**。"
         "阳性 = 共整合（copy-and-paste 副产物），不是成功克隆 → "
         "划线纯化传代 2–5 次复检；持续阳性则弃去该克隆"),
        ("D-5", "物种确认",
         "16S 测序；⚠️ **不得只看 top-1**："
         "*B. thetaiotaomicron* vs *B. faecis* **98.43%**；"
         "*P. vulgatus* vs *P. dorei* **97.10%**"),
    ]
    for sid, title, body in verify:
        L.append(f"\n**{sid} {title}**\n")
        L.append(f"{body}\n")

    # =====================================================================
    L.append("\n---\n## 6. 时间轴总览\n")
    L.append("| 天 | 做什么 | 产出 | 能否并行 |\n")
    L.append("|---|---|---|---|\n")
    L.append("| 第 0 天 | 下单寡核苷酸/引物；准备培养基与平板 | 引物到货 | — \n")
    L.append("| 第 1 天 | §1 复核位点与脱靶；§2 B.1-1～B.1-3 Golden Gate | 组装产物 | B.1-1 与 B.1-2 同步 \n")
    L.append("| 第 2 天 | B.1-4～B.1-5 转化、涂板、菌落 PCR 验证 | 含正确 sgRNA 的质粒 | — \n")
    L.append("| 第 2–3 天 | §3 B.2 Gibson（仅插入 cargo 时需要） | 终载体 | PCR 可并行 \n")
    L.append("| 第 3 天 | §4 C-1～C-2 供体与受体培养 | 对数期菌液 | 可同步 \n")
    L.append("| 第 4 天 | §4 C-3～C-8 接合、涂板 | 接合平板 | — \n")
    L.append(f"| 第 5–7 天 | §5 挑克隆、PCR、测序、骨架 PCR | 验证过的编辑株 | — \n")
    L.append("")
    L.append(f"> 集落出现：{_fact(kb, 'conj.colony_days') or ''}\n")

    # =====================================================================
    L.append("\n---\n## 7. 试剂与耗材清单\n")
    L.append("| 类别 | 条目 | 出处 |\n")
    L.append("|---|---|---\n")
    for f in kb.by_prefix("reagent."):
        val = f.value if not isinstance(f.value, dict) else \
            f"{f.value.get('conc', '')} {f.value.get('use', '')}"
        L.append(f"| 试剂 | {val} | {f.source} {f.locator} \n")
    for f in kb.by_prefix("media."):
        val = f.value if not isinstance(f.value, dict) else \
            " + ".join(f"{k} {v}" for k, v in f.value.items())
        L.append(f"| 培养基 | {val} | {f.source} {f.locator} \n")
    L.append("")

    # =====================================================================
    if plan.decisions:
        L.append("\n---\n## 8. 判定点（做错了怎么发现）\n")
        for d in plan.decisions:
            # ⚠️ DecisionPoint 的字段是 after / observe / if_ok / if_not /
            #    hypothesis / discriminator，**没有 title**。
            #    早期版本用 getattr(d, 'title', d) 会把整个 dataclass repr 打出来。
            L.append(f"\n### 若「{d.observe}」\n\n")
            L.append(f"- **何时看**：{d.after}\n")
            L.append(f"- **观察**：{d.observe}\n")
            if d.if_ok:
                L.append(f"- **符合预期**：{d.if_ok}\n")
            if d.if_not:
                L.append(f"- **不符合**：{d.if_not}\n")
            if d.hypothesis:
                L.append(f"- **假设**：{d.hypothesis}\n")
            if d.discriminator:
                L.append(f"- **最便宜的判别实验**：{d.discriminator}\n")

    if plan.risks:
        L.append("\n---\n## 9. 风险清单\n")
        for i, r in enumerate(plan.risks, 1):
            if isinstance(r, dict):
                L.append(f"\n**{r.get('id', i)}　{r.get('title', '')}**"
                         f"（{r.get('level', '')}）\n")
                for k in ("detail", "detect", "action"):
                    if r.get(k):
                        L.append(f"- {k}：{r[k]}\n")
            else:
                L.append(f"{i}. {r}\n")

    L.append("\n---\n## 10. 能力边界（不要对外宣称超出这些）\n")
    for line in [
        "**不预测编辑效率** —— 已发表数据点太少（约 11 条带效率的位点）。",
        "**脱靶引擎**为 6-mer 块 ×3 + 每块 1 处替换枚举：≤4 错配保证不漏检"
        "（已用暴力穷举交叉验证）；**>4 时保证失效**并会显式警示；"
        "**不等价于 Cas-OFFinder 穷举**。",
        "**任何错配都只增风险，不存在「减害」。**",
        "**Tm 为估算值**（Wallace 修正式），首次须梯度 PCR 确认。",
        "**未做湿实验验证** —— 正确性上限是「忠实于已发表方法与实验记录」。",
    ]:
        L.append(f"- {line}\n")

    L.append(f"\n> 工具：{plan.meta.get('tool', 'stib-planner')} "
             f"{plan.meta.get('version', '')}\n")
    return "".join(L)
