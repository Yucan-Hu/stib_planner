"""
STIB Planner — 知识库层。

设计原则（README §1 核心判据）:
  **所有条件值必须带来源 locator，能点回论文页码 / 补充表 / 记录库日期。**
  不可回溯的值一律不进入知识库。

本模块是数据的**唯一权威来源**；CLI 与插件都从它读。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any

# ---------------------------------------------------------------------------
# 出处（provenance）
# ---------------------------------------------------------------------------

#: 来源标识 → 人类可读定位。
#: 每条知识都必须引用其中之一。
#: 本工具所依据论文的完整出处。
#: 只在这里写一次，下面各条来源共享它 —— 避免同一篇文献出现两种写法。
#: 课题组 2026-10 确认 DOI 为 10.1016/j.cels.2026.101650。
PAPER_CITATION = (
    "Yucan Hu, Qingmei Li, Yingying Li, Yu Zeng, Linggang Zheng, Juntao Shen, "
    "Xiang Gao, Guo-Ping Zhao, Wei Zhao, Lei Dai. "
    "Targeted genomic editing of human gut Bacteroides species based on "
    "CRISPR-associated transposases. Cell Systems, 2026, 101650. "
    "doi:10.1016/j.cels.2026.101650"
)
PAPER_DOI = "10.1016/j.cels.2026.101650"
PAPER_URL = "https://doi.org/10.1016/j.cels.2026.101650"

SOURCES: dict[str, str] = {
    "PAPER": f"论文：{PAPER_CITATION}",
    "PAPER_METHODS": f"论文 STAR Methods（{PAPER_DOI}）",
    "PAPER_RESULTS": f"论文 RESULTS（{PAPER_DOI}）",
    "PAPER_FIG2": f"论文 Figure 2（插入窗口/定向性；{PAPER_DOI}）",
    "PAPER_FIG4": f"论文 Figure 4（菌群编辑，56-63 bp 窗口；{PAPER_DOI}）",
    "TABLE_S2": "Table S2. sgRNA used in this study.xlsx",
    "TABLE_S3": "Table S3. Strains and culture medium.xlsx",
    "TABLE_S4": "Table S4. 16S rDNA verification of isolates (BLASTn)",
    "TABLE_S7": "Table S7. Primers and probes used in this study.xlsx",
    "SOP_MAIN": "039_STIB 使用protocol.md（中文 SOP，442 行）",
    "SOP_EN": "039_STIB protocol (EN).md（英文 SOP，441 行）",
    "SOP_BUILD": "039_STIB 载体构建.md",
    "NOTEBOOK": "实验记录库 06__experiment/061__CAST（269 条，2023-06 → 2026-02）",
    "PEER_REVIEW": "透明同行评审记录（comments.txt / email.txt）",
}

#: 证据等级。abstract_only 等低等级条目在输出时必须带警示。
EVIDENCE_LEVELS = ("full_text", "supplement_table", "lab_notebook", "peer_review", "abstract_only")


@dataclass
class Fact:
    """一条带出处的知识。"""

    key: str
    value: Any
    source: str                 # SOURCES 的键
    locator: str = ""           # 章节/表号/日期
    note: str = ""
    evidence: str = "full_text"

    def ref(self) -> str:
        return f"{SOURCES.get(self.source, self.source)}" + (f" — {self.locator}" if self.locator else "")

    def as_dict(self) -> dict:
        d = asdict(self)
        d["ref"] = self.ref()
        return d


@dataclass
class Conflict:
    """
    一处**可判定的矛盾**：不同来源对同一事实给出不同说法。

    这类条目**禁止静默合并**（README §4.2 硬约束 2）。
    """

    key: str
    claims: list[dict]          # [{source, locator, statement}]
    kind: str                   # 'contradiction' | 'open_question'
    severity: str               # 'high' | 'medium' | 'low'
    guidance: str = ""          # 建议动作
    discriminator: str = ""     # 若是 open_question：最便宜的判别实验

    def as_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# 知识库
# ---------------------------------------------------------------------------

class KnowledgeBase:
    """STIB 领域知识，全部带出处。"""

    def __init__(self) -> None:
        self.facts: list[Fact] = []
        self.conflicts: list[Conflict] = []
        self._build()

    # -- 查询 ---------------------------------------------------------------

    def get(self, key: str) -> Fact | None:
        for f in self.facts:
            if f.key == key:
                return f
        return None

    def by_prefix(self, prefix: str) -> list[Fact]:
        return [f for f in self.facts if f.key.startswith(prefix)]

    def conflicts_of(self, kind: str | None = None) -> list[Conflict]:
        return [c for c in self.conflicts if kind is None or c.kind == kind]

    # -- 构建 ---------------------------------------------------------------

    def _add(self, key: str, value: Any, source: str, locator: str = "",
             note: str = "", evidence: str = "full_text") -> None:
        if source not in SOURCES:
            raise ValueError(f"未知来源 {source!r}：知识条目必须引用 SOURCES 中的出处")
        self.facts.append(Fact(key, value, source, locator, note, evidence))

    def _build(self) -> None:
        self._build_system()
        self._build_reagents()
        self._build_sites()
        self._build_loci()
        self._build_loci_inputs()
        self._build_risks()
        self._build_conflicts()

    # 系统机制
    def _build_system(self) -> None:
        self._add("system.name", "STIB (ShCAST-based Transient Insertion system for Bacteroides)",
                  "PAPER_RESULTS", "SUMMARY")
        self._add("system.origin", "ShCAST / Type V-K CRISPR-associated transposase (Scytonema hofmanni)",
                  "SOP_MAIN", "§1 背景")
        self._add("system.pam", "GTN（论文与 SOP 写 NGTN，但 Table S2 全部 25 条实测为 GTN）",
                  "TABLE_S2", "PAM 列",
                  note="扫描时按 GTN 处理；报告会同时提示 NGTN 的表述差异", evidence="supplement_table")
        self._add("system.spacer_len", 23, "SOP_MAIN", "§4 阶段A 第3条")
        self._add("system.insert_window", (56, 66), "PAPER_FIG2",
                  "Figure 2E/2F 与 Figure 4I-K 的并集（v0.6 口径更正）",
                  note="**这是两次独立观察到的分布，不是验收标准。** "
                       "Figure 2E/2F（纯培养，Bt01 98.2% / Bt02 77.8%）落在 60-66 bp；"
                       "Figure 4I-K（菌群编辑）落在 56-63 bp。"
                       "**区间内偏移不构成风险** —— 仅方向敏感/框内插入（启动子捕获、融合蛋白）"
                       "才需精确读坐标。v0.6 前以 (60,66) 为单一主值，"
                       "会让使用者为了挤进 66 bp 反复更换 spacer —— 那是白做工。")
        self._add("system.tsd_len", 5, "SOP_MAIN", "§1 背景",
                  note="5 bp target site duplication")
        self._add("system.orientation", "单向：LE 靠近靶位点（T-LR >99%）",
                  "PAPER_FIG2", "Figure 2C")
        self._add("system.vector", "R6K 复制子自杀载体；在 pir- 拟杆菌内不能复制 → 瞬时递送，通常无需消除",
                  "SOP_MAIN", "§4 阶段D 质粒消除")

    # 试剂
    def _build_reagents(self) -> None:
        self._add("reagent.amp", {"conc": "100 µg/mL", "use": "供体 E. coli 筛选（LB）"},
                  "PAPER_METHODS", "Transposition experiment")
        self._add("reagent.genta", {"conc": "200 µg/mL", "use": "杀供体 / 反筛"},
                  "PAPER_METHODS", "Transposition experiment")
        self._add("reagent.cm_bt_bo", {"conc": "25 µg/mL", "use": "B. thetaiotaomicron / B. ovatus 筛选"},
                  "PAPER_METHODS", "Transposition experiment")
        self._add("reagent.cm_other", {"conc": "12.5 µg/mL", "use": "其他 Bacteroides spp. 筛选"},
                  "PAPER_METHODS", "Transposition experiment")
        self._add("reagent.erm", {"conc": "25 µg/mL", "use": "ermG 标记"},
                  "TABLE_S3", "Culture agar 列", evidence="supplement_table")
        self._add("reagent.dap", {"conc": "60 µg/mL (0.06 g/L)", "use": "ΔdapA 缺陷型供体（S17DMK/WM3064）必需"},
                  "PAPER_METHODS", "Bacterial strains")
        self._add("reagent.bsai", "BsaI-HFv2 (NEB R3733)", "SOP_MAIN", "§2.2")
        self._add("reagent.gibson", "NEBuilder HiFi Assembly (NEB E2621S)", "SOP_MAIN", "§2.2")
        self._add("reagent.blood", "8% (v/v) 脱纤维血（马血或绵羊血均可）", "SOP_MAIN", "§9 配方2",
                  note="**v0.6 收敛：种类不关键，关键是无菌、脱纤维。** 三条来源只是措辞不同："
                       "PAPER_METHODS（Transposition experiment）写 '8% horse blood'；"
                       "SOP_MAIN §9 配方2 写脱纤维马血；TABLE_S3 'Culture agar' 列写 "
                       "'8% sheep blod'（羊血，原表拼写有误）。"
                       "**血平板污染会直接毁掉整批接合**，这是唯一需要核实的点。")

        # ── B.1 sgRNA 更换：BsaI Golden Gate 的三个体系（v0.11）──
        #    全部来自 SOP_MAIN（039_STIB 使用protocol.md）§B.1 的表 2/3/4。
        #    这些是**精确体积**，不是示意值 —— 协议里要能照着加样。
        self._add("gg.phospho_anneal", {
            "total_µL": 20,
            "reagents": [
                ("sgRNA-oligoF (10 µM)", 2),
                ("sgRNA-oligoR (10 µM)", 2),
                ("T4 PNK", 1),
                ("10× T4 连接酶缓冲液", 2),
                ("ddH2O", 13),
            ],
            "program": [
                ("磷酸化", "37 °C", "30 min"),
                ("退火", "95 °C", "5 min"),
                ("Touch down", "95 → 16 °C，每循环 −1 °C", "15 s × 79"),
                ("保温", "16 °C", "∞"),
            ],
            "total_time": "约 55 min",
        }, "SOP_MAIN", "§B.1 步骤1 / 表2",
            note="与步骤2（载体酶切）**同步进行**。Touch down 共 79 循环，"
                 "自 95 °C 每循环降 1 °C 至 16 °C。")

        self._add("gg.vector_digest", {
            "total_µL": 10,
            "reagents": [
                ("质粒 (0.015–0.05 pmol)", 2),
                ("BsaI-HFv2 (20 U)", 1),
                ("10× T4 连接酶缓冲液", 1),
                ("ddH2O", 6),
            ],
            "program": [("酶切", "37 °C", "1 h"), ("ramp", "→ 4 °C", "")],
        }, "SOP_MAIN", "§B.1 步骤2 / 表3",
            note="与步骤1同步进行。质粒用量 0.015–0.05 pmol。")

        self._add("gg.assembly", {
            "total_µL": 20,
            "reagents": [
                ("步骤2 酶切产物", 10),
                ("退火寡核苷酸（步骤1）", 1),
                ("Hi-T4 连接酶", 1),
                ("10× T4 连接酶缓冲液", 1),
                ("ddH2O", 7),
            ],
            "program": [
                ("酶切", "37 °C", "5 min", "15 循环"),
                ("连接", "25 °C", "5 min", "返回上一步"),
                ("酶切", "37 °C", "30 min", "1"),
                ("失活", "80 °C", "10 min", "1"),
                ("保温", "4 °C", "∞", "1"),
            ],
            "total_time": "约 3 h",
        }, "SOP_MAIN", "§B.1 步骤3 / 表4-1、4-2",
            note="步骤 1–2 循环 15 次。若用**不耐热** T4 连接酶，"
                 "连接步骤 25 °C 改为 **16 °C**；时间紧张时可 37 °C 恒温 1–1.5 h。")

        # ── B.2 cargo 更换：高保真 PCR 与 Gibson（v0.11）──
        self._add("pcr.phanta_max", {
            "total_µL": 50,
            "reagents": [
                ("2× PhantaMax", 25),
                ("引物 F", 2),
                ("引物 R", 2),
                ("模板", 0.5),
                ("ddH2O", 20.5),
            ],
            "program": [
                ("预变性", "95 °C", "40 s（质粒）/ 3 min（基因组·长片段）", "1"),
                ("变性", "95 °C", "15 s", "29–30"),
                ("退火", "Tm", "15 s", ""),
                ("延伸", "72 °C", "30–60 s/kb", ""),
                ("终延伸", "72 °C", "3 min", "1"),
                ("保温", "16 °C", "∞", "1"),
            ],
        }, "SOP_MAIN", "§B.2 步骤2 / 表6",
            note="预变性时间依模板调整：质粒 40 s，基因组或长片段 3 min。")

        self._add("pcr.product_cleanup",
                  "跑胶确认条带单一 → DMT 酶 37 °C 消化 1 h（去除甲基化模板）→ 纯化回收，"
                  "浓度宜 >20 ng/µL；条带不单一者切胶回收",
                  "SOP_MAIN", "§B.2 步骤3")

        self._add("gibson.reaction", {
            "total_µL": 10,
            "reagents": [
                ("2× Basic Assembly Mix", 5),
                ("Frag1（载体片段）", "适量"),
                ("Frag2/Frag3（插入片段）", "适量"),
                ("ddH2O", "补至 10"),
            ],
            "program": [
                ("孵育", "50 °C", "两片段 15–20 min；多片段 40–60 min"),
            ],
            "pmol": {"载体片段": 0.08, "插入片段": 0.16, "允许范围": "0.01–0.25 pmol"},
            "formula": "pmol = weight(ng) / (length(bp) × 0.65)",
        }, "SOP_MAIN", "§B.2 步骤4 / 表7")

        self._add("transform.heat_shock", [
            "冰上解冻 100 µL 感受态 5–10 min",
            "加 10 µL 组装产物，轻弹混匀 4–5 次",
            "冰浴 30 min",
            "42 °C 热激 60–90 s",
            "冰浴 2 min",
            "加 900 µL 室温 SOC（或 LB）",
            "37 °C、250 rpm 振荡 60 min",
            "100 µL 涂板；长片段大质粒可 30 °C 培养 24 h",
        ], "SOP_MAIN", "§B.2 步骤5（NEB 热激法）",
            note="第 8 步涂布 LB + Amp 100 µg/mL。Golden Gate 产物同样用此法（取 5–10 µL）。")

        self._add("verify.colony_pcr",
                  "挑克隆接种 3 mL LB + Amp 过夜 → 菌落 PCR + 测序，"
                  "确认 LE / cargo / RE 各 junction 正确；"
                  "换 sgRNA 时扩增 sgRNA 表达盒测序确认 spacer",
                  "SOP_MAIN", "§B.1 步骤6、§B.2 步骤6")

        # 培养基
        self._add("media.BHIS", {
            "BHI": "37 g/L",
            "L-cysteine": "0.5 g/L",
            "VK3": "0.5 mg/L（1 mg/mL 母液 0.2 mL/400 mL）",
            "Hemin": "5 mg/L（0.5 mg/mL 母液 4 mL/400 mL）",
            "Agar(固体)": "15 g/L",
        }, "SOP_MAIN", "§9 配方1")
        self._add("media.blood_plate", "BHIS-琼脂冷却至 50 °C 后加 8%(v/v) 脱纤维马血（厌氧箱内操作）",
                  "SOP_MAIN", "§9 配方2")
        self._add("media.M9S", "M9 基础 + 50 mg/L L-半胱氨酸 + 5 mg/L 氯化血红素 + 2.5 µg/L VK1 + 2 mg/L FeSO4·7H2O + 5 µg/L 维生素B12 + 0.7% 琼脂糖",
                  "PAPER_METHODS", "Bacterial strains")
        self._add("media.anaerobic_atm", "85% N2 / 10% H2 / 5% CO2 (Coy 厌氧箱)", "SOP_MAIN", "§3 仪器")

    # 接合与筛选条件
    def _build_sites(self) -> None:
        self._add("conj.donor_od", "供体 1:100 转接后 2.5-3 h，OD600 ≈ 0.4（对数前期）",
                  "SOP_MAIN", "§4 阶段C Day2 第1条")
        self._add("conj.recipient_od", "受体 1:100 转接 BHIS，37 °C 厌氧静置约 6 h，OD600 ≈ 0.8",
                  "SOP_MAIN", "§4 阶段C Day2 第2-3条")
        self._add("conj.ratio", "供体:受体 1:1（各 1 mL）；也可 1:5（200 µL : 1000 µL）",
                  "SOP_MAIN", "§4 阶段C Day2 第4条")
        self._add("conj.spin", "4000 × g，20 min，弃上清，约 100 µL BHIS 重悬",
                  "SOP_MAIN", "§4 阶段C Day2 第4条")
        self._add("conj.spot", "滴加至无抗 BHIB 血平板中心（硬币大小），约 25 mL/板",
                  "SOP_MAIN", "§4 阶段C Day2 第5条")
        self._add("conj.mating", "37 °C 静置 20-24 h，好氧（受体以菌苔形式短期耐氧）",
                  "SOP_MAIN", "§4 阶段C Day2 第5条",
                  note="⚠️ 记录库显示该条件从未收敛，见 open_questions")
        self._add("conj.recovery", "棉签刮取菌苔 → 1 mL BHI 重悬 → 梯度稀释",
                  "SOP_MAIN", "§4 阶段C Day3 第6条")
        self._add("conj.plating", "原液与 10⁻² 各 100 µL 涂 BHIS + Genta 200 + Cm 25；37 °C 厌氧 2-3 天",
                  "SOP_MAIN", "§4 阶段C Day3 第7条")
        self._add("conj.colony_days", {"BHIS/M9S-glucose": "2-3 天", "M9S-inulin": "3-4 天", "B. vulgatus": "可延至 7 天"},
                  "SOP_MAIN", "§4 阶段C Day3 第7条",
                  note="B. vulgatus 7 天数据来自记录库 20240111", evidence="lab_notebook")

        # 鉴定
        self._add("verify.sanger", "LE/RE junction 测序须确认三点：(i) 整合落在 PAM 下游 56–66 bp 观察区间内（区间内均可接受）；(ii) 存在 5 bp TSD；(iii) 边界与载体 LE/RE 一致",
                  "SOP_MAIN", "§4 阶段D D.2")
        self._add("verify.backbone_pcr", "骨架特异引物（ampR/ori）PCR，**应为阴性**；阳性者划线纯化复检，仍阳性则弃去",
                  "SOP_MAIN", "§4 阶段D 表8 / §8")
        self._add("verify.16s", "27F/1492R 全长 16S，Sanger 双向，EzBioCloud 比对",
                  "SOP_MAIN", "§4 阶段D D.3")
        self._add("verify.16s_lysis", "菌液 + ALP 裂解液 1:9（2 µL + 18 µL），90 °C 10 min，离心取上清",
                  "SOP_MAIN", "§4 阶段D D.3 第1条")

    # 位点（gene → locus）映射所需的最小结构
    def _build_loci(self) -> None:
        # Bt01 对应的论文位点定义
        self._add("locus.Bt01.spacer", "ATTCGGAAGAAGACGTGGTATCA", "TABLE_S2", "sgRNA26 行",
                  note="靶向 tdk；论文主用位点", evidence="supplement_table")
        self._add("locus.Bt01.pam", "GTT", "TABLE_S2", "sgRNA26 行", evidence="supplement_table")
        self._add("locus.Bt01.gene", "tdk", "TABLE_S2", "spacer_location 列", evidence="supplement_table")

        # 论文报告的多位点效率（用于给用户预期）
        self._add("eff.bt_sites", "Bt 位点编辑效率 29%-91%（Figure 2A，p=0.0188）",
                  "PAPER_RESULTS", "Systematic characterization")
        self._add("eff.cargo_size", "载货 2.5 / 5.5 / 8.4 kb 时效率均 >85%，无显著差异（p=0.8343）",
                  "PAPER_RESULTS", "Figure 2B")
        self._add("eff.species", {
            "B. ovatus": ">85%（Bo02 例外，30%）",
            "B. salyersiae DA1247": "6%-100%",
            "B. vulgatus": "Bv02/Bv03 为 0%，其余 >93%",
            "B. fragilis / B. uniformis（靶向 16S）": ">90%",
        }, "PAPER_RESULTS", "Genomic editing of different Bacteroides species")
        self._add("eff.ontarget", ">97% on-target（摘要）", "PAPER_RESULTS", "SUMMARY",
                  note="⚠️ 记录库显示 RBS 组在 2023-11 重复失败，见 conflicts")

    # 风险
    def _build_risks(self) -> None:
        R = [
            dict(id="RK1", level="high", title="共整合（copy-and-paste 副产物）",
                 detail="Type V-K 缺 TnsA，天然经复制型通路转座，产生夹带载体骨架的共整合。STIB 用切口归巢内切酶融合 TnsB + ATPase 融合 Cas12k 抑制，但仍需筛查。",
                 detect="骨架特异引物（ampR / ori 区域）PCR，**阴性为正确**",
                 action="阳性者划线纯化传代后复检；持续阳性则弃去该克隆并重挑。"
                        "⚠️ **共整合与位点选择是否存在因果关系尚不确定** —— 不要据此断定位点有问题而盲目换位点",
                 source="SOP_MAIN", locator="§8 Box 1 第1条"),
            dict(id="RK2", level="high", title="脱靶插入（含 RNA 引导与 RNA 非依赖两类）",
                 detail="**关键区在 PAM 近端**：ShCAST 的 PAM 位于 spacer 上游，Cas12k 与 RNA:DNA 异源双链的接触集中在 **前 17 nt（helix）**，3' 端 6 nt（tail）为延伸尾部。因此同样错配数下，**错配落在前 17 nt 的位点更危险**。（依据 *Structural basis for the assembly of the type V CRISPR-associated transposon complex*；v0.6 更正 —— 此前误套 Cas9 的 3' 端种子区惯例）",
                 detect="用本工具的 offtarget 做全基因组 ≤4 错配枚举（已交叉验证）；发表级结果再用 Cas-OFFinder 复核",
                 action="排除 helix 段（前 17 nt）存在近同源位点的 spacer；重要克隆以 WGS 复核",
                 source="SOP_MAIN", locator="§4 阶段A 第5条 / §8 第3条"),
            dict(id="RK3", level="high", title="插入方向与位置固定性",
                 detail="ShCAST 单向插入，落点集中在 PAM 下游 56–66 bp。"
                        "**该区间是观察分布，不是验收标准** —— Figure 2（纯培养）60–66 bp、"
                        "Figure 4（菌群）56–63 bp。区间内偏移不构成风险。",
                 detect="设计阶段即确认窗口位置与方向（启动子捕获、框内插入等方向敏感应用尤其关键）",
                 action="窗口若落在必需基因编码区或上游调控元件，必须换位点",
                 source="SOP_MAIN", locator="§4 阶段A 第4条 / §8 第4条"),
            dict(id="RK4", level="high", title="ΔdapA 缺陷型供体的血平板漏加 DAP → 转化子为零",
                 detail="记录库 2024-03 实测：S17-ΔdapA 与 WM3064 在 3 生物学 × 3 技术重复下**全部为 0 CFU**，事后发现血平板未加 DAP。",
                 detect="接合前核对血平板配方含 DAP 0.06 g/L",
                 action="漏加则整批作废重做；建议在平板标签上写明 +DAP",
                 source="NOTEBOOK", locator="20240322-确定donor-接合",
                 note="这是论文与 SOP 都未提示的操作陷阱"),
            dict(id="RK5", level="high", title="16S 无法区分近缘种 → 物种误判",
                 detail="Table S4 显示 B. thetaiotaomicron 的 top hits 中 B. faecis MAJ27 相似度 **98.43%**；P. vulgatus 之后 P. dorei **97.10%**。记录库中确有把 bf 误判为 bt 的实例。",
                 detect="16S 比对须看 top-5 hits 及相似度差值，不可只看第一条",
                 action="必要时补做 species-specific PCR 或 WGS；交接菌株时双人核对",
                 source="TABLE_S4", locator="B.Theta / B.vulgatus 比对结果"),
            dict(id="RK6", level="medium", title="菊粉筛选非杀菌 → 背景生长",
                 detail="菊粉选择不杀菌，必然出现背景生长，integrants 必须重新划线纯化。",
                 detect="观察平板是否有弥散背景",
                 action="always restreak for isolation；记录库 2024-07 另有 RBS40 突变导致的构建失败",
                 source="SOP_MAIN", locator="§7"),
            dict(id="RK7", level="medium", title="质粒骨架自发突变（RBS 漂移）",
                 detail="记录库 2026-02 记录 psh4300C 的 RBS 发生突变需要回补，且任务目的是样品外送。",
                 detect="每次从冻存重新划线制备；外送/交接前对 RBS 与关键元件测序核对",
                 action="避免长期传代累积突变；建立骨架测序基线",
                 source="NOTEBOOK", locator="20260227-CAST-恢复突变psh4300C"),
            dict(id="RK8", level="medium", title="抗生素失效 / 浓度不足 → 假阳性",
                 detail="记录库 2025-08 实测 12.5 µg/mL 时 bv 组与 NC 组同时生长；2024-11 记录「5FudR 失效」。",
                 detect="每次实验设无抗对照与 NC；先做 MIC 梯度（25/12.5/6.25/3.125 µg/mL）",
                 action="NC 生长即实验作废；按目标物种调整 Cm（bt/bo 25，其他 12.5）",
                 source="NOTEBOOK", locator="20250819 / 20241124"),
            dict(id="RK9", level="medium", title="慢生长物种假阴性",
                 detail="记录库 2024-01 记录 bv 集落延迟至 **7 天** 才出现（常规读数 2-3 天），当时 bv 菌况不佳。",
                 detect="对 B. vulgatus 等慢生长种，读数延长至 7 天",
                 action="勿在 2-3 天未见克隆时判定为阴性",
                 source="NOTEBOOK", locator="20240111-CAST-bobs重新菌落-高通量-PCR鉴定"),
            dict(id="RK10", level="medium", title="近缘种 16S 通用靶点可能误伤非靶标物种",
                 detail="sgRNA81 靶向 16S 保守区，用于 B. fragilis / B. uniformis 通用编辑；但 16S 保守区在其他物种中亦存在。",
                 detect="对目标菌群做 PAM+spacer 的基因组/宏基因组比对",
                 action="菌群编辑中若需物种特异，应在鉴定阶段用 16S 全长确认物种身份",
                 source="PAPER_RESULTS", locator="Targeted insertion of gene-of-interest"),
            # ── v0.6：本领域**真实**的实验室安全边界（而非通用生物安全审查）──
            dict(id="RK11", level="high", title="筛选标记应为非临床相关抗生素",
                 detail="**这是论文自己提出的要求。** 论文原文指出：氯霉素与头孢西丁抗性基因曾被"
                        "用于筛选拟杆菌复制型质粒，但**头孢西丁与碳青霉烯类是临床相关抗生素，"
                        "实验室应尽量避免使用，以防抗性基因扩散**。"
                        "本工具因此不建议把临床关键抗生素作为筛选标记。",
                 detect="检查设计中的筛选标记是否为临床关键抗生素（碳青霉烯类、头孢西丁、"
                        "黏菌素、替加环素、万古霉素等）",
                 action="**改用非抗生素筛选**（如论文的菊粉利用盒 pLGB28/29/31）"
                        "或非临床相关标记（红霉素/四环素用于体外研究时需评估）",
                 source="PAPER_RESULTS", locator="DISCUSSION（抗生素选择）",
                 note="论文明确点名的两个：cefoxitin（头孢西丁）、carbapenems（碳青霉烯类）"),
            dict(id="RK12", level="medium", title="工程菌株的废弃物与灭活",
                 detail="编辑后的拟杆菌若携带抗生素抗性标记，其培养物与废弃物应按机构规定灭活，"
                        "不得直接排入下水道；菌株保藏与转移须符合机构生物安全等级。",
                 detect="实验前确认：本机构对基因工程微生物废弃物的灭活要求、菌株转移审批流程",
                 action="高压灭菌后废弃；菌株外送/交接前核对审批与包装要求",
                 source="SOP_MAIN", locator="（机构规定，非论文内容）",
                 note="这是实验室合规要求，不是方法学风险"),
        ]
        for r in R:
            self._add(f"risk.{r['id']}", r, r["source"], r["locator"],
                      note=r.get("note", ""),
                      evidence="lab_notebook" if r["source"] == "NOTEBOOK" else "full_text")

    # 冲突与未收敛变量
    def _build_conflicts(self) -> None:
        """
        冲突台账（v0.6 依据课题组反馈最终收敛）。

        已删除（用户明确指示）：
          - final_construct_choice → **删除**。只做 psh43 系列。
          - table_s5_missing       → **删除**。
          - syncom_size            → **删除**（按 Bacteroides 内按种记）。
          - mating_atmosphere      → 改为**按菌株的规则**（见 aerobic_capable）。

        已降级为 info（不再是待解决矛盾）：
          - blood_type    → 课题组已裁定：马血/绵羊血**都可以**，关键是无菌、脱纤维。
                           三条来源只是措辞不同，**不存在分歧** → 降为配方注记。
          - insert_window → 两个数值是**不同实验条件下观察到的分布**，既不是互斥主张，
                           也不是验收标准 → 降为观察记录，可接受窗口取并集 56-66 bp。

        ⚠️ 降级不等于丢弃：每条来源与出处仍完整保留（见 system.insert_window 与
        reagent.blood_plate 的 note），只是不再进入「待解决矛盾」清单 ——
        否则 open_questions 会让人以为还有悬而未决的分歧，从而去做无谓的确认实验。
        """
        self.conflicts: list[Conflict] = []

    def _build_loci_inputs(self) -> None:
        """位点设计输入（v0.3 新增，依据课题组反馈）。"""
        # 耐氧能力：决定接合气氛
        self._add("o2.aerobic_capable", {
            "species": ["B. thetaiotaomicron", "B. fragilis", "B. vulgatus",
                        "B. uniformis", "B. ovatus", "B. salyersiae"],
            "rule": "论文中测试过的这些物种**耐氧能力较强，可直接好氧接合**",
        }, "PAPER_RESULTS", "Genomic editing of different Bacteroides species",
            note="bf/bv/bu/bo/bs 在论文中均测试过；bt 为主力物种")
        self._add("o2.anaerobic_required", {
            "cases": [
                "客户菌株耐氧能力不好（未在论文中测试过的非模式种）",
                "菌群编辑（SynCom / 粪便菌群）",
            ],
            "rule": "这两类情况**必须厌氧接合**（或先做好氧/厌氧平行对照）",
        }, "SOP_MAIN", "§4 阶段C Day2 第5条",
            note="SOP 写死好氧，仅适用于耐氧菌株的单菌接合")

        # 构建选择：只做 psh43 系列
        self._add("construct.psh43_only", {
            "rule": "**只设计 psh43 系列**（Cas12k–TnsC 融合）",
            "reason": "课题组明确指示：只做 psh43 系列，其他构建（psh33xx 低-TnsC、psh53xx TniQ 融合）不做",
        }, "PAPER_RESULTS", "Figure 1C/1D",
            note="论文推 psh43 为最终方案，特异性 97.1%、效率 >80%")

        # PAM 偏好
        self._add("pam.preference", [
            {"pam": "GTT", "tier": 1, "reason": "课题组做过最多，经验最足"},
            {"pam": "GTA", "tier": 2, "reason": "GTN 中的次优选"},
            {"pam": "其他 GTN", "tier": 3, "reason": "可用但经验较少；论文与 SOP 写 NGTN，实测全部为 GTN"},
        ], "TABLE_S2", "PAM 列",
            note="Table S2 全部 25 条实测为 GTN（24×GTT + GTC + GTA）",
            evidence="supplement_table")

        # 插入位点放置规则
        self._add("locus.placement_gene_inactivation", {
            "general": "插入位点置于目标基因 **CDS 中央 80% 区域**（避开两端）",
            "stricter": "论文采用更严格的做法：**CDS 前 1/3 到 2/3 之间**",
            "reason": "避免插入在 CDS 两端而未能有效破坏功能；同时避开上游调控元件",
        }, "PAPER_METHODS", "Guide RNA design",
            note="论文：'gRNAs were designed to direct transposon insertion into the first third to two-thirds of target genes'")

        # ShCAST 插入的序列偏好（来自 TnsC/TnsB 的序列接触）
        self._add("shcast.sequence_features", {
            "tsd_preference": "TSD 内含 **GCWGC** 偏好（TnsB 偏好在此插入）",
            "tsd_flank_symmetry": "TSD 中心 **±5 bp** 处偏好 **(A/T)**（由两个 TnsB 单体的 K290 接触）",
            "bp_12_flank": "TSD **±12 bp** 处有核苷酸偏好（TnsB IIb 结构域 R416/T417/Q425/N428 接触）",
            "at_rich_upstream": {
                "length": "~25 bp",
                "role": "可容纳两圈 dsDNA 结合的 TnsC 丝状体",
                "note": "位于整合位点的**单侧**；TnsC 丝结合 AT 富集 DNA 后定向招募 TnsB，"
                        "决定 TnsC-LR 方向（与 RNA 引导位点同向）",
            },
            "orientation": "所有整合位点共享 **TnsC-LR** 取向",
        }, "PAPER_METHODS", "Guide RNA design / Tn-seq 分析",
            note="⚠️ 功能描述来自课题组提供的 ShCAST 机制原文（Science 2019 系统的后续分析）。"
                 "**相对 PAM 的精确坐标偏移尚未从原文核实**，故本工具只做定性偏好打分，"
                 "不硬编码偏移量。")



# ---------------------------------------------------------------------------
# 校验
# ---------------------------------------------------------------------------

def validate(kb: KnowledgeBase) -> list[str]:
    """
    校验知识库，返回问题列表（空 = 通过）。

    硬约束（README §1）:
      1. 来源锚定 —— 每条 fact 必须引用已知 SOURCES 且 evidence 合法
      2. 冲突显式 —— 每个 conflict 至少两个 claim，且 severity / kind 合法
    """
    problems: list[str] = []

    for f in kb.facts:
        if f.source not in SOURCES:
            problems.append(f"[来源未锚定] {f.key} 引用了未知来源 {f.source!r}")
        if f.evidence not in EVIDENCE_LEVELS:
            problems.append(f"[证据等级非法] {f.key} evidence={f.evidence!r}")
        if f.value in (None, "", [], {}):
            problems.append(f"[空值] {f.key} 的 value 为空")

    for c in kb.conflicts:
        if len(c.claims) < 2:
            problems.append(f"[冲突不成立] {c.key} 只有 {len(c.claims)} 条 claim，无法构成矛盾")
        if c.kind not in ("contradiction", "open_question"):
            problems.append(f"[kind 非法] {c.key} kind={c.kind!r}")
        if c.severity not in ("high", "medium", "low"):
            problems.append(f"[severity 非法] {c.key} severity={c.severity!r}")
        if c.kind == "open_question" and not c.discriminator:
            problems.append(f"[缺判别实验] 未收敛变量 {c.key} 必须给出 discriminator")

    return problems


if __name__ == "__main__":
    kb = KnowledgeBase()
    probs = validate(kb)
    print(f"facts:     {len(kb.facts)}")
    print(f"conflicts: {len(kb.conflicts)}")
    print(f"problems:  {len(probs)}")
    for p in probs:
        print("  -", p)
    raise SystemExit(1 if probs else 0)
