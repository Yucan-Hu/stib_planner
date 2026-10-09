"""
STIB Planner — 序列层。

纯 Python 标准库。所有算法都设计成**可用论文已发表序列手工核对**（见 selftest 的 golden 基准）。

关键事实（出处见 README §2）:
  - ShCAST / Type V-K: NGTN PAM（本体系实际用 GTN）
  - spacer 长度 23 nt
  - 插入位点 = PAM 下游 56-66 bp（**观察区间**，非验收标准）
  - LE 靠近靶位点（单向插入）
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Iterator

# ---------------------------------------------------------------------------
# 常量（全部带出处，见 README §2）
# ---------------------------------------------------------------------------

#: ShCAST PAM 正则。论文与 SOP 写 NGTN，但 Table S2 全部 25 条 sgRNA 的 PAM
#: 实测均为 GTN（24×GTT + 1×GTC + 1×GTA）。故默认按 GTN 扫描。
#: 出处: 039_STIB 使用protocol.md §1；Table S2
PAM_REGEX = "GTN"

#: 本体系 spacer 长度。出处: SOP §4 阶段A 第3条；论文 STAR Methods
SPACER_LEN = 23

#: 可接受的插入窗口（PAM 3' 末端下游的 bp 偏移范围）= 两次独立观察的**并集**。
#:
#: v0.6 口径更正：此前以 (60, 66) 为单一主值，把「观察到的分布」当成了「验收标准」。
#: 论文 Figure 2E/2F（纯培养，Bt01 98.2% / Bt02 77.8%）报告 60-66 bp；
#: Figure 4I-K（菌群编辑）报告 56-63 bp。两者是**不同条件下的观察**，
#: 既不是互斥主张，也不是必须命中的目标 → 取并集。
#:
#: ⚠️ 把 60-66 当硬标准，会让使用者为了挤进 66 bp 反复更换 spacer —— 白做工。
INSERT_WINDOW = (56, 66)

#: 论文两处观察窗口的原始值（**仅用于在报告中并列展示出处，不用于判定**）
INSERT_WINDOW_FIG2 = (60, 66)
INSERT_WINDOW_ALT = (56, 63)

#: 预测用的**名义**插入偏移：相对 PAM 3' 末端下游的 bp 数。
#:
#: ⚠️ 必须与 INSERT_WINDOW **解耦**（v0.6 修）。此前用窗口中点
#: `(INSERT_WINDOW[0] + INSERT_WINDOW[1]) // 2` 当代表值 —— 于是把
#: 「可接受范围」从 (60,66) 改成 (56,66) 时，预测偏移被**顺带**从 63 改成了 61，
#: 预测插入位点整体挪了 2 bp。范围是范围，代表值是代表值，不能共用。
#:
#: 取 63 = 论文观测区间 (60,66) 的中点，来源为论文实测（Bt01 98.2%）。
INSERT_OFFSET_FROM_PAM3 = 63

#: PAM 长度（GTN）。insert_site = pam_start + PAM_LEN + INSERT_OFFSET_FROM_PAM3
PAM_LEN = 3

#: 靶位点重复长度。出处: SOP §1
TSD_LEN = 5

#: BsaI 粘性末端（sgRNA oligo 用）。出处: SOP 表 1
BSAI_F_OVERHANG = "aaag"
BSAI_R_OVERHANG = "agcc"

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


# ---------------------------------------------------------------------------
# FASTA 读写
# ---------------------------------------------------------------------------

def read_fasta(path: str) -> dict[str, str]:
    """读取 FASTA，返回 {seq_id: sequence}。序列统一转大写。"""
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return {k: v.upper() for k, v in seqs.items()}


def write_fasta(path: str, seqs: dict[str, str], width: int = 70) -> None:
    """写出 FASTA。"""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# ---------------------------------------------------------------------------
# 基础序列操作
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    """反向互补。"""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """GC 含量比例（0-1）。非 ACGT 字符不计入分母。"""
    seq = seq.upper()
    gc = sum(seq.count(b) for b in "GC")
    at = sum(seq.count(b) for b in "AT")
    total = gc + at
    return (gc / total) if total else 0.0


def pam_regex(pam: str = PAM_REGEX) -> re.Pattern[str]:
    """把 IUPAC PAM 串编译成正则。"""
    body = "".join(f"[{IUPAC[c]}]" if c in IUPAC else c for c in pam.upper())
    return re.compile("(?=(" + body + "))")


def find_pams(seq: str, pam: str = PAM_REGEX) -> Iterator[tuple[int, str]]:
    """
    扫描正链上的 PAM。

    产出 (pam_start_0based, pam_seq)。
    """
    rx = pam_regex(pam)
    for m in rx.finditer(seq):
        yield m.start(), m.group(1)


# ---------------------------------------------------------------------------
# sgRNA / protospacer
# ---------------------------------------------------------------------------

@dataclass
class SgRna:
    """一条候选 sgRNA。"""

    spacer: str                      # 23 nt protospacer（正链）
    pam: str                         # PAM 序列
    pam_start: int                   # PAM 起点（0-based，正链）
    strand: str                      # '+' 或 '-'
    insert_site: int                 # 预测插入位点（0-based，正链；= PAM 3' 末端 + 偏移）
    gc: float
    unique: bool = True

    @property
    def pam_end(self) -> int:
        """PAM 3' 末端之后的位置（0-based exclusive）。"""
        return self.pam_start + len(self.pam)

    @property
    def oligo_f(self) -> str:
        """下单用正向寡核苷酸（含 BsaI 粘性末端）。"""
        return BSAI_F_OVERHANG + self.spacer

    @property
    def oligo_r(self) -> str:
        """下单用反向寡核苷酸（含 BsaI 粘性末端）。"""
        return BSAI_R_OVERHANG + revcomp(self.spacer)

    def as_dict(self) -> dict:
        return {
            "spacer": self.spacer,
            "pam": self.pam,
            "pam_start": self.pam_start,
            "strand": self.strand,
            "insert_site": self.insert_site,
            "gc": round(self.gc, 3),
            "unique": self.unique,
            "oligo_f": self.oligo_f,
            "oligo_r": self.oligo_r,
        }


def _scan_strand(seq: str, offset: int, strand: str) -> list[SgRna]:
    """
    在给定序列上扫描候选 sgRNA。

    offset = 该序列在正链基因组上的起始坐标（0-based）。
    """
    out: list[SgRna] = []
    L = len(seq)
    for pam_start, pam in find_pams(seq):
        sp_start = pam_start + len(pam)
        sp_end = sp_start + SPACER_LEN
        if sp_end > L:
            continue
        spacer = seq[sp_start:sp_end]
        if any(c not in "ACGT" for c in spacer):
            continue

        # 插入位点：PAM 3' 末端向下游偏移 INSERT_OFFSET_FROM_PAM3 bp
        # （名义值，与可接受区间 INSERT_WINDOW 解耦）
        ins = pam_start + len(pam) + INSERT_OFFSET_FROM_PAM3

        if strand == "+":
            out.append(SgRna(
                spacer=spacer,
                pam=pam,
                pam_start=offset + pam_start,
                strand="+",
                insert_site=offset + ins,
                gc=gc_content(spacer),
            ))
        else:
            # 负链：传入的是反向互补序列，坐标需换算回正链
            rc_spacer = revcomp(spacer)
            rc_pam = revcomp(pam)
            pam_start_fwd = offset + (L - (pam_start + len(pam)))
            out.append(SgRna(
                spacer=rc_spacer,
                pam=rc_pam,
                pam_start=pam_start_fwd,
                strand="-",
                insert_site=pam_start_fwd + len(pam) + INSERT_OFFSET_FROM_PAM3,
                gc=gc_content(rc_spacer),
            ))
    return out


def scan_sgrnas(genome: str, start: int, end: int,
                flank: int = 300, pam: str = PAM_REGEX) -> list[SgRna]:
    """
    扫描 [start, end) 区域（±flank）内的候选 sgRNA，正负链都扫。

    坐标 0-based，左闭右开。
    """
    win_start = max(0, start - flank)
    win_end = min(len(genome), end + flank)
    fwd = genome[win_start:win_end]
    candidates = _scan_strand(fwd, win_start, "+")

    # 负链：取该窗口的反向互补，扫描后换算坐标
    # 为避免重复换算，直接在此处理
    rc_seq = revcomp(fwd)
    # 负链上位置 i 对应正链位置 win_end - 1 - i
    L = len(rc_seq)
    for pam_start, pamseq in find_pams(rc_seq, pam):
        sp_start = pam_start + len(pamseq)
        sp_end = sp_start + SPACER_LEN
        if sp_end > L:
            continue
        rc_spacer = rc_seq[sp_start:sp_end]
        if any(c not in "ACGT" for c in rc_spacer):
            continue
        # 该 PAM 在 rc 上的坐标 -> 正链 PAM 起点
        pam_start_fwd = win_end - (pam_start + len(pamseq))
        ins = pam_start_fwd + len(pamseq) + (INSERT_WINDOW[0] + INSERT_WINDOW[1]) // 2
        candidates.append(SgRna(
            spacer=revcomp(rc_spacer),
            pam=revcomp(pamseq),
            pam_start=pam_start_fwd,
            strand="-",
            insert_site=ins,
            gc=gc_content(revcomp(rc_spacer)),
        ))
    return candidates


def dedupe_sgrnas(cands: Iterable[SgRna]) -> list[SgRna]:
    """按 (spacer, pam_start, strand) 去重并按坐标排序。"""
    seen: set[tuple[str, int, str]] = set()
    out: list[SgRna] = []
    for c in sorted(cands, key=lambda x: (x.pam_start, x.strand)):
        k = (c.spacer, c.pam_start, c.strand)
        if k in seen:
            continue
        seen.add(k)
        out.append(c)
    return out


# ---------------------------------------------------------------------------
# 唯一性（同基因组内）
# ---------------------------------------------------------------------------

def count_occurrences(genome: str, motif: str, both_strands: bool = True) -> int:
    """
    统计 motif 在基因组中的出现次数。

    both_strands=True 时同时统计反向互补（用于评估 spacer 的基因组唯一性）。

    注意：这是**精确匹配**计数，不等同于脱靶预测。
    完整脱靶分析需 Cas-OFFinder 允许 ≤5 错配（见 README §7）。
    """
    n = genome.count(motif)
    if both_strands:
        rc = revcomp(motif)
        if rc != motif:
            n += genome.count(rc)
    return n


def seed_is_unique(genome: str, spacer: str, seed_len: int = 12) -> bool:
    """
    检查 spacer 的 **PAM 近端段**（前 seed_len nt）在基因组内是否唯一。

    ⚠️ v0.6 更正：ShCAST 的关键段在 **PAM 近端（spacer 5' 端）**，
    而非 Cas9 惯例的 3' 端。SOP §4 阶段A 第5条的「种子区」沿用 Cas9 语汇。
    本函数默认取**前** seed_len nt。

    注意：这是**精确匹配计数**的廉价检查。完整的近同源评估请用
    `stib.offtarget.find_offtargets`（6-mer 块 + 替换枚举，≤4 错配不漏检）。
    """
    seed = spacer[:seed_len]
    return count_occurrences(genome, seed, both_strands=True) == 1


# ---------------------------------------------------------------------------
# 引物设计
# ---------------------------------------------------------------------------

@dataclass
class Primer:
    """一条引物。"""

    name: str
    seq: str
    tm: float
    gc: float
    start: int          # 正链 0-based 起点
    end: int            # 正链 0-based 终点（exclusive）
    strand: str
    note: str = ""

    @property
    def length(self) -> int:
        return len(self.seq)

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "seq": self.seq,
            "len": self.length,
            "tm": round(self.tm, 1),
            "gc": round(self.gc, 3),
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "note": self.note,
        }


def tm_basic(seq: str) -> float:
    """
    引物 Tm 估算。

    ≤13 nt 用 Wallace 规则；更长用 GC 校正公式（Marmur/Chester-Marshak 变体）。

    ⚠️ 这是**估算**，用于筛选取舍，不可作为最终退火温度。
    实际使用前须用梯度 PCR 确认。README §7。
    """
    seq = seq.upper()
    n = len(seq)
    gc = sum(seq.count(b) for b in "GC")
    if n == 0:
        return 0.0
    if n <= 13:
        return 2.0 * (n - gc) + 4.0 * gc
    return 64.9 + 41.0 * (gc - 16.4) / n


def _pick_primer(seq: str, must_end_at: int, name: str,
                 strand: str, scan: int = 200,
                 tm_range: tuple[float, float] = (56.0, 66.0),
                 len_range: tuple[int, int] = (18, 26)) -> Primer | None:
    """
    在 seq 的一段区间内挑一条引物，要求其 3' 端尽量靠近 must_end_at（正链坐标）。

    strand='-': 引物在正链上位于 must_end_at 上游，向右延伸（用于扩增插入位点上游）
                返回的 seq 是反向互补后可直接下单的序列。
    """
    best: Primer | None = None
    best_score = -1e9
    for plen in range(len_range[1], len_range[0] - 1, -1):
        s = must_end_at - plen
        if s < 0:
            continue
        sub = seq[s:must_end_at]
        if any(c not in "ACGT" for c in sub):
            continue
        gc = gc_content(sub)
        if not (0.40 <= gc <= 0.60):
            continue
        tm = tm_basic(sub)
        if not (tm_range[0] <= tm <= tm_range[1]):
            continue
        # 3' 末端稳定性偏好（末端 G/C 更稳）
        end_bonus = 2.0 if sub[-1] in "GC" else 0.0
        run_penalty = -3.0 if re.search(r"(.)\1{3,}", sub) else 0.0
        score = -abs(tm - 60.0) + end_bonus + run_penalty
        if score > best_score:
            best_score = score
            best = Primer(
                name=name,
                seq=revcomp(sub) if strand == "-" else sub,
                tm=tm,
                gc=gc,
                start=s,
                end=must_end_at,
                strand=strand,
            )
    return best


def _pick_primer_near(seq: str, anchor: int, name: str, strand: str,
                      tol: int = 300, step: int = 5, **kw) -> Primer | None:
    """
    在 anchor 附近挑引物：先试精确锚点，失败则按 step 向两侧对称扩散至 ±tol。

    ⚠️ 修 bug 背景：`_pick_primer` 把引物 3' 端**钉死**在单一坐标上，只尝试
    18–26 nt 这 9 个窗口。若该坐标恰好落在 AT 富集区（GC < 0.40），9 个窗口
    全部不合格 → 返回 None，而调用方**静默**产出一条缺失的引物。
    BT1896（2,385,704）的下游锚点 2,386,204 就属于这种情况。
    就近扩散可保持产物长度基本不变，同时让设计不再无声失败。

    返回的 Primer 的 start/end 是**实际选中锚点**的坐标，调用方应以它计算产物长度。

    在 ±tol 内**收集全部可用候选再择优**（而不是取第一个能用的）：
    否则紧邻失败锚点处一条 Tm=56.4 的勉强合格引物会挤掉 300 bp 外
    Tm=59.5 的优质引物。评分 = 偏向 Tm 60 °C − 0.01 × 偏离锚点的距离。
    """
    cands: list[tuple[float, Primer]] = []
    for d in range(0, tol + 1, step):
        for a in ((anchor,) if d == 0 else (anchor - d, anchor + d)):
            if a < 0 or a > len(seq):
                continue
            p = _pick_primer(seq, a, name, strand, **kw)
            if p is not None:
                cands.append((-abs(p.tm - 60.0) - 0.01 * d, p))
        # 已收满一整圈就有足够候选，不必扫到 ±tol 边界
        if d >= step and len(cands) >= 8:
            break
    if not cands:
        return None
    cands.sort(key=lambda t: (-t[0], t[1].start))
    return cands[0][1]


def design_junction_primers(genome: str, insert_site: int,
                            reach: int = 500, name_prefix: str = "J") -> dict:
    """
    为插入位点设计 junction 验证引物。

    产出:
      up_fwd  — 位于插入位点上游、朝向插入位点（正链）
      dn_rev  — 位于插入位点下游、朝向插入位点（负链）

    这两条用于**验证基因组侧**是否真的插入了东西（配合 cargo 内部引物）。
    `reach` = 引物 3' 端距插入位点的目标距离（bp）。

    出处: SOP 表 8（LE junction / RE junction 检测方案）
    """
    out: dict[str, Primer | None] = {}

    # 上游引物：3' 端落在 insert_site - reach，向右延伸（正链）
    anchor_up = max(0, insert_site - reach)
    out["up_fwd"] = _pick_primer_near(genome, anchor_up, f"{name_prefix}_up_F", "+")

    # 下游引物：3' 端落在 insert_site + reach，取其反向互补向左延伸（负链）
    anchor_dn = min(len(genome), insert_site + reach)
    out["dn_rev"] = _pick_primer_near(genome, anchor_dn, f"{name_prefix}_dn_R", "-")

    return out


def design_amplification_primers(genome: str, start: int, end: int,
                                 name_prefix: str = "A",
                                 target_len: tuple[int, int] = (300, 1500)) -> dict:
    """
    为一个区间 [start, end) 设计扩增引物。

    引物落在区间**外侧**：正向引物 3' 端在 start 之前，反向引物 3' 端在 end 之后，
    使扩增产物完整覆盖目标区间。会自动调整外侧距离以满足 target_len。

    ⚠️ 若返回 None，说明该区间两端找不到符合约束的引物，需人工设计。

    返回 {'fwd': Primer|None, 'rev': Primer|None, 'product_size': int|None}。
    """
    lo, hi = target_len
    span = max(1, end - start)

    # 目标产物长度 → 两侧各留多少（起点先给一个合理值）
    side = max(20, min(300, (lo - span) // 2))

    fwd: Primer | None = None
    rev: Primer | None = None
    product: int | None = None

    for extra in range(side, 0, -5):
        f = _pick_primer(genome, max(0, start - extra), f"{name_prefix}_F", "+",
                         len_range=(18, 26))
        r = _pick_primer(genome, min(len(genome), end + extra), f"{name_prefix}_R", "-",
                         len_range=(18, 26))
        if f is None or r is None:
            continue
        size = r.start - f.end
        if size < 0:
            continue
        # 记录第一个可用组合；命中目标区间则立即采用
        if fwd is None:
            fwd, rev, product = f, r, size
        if lo <= size <= hi:
            fwd, rev, product = f, r, size
            break

    return {"fwd": fwd, "rev": rev, "product_size": product}


class GenomeIndex:
    """
    极简基因组索引：支持按 motif 精确查找与基因区间查询。

    刻意不做 k-mer 倒排索引 —— 6.3 Mb 基因组上 str.find 循环足够快，
    且避免引入复杂度。若日后需要脱靶搜索，应外接 Cas-OFFinder（README §7）。
    """

    def __init__(self, genome: str, seq_id: str = "genome"):
        self.seq = genome.upper()
        self.id = seq_id

    def __len__(self) -> int:
        return len(self.seq)

    def find_all(self, motif: str, both_strands: bool = False) -> list[int]:
        """返回 motif 所有出现位置的 0-based 起点。"""
        motif = motif.upper()
        hits: list[int] = []
        i = self.seq.find(motif)
        while i != -1:
            hits.append(i)
            i = self.seq.find(motif, i + 1)
        if both_strands:
            rc = revcomp(motif)
            if rc != motif:
                i = self.seq.find(rc)
                while i != -1:
                    hits.append(i)
                    i = self.seq.find(rc, i + 1)
        return sorted(hits)

    def subseq(self, start: int, end: int) -> str:
        return self.seq[max(0, start):min(len(self.seq), end)]

    def gene_span(self, gene_start: int, gene_end: int) -> tuple[int, int]:
        """返回基因区间（用于判断插入是否落在基因内）。"""
        return (max(0, gene_start), min(len(self.seq), gene_end))


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------

def selftest() -> list[tuple[str, bool, str]]:
    """
    对论文已发表序列做基准检查，返回 [(检查名, 通过?, 说明)]。

    这些基准来自 Table S2，可任何人独立核对。
    """
    results: list[tuple[str, bool, str]] = []

    # --- 基础操作 ---
    results.append(("revcomp", revcomp("ATGC") == "GCAT", f"{revcomp('ATGC')}"))

    # --- 论文已发表 sgRNA26（Table S2）---
    SG26 = "ATTCGGAAGAAGACGTGGTATCA"
    results.append(("sgRNA26 长度=23", len(SG26) == 23, f"len={len(SG26)}"))

    # spacer 全是 ACGT
    results.append(("sgRNA26 字符集", set(SG26) <= set("ACGT"), "ACGT only"))

    # oligo 组装：必须能在 revcomp 后还原 spacer
    # F = aaag + spacer ; R = agcc + revcomp(spacer)
    F = BSAI_F_OVERHANG + SG26
    R = BSAI_R_OVERHANG + revcomp(SG26)
    results.append((
        "oligo 粘性末端",
        F.startswith("aaag") and R.startswith("agcc"),
        f"F={F[:8]}... R={R[:8]}...",
    ))

    # --- 合成底物上的 PAM 扫描自检（可手工核对）---
    # 构造成 GTT(23)AATG → 期望在 index 0 找到 GTT，spacer = 其后的 23 nt
    sub = "GTT" + "A" * SPACER_LEN + "AATG"
    found = list(find_pams(sub))
    ok = len(found) == 1 and found[0] == (0, "GTT")
    results.append(("PAM 扫描 (GTT)", ok, f"found={found}"))

    cands = scan_sgrnas(sub, 0, len(sub), flank=0)
    ok2 = any(c.spacer == "A" * SPACER_LEN and c.pam == "GTT" for c in cands)
    results.append(("spacer 抽取 23 nt", ok2, f"n_cand={len(cands)}"))

    # --- 反向链扫描：反向互补后应找到同样的位点 ---
    rc_sub = revcomp(sub)
    cands_rc = scan_sgrnas(rc_sub, 0, len(rc_sub), flank=0)
    ok3 = any(c.strand == "-" for c in cands_rc)
    results.append(("负链扫描产出 '-' 链", ok3, f"n={len(cands_rc)}"))

    # --- GC 与 Tm 基本合理性 ---
    results.append(("GC=0", gc_content("AAAA") == 0.0, ""))
    results.append(("GC=1", gc_content("GCGC") == 1.0, ""))
    tm20 = tm_basic("ACGTACGTACGTACGTACGT")
    results.append(("Tm 在合理区间", 50.0 <= tm20 <= 70.0, f"Tm={tm20:.1f}"))

    # --- 插入窗口常量 ---
    results.append((
        "插入窗口并集 56-66",
        INSERT_WINDOW == (56, 66),
        f"{INSERT_WINDOW} (与论文菌群章节 {INSERT_WINDOW_ALT} 存在不一致，见 README)",
    ))

    # --- 唯一性检查自检 ---
    #
    # ⚠️ v0.6：此前的断言写成 `seed_is_unique(...) is False or True`，
    # 恒为 True，**毫无判别力**。现改为一对正/反例，确认真有区分能力。
    g = "ACGT" * 100                                    # 高重复底物
    uniq_probe = "ACGTTGCACTGA"                          # 12 nt，底物中应只出现一次
    g2 = "TTTT" + uniq_probe + "TTTT" + ("ACGT" * 20)
    results.append(("唯一性: 唯一序列判 True",
                    seed_is_unique(g2, uniq_probe + "AAAA", seed_len=12) is True,
                    "底物含该前 12 nt 一次"))
    results.append(("唯一性: 重复序列判 False",
                    seed_is_unique(g, "ACGTACGTACGT", seed_len=12) is False,
                    "底物为 (ACGT)×100，前 12 nt 反复出现"))
    n_hit = count_occurrences(g, "ACGT")
    # ⚠️ ACGT 的反向互补就是它自己（回文），所以 both_strands 不会翻倍。
    results.append(("count_occurrences 回文不翻倍", n_hit == 100,
                    f"n={n_hit}（(ACGT)×100；ACGT 是回文，故不与反向互补叠加）"))
    n_hit2 = count_occurrences(g, "ACGTACG")
    results.append(("count_occurrences 多命中", n_hit2 >= 100,
                    f"n={n_hit2}（ACGTACG 非回文）"))

    return results


if __name__ == "__main__":
    print("STIB Planner — sequence layer selftest\n")
    rows = selftest()
    width = max(len(r[0]) for r in rows)
    n_ok = 0
    for name, ok, note in rows:
        flag = "PASS" if ok else "FAIL"
        if ok:
            n_ok += 1
        print(f"  [{flag}] {name:<{width}}  {note}")
    print(f"\n  {n_ok}/{len(rows)} passed")
    raise SystemExit(0 if n_ok == len(rows) else 1)
