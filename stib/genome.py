"""
STIB Planner — 基因组区域解析。

═══ 用户的诉求 ═══

「针对这个拟杆菌基因组文件的 A 基因设计插入失活的引物」——
一句话里只给了**基因组文件 + 基因名**，工具必须自己找到坐标。

═══ 为什么不能只支持 FASTA ═══

**纯 FASTA 只有序列、没有注释**，从序列本身无法知道哪一段是 A 基因。
因此：

  GenBank（.gb/.gbk/.gbff）  → 带 FEATURES 表，**可以按基因名定位** ✅
  FASTA（.fa/.fasta/.fna）   → 只有序列，必须另给坐标或序列片段 ⚠️

这不是实现偷懒，是**信息论限制**：基因名不在碱基序列里。

═══ 本模块提供 ═══

  parse_genbank(path)          → GenBankRecord（序列 + 特征表）
  load_genome_any(path)        → 自动识别格式
  GenBankRecord.find_gene(name) → 按基因名/位点名/locus_tag 查坐标
  resolve_region(...)          → 统一的区域解析入口（基因名 / 坐标 / 序列片段）
"""

from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass, field
from typing import Any

# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class Feature:
    """GenBank FEATURES 表中的一个特征。"""

    key: str                     # CDS / tRNA / rRNA / ncRNA / …
    start: int                   # 0-based，闭区间起点
    end: int                     # 0-based，**开区间**终点（= GenBank 的 end）
    strand: str                  # '+' / '-'
    qualifiers: dict[str, str] = field(default_factory=dict)
    #: 一个 CDS 可以有多个位置区间（join）；这里保存全部，start/end 是总的起止
    parts: list[tuple[int, int]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    def name(self) -> str:
        """尽量取一个人类可读的名字。"""
        for k in ("gene", "locus_tag", "product", "note"):
            v = self.qualifiers.get(k)
            if v:
                return v
        return f"{self.key}@{self.start}-{self.end}"

    def aliases(self) -> set[str]:
        """这个特征所有可用于查询的名字（小写）。"""
        out: set[str] = set()
        for k in ("gene", "locus_tag", "old_locus_tag", "product"):
            v = self.qualifiers.get(k)
            if v:
                out.add(v.strip().lower())
                # product 常写成 "thymidine kinase" → 也允许按整个串查
        return out

    def as_dict(self) -> dict:
        return {
            "key": self.key, "start": self.start, "end": self.end,
            "strand": self.strand, "name": self.name(),
            "qualifiers": self.qualifiers, "parts": self.parts,
        }


@dataclass
class GenBankRecord:
    """一个 GenBank 记录：序列 + 特征表。"""

    locus: str
    seq: str
    features: list[Feature]
    definition: str = ""
    accession: str = ""
    source: str = ""

    def __len__(self) -> int:
        return len(self.seq)

    # -- 查询 ---------------------------------------------------------------

    def coding_features(self) -> list[Feature]:
        return [f for f in self.features if f.key in ("CDS", "tRNA", "rRNA", "ncRNA")]

    def find_gene(self, name: str, exact: bool = True) -> list[Feature]:
        """
        按名字查特征。不区分大小写。

        查询顺序：gene → locus_tag → old_locus_tag → product。
        exact=False 时退化为子串匹配（用于名字记不全的情况）。
        """
        q = name.strip().lower()
        if not q:
            return []
        feats = self.coding_features() or self.features

        hits: list[Feature] = []
        # 1) 精确匹配 gene / locus_tag / product
        for f in feats:
            vals = {f.qualifiers.get(k, "").strip().lower()
                    for k in ("gene", "locus_tag", "old_locus_tag", "product")}
            vals.discard("")
            if q in vals:
                hits.append(f)
        if hits:
            return hits
        # 2) 子串匹配
        if not exact:
            for f in feats:
                for k in ("gene", "locus_tag", "old_locus_tag", "product"):
                    v = f.qualifiers.get(k, "").strip().lower()
                    if v and q in v:
                        hits.append(f)
                        break
        return hits

    def search_genes(self, fragment: str, limit: int = 20) -> list[Feature]:
        """列出名字中包含该片段的基因（用于「名字记不清」时给候选）。"""
        q = fragment.strip().lower()
        out: list[Feature] = []
        for f in self.coding_features():
            for k in ("gene", "locus_tag", "product"):
                v = f.qualifiers.get(k, "").strip().lower()
                if v and q in v:
                    out.append(f)
                    break
            if len(out) >= limit:
                break
        return out

    def summary(self) -> dict:
        from collections import Counter
        c = Counter(f.key for f in self.features)
        return {
            "locus": self.locus,
            "accession": self.accession,
            "length": len(self.seq),
            "definition": self.definition,
            "feature_counts": dict(c),
            "n_genes_with_name": sum(
                1 for f in self.coding_features() if f.qualifiers.get("gene")),
        }


# ---------------------------------------------------------------------------
# GenBank 解析
# ---------------------------------------------------------------------------

_LOCUS_RE = re.compile(r"^LOCUS\s+(\S+)")
_DEF_RE = re.compile(r"^DEFINITION\s+(.*)$")
_ACC_RE = re.compile(r"^ACCESSION\s+(\S+)")
#: 位置串，如 1234..5678 / complement(1234..5678) / join(1..9,20..30)
_POS_RE = re.compile(r"(\d+)\s*\.\.\s*([<>]?\d+)")


def _open_text(path: str):
    """
    打开文本文件。

    用 **utf-8-sig**：有 BOM 时剥掉、无 BOM 时等同于 utf-8。
    Windows 工具导出的 GenBank/FASTA 常带 BOM；若按 utf-8 读，
    首行会变成 '\\ufeffLOCUS...'，导致 startswith("LOCUS") 失败、
    解析**静默返回空结果**（实测踩到）。
    """
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return open(path, "r", encoding="utf-8-sig", errors="replace")


def _parse_location(loc: str, offset: int = 0) -> tuple[list[tuple[int, int]], str]:
    """
    解析 GenBank 位置串 → ([(start, end), ...] 0-based 半开区间, strand)。

    支持 complement(...)、join(...)、多行续行、< > 部分边界。
    """
    strand = "-" if "complement" in loc else "+"
    parts: list[tuple[int, int]] = []
    for m in _POS_RE.finditer(loc):
        a = int(m.group(1).lstrip("<>"))
        b = int(m.group(2).lstrip("<>"))
        if b < a:
            a, b = b, a
        # GenBank 是 1-based 闭区间 → 0-based 半开
        parts.append((a - 1 - offset, b - offset))
    if not parts:
        return [], strand
    parts.sort()
    return parts, strand


def _parse_qualifiers(lines: list[str]) -> dict[str, str]:
    """
    把 FEATURES 块的续行解析成 qualifiers。

    GenBank 的 qualifier 形如：
        /gene="tdk"
        /product="thymidine kinase"
        /note="multi
               line"
        /translation="MKK..."
    """
    text = "\n".join(lines)
    quals: dict[str, str] = {}
    # /key="value"（value 可跨行，到下一个 /key 或结束）
    for m in re.finditer(r'/(\w+)\s*=\s*"((?:[^"\\]|\\.)*)"', text, re.S):
        key = m.group(1)
        val = m.group(2).replace("\\\n", "").replace("\n", " ")
        val = re.sub(r"\s+", " ", val).strip()
        quals.setdefault(key, val)
    # /key（无值，如 /pseudo）
    for m in re.finditer(r"/(\w+)(?!\s*=)", text):
        quals.setdefault(m.group(1), "")
    # translation 单独处理：需要去掉所有空白
    if "translation" in quals:
        quals["translation"] = re.sub(r"\s+", "", quals["translation"])
    return quals


def parse_genbank(path: str, want_features: bool = True) -> list[GenBankRecord]:
    """
    解析 GenBank 文件为若干记录。

    只做我们需要的部分：LOCUS / ACCESSION / DEFINITION / FEATURES / ORIGIN。
    （不解析 CONTIG 组装式记录 —— 那类没有 ORIGIN 序列。）

    ⚠️ 必须容忍 **UTF-8 BOM**。实测：从 Windows 工具导出的 GenBank 文件首行
    常带 `\\ufeff`，而 `"\\ufeffLOCUS...".startswith("LOCUS")` 为 False ——
    整个解析会静默返回 0 条记录（不报错），用户看到的是"文件读不出来"。
    """
    with _open_text(path) as fh:
        raw_lines = fh.read().splitlines()

    # 剥掉 BOM（只在首行可能出现）
    if raw_lines and raw_lines[0].startswith("\ufeff"):
        raw_lines[0] = raw_lines[0].lstrip("\ufeff")

    records: list[GenBankRecord] = []
    i = 0
    n = len(raw_lines)

    while i < n:
        line = raw_lines[i]
        if not line.startswith("LOCUS"):
            i += 1
            continue

        m = _LOCUS_RE.match(line)
        locus = m.group(1) if m else "unknown"
        accession = ""
        definition = ""
        seq_chunks: list[str] = []
        features: list[Feature] = []

        i += 1
        # ---- 头部 ----
        while i < n and not raw_lines[i].startswith("FEATURES") \
                and not raw_lines[i].startswith("ORIGIN") \
                and not raw_lines[i].startswith("LOCUS"):
            l = raw_lines[i]
            if l.startswith("ACCESSION"):
                mm = _ACC_RE.match(l)
                if mm:
                    accession = mm.group(1)
            elif l.startswith("DEFINITION"):
                mm = _DEF_RE.match(l)
                if mm:
                    definition = mm.group(1).strip()
            elif l.startswith(" " * 5) and definition and not definition.endswith("."):
                definition += " " + l.strip()          # DEFINITION 续行
            i += 1

        # ---- FEATURES ----
        if want_features and i < n and raw_lines[i].startswith("FEATURES"):
            i += 1
            cur_key: str | None = None
            cur_loc: list[str] = []
            cur_qual_lines: list[str] = []

            def flush() -> None:
                if cur_key is None:
                    return
                loc = " ".join(cur_loc)
                parts, strand = _parse_location(loc)
                if parts:
                    features.append(Feature(
                        key=cur_key,
                        start=min(p[0] for p in parts),
                        end=max(p[1] for p in parts),
                        strand=strand,
                        qualifiers=_parse_qualifiers(cur_qual_lines),
                        parts=parts,
                    ))

            while i < n:
                l = raw_lines[i]
                if l.startswith("ORIGIN") or l.startswith("//") or l.startswith("LOCUS"):
                    break
                if not l.strip():
                    i += 1
                    continue
                # 特征键在 5..20 列
                m2 = re.match(r"^ {5}(\S+)\s+(.*)$", l)
                if m2:
                    flush()
                    cur_key = m2.group(1)
                    cur_loc = [m2.group(2)]
                    cur_qual_lines = []
                else:
                    # 续行：位置串续行 or qualifier 续行
                    body = l[21:] if len(l) > 21 else l.strip()
                    if cur_key is not None:
                        if cur_qual_lines or body.lstrip().startswith("/"):
                            cur_qual_lines.append(body)
                        else:
                            cur_loc.append(l.strip())
                i += 1
            flush()

        # ---- ORIGIN ----
        if i < n and raw_lines[i].startswith("ORIGIN"):
            i += 1
            while i < n and not raw_lines[i].startswith("//"):
                seq_chunks.append(re.sub(r"[^A-Za-z]", "", raw_lines[i]))
                i += 1
        # 跳过到记录结尾
        while i < n and not raw_lines[i].startswith("//"):
            i += 1
        i += 1

        seq = "".join(seq_chunks).upper()
        if seq:
            records.append(GenBankRecord(
                locus=locus, seq=seq, features=features,
                definition=definition, accession=accession,
            ))
    return records


# ---------------------------------------------------------------------------
# 统一加载
# ---------------------------------------------------------------------------

#: 识别为 GenBank 的扩展名
GENBANK_EXTS = (".gb", ".gbk", ".gbff", ".genbank")
#: 识别为 FASTA 的扩展名
FASTA_EXTS = (".fa", ".fasta", ".fna", ".fas", ".ffn", ".faa")
#: 不带扩展名的 accession 默认按 FASTA 处理（NC_*.fasta 是最常见的缓存形式）
DEFAULT_EXT = ".fasta"


def detect_format(path: str) -> str:
    """
    返回 "genbank" / "fasta" / "unknown"。

    判定顺序：
      1. 扩展名（.gb/.gbk/.gbff/.genbank → genbank；.fa/.fasta/… → fasta）
      2. 扩展名不认识**或没有** → **按内容嗅探**（首行 LOCUS / >）

    ⚠️ 早期版本此处有一段死代码（第一行 `ext = ...` 被下一行覆盖），
    且扩展名分支命中后就直接返回、没有回退 —— 于是一个被改名成 `.fasta`
    的 GenBank 文件会被判成 FASTA，再按 FASTA 解析就读不出序列。
    """
    base = path
    if base.lower().endswith(".gz"):
        base = os.path.splitext(base)[0]
    ext = os.path.splitext(base)[1].lower()

    if ext in GENBANK_EXTS:
        return "genbank"
    if ext in FASTA_EXTS:
        return "fasta"

    # 未知（或无）扩展名 → 嗅探内容
    try:
        with _open_text(path) as fh:
            for _ in range(500):
                line = fh.readline()
                if not line:
                    break
                if line.startswith("LOCUS"):
                    return "genbank"
                if line.startswith(">"):
                    return "fasta"
                if line.strip():
                    break
    except OSError:
        pass
    return "unknown"


def load_genome_any(path: str) -> tuple[str, str, list[Feature], dict]:
    """
    加载基因组，自动识别 FASTA / GenBank。

    返回 (seq_id, sequence, features, meta)。

    FASTA 的 features 为空 —— **这是格式限制，不是缺陷**：
    纯序列里没有基因名信息。调用方应据此提示用户改用 GenBank
    或显式提供坐标/序列片段。
    """
    fmt = detect_format(path)

    # ── GenBank ──
    if fmt == "genbank":
        recs = parse_genbank(path)
        if not recs:
            # 扩展名说是 GenBank，但解析不出记录 —— 可能其实是 FASTA
            recs_f = _try_fasta(path)
            if recs_f:
                return recs_f
            raise SystemExit(f"ERROR: GenBank 文件里没有可用的 ORIGIN 序列: {path}")
        if len(recs) > 1:
            print(f"提示: 文件含 {len(recs)} 条记录，使用第一条 "
                  f"({recs[0].locus}, {len(recs[0].seq):,} bp)")
        r = recs[0]
        meta = {
            "format": "genbank",
            "locus": r.locus,
            "accession": r.accession,
            "definition": r.definition,
            "n_features": len(r.features),
            "n_cds": len(r.coding_features()),
            "has_annotations": True,
        }
        return r.accession or r.locus, r.seq, r.features, meta

    # ── FASTA（或未知格式）──
    got = _try_fasta(path)
    if got is not None:
        sid, seq, meta = got
        return sid, seq, [], meta

    # FASTA 读不出来 —— **兜底再试一次 GenBank**。
    # 场景：文件被改成 .fasta 扩展名（实为 GenBank 文本）。扩展名优先的判定
    # 会把它当 FASTA，read_fasta 读不到任何 ">" 开头的记录 → 走这里自愈。
    recs = parse_genbank(path)
    if recs:
        print(f"提示: 扩展名是 FASTA 但内容像 GenBank，已按 GenBank 解析（{path}）")
        if len(recs) > 1:
            print(f"      文件含 {len(recs)} 条记录，使用第一条 ({recs[0].locus})")
        r = recs[0]
        return r.accession or r.locus, r.seq, r.features, {
            "format": "genbank", "locus": r.locus, "accession": r.accession,
            "definition": r.definition, "n_features": len(r.features),
            "n_cds": len(r.coding_features()), "has_annotations": True,
            "note": "扩展名与内容不符，已按内容解析",
        }

    raise SystemExit(
        f"ERROR: 无法从文件中读出序列: {path}\n"
        "    已尝试 FASTA 与 GenBank 两种解析。请确认文件格式"
        "（GenBank 应以 LOCUS 开头，FASTA 应以 > 开头）。")


def _try_fasta(path: str):
    """尝试按 FASTA 解析，失败返回 None。"""
    from .sequence import read_fasta
    try:
        seqs = read_fasta(path)
    except OSError:
        return None
    if not seqs:
        return None
    seqs = {k: v for k, v in seqs.items() if v}      # 丢掉空序列
    if not seqs:
        return None
    if len(seqs) > 1:
        print(f"提示: 文件含 {len(seqs)} 条序列，使用第一条")
    sid = next(iter(seqs))
    return sid, seqs[sid], {
        "format": "fasta", "locus": sid,
        "n_features": 0, "n_cds": 0, "has_annotations": False,
    }


# ---------------------------------------------------------------------------
# 统一的区域解析
# ---------------------------------------------------------------------------

class RegionError(ValueError):
    """区域无法解析时抛出，附带可操作的提示。"""


def resolve_region(seq: str, features: list[Feature], *,
                   gene: str | None = None,
                   start: int | None = None,
                   end: int | None = None,
                   motif: str | None = None,
                   pad: int = 500,
                   ) -> tuple[int, int, dict]:
    """
    把用户的输入解析为一个 (start, end) 目标区间。

    优先级：坐标 > 序列片段 > 基因名。

    返回 (start, end, resolver_info)，resolver_info 记录**这次是用什么定位的**，
    以便在报告里如实标注（用户需要知道区间是怎么来的）。

    失败时抛 RegionError，消息里给出**下一步该做什么**，而不是只说"没找到"。
    """
    n = len(seq)

    # ---- 1) 显式坐标 ----
    if start is not None and end is not None:
        s, e = int(start), int(end)
        if s < 0 or e > n or s >= e:
            raise RegionError(
                f"坐标区间非法: [{s:,}, {e:,})；基因组长度 {n:,}。"
                "起点须 ≥0、终点须 ≤基因组长度、且起点 < 终点。")
        return s, e, {"method": "coordinates", "note": f"显式坐标 [{s:,}, {e:,})"}

    # ---- 2) 序列片段 ----
    if motif:
        from .sequence import revcomp
        m = motif.upper().strip()
        hits = [i for i in range(n) if seq.startswith(m, i)]
        strand_note = "正链"
        if not hits:
            rc = revcomp(m)
            hits = [i for i in range(n) if seq.startswith(rc, i)]
            strand_note = "反向互补链"
        if not hits:
            raise RegionError(
                f"序列片段在基因组中未找到（正链与反向互补链都没有）: {m[:40]}…\n"
                "可能原因：(a) 片段属于别的菌株/版本；(b) 复制时带了多余字符；"
                "(c) 片段跨了拼接边界。\n"
                "替代做法：改用 --gene（需 GenBank 注释）或 --start/--end 坐标。")
        if len(hits) > 1:
            print(f"警告: 片段命中 {len(hits)} 处，使用第一处 @ {hits[0]:,}"
                  f"（其余: {', '.join(f'{h:,}' for h in hits[1:6])}）")
        s = hits[0]
        return (max(0, s - pad), min(n, s + len(m) + pad),
                {"method": "motif", "note": f"序列片段定位 @ {s:,}（{strand_note}），"
                                             f"两侧各扩 {pad} bp"})

    # ---- 3) 基因名 ----
    if gene:
        if not features:
            raise RegionError(
                f"无法按基因名「{gene}」定位：当前基因组文件**没有注释**。\n"
                "纯 FASTA 只有碱基序列，基因名不在序列里，因此无法从 FASTA 查到坐标。\n"
                "三个可选做法（任选其一）：\n"
                "  1) 改用 GenBank 格式的基因组（.gb/.gbk/.gbff，自带 FEATURES 注释）"
                "—— 最省事，之后就能直接按基因名设计；\n"
                "  2) 提供基因中的一段序列（--motif），工具会在基因组里定位它；\n"
                "  3) 直接给坐标（--start/--end）。")
        q = gene.strip().lower()
        feats = [f for f in features if f.key in ("CDS", "tRNA", "rRNA", "ncRNA")]
        hits = []
        for f in feats:
            vals = {f.qualifiers.get(k, "").strip().lower()
                    for k in ("gene", "locus_tag", "old_locus_tag")}
            vals.discard("")
            if q in vals:
                hits.append(f)
        info_extra = ""
        if not hits:
            # 退化：任意 qualifier 子串匹配
            soft = [f for f in feats
                    if any(q in f.qualifiers.get(k, "").lower()
                           for k in ("gene", "locus_tag", "old_locus_tag", "product"))]
            if soft:
                names = "、".join(sorted({f.name() for f in soft})[:8])
                raise RegionError(
                    f"未找到名为「{gene}」的基因，但有相近的：{names}\n"
                    "请用上面列出的准确名字重试。")
            raise RegionError(
                f"未找到名为「{gene}」的基因（已在 {len(feats)} 个注释特征中查找 "
                "gene / locus_tag / old_locus_tag）。\n"
                "可能该基因在这个注释里没有 gene 名，只有 locus_tag。"
                "可先列出候选：stib.cli genes --search <关键词>")

        if len(hits) > 1:
            # 同名多拷贝（如多个 rRNA 或重复基因）—— 必须让用户选，不能默默取第一个
            spots = "、".join(f"{f.start:,}–{f.end:,}({f.strand})" for f in hits[:8])
            raise RegionError(
                f"基因名「{gene}」在该基因组中命中 **{len(hits)} 处**：{spots}\n"
                "这通常意味着它是多拷贝基因（如 rRNA、转座酶、重复家族）。\n"
                "**不能自动选第一个** —— 多拷贝基因的插入可能落在任一份拷贝上，"
                "实验结果会不可解释。\n"
                f"请用 --start/--end 指定具体那一处，或先用 "
                f"stib.cli genes --search {gene} 查看全部拷贝坐标。")

        f = hits[0]
        s, e = f.start, f.end
        tag = f.qualifiers.get("gene") or f.qualifiers.get("locus_tag") or gene
        prod = f.qualifiers.get("product", "")
        info_extra = (f"基因 {tag}" + (f"（{prod}）" if prod else "")
                      + f"，CDS {s:,}–{e:,}（{f.strand} 链，{e - s:,} bp）")
        return s, e, {"method": "gene", "gene": tag, "product": prod,
                      "strand": f.strand, "note": info_extra}

    raise RegionError(
        "没有指定目标区间。请提供以下任一项：\n"
        "  --gene NAME      基因名（需 GenBank 格式基因组）\n"
        "  --motif SEQ      基因中的一段序列（FASTA 也可用）\n"
        "  --start/--end    直接坐标")
