"""
STIB Planner — 构建设计输出（SnapGene .dna / GenBank / FASTA）。

═══ 为什么要有这个模块 ═══

用户拿到一堆引物和坐标后，还需要在 SnapGene 里**看到**整个构建：
cargo 在哪、LE/RE 边界在哪、引物结合在哪。手工搭这个图很费时且容易搭错。

本模块把组装方案直接输出成可打开的图谱文件。

═══ 关于 .dna 格式的诚实说明 ═══

SnapGene 的 `.dna` 是**私有二进制格式**。本实现依据 Biopython `Bio.SeqIO.SnapGeneIO`
（Damien Goutte-Gattat，2017-2019）的逆向工程结果：

    文件 = 若干 packet 顺序拼接
    packet = [1 字节类型][4 字节大端长度][数据]

    0x09 cookie    "SnapGene"(8B) + seq_type(2B) + exp_version(2B) + imp_version(2B)
    0x00 DNA       flags(1B，bit0=环状) + ASCII 序列（**明文，不压缩**）
    0x06 Notes     XML: LastModified / AccessionNumber / Comments
    0x0A Features  XML: Feature(type,directionality,name) > Segment(range) / Q > V
    0x05 Primers   XML: Primer(name) > BindingSite(range)

⚠️ **本实现做了往返校验，但未在真实 SnapGene / 其它阅读器中打开验证过。**
   校验方式：用独立的解析器（复刻 Biopython 的读取逻辑）读回自己写的文件，
   逐项比对序列、特征数、引物数、坐标与链向。这能排除编码错误，
   但**不能证明 SnapGene 一定接受**。
   因此**同时输出 `.gb`（GenBank）作为保底** —— SnapGene 能原生打开 .gb，
   且任何其它工具都能读。
"""

from __future__ import annotations

import io
import struct
import xml.sax.saxutils as saxutils
from dataclasses import dataclass, field
from datetime import date
from typing import Any

#: packet 类型
PACKET_DNA = 0x00
PACKET_PRIMERS = 0x05
PACKET_NOTES = 0x06
PACKET_COOKIE = 0x09
PACKET_FEATURES = 0x0A

#: cookie 里的 seq_type：1 = DNA
SEQ_TYPE_DNA = 1
#: cookie 里的版本号。Biopython 的解析器不校验具体数值，只校验 magic；
#: 这里填常见值（SnapGene 自身会做宽松兼容）。
EXP_VERSION = 1
IMP_VERSION = 14

#: 引物 BindingSite 的存储偏移。Biopython 读取引物时对 start 与 end **各** +1
#: （见 SnapGeneIO._parse_location 的 is_primer 分支）。
#: 写出端据此取 s = start、e = end - 1（推导见 _primers_xml）。
PRIMER_COORD_SHIFT = -1


def _xml_escape(s: str) -> str:
    return saxutils.escape(str(s), {'"': "&quot;"})


# ---------------------------------------------------------------------------
# 输出模型
# ---------------------------------------------------------------------------

@dataclass
class OutFeature:
    """要写进图谱文件的一个特征。"""

    name: str
    start: int                  # 0-based 半开
    end: int
    strand: int = 1             # +1 / -1
    kind: str = "misc_feature"
    qualifiers: dict[str, str] = field(default_factory=dict)
    color: str = ""             # SnapGene 颜色（可空）

    def clamped(self, n: int) -> "OutFeature":
        return OutFeature(
            name=self.name, start=max(0, min(self.start, n)),
            end=max(0, min(self.end, n)), strand=self.strand,
            kind=self.kind, qualifiers=dict(self.qualifiers), color=self.color,
        )


@dataclass
class OutPrimer:
    """要写进图谱文件的一条引物。"""

    name: str
    seq: str
    start: int                  # 0-based，引物 5' 端在模板上的位置
    end: int
    strand: int = 1
    tm: float = 0.0

    def clamped(self, n: int) -> "OutPrimer":
        return OutPrimer(
            name=self.name, seq=self.seq,
            start=max(0, min(self.start, n)), end=max(0, min(self.end, n)),
            strand=self.strand, tm=self.tm,
        )


@dataclass
class ConstructMap:
    """一个构建的图谱内容。"""

    name: str
    seq: str
    topology: str = "linear"        # "linear" / "circular"
    features: list[OutFeature] = field(default_factory=list)
    primers: list[OutPrimer] = field(default_factory=list)
    accession: str = ""
    comments: str = ""

    def validate(self) -> list[str]:
        """返回问题清单（空 = 干净）。"""
        problems: list[str] = []
        n = len(self.seq)
        if n == 0:
            problems.append("序列为空")
        bad = sorted({c for c in self.seq.upper() if c not in "ACGT"})
        if bad:
            problems.append(f"序列含非 ACGT 字符: {''.join(bad)}")
        for f in self.features:
            if not (0 <= f.start < f.end <= n):
                problems.append(
                    f"特征 {f.name!r} 坐标越界: [{f.start}, {f.end}) vs 长度 {n}")
        for p in self.primers:
            if not (0 <= p.start < p.end <= n):
                problems.append(
                    f"引物 {p.name!r} 坐标越界: [{p.start}, {p.end}) vs 长度 {n}")
        return problems


# ---------------------------------------------------------------------------
# packet 组装
# ---------------------------------------------------------------------------

def _packet(ptype: int, payload: bytes) -> bytes:
    return struct.pack(">BI", ptype, len(payload)) + payload


def _cookie_payload() -> bytes:
    return struct.pack(">8sHHH", b"SnapGene", SEQ_TYPE_DNA, EXP_VERSION, IMP_VERSION)


def _dna_payload(seq: str, topology: str) -> bytes:
    flags = 0x01 if topology == "circular" else 0x00
    return struct.pack(">B", flags) + seq.upper().encode("ASCII")


def _notes_xml(c: ConstructMap) -> str:
    today = date.today().strftime("%Y.%m.%d")
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        "<Notes>",
        "<Type>Synthetic</Type>",
        f"<LastModified>{today}</LastModified>",
    ]
    if c.accession:
        parts.append(f"<AccessionNumber>{_xml_escape(c.accession)}</AccessionNumber>")
    comment = c.comments or c.name
    parts.append(f"<Comments>{_xml_escape(comment)}</Comments>")
    parts.append("</Notes>")
    return "".join(parts)


def _features_xml(c: ConstructMap) -> str:
    n = len(c.seq)
    buf = io.StringIO()
    buf.write('<?xml version="1.0" encoding="UTF-8"?>')
    buf.write("<Features>")
    for f in c.features:
        f = f.clamped(n)
        directionality = 2 if f.strand < 0 else 1
        buf.write(
            f'<Feature name="{_xml_escape(f.name)}" type="{_xml_escape(f.kind)}" '
            f'directionality="{directionality}">')
        # SnapGene 坐标是 1-based 闭区间：0-based [s, e) → 1-based s+1 .. e
        buf.write(f'<Segment range="{f.start + 1}-{f.end}" type="standard" />')
        for k, v in f.qualifiers.items():
            if v is None or v == "":
                continue
            buf.write(f'<Q name="{_xml_escape(k)}">'
                      f'<V text="{_xml_escape(v)}" /></Q>')
        buf.write("</Feature>")
    buf.write("</Features>")
    return buf.getvalue()


def _primers_xml(c: ConstructMap) -> str:
    n = len(c.seq)
    buf = io.StringIO()
    buf.write('<?xml version="1.0" encoding="UTF-8"?>')
    buf.write("<Primers>")
    # Biopython 说明：Primers 包里未被 HybridizationParams 满足的结合位点会被
    # SnapGene 丢弃。给出宽松参数，确保我们写入的位点都会显示。
    buf.write('<HybridizationParams minContinuousMatchLen="10" '
              'minMeltingTemperature="0" />')
    for p in c.primers:
        p = p.clamped(n)
        buf.write(f'<Primer name="{_xml_escape(p.name)}" '
                  f'sequence="{_xml_escape(p.seq)}">')
        # ⚠️ 引物 BindingSite 坐标换算是本模块最容易错的一处。
        #    推导（而非试错）如下：
        #
        # 内部约定：引物占据 **半开区间 [start, end)**（0-based）。
        # SnapGene 存储：`range="s-e"`，1-based **闭区间**。
        # Biopython 读取（SnapGeneIO._parse_location，is_primer=True）：
        #
        #     start = s - 1        # 1-based 闭 → 0-based
        #     start += 1           # 引物专属偏移
        #     end   += 1
        #  ⇒ start_read = s
        #  ⇒ end_read   = e + 1
        #
        # 要求读回与内部一致（start_read == start、end_read == end）：
        #
        #     s = start            e = end - 1
        #
        # 例：内部 [41, 85) → 写出 range="41-84"。
        # 注：这与「基因特征的 1-based 闭区间」（42..85）不同 ——
        # SnapGene 对引物就是再小 1，这是 Biopython 注释所称的
        # "shifted by -1 for some reasons"。
        s = p.start
        e = p.end - 1
        buf.write(f'<BindingSite range="{s}-{e}" '
                  f'bound="1" strand="{1 if p.strand > 0 else 2}" />')
        buf.write("</Primer>")
    buf.write("</Primers>")
    return buf.getvalue()


def write_snapgene(c: ConstructMap) -> bytes:
    """
    生成 SnapGene `.dna` 字节流。

    包顺序：cookie(0x09) → DNA(0x00) → Notes(0x06) → Features(0x0A) → Primers(0x05)。
    cookie 必须是第一个包，否则 SnapGene 不认。
    """
    problems = c.validate()
    if problems:
        raise ValueError("构建图谱有误，拒绝写出：\n  - " + "\n  - ".join(problems))

    out = bytearray()
    out += _packet(PACKET_COOKIE, _cookie_payload())
    out += _packet(PACKET_DNA, _dna_payload(c.seq, c.topology))
    out += _packet(PACKET_NOTES, _notes_xml(c).encode("UTF-8"))
    if c.features:
        out += _packet(PACKET_FEATURES, _features_xml(c).encode("UTF-8"))
    if c.primers:
        out += _packet(PACKET_PRIMERS, _primers_xml(c).encode("UTF-8"))
    return bytes(out)


def write_snapgene_file(path: str, c: ConstructMap) -> None:
    with open(path, "wb") as fh:
        fh.write(write_snapgene(c))


# ---------------------------------------------------------------------------
# 独立解析器（用于往返校验；刻意不复用写出代码逻辑）
# ---------------------------------------------------------------------------

def read_snapgene(data: bytes) -> dict:
    """
    独立解析 SnapGene 字节流 —— 复刻 Biopython SnapGeneIO 的读取逻辑。

    用途：**往返校验**。写出后读回，逐项比对。
    刻意与写出代码分开实现，避免"用同一套错误逻辑自证"。
    """
    import xml.dom.minidom as minidom

    pos = 0
    n = len(data)
    cookie = None
    seq = None
    topology = None
    features: list[dict] = []
    primers: list[dict] = []
    notes: dict = {}

    first = True
    while pos < n:
        if pos + 5 > n:
            raise ValueError(f"packet 头不完整 @ {pos}")
        ptype = data[pos]
        length = struct.unpack(">I", data[pos + 1:pos + 5])[0]
        payload = data[pos + 5:pos + 5 + length]
        if len(payload) != length:
            raise ValueError(f"packet 数据不完整 @ {pos}")
        pos += 5 + length

        if first and ptype != PACKET_COOKIE:
            raise ValueError(f"第一个 packet 不是 cookie(0x09) 而是 0x{ptype:02x}")
        first = False

        if ptype == PACKET_COOKIE:
            magic, st, ev, iv = struct.unpack(">8sHHH", payload)
            cookie = {"magic": magic.decode("ASCII"), "seq_type": st,
                      "exp_version": ev, "imp_version": iv}
            if cookie["magic"] != "SnapGene":
                raise ValueError("cookie 不是 'SnapGene'")

        elif ptype == PACKET_DNA:
            flags = payload[0]
            seq = payload[1:].decode("ASCII")
            topology = "circular" if (flags & 0x01) else "linear"

        elif ptype == PACKET_NOTES:
            xml = minidom.parseString(payload.decode("UTF-8"))
            for tag in ("Type", "LastModified", "AccessionNumber", "Comments"):
                els = xml.getElementsByTagName(tag)
                if els and els[0].childNodes:
                    notes[tag] = els[0].childNodes[0].nodeValue

        elif ptype == PACKET_FEATURES:
            xml = minidom.parseString(payload.decode("UTF-8"))
            for f in xml.getElementsByTagName("Feature"):
                seg = f.getElementsByTagName("Segment")
                rng = seg[0].getAttribute("range") if seg else ""
                s, e = (int(x) for x in rng.split("-"))
                directionality = int(f.getAttribute("directionality") or "1")
                quals = {}
                for q in f.getElementsByTagName("Q"):
                    vs = [v.getAttribute("text") for v in q.getElementsByTagName("V")]
                    quals[q.getAttribute("name")] = vs
                features.append({
                    "name": f.getAttribute("name"),
                    "type": f.getAttribute("type"),
                    "strand": -1 if directionality == 2 else 1,
                    # 1-based 闭 → 0-based 半开
                    "start": s - 1, "end": e, "qualifiers": quals,
                })

        elif ptype == PACKET_PRIMERS:
            xml = minidom.parseString(payload.decode("UTF-8"))
            for p in xml.getElementsByTagName("Primer"):
                bs = p.getElementsByTagName("BindingSite")
                rng = bs[0].getAttribute("range") if bs else ""
                s, e = (int(x) for x in rng.split("-"))
                primers.append({
                    "name": p.getAttribute("name"),
                    "sequence": p.getAttribute("sequence"),
                    # 严格复刻 Biopython _parse_location(is_primer=True)：
                    #   start = s - 1        (1-based → 0-based)
                    #   start += 1 / end += 1  (引物专属偏移，各一次)
                    # ⇒ start_read = s ; end_read = e + 1
                    "start": s,
                    "end": e + 1,
                    "strand": 1 if bs and bs[0].getAttribute("strand") == "1" else -1,
                })

    return {
        "cookie": cookie, "seq": seq, "topology": topology,
        "features": features, "primers": primers, "notes": notes,
    }


def roundtrip_check(c: ConstructMap, verbose: bool = False) -> tuple[bool, list[str]]:
    """
    写出 → 独立读回 → 逐项比对。

    返回 (是否一致, 差异清单)。
    """
    diffs: list[str] = []
    try:
        blob = write_snapgene(c)
    except Exception as exc:  # noqa: BLE001
        return False, [f"写出失败: {exc!r}"]

    try:
        back = read_snapgene(blob)
    except Exception as exc:  # noqa: BLE001
        return False, [f"读回失败: {exc!r}"]

    if back["cookie"] is None or back["cookie"]["magic"] != "SnapGene":
        diffs.append("cookie 缺失或错误")
    if back["seq"] != c.seq.upper():
        diffs.append(f"序列不一致（写出 {len(c.seq)} bp，读回 "
                     f"{len(back['seq']) if back['seq'] else 0} bp）")
    if back["topology"] != c.topology:
        diffs.append(f"拓扑不一致: {back['topology']} != {c.topology}")
    if len(back["features"]) != len(c.features):
        diffs.append(f"特征数不一致: {len(back['features'])} != {len(c.features)}")

    n = len(c.seq)
    for want, got in zip(c.features, back["features"]):
        w = want.clamped(n)
        if (w.start, w.end) != (got["start"], got["end"]):
            diffs.append(f"特征 {w.name!r} 坐标不一致: "
                         f"({w.start},{w.end}) vs ({got['start']},{got['end']})")
        if w.name != got["name"]:
            diffs.append(f"特征名不一致: {w.name!r} vs {got['name']!r}")
        if w.strand != got["strand"]:
            diffs.append(f"特征 {w.name!r} 链向不一致: {w.strand} vs {got['strand']}")

    if len(back["primers"]) != len(c.primers):
        diffs.append(f"引物数不一致: {len(back['primers'])} != {len(c.primers)}")
    for want, got in zip(c.primers, back["primers"]):
        w = want.clamped(n)
        # 写出端已按读取逻辑的逆变换处理，故此处直接比对
        if (w.start, w.end) != (got["start"], got["end"]):
            diffs.append(f"引物 {w.name!r} 坐标不一致: "
                         f"({w.start},{w.end}) vs ({got['start']},{got['end']})")
        if w.seq != got["sequence"]:
            diffs.append(f"引物 {w.name!r} 序列不一致")

    if verbose:
        print(f"  往返：{len(back['features'])} 特征 / {len(back['primers'])} 引物 / "
              f"{len(back['seq'] or '')} bp")
    return (not diffs), diffs


# ---------------------------------------------------------------------------
# GenBank 输出（保底格式：SnapGene 与所有工具都能打开）
# ---------------------------------------------------------------------------

def write_genbank(c: ConstructMap) -> str:
    """把构建图谱写成 GenBank 文本（SnapGene 原生支持导入）。"""
    problems = c.validate()
    if problems:
        raise ValueError("构建图谱有误，拒绝写出：\n  - " + "\n  - ".join(problems))

    n = len(c.seq)
    today = date.today().strftime("%d-%b-%Y").upper()
    loc = f"{c.name[:16]:<16}"
    out: list[str] = []
    out.append(f"LOCUS       {loc} {n:>11} bp    DNA     "
               f"{'circular' if c.topology == 'circular' else 'linear'}   "
               f"SYN {today}")
    out.append(f"DEFINITION  {c.comments or c.name}.")
    if c.accession:
        out.append(f"ACCESSION   {c.accession}")
    out.append("FEATURES             Location/Qualifiers")
    out.append("     source          1..%d" % n)
    out.append('                     /organism="synthetic DNA construct"')
    out.append('                     /mol_type="other DNA"')

    def q(key: str, val: str) -> None:
        val = val.replace("\n", " ")
        out.append(f'                     /{key}="{val}"')

    for f in c.features:
        f = f.clamped(n)
        rng = f"{f.start + 1}..{f.end}"
        if f.strand < 0:
            rng = f"complement({rng})"
        out.append(f"     {f.kind[:15]:<15} {rng}")
        if f.name:
            q("label", f.name)
        for k, v in f.qualifiers.items():
            if v:
                q(k, v)
    for p in c.primers:
        p = p.clamped(n)
        rng = f"{p.start + 1}..{p.end}"
        if p.strand < 0:
            rng = f"complement({rng})"
        out.append(f"     primer_bind     {rng}")
        q("label", p.name)
        q("sequence", p.seq)

    out.append("ORIGIN")
    s = c.seq.lower()
    for i in range(0, n, 60):
        chunk = s[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        out.append(f"{i + 1:>9} {groups}")
    out.append("//")
    return "\n".join(out) + "\n"


def write_fasta(c: ConstructMap, width: int = 70) -> str:
    lines = [f">{c.name} {c.comments}".rstrip()]
    s = c.seq.upper()
    for i in range(0, len(s), width):
        lines.append(s[i:i + width])
    return "\n".join(lines) + "\n"
