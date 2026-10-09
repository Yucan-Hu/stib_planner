"""
STIB Planner — 构建规范的读写。

把 `stib/data/construct_psh43.yml` 里的载体内涵（元件坐标、臂边界、引物、
cargo 槽位）读成可用的对象，并**从载体文件本身取全部注释** ——
输出构建图谱时必须保留骨架原有的一切（sgRNA 盒、酶、启动子、抗性、ori…），
否则用户拿到的图谱就"丢信息"了。

⚠️ 特征平移规则：
   在载体坐标 `insert_at` 处插入长度 `delta` 的序列（delta 可正可负），
   则该位置**之后**的所有特征坐标 +delta；
   跨越插入点的特征会被截断/拆分（本工具会显式告警，不静默处理）。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import yamlmini


@dataclass
class Arm:
    core_start: int
    core_end: int
    extended_start: int
    extended_end: int
    rule: str = ""


@dataclass
class CargoSlot:
    start: int
    end: int
    length: int
    note: str = ""
    current_contents: list[dict] = field(default_factory=list)
    assembly_strategy: str = "UNCONFIRMED"
    #: 「保留 CmR 筛选」时的操作区间 —— 新 cargo 插在 LE 之后、CmR 之前
    insert_between: list[int] | None = None
    #: 可选插入位点（v0.8，依据课题组手绘图谱）
    sites: list["InsertSite"] = field(default_factory=list)


@dataclass
class InsertSite:
    """一个可选的 cargo 插入位点。"""

    id: str
    label: str
    region: tuple[int, int]        # (start, end) 0-based 半开
    #: cargo 应处的链：1 = 正链，-1 = 负链（与 CmR 同向）
    cargo_strand: int = 1
    keep_elements: list[str] = field(default_factory=list)
    #: 整段替换语义：cargo 长度不受区间长度限制（区间内元件全部被替换）
    replace_all: bool = False
    note: str = ""
    orientation_rule: str = ""


@dataclass
class ConstructSpec:
    """psh43 载体的构建规范。"""

    path: str
    meta: dict
    elements: list[dict]
    cargo_slot: CargoSlot
    arms: dict[str, Arm]
    primers: dict[str, dict]
    gibson: dict
    cargo_checks: dict
    open_questions: list[dict]
    raw: dict

    # -- 便利访问 -----------------------------------------------------------
    @property
    def le(self) -> Arm:
        return self.arms["LE"]

    @property
    def re(self) -> Arm:
        return self.arms["RE"]

    def element_named(self, name: str) -> dict | None:
        for e in self.elements:
            if e.get("name") == name:
                return e
        return None

    def cmr_region(self) -> tuple[int, int] | None:
        """CmR + PcfxA 的合并区间（用于「是否保留筛选标记」的判断）。"""
        cmr = self.element_named("CmR")
        pcf = self.element_named("PcfxA")
        if not cmr:
            return None
        a = min(cmr["start"], (pcf or cmr)["start"])
        b = max(cmr["end"], (pcf or cmr)["end"])
        return a, b

    def gap_between_le_and_selection(self) -> tuple[int, int]:
        """
        LE 末端 → 筛选盒起点之间可放新 cargo 的区间。

        这是「保留 CmR」策略下的操作区间（v0.7，依据课题组澄清）。
        """
        if self.cargo_slot.insert_between:
            return self.cargo_slot.insert_between[0], self.cargo_slot.insert_between[1]
        cm = self.cmr_region()
        end = cm[0] if cm else self.cargo_slot.end
        return self.le.core_end, end

    def site_named(self, site_id: str) -> "InsertSite | None":
        for s_ in self.cargo_slot.sites:
            if s_.id == site_id:
                return s_
        return None

    def blockers_in(self, start: int, end: int) -> list[dict]:
        """
        返回 [start, end) 内**不建议删除**的已注释元件（终止子等）。

        用途：这段区间虽然「没有 cargo」，但里面常有转录终止子 ——
        直接替换会失去转录隔离、造成读通（read-through）。
        工具必须把它算出来并提示，而不是默默当成可用空间。
        """
        out: list[dict] = []
        for e in self.elements:
            nm = str(e.get("name", ""))
            if ("terminator" in nm.lower()) or ("终止子" in nm):
                if e.get("start") is not None and e.get("end") is not None \
                        and int(e["start"]) >= start and int(e["end"]) <= end:
                    out.append(e)
        return out

    def free_gaps(self, start: int, end: int) -> list[tuple[int, int]]:
        """[start, end) 内扣除 blockers 后**真正可用**的连续空档。"""
        bl = sorted((int(e["start"]), int(e["end"])) for e in self.blockers_in(start, end))
        gaps: list[tuple[int, int]] = []
        cur = start
        for a, b in bl:
            if a > cur:
                gaps.append((cur, a))
            cur = max(cur, b)
        if cur < end:
            gaps.append((cur, end))
        return gaps


def load_spec(path: str) -> ConstructSpec:
    raw = yamlmini.load(path)
    slot_raw = raw.get("cargo_slot") or {}
    arms_raw = raw.get("arms") or {}

    arms: dict[str, Arm] = {}
    for k, v in arms_raw.items():
        arms[k] = Arm(
            core_start=int(v.get("core_start", 0)),
            core_end=int(v.get("core_end", 0)),
            extended_start=int(v.get("extended_start", v.get("core_start", 0))),
            extended_end=int(v.get("extended_end", v.get("core_end", 0))),
            rule=str(v.get("rule", "")),
        )

    sites: list[InsertSite] = []
    for raw_site in slot_raw.get("sites") or []:
        rg = raw_site.get("region") or [0, 0]
        sites.append(InsertSite(
            id=str(raw_site.get("id", "?")),
            label=str(raw_site.get("label", "")),
            region=(int(rg[0]), int(rg[1])),
            cargo_strand=int(raw_site.get("cargo_strand", 1)),
            keep_elements=list(raw_site.get("keep_elements") or []),
            replace_all=bool(raw_site.get("replace_all", False)),
            note=str(raw_site.get("note", "")),
            orientation_rule=str(raw_site.get("orientation_rule", "")),
        ))

    ib = slot_raw.get("insert_between")
    slot = CargoSlot(
        start=int(slot_raw.get("start", 0)),
        end=int(slot_raw.get("end", 0)),
        length=int(slot_raw.get("length", 0)),
        note=str(slot_raw.get("note", "")),
        current_contents=list(slot_raw.get("current_contents") or []),
        assembly_strategy=str(slot_raw.get("assembly_strategy", "UNCONFIRMED")),
        insert_between=[int(ib[0]), int(ib[1])] if ib else None,
        sites=sites,
    )

    return ConstructSpec(
        path=path,
        meta=raw.get("meta") or {},
        elements=list(raw.get("elements") or []),
        cargo_slot=slot,
        arms=arms,
        primers=dict(raw.get("primers") or {}),
        gibson=dict(raw.get("gibson") or {}),
        cargo_checks=dict(raw.get("cargo_checks") or {}),
        open_questions=list(raw.get("open_questions") or []),
        raw=raw,
    )


# ---------------------------------------------------------------------------
# 从载体文件取全部注释 + 按插入平移
# ---------------------------------------------------------------------------

def features_from_snapgene(bb: dict) -> list[dict]:
    """
    从 SnapGene 解析结果取出**全部**特征（除 source）。

    ⚠️ 这是「输出图谱必须保留骨架全部信息」的实现点。
    """
    out: list[dict] = []
    for f in bb.get("features") or []:
        if f.get("type") == "source":
            continue
        label = (f.get("qualifiers", {}).get("label") or [None])[0]
        note = (f.get("qualifiers", {}).get("note") or [None])[0]
        out.append({
            "name": label or f.get("name") or f.get("type") or "feature",
            "start": int(f["start"]),
            "end": int(f["end"]),
            "strand": int(f.get("strand", 1)),
            "kind": f.get("type") or "misc_feature",
            "note": note or "",
            "label": label or "",
        })
    return out


def shift_features(features: list[dict], *, edit_start: int, edit_end: int,
                   delta: int) -> tuple[list[dict], list[str]]:
    """
    把「用 edit_start..edit_end 替换为 delta 净长度变化」反映到特征坐标上。

    - 完全在 edit_start 之前的特征：不变
    - 完全在 edit_end 之后的特征：+delta
    - 与 [edit_start, edit_end) 相交的特征：**被移除**（已被这次编辑取代），
      并记录一条告警 —— 使用者必须知道哪些元件没了。

    ⚠️ v0.9 修正：早期版本对相交元件「保留原坐标」，导致坐标越界 ——
    整段替换 LE–RE（8,590–9,784）时，CmR 仍指向 [8932, 9592)，
    而新序列只剩 9,089 bp，写出图谱时直接失败。
    **相交即被移除**才是正确语义。

    返回 (新特征列表, 告警列表)。
    """
    warnings: list[str] = []
    out: list[dict] = []
    for f in features:
        s, e = f["start"], f["end"]
        if e <= edit_start:
            out.append(dict(f))
        elif s >= edit_end:
            g = dict(f)
            g["start"] = s + delta
            g["end"] = e + delta
            out.append(g)
        else:
            warnings.append(
                f"编辑区间 [{edit_start:,}, {edit_end:,}) 覆盖了元件 "
                f"「{f['name']}」[{s:,}, {e:,}) —— 该元件**已被移除**，"
                "不会出现在输出图谱中")
    return out, warnings


def features_to_construct_map(features: list[dict], n: int):
    """把特征 dict 列表转成 buildout.OutFeature 列表（越界会被 clamp）。"""
    from .buildout import OutFeature

    out = []
    for f in features:
        quals: dict[str, str] = {}
        if f.get("note"):
            quals["note"] = f["note"]
        if f.get("label") and f.get("label") != f.get("name"):
            quals["label"] = f["label"]
        out.append(OutFeature(
            name=f["name"], start=f["start"], end=f["end"],
            strand=f.get("strand", 1), kind=f.get("kind", "misc_feature"),
            qualifiers=quals,
        ))
    return out
