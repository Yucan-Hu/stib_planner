"""
STIB Planner — ShCAST 序列偏好评分。

背景（依据课题组提供的 ShCAST 机制原文）:
  ShCAST 的整合位点并非只由 sgRNA 决定，DNA 序列本身也有偏好：

  1. **TSD（靶位点重复，5 bp）内偏好 GCWGC** —— TnsB 优先在此插入
  2. **TSD 中心 ±5 bp 处偏好 (A/T)** —— 由两个 TnsB 单体的 K290 残基接触
  3. **TSD ±12 bp 处的核苷酸偏好** —— TnsB IIb 结构域（R416/T417/Q425/N428）接触
  4. **整合位点上游 ~25 bp 的 AT 富集区** —— 可容纳两圈 dsDNA 结合的 TnsC 丝状体；
     该区仅在整合位点**单侧**，TnsC 丝结合后定向招募 TnsB，决定 **TnsC-LR** 取向
     （与 RNA 引导位点同向）

═══════════════════════════════════════════════════════════════════════
⚠️ 坐标约定的不确定性（必须显式声明，不得掩盖）
═══════════════════════════════════════════════════════════════════════

我们**知道**（来自 STIB 论文与 SOP）:
    整合位点 = PAM 下游 **56–66 bp**（观察区间；名义偏移取中点）
    插入时在整合位点处产生 **5 bp TSD**

我们**不知道**（课题组提供的机制原文未给出可核实的相对坐标）:
    AT 富集 25 bp 区**相对**整合位点/TSD 的精确偏移

因此本模块的立场是:
  - 「TSD 内容偏好」（规则 1）与「±5 / ±12 侧翼偏好」（规则 2、3）
    都以**整合位点自身**为锚 —— 这三个是本模块给分的依据，坐标定义明确。
  - AT 富集区只报**定性统计量**（不硬编码偏移、不参与硬过滤），
    并在输出中标注「偏移待核实」。

换句话说: **能算的算，不能确定的不装作确定。**
"""

from __future__ import annotations

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# 可调参数（全部集中于此，便于按实测校准）
# ---------------------------------------------------------------------------

#: 整合位点相对 PAM 3' 末端的**名义**偏移。
#: 与可接受区间 INSERT_WINDOW 解耦 —— 范围是范围，代表值是代表值。
#: 权威定义在 sequence.INSERT_OFFSET_FROM_PAM3，此处仅为兼容别名。
INSERT_OFFSET_DEFAULT = 63

#: TSD 长度（论文与 SOP 均为 5 bp）
TSD_LEN = 5

#: TnsB 偏好的 TSD 模体（IUPAC）。出处：机制原文「GCWGC within the TSD」
TSD_PREFERRED_MOTIF = "GCWGC"

#: ±5 bp 处偏好 (A/T)。出处：机制原文「(A/T) at symmetric positions located ±5 bp
#: from the TSD center」
AT_PREFERRED_OFFSET = 5

#: ±12 bp 处的偏好（机制原文给的是接触残基，未给具体偏好碱基）→ 仅统计，不给分
CONTACT_OFFSET_12 = 12

#: AT 富集区长度（机制原文「~25 bp」）。⚠️ 相对偏移未核实，故只做定性统计。
AT_RICH_LEN = 25

_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _matches_iupac(seq: str, pattern: str) -> bool:
    """判断 seq 是否匹配 IUPAC 模体（长度须相同）。"""
    if len(seq) != len(pattern):
        return False
    for s, p in zip(seq.upper(), pattern.upper()):
        allowed = _IUPAC.get(p)
        if allowed is None or s not in allowed:
            return False
    return True


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class SeqFeatures:
    """一个整合位点的 ShCAST 序列特征。"""

    insert_site: int
    tsd: str = ""                    # 5 bp 靶位点重复（= 整合位点处的基因组序列）
    tsd_preferred: bool = False      # 是否含 TnsB 偏好的 GCWGC
    at_at_pm5: tuple[str, str] = ("", "")   # TSD 中心 ±5 bp 的碱基
    at_pm5_ok: bool = False          # 两者是否都是 A/T
    bases_pm12: tuple[str, str] = ("", "")
    at_rich_upstream: dict = field(default_factory=dict)
    score: float = 0.0
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "insert_site": self.insert_site,
            "tsd": self.tsd,
            "tsd_preferred_GCWGC": self.tsd_preferred,
            "bases_at_center_pm5": list(self.at_at_pm5),
            "at_pm5_ok": self.at_pm5_ok,
            "bases_pm12": list(self.bases_pm12),
            "at_rich_upstream": self.at_rich_upstream,
            "seq_feature_score": round(self.score, 2),
            "notes": self.notes,
        }


# ---------------------------------------------------------------------------
# 主计算
# ---------------------------------------------------------------------------

def compute_features(genome: str,
                     insert_site: int,
                     tsd_len: int = TSD_LEN,
                     at_rich_len: int = AT_RICH_LEN) -> SeqFeatures:
    """
    计算一个整合位点的 ShCAST 序列特征。

    坐标: insert_site 为 0-based（整合位点在正链上的位置，即 TSD 起点）。
    """
    n = len(genome)
    f = SeqFeatures(insert_site=insert_site)

    # --- 边界检查 ---
    need_left = max(CONTACT_OFFSET_12 + 2, at_rich_len + tsd_len)
    need_right = at_rich_len + 2
    if insert_site - need_left < 0 or insert_site + need_right > n:
        f.notes.append("整合位点距序列边界过近，部分特征无法计算")
        return f

    # --- 1. TSD 内容偏好（GCWGC）---
    half = tsd_len // 2                      # 5 bp → half = 2
    center = insert_site + half
    # TSD 模体可落在 TSD 窗口内的任何相位：取 5 bp 窗口本身
    f.tsd = genome[insert_site:insert_site + tsd_len].upper()
    if len(f.tsd) == tsd_len:
        f.tsd_preferred = _matches_iupac(f.tsd, TSD_PREFERRED_MOTIF)
        if f.tsd_preferred:
            f.score += 2.0
        else:
            f.notes.append(f"TSD '{f.tsd}' 不含 TnsB 偏好模体 {TSD_PREFERRED_MOTIF}")

    # --- 2. TSD 中心 ±5 bp 处偏好 (A/T) ---
    l5 = genome[center - AT_PREFERRED_OFFSET]
    r5 = genome[center + AT_PREFERRED_OFFSET]
    f.at_at_pm5 = (l5, r5)
    f.at_pm5_ok = (l5 in "AT") and (r5 in "AT")
    if f.at_pm5_ok:
        f.score += 1.5
    else:
        f.notes.append(f"TSD 中心 ±5 bp 为 {l5}/{r5}，非双 (A/T)（偏好项，非硬性）")

    # --- 3. ±12 bp 处（仅统计，不给分：机制原文未给偏好碱基）---
    f.bases_pm12 = (genome[center - CONTACT_OFFSET_12],
                    genome[center + CONTACT_OFFSET_12])
    f.notes.append("±12 bp 由 TnsB IIb 结构域接触；原文未给偏好碱基，故仅统计不计分")

    # --- 4. 上游 AT 富集区（定性统计，⚠️ 相对偏移未核实）---
    # 取整合位点上游 at_rich_len 的窗口做统计，但明确标注偏移未核实
    up = genome[max(0, insert_site - at_rich_len):insert_site].upper()
    if up:
        at = sum(1 for c in up if c in "AT")
        f.at_rich_upstream = {
            "window": f"整合位点上游 {len(up)} bp（⚠️ 相对偏移未从原文核实）",
            "length": len(up),
            "AT_count": at,
            "AT_fraction": round(at / len(up), 3),
            "is_at_rich": (at / len(up)) >= 0.60,
        }

    return f


def score_of(features: SeqFeatures) -> float:
    """返回序列特征分（越大越好）。"""
    return features.score


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------

def selftest() -> list[tuple[str, bool, str]]:
    """合成底物自检：TSD 与 ±5 位置的偏好应可手工核对。"""
    rows: list[tuple[str, bool, str]] = []

    # 构造：让 TSD = GCAGC（匹配 GCWGC，第 3 位 A ∈ {A,T}），且中心 ±5 为 A/T
    # 布局:  [上游 40]  [TSD 起点]
    # TSD = G C A G C  → half = 2, center = TSD 起点 + 2
    # 需要 center-5 与 center+5 都是 A/T
    tsd = "GCAGC"
    left_pad = "T" * 40
    # 让 center 在 TSD 内的位置 +2；center±5 落在 TSD 之外或之内均可
    # 简化：把 TSD 放在左 pad 之后，然后在 center±5 强制放 A/T
    seq = list(left_pad + tsd + "T" * 40)
    insert_site = len(left_pad)
    center = insert_site + 2
    seq[center - 5] = "A"
    seq[center + 5] = "T"
    genome = "".join(seq)

    f = compute_features(genome, insert_site)
    rows.append(("TSD 抽取正确", f.tsd == tsd, f"tsd={f.tsd} 期望 {tsd}"))
    rows.append(("GCWGC 偏好识别（正例）", f.tsd_preferred is True,
                 f"tsd={f.tsd} preferred={f.tsd_preferred}"))
    rows.append(("±5 位置 A/T 识别", f.at_pm5_ok is True,
                 f"±5 = {f.at_at_pm5}"))
    rows.append(("序列特征分为正", f.score >= 2.0, f"score={f.score}（期望 ≥2.0）"))

    # 反面 1：GCGGC —— 第 3 位 G 不满足 W（W = A/T）
    f2 = compute_features("T" * 40 + "GCGGC" + "T" * 40, 40)
    rows.append(("非 GCWGC 不判偏好（相位）", f2.tsd_preferred is False,
                 f"tsd={f2.tsd} preferred={f2.tsd_preferred}"))

    # 反面 2：GCATC —— 第 3 位 A 满足 W，但第 4 位 T 不满足 G
    f3 = compute_features("T" * 40 + "GCATC" + "T" * 40, 40)
    rows.append(("非 GCWGC 不判偏好（碱基）", f3.tsd_preferred is False,
                 f"tsd={f3.tsd} preferred={f3.tsd_preferred}"))

    # 边界：接近序列起点不应崩溃
    try:
        f4 = compute_features("ACGT" * 50, 2)
        rows.append(("近边界不崩溃", True, (f4.notes[0] if f4.notes else "ok")))
    except Exception as exc:  # noqa: BLE001
        rows.append(("近边界不崩溃", False, repr(exc)))

    return rows


if __name__ == "__main__":
    out = selftest()
    w = max(len(r[0]) for r in out)
    ok = 0
    for name, passed, note in out:
        print(f"[{'PASS' if passed else 'FAIL'}] {name:<{w}}  {note}")
        ok += 1 if passed else 0
    print(f"\n{ok}/{len(out)} passed")
    raise SystemExit(0 if ok == len(out) else 1)
