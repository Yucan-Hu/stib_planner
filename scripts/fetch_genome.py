#!/usr/bin/env python3
"""
取回并拼接完整基因组（分块请求，避免单次超大响应）。

用法:
  python scripts/fetch_genome.py NC_004663.1 [--out data] [--chunk 500000]

为什么分块: NCBI efetch 对超长序列单次返回可能截断或被网关限制。
分块 + 拼接 + 长度校验是更稳的做法。
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def _get(url: str, retries: int = 3) -> str:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "stib-planner/0.1 (research; contact: lab)"}
            )
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise SystemExit(f"ERROR: 请求失败 {url}\n  {last}")


def get_length(accession: str) -> int:
    """从 esummary 取序列长度。"""
    q = urllib.parse.urlencode({"db": "nuccore", "id": accession, "retmode": "json"})
    text = _get(f"{ESUMMARY}?{q}")
    import json

    data = json.loads(text)
    result = data.get("result", {})
    uid = result.get("uids", [None])[0]
    if uid is None:
        raise SystemExit(f"ERROR: 无法解析 esummary: {text[:300]}")
    slen = result[uid].get("slen")
    if not slen:
        raise SystemExit(f"ERROR: esummary 中无 slen 字段")
    return int(slen)


def fetch_chunk(accession: str, start: int, stop: int) -> str:
    """取回 [start, stop]（1-based，闭区间）的 FASTA 序列部分（不含 header）。"""
    q = urllib.parse.urlencode({
        "db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text",
        "seq_start": start, "seq_stop": stop,
    })
    text = _get(f"{EFETCH}?{q}")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith(">"):
        raise SystemExit(f"ERROR: 非 FASTA 响应（{start}-{stop}）:\n{text[:300]}")
    return "".join(lines[1:])


def main() -> int:
    ap = argparse.ArgumentParser(description="分块取回完整基因组")
    ap.add_argument("accession")
    ap.add_argument("--out", default="data")
    ap.add_argument("--chunk", type=int, default=500_000)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    dest = os.path.join(args.out, f"{args.accession}.fasta")

    print(f"查询长度: {args.accession}")
    total = get_length(args.accession)
    print(f"  长度 = {total:,} bp")

    parts: list[str] = []
    pos = 1
    while pos <= total:
        stop = min(pos + args.chunk - 1, total)
        print(f"  取回 {pos:>10,} - {stop:>10,} ...", end="", flush=True)
        seq = fetch_chunk(args.accession, pos, stop)
        expect = stop - pos + 1
        status = "OK" if len(seq) == expect else f"长度不符 (得 {len(seq)}, 期望 {expect})"
        print(f" {status}")
        parts.append(seq)
        pos = stop + 1
        time.sleep(0.4)  # 遵守 NCBI 建议的请求间隔

    genome = "".join(parts)
    if len(genome) != total:
        print(f"WARNING: 拼接长度 {len(genome):,} != 期望 {total:,}", file=sys.stderr)

    with open(dest, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f">{args.accession} {os.path.basename(dest)}\n")
        for i in range(0, len(genome), 70):
            fh.write(genome[i:i + 70] + "\n")

    print(f"\n完成: {dest}")
    print(f"  序列长度 = {len(genome):,} bp")
    return 0


if __name__ == "__main__":
    sys.exit(main())
