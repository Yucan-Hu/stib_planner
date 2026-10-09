"""
STIB Planner — 最小 YAML 读取器。

═══ 为什么不用 PyYAML ═══

本项目**刻意零第三方依赖**（用户要在离线实验室机器上直接跑）。
这里实现一个只覆盖规范文件实际用法的 YAML 子集。

支持：
  - 嵌套映射（按缩进）
  - 列表（`- item`、`- key: value`、`- {inline: map}`）
  - 块标量  `key: >` / `key: |`（后续**更深缩进**的行拼成一段文本）
  - 行内数组 `[a, b]`、行内映射 `{k: v}`
  - 整行注释

**不支持**（遇到抛 YamlError，而不是静默给错值）：
  锚点/别名、多文档、制表符缩进、引号内转义。

设计取向：宁可抛错并说清不支持什么，也不静默解析出错误坐标 ——
这份文件里的坐标决定 cargo 放在哪里。
"""

from __future__ import annotations

from typing import Any


class YamlError(ValueError):
    """YAML 子集不支持的语法。"""


def _strip_comment(line: str) -> str:
    """去掉行尾注释（引号内的 # 不算）。"""
    out: list[str] = []
    quote: str | None = None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
            out.append(ch)
            continue
        if ch == "#":
            break
        out.append(ch)
    return "".join(out).rstrip()


def _split_flow(s: str) -> list[str]:
    """按顶层逗号切分（忽略引号与嵌套括号内的逗号）。"""
    parts: list[str] = []
    depth = 0
    quote: str | None = None
    cur: list[str] = []
    for ch in s:
        if quote:
            cur.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
            cur.append(ch)
            continue
        if ch in "[{":
            depth += 1
        elif ch in "]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
            continue
        cur.append(ch)
    if cur:
        parts.append("".join(cur).strip())
    return [p for p in parts if p]


def _scalar(tok: str) -> Any:
    t = tok.strip()
    if t == "" or t == "~" or t.lower() == "null":
        return None
    if t.lower() == "true":
        return True
    if t.lower() == "false":
        return False
    if len(t) >= 2 and t[0] == t[-1] and t[0] in "\"'":
        return t[1:-1]
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t)
    except ValueError:
        pass
    return t


def _flow(tok: str) -> Any:
    """解析行内 [ ... ] / { ... }（一层，元素为标量或行内映射）。"""
    t = tok.strip()
    if t.startswith("[") and t.endswith("]"):
        inner = t[1:-1].strip()
        return [] if not inner else [_flow(x) for x in _split_flow(inner)]
    if t.startswith("{") and t.endswith("}"):
        inner = t[1:-1].strip()
        if not inner:
            return {}
        out: dict[str, Any] = {}
        for part in _split_flow(inner):
            if ":" not in part:
                raise YamlError(f"行内映射元素缺少 ':': {part!r}")
            k, v = part.split(":", 1)
            out[k.strip()] = _flow(v)
        return out
    return _scalar(t)


class _Parser:
    def __init__(self, text: str) -> None:
        raw = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
        self.lines: list[tuple[int, str]] = []
        for ln in raw:
            indent_part = ln[: len(ln) - len(ln.lstrip(" \t"))]
            if "\t" in indent_part:
                raise YamlError("不支持制表符缩进；请用空格")
            body = _strip_comment(ln)
            if not body.strip():
                continue
            self.lines.append((len(body) - len(body.lstrip(" ")), body.strip()))
        self.pos = 0

    def peek(self) -> tuple[int, str] | None:
        return self.lines[self.pos] if self.pos < len(self.lines) else None

    def _block_scalar(self, parent_indent: int) -> str:
        """
        收集块标量内容：必须**比父键更深缩进**。
        遇到同级或更浅的键即停止 —— 上一版正是这里没停，把后续内容全吃掉了。
        """
        chunks: list[str] = []
        while self.pos < len(self.lines):
            ind, body = self.lines[self.pos]
            if ind <= parent_indent:
                break
            chunks.append(body)
            self.pos += 1
        return " ".join(chunks)

    def parse_block(self) -> Any:
        head = self.peek()
        if head is None:
            return None
        indent0, body0 = head
        if body0.startswith("- ") or body0 == "-":
            return self._parse_list(indent0)
        return self._parse_map(indent0)

    def _parse_list(self, indent0: int) -> list[Any]:
        items: list[Any] = []
        while self.pos < len(self.lines):
            ind, body = self.lines[self.pos]
            if ind != indent0 or not (body.startswith("- ") or body == "-"):
                break
            rest = body[1:].strip()
            self.pos += 1
            if rest == "":
                items.append(self.parse_nested(indent0))
                continue
            if ":" in rest and not rest.startswith(("{", "[", "\"", "'")):
                items.append(self._parse_item_map(rest, indent0))
                continue
            items.append(_flow(rest))
        return items

    def _parse_item_map(self, first: str, item_indent: int) -> dict[str, Any]:
        """解析 `- k: v` 及后续同项键（缩进 > item_indent）。"""
        m: dict[str, Any] = {}
        k, v = first.split(":", 1)
        key, val = k.strip(), v.strip()
        if val.split()[0] in (">", "|", ">-", "|-", ">+", "|+") if val else False:
            m[key] = self._block_scalar(item_indent)
        elif val == "":
            m[key] = self.parse_nested(item_indent)
        else:
            m[key] = _flow(val)
        while self.pos < len(self.lines):
            ind, body = self.lines[self.pos]
            if ind <= item_indent or body.startswith("- "):
                break
            if ":" not in body:
                raise YamlError(f"无法解析的行: {body!r}")
            kk, vv = body.split(":", 1)
            self.pos += 1
            key2, val2 = kk.strip(), vv.strip()
            if val2.split()[0] in (">", "|", ">-", "|-", ">+", "|+") if val2 else False:
                m[key2] = self._block_scalar(ind)
            elif val2 == "":
                m[key2] = self.parse_nested(ind)
            else:
                m[key2] = _flow(val2)
        return m

    def _parse_map(self, indent0: int) -> dict[str, Any]:
        out: dict[str, Any] = {}
        while self.pos < len(self.lines):
            ind, body = self.lines[self.pos]
            if ind != indent0 or body.startswith("- "):
                break
            if ":" not in body:
                raise YamlError(f"无法解析的行: {body!r}")
            k, v = body.split(":", 1)
            key, val = k.strip(), v.strip()
            self.pos += 1
            if val.split()[0] in (">", "|", ">-", "|-", ">+", "|+") if val else False:
                out[key] = self._block_scalar(ind)
            elif val == "":
                out[key] = self.parse_nested(ind)
            else:
                out[key] = _flow(val)
        return out

    def parse_nested(self, parent_indent: int) -> Any:
        """解析作为某键值的嵌套块（缩进必须 > parent_indent）。"""
        head = self.peek()
        if head is None or head[0] <= parent_indent:
            return None
        return self._parse_list(head[0]) if (head[1].startswith("- ") or head[1] == "-") \
            else self._parse_map(head[0])


def loads(text: str) -> Any:
    return _Parser(text).parse_block()


def load(path: str) -> Any:
    with open(path, "r", encoding="utf-8-sig") as fh:
        return loads(fh.read())
