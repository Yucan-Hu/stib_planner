# `out/` — 运行产出目录

这个目录存放工具生成的方案、图谱与 JSON。**目录内容是运行结果，已在
`.gitignore` 中排除**（本 `README.md` 是唯一被保留的文件，用于说明用途）。

## 会在这里生成什么

```bash
python -m stib.cli plan --genome testdata/mini.gb --gene thrA \
    --out out/plan.md --html out/plan.html --json out/plan.json
```

| 文件 | 内容 |
|---|---|
| `plan.md` | 方案摘要（或加 `--detailed` 时的**完整实验协议**） |
| `plan.html` | **自包含 HTML 可视化** —— 双击即可打开，内联 SVG，无需联网 |
| `plan.json` | 机器可读结果 |

`cargo` 子命令还会在这里生成 SnapGene 图谱：

```bash
python -m stib.cli cargo --cargo testdata/demo_cargo.fa \
    --backbone testdata/psh4300C.dna --site after_selection \
    --cargo-strand -1 --out out
# → out/<载体名>_<cargo名>.dna / .gb / .fa
```

## 想放别处

`--out` 接受任意路径：

```bash
python -m stib.cli plan ... --out D:/lab/2026-10/plan.md --html D:/lab/2026-10/plan.html
```

**建议按实验日期分目录**，便于与实验记录本对应。
