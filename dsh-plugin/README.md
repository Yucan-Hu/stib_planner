# dsh-stib-planner

DeepSeek Harness 插件：**拟杆菌 STIB 基因编辑方案生成器**。

面向**第一次做拟杆菌编辑的课题组** —— 补 SOP 的入口（位点怎么选）与诊断（出问题怎么判），
而不只是复述流程。

---

## 装什么

本插件是**薄壳**：领域逻辑全部在独立的 Python 包 `stib-planner` 里，插件只做两件事：

1. 把 `stib-planner` CLI 暴露成 9 个模型可调用的确定性工具
2. 随包提供一个 Skill（Markdown），承载流程与判定规则

**依赖**：

| 依赖 | 说明 |
|---|---|
| Node.js `^22.19.0 \|\| >=24.0.0` | DSH 本身的要求 |
| DSH `0.2.0-rc.2` | 插件 peer 锁此版本，见下方「版本兼容」 |
| Python 3.10+ | 运行 `stib-planner`（插件通过 `ctx.subprocess` 调用） |
| `stib-planner` 包 | 需与本插件一并部署，路径由配置项 `plannerDir` 指定 |

`stib-planner` 是**纯标准库**实现，无需 pip 安装任何第三方包。

---

## 提供的 9 个工具

| 工具 | 用途 |
|---|---|
| `stib_plan` | 单基因完整方案：候选排序 + 脱靶评估 + 寡核苷酸 + 验证引物 + 流程 + 判定点 + 风险清单 |
| `stib_find_sgrnas` | 只列某区域的候选 sgRNA |
| `stib_offtarget` | 单独的全基因组脱靶评估（支持内置已验证面板） |
| `stib_batch` | 批量设计：宽区间内一次给多条候选 |
| `stib_open_questions` | 列出已核实的材料矛盾 |
| `stib_selftest` | 自检（论文序列 + 知识库 + 脱靶引擎） |
| `stib_find_gene` | 列出/搜索基因（**基因名不确定时先查这个**；需 GenBank 注释） |
| `stib_neutral_sites` | 中性位点查找 + 位点内 sgRNA 设计（插入大片段 cargo 用） |
| `stib_cargo_assembly` | cargo 的 Gibson 组装引物 + 输出 SnapGene `.dna` 图谱 |

---

## 配置

在 `cordis.patch.yml` 里改这三项（**必须是绝对路径**）：

```yaml
- insert:
    - id: stib-planner
      name: dsh-stib-planner
      config:
        python: /absolute/path/to/python        # 或 PATH 上的 python
        plannerDir: /absolute/path/to/stib_planner
        dataDir: /absolute/path/to/stib_planner/data
        timeoutMs: 180000
        maxOutputBytes: 1048576
```

**工具里的 `accession` 参数支持三种写法**：

| 写法 | 例子 | 解析方式 |
|---|---|---|
| 完整路径 | `E:\lab\genomes\BT.gb` | 原样使用 |
| 相对路径 | `testdata/mini.gb` | 相对 `plannerDir` |
| 裸 accession | `NC_004663.1` | 在 `dataDir` 下按 `.gb/.gbk/.gbff/.fasta/.fa/.fna/.fas` 依次探测 |

⚠️ 早期版本对任何输入都硬拼 `.fasta`，导致 `.gb` 路径完全不可达
（`stib_find_gene` 直接 `FileNotFoundError`）。现已修正，并在文件扩展名
与实际内容不符时按内容自动纠正。

---

## ⚠️ 版本兼容：这一条最容易踩坑

DSH 的兼容性检查**只看 `peerDependencies`**（`engines.dsh` 不生效）。
运行时是**预发布版**，普通 caret 范围**匹配不上预发布版** ——
不匹配时插件行会被**静默拒绝**，表现为「装了但什么都没发生」。

因此本包把 peer **精确锁死**：

```json
"peerDependencies": {
  "@deepseek-ai/cordis": "4.0.4",
  "@deepseek-ai/dsh-tools": "0.2.0-rc.2",
  "@deepseek-ai/dsh-subprocess": "0.2.0-rc.2",
  "@deepseek-ai/schemastery": "3.18.4"
}
```

**升级 DSH 后必须重新核对这几个版本号。**
若被拒绝，用 `dsh plugin --profile <p> version-exemptions` 查看，
再决定是否豁免（豁免有风险，需谨慎）。

---

## 开发

```bash
npm install
npm run build      # tsc -> lib/
```

编译产物：`lib/index.js` 与 `lib/types/index.d.ts`。

**冒烟测试**（确认导出符合 DSH 契约）：

```bash
node --input-type=module -e "
const m = await import('./lib/index.js');
console.log(Object.keys(m).sort().join(', '));
"
# 期望：Config, apply, inject, name
```

⚠️ 本包**不得**有 `default` 导出 —— 官方 AGENTS.md 明确要求
对无 `inject` 的 bundle 加显式守卫 `expect('default' in mod).toBe(false)`。

---

## 安装到桌面版

```bash
# 1) 完全退出 DSH 桌面版（profile 被 Desktop 独占）
# 2) 用官方 CLI 安装
dsh plugin --profile desktop add /absolute/path/to/dsh-plugin
# 3) 重开桌面版
```

装好后应出现上述 9 个工具。

---

## 能力边界（不要对外宣称超出这些）

- **不预测编辑效率** —— 已发表数据点太少（约 11 条带效率的位点），
  抽不出具统计意义的模型。排序只反映「是否可用」与可解释偏好。
- **脱靶引擎为 6-mer 块 ×3 + 每块 1 处替换枚举**：≤4 错配保证不漏检
  （鸽巢上限 5；已用暴力穷举交叉验证）。**>4 时保证失效**且会显式警示。
  **不等价于 Cas-OFFinder 穷举。**
- **任何错配都只增风险，不存在「减害」。**
- **不设计 cargo 内部序列与骨架引物**（需用户提供序列）。
- **未做湿实验验证。** 正确性上限是「忠实于已发表方法与实验记录」。

---

## 许可与安全

MIT，见 [LICENSE](LICENSE)。

**本项目不提供受管制生物、致病性增强、抗性播散类编辑目标的方案。**
涉及此类目标时，插件与底层引擎均应在知识库入口拒绝。
