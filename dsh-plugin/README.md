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
| Node.js `>=24.9.0` | DSH 0.2.1-rc.2 起不再支持 Node 22 |
| DSH `^0.2.2-rc.3` | 插件 peer 声明；见下方「版本兼容」 |
| Python 3.10+ | 运行 stib-planner。**无需配置** —— 插件自动探测 |
| `stib/` Python 包 | **随插件发布**，无需另行部署 |

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

**默认零配置。** `cordis.patch.yml` 里只需保留超时上限：

```yaml
- insert:
    - id: stib-planner
      name: dsh-stib-planner
      config:
        timeoutMs: 300000
        maxOutputBytes: 2097152
```

以下三项**都是可选的覆盖项**，只有特殊情形才需要填：

| 配置项 | 默认行为 | 什么时候要填 |
|---|---|---|
| `python` | 自动探测（见下），取第一个 **3.10+** 的 | 机器上有多个 Python，想指定其中一个 |
| `plannerDir` | **插件包根目录**（其下就是随包发布的 `stib/`） | 想改用仓库里另一份 Python 代码 |
| `dataDir` | `<插件目录>/data` | 想把基因组缓存放到别处 |

### Python 解释器是怎么找到的

按顺序探测，第一个满足 `>= 3.10` 的胜出（结果缓存，只探一次）：

1. `python3`、`python`（PATH 上）
2. **DSH 自带运行时**：`<home>/dsh-runtimes/*/dependencies/python/{python.exe, bin/python3}`
   - `<home>` 依次取 `DSH_HOME`、`~/.dsh`、`DSH_PROFILE_DIR` 的上两级
   - ⚠️ **只靠 `DSH_HOME` 不够**：实测插件运行在 Host 进程里时看不到该变量
     （它只注入给终端/工具会话），所以必须有 `~/.dsh` 这个与进程环境无关的兜底
3. macOS：`/opt/homebrew/bin/python3`（Apple Silicon）、`/usr/local/bin/python3`
4. 带小版本的名字：`python3.13` → `python3.10`

> **为什么必须探测而不是直接用 `python`**：
> Windows 上 PATH 里的 `python` 常见是 **3.9**（本机实测 3.9.2）；
> macOS 上 Xcode 命令行工具自带的 `python3` 是 **3.9.6**。
> 两者都低于要求的 3.10，直接用会以看不懂的方式失败。
> DSH 自带运行时是 3.12，所以把它列为主要候选。

### 平台支持

| 平台 | 状态 | 说明 |
|---|---|---|
| Windows | ✅ 已实测 | Python 3.12（DSH 自带） |
| macOS | ⚠️ 代码已按 POSIX 路径写，**未实机验证** | 探测候选含 Homebrew 与 `bin/python3` 布局 |
| Linux | ⚠️ 同上 | 依赖 PATH 上的 python3 或 DSH 自带运行时 |

跨平台相关的实现要点（已逐项核对）：

- 路径全部走 `path.join` / `os.path.join`，**没有**字符串拼 `/`
- Python 侧所有文本写入都显式 `encoding="utf-8", newline="\n"`（两平台都出 LF）
- `.dna` 用**二进制**写入，不受换行转换影响
- 读文件用 `utf-8-sig`，**带 BOM 也不会污染首条记录**
- 无 `sys.platform` 分支、无 `winreg`/`msvcrt`、无外部命令调用

**工具里的 `accession` 参数支持三种写法**：

| 写法 | 例子 | 解析方式 |
|---|---|---|
| 完整路径 | `E:\lab\genomes\BT.gb` | 原样使用 |
| 相对路径 | `testdata/mini.gb` | 相对插件包目录 |
| 裸 accession | `NC_004663.1` | 在 `dataDir` 下按 `.gb/.gbk/.gbff/.fasta/.fa/.fna/.fas` 依次探测 |

⚠️ 早期版本对任何输入都硬拼 `.fasta`，导致 `.gb` 路径完全不可达
（`stib_find_gene` 直接 `FileNotFoundError`）。现已修正，并在文件扩展名
与实际内容不符时按内容自动纠正。

---

## ⚠️ 版本兼容：这一条最容易踩坑

从 DSH `0.2.1-alpha.3` 起，兼容性判定规则是：
**先读包自身的 `engines.dsh`，没有才回退到 `@deepseek-ai/dsh-*` 的
`peerDependencies`**。本插件不声明 `engines.dsh`，走 peer 规则。

当前声明（对应 DSH `0.2.2-rc.3`）：

```json
"peerDependencies": {
  "@deepseek-ai/cordis": "~4.1.0",
  "@deepseek-ai/dsh-tools": "^0.2.2-rc.3",
  "@deepseek-ai/dsh-subprocess": "^0.2.2-rc.3",
  "@deepseek-ai/schemastery": "~3.19.0"
}
```

`^0.2.2-rc.3` 覆盖同一 minor（`0.2.x`）内的后续预发布与正式版，
所以 `0.2.2` 正式版、`0.2.3-rc.1` 都仍可用。

**若被拒绝**（表现为 `plugin list` 出现 `is incompatible with dsh ...`
且插件不工作）：

```bash
dsh plugin --profile <p> list          # 看拒绝信息里的 requirements
```

然后按需二选一：
- **推荐**：改 `peerDependencies` 到你的 DSH 版本后重新安装；
- 临时：`dsh plugin --profile <p> allow-version` 豁免（**有崩溃与数据损失风险**）。

升级 DSH 跨 minor（如 `0.2.x` → `0.3.x`）时，**必须重新核对 API**，
不能只改版本号：`ToolDefinition` 与 `SubprocessSpawnSpec` 都是会变的契约。

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
