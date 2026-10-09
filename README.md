# STIB Planner — 拟杆菌 STIB 基因编辑方案生成器

**v0.0.1（预发布）** ｜ 纯 Python 标准库，**零第三方依赖** ｜ 可离线运行

> ⚠️ **未做湿实验验证。** 本工具的全部领域知识、条件值与判定规则均来自下方论文
> 及其配套 SOP，但它**不替代**原始来源，也不保证在你的菌株与实验条件下成立。
> 发表级结论前请人工复核。

---

## 引用

使用本工具产出的方案时，**请引用它所依据的论文** —— 工具里的每一个条件值、
判定规则与风险提示都来自这篇工作：

> Hu Y, Li Q, Li Y, Zeng Y, Zheng L, Shen J, Gao X, Zhao G-P, Zhao W, Dai L.
> **Targeted genomic editing of human gut *Bacteroides* species based on
> CRISPR-associated transposases.** *Cell Systems*, 2026, 101650.
> doi:[10.1016/j.cels.2026.101650](https://doi.org/10.1016/j.cels.2026.101650)

机器可读的引用信息见 [`CITATION.cff`](CITATION.cff)。

---

## 解决什么问题

课题组已有 STIB 操作 SOP，覆盖"怎么做"。本工具补 SOP 的两端：

| | SOP 已有 | 本工具补的 |
|---|---|---|
| 步骤与参数 | ✅ | — |
| **位点怎么选** | ❌ | 给定基因 → 候选 sgRNA + 全基因组脱靶评估 |
| **出问题怎么查** | 部分 | 症状 → 假设 → **最便宜的判别实验** |

具体产出：

- **sgRNA 与引物**：可直接下单的寡核苷酸、跨 junction 的验证引物
- **全基因组脱靶评估**：正负链，≤4 错配枚举，按错配位置分级
- **定制实验协议**：Golden Gate 与 Gibson 的精确体系体积与程序、接合分步、时间轴
- **中性位点查找**：插入大片段 cargo 时找不破坏宿主功能的落点
- **cargo 组装**：Gibson 引物 + 可打开的 SnapGene 图谱
- **HTML 可视化报告**：候选位置图、评分条形图、脱靶明细

---

## 安装

### 推荐：装 DeepSeek Harness 桌面版，然后在对话框里让它自己装

**不需要命令行、不需要 git、不需要懂 Python。** 全程在对话框里完成。

#### 第 1 步：装 DeepSeek Harness 桌面版

到官方仓库下载并安装（Windows / macOS / Linux 都有）：

> **<https://github.com/deepseek-ai/deepseek-harness>**

装好后打开，用你的方式登录。**先随便新建一个会话确认能用**。

#### 第 2 步：在对话框里说这句话

> **"帮我安装 dsh-stib-planner 插件，这是一个拟杆菌 STIB 基因编辑方案生成器，
> 仓库在 https://github.com/Yucan-Hu/stib_planner"**

agent 会：

1. 把仓库克隆到本地
2. 调用 `dsh plugin add` 装进你的 profile
3. **提示你完全退出并重开 DSH**（这一步不能省，插件在启动时加载）

#### 第 3 步：重开后验证

新开一个会话，说：

> **"跑一下 stib_selftest"**

看到 `RESULT: ALL PASS` 就说明装好了。

#### 第 4 步：开始用

> **"用 E:\path\to\BT.gb，给 tdk 基因设计插入失活的方案。"**

只需要说清**基因组文件 + 基因名**。想更细就说：

> **"给我一份详细到能照着做实验的 protocol，输出 HTML 和 md。"**

> ⚠️ **如果 agent 说装不上**：最常见原因是 DSH 版本与插件的
> `peerDependencies` 不匹配（本插件锁定了开发时的版本）。
> 把 agent 报的原文贴给我，可以快速定位。

---

### 方式二：直接用 CLI（不需要 DSH）

适合想先看看产出、或要集成进自己流程的人。

```bash
git clone https://github.com/Yucan-Hu/stib_planner.git
cd stib_planner
python -m stib.cli --help
```

**只需要 Python 3.10+。** 工具是纯标准库实现，**无需 pip 安装任何东西**。

---

## 快速使用

### 先验证环境没问题

```bash
python -m stib.cli selftest
```

期望输出 `RESULT: ALL PASS`（3 层共 27 项：序列层 / 知识库层 / 脱靶引擎）。
**这一步不需要任何外部文件。**

### 拿一个真实基因组试

```bash
# 取回参考基因组（6.1 MB，只需一次）
python -m stib.cli fetch NC_004663.1

# 按基因名生成方案（需要 GenBank 格式才能按基因名定位）
python -m stib.cli plan --genome testdata/mini.gb --gene thrA \
    --out plan.md --html plan.html

# 要"照着做"的完整实验协议，加 --detailed
python -m stib.cli plan --genome testdata/mini.gb --gene thrA --detailed \
    --out protocol.md --html protocol.html
```

产出 `plan.md`（方案）+ `plan.html`（可视化，**双击即可打开，无需联网**）。

### 全部命令

```
selftest  check  facts  conflicts  fetch  scan  plan  batch
genes  neutral  cargo  offtarget
```

---

## 三条不可协商的设计原则

### 1. 不做效率预测

已发表数据点太少（约 11 条带效率的位点），抽不出具统计意义的评估标准。
硬做只会产出**假精度** —— 看起来科学，实际无法校准。本工具只做可解释的
排序（PAM 优选、落点位置、序列偏好、脱靶风险），并明确说明排序不代表效率。

### 2. 脱靶评估是硬判据，不是打分项

高风险的脱靶位点会让候选**直接被排除**，而不是扣几分了事。

- **引擎**：6-mer 块 ×3 + 每块允许 1 处替换，全基因组正负链枚举
- **保证**：**≤4 错配不漏检**（鸽巢：需 ≥6 个错配才能弄脏全部 3 块）。
  **已用暴力穷举交叉验证**（4 组随机种子 × 80 条 spacer，完全一致）
- **边界**：>4 错配时保证失效、**可能少报** → 输出中显式警示
- **不等价声明**：与 Cas-OFFinder 的穷举**不等价**。发表级结果请用它复核

> **分区与 Cas9 相反**：ShCAST 的 PAM 位于 spacer **上游**，因此关键区是
> **PAM 近端前 17 nt（helix）**，3' 端 6 nt（tail）为延伸尾部。
> 错配落在前 17 nt 的位点**更危险**。

### 3. 未收敛变量必须显式呈现

来源之间**有冲突**的数值（如插入窗口在不同章节写 60–66 与 56–63 bp），
**并列列出并标注**，绝不悄悄取一个"看起来合理"的值。

每个条件值都带出处（论文 / 补充表 / SOP 章节 / 实验记录日期），
使用者能核对回原始来源。

---

## 技术上不做的事

- ❌ **不预测编辑效率** —— 数据点太少，无可用校准模型
- ❌ **不做 Cas-OFFinder 等价的穷举** —— 本引擎是 6-mer 种子 + 替换枚举；
  ≤4 不漏检，>4 保证失效（会警示）
- ❌ **不设计 cargo 内部序列** —— 需用户提供目的片段
- ❌ **不设计骨架 / LE / RE 引物** —— 需用户提供骨架与 cargo 序列
- ❌ **不验证引物二聚体 / 发夹结构** —— Tm 为估算值，须梯度 PCR 确认
- ❌ **不替代课题组 SOP** —— 流程部分引用既有 SOP，不重写
- ❌ **不处理区域抽取的 GenBank 坐标偏移** —— 特征坐标是绝对的、
  序列是相对的，目前未处理这种情形
- ❌ **不做湿实验验证** —— 正确性上限是"忠实于已发表方法与实验记录"

---

## 数据说明

仓库**不含基因组**（体积大且可随时取回），需要自己 `fetch`。
但包含以下用于自检与演示的小文件：

| 路径 | 内容 | 说明 |
|---|---|---|
| `stib/data/construct_psh43.yml` | psh43 构建规范 | **只有坐标与元件名，无序列** |
| `testdata/psh4300C.dna` | psh43 载体 | 10,040 bp 完整序列 + 29 个注释 |
| `testdata/mini.gb` | 小型 GenBank 示例 | 演示按基因名定位 |
| `testdata/demo_cargo.fa` | 示例 cargo | 演示组装流程 |

---

## 更多文档

| 文档 | 内容 |
|---|---|
| [USER_GUIDE.md](USER_GUIDE.md) | 面向第一次使用者的对话指导、可直接抄的问法、常见踩坑 |
| [LICENSE](LICENSE) | MIT |
