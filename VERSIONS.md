# 版本依赖与兼容性

> 对应插件版本 **v0.0.9**（2026-11-16）
> 这份文档回答三个问题：**它需要什么**、**兼容性是怎么判的**、**升级时该核对什么**。

---

## 1. 实测环境

下表是**真正跑过**的组合，不是"应该能用"：

| 组件 | 版本 | 说明 |
|---|---|---|
| **DeepSeek Harness** | `0.2.2-rc.3` | 插件 `peerDependencies` 的目标版本 |
| **Node.js** | `24.21.0` | DSH 自带（Electron `44.7.0` 的 Node） |
| **Python** | `3.12.14` | DSH 自带运行时；也可用 PATH 上的 3.10+ |
| **操作系统** | Windows 11 | ✅ 全流程实测 |
| macOS / Linux | — | ⚠️ **仅代码审计，未实机验证**（见 ROADMAP P1-1） |

Python 侧**不使用任何第三方包**（纯标准库），所以不存在 `pip install` 这一步。

---

## 2. 插件声明的依赖

`dsh-plugin/package.json`：

| 包 | 声明 | 含义 |
|---|---|---|
| `@deepseek-ai/cordis` | `~4.1.0` | 允许 `4.1.x`；cordis 的次版本升级曾改过 API |
| `@deepseek-ai/dsh-tools` | `^0.2.2-rc.3` | `>=0.2.2-rc.3 <0.3.0` |
| `@deepseek-ai/dsh-subprocess` | `^0.2.2-rc.3` | 同上 |
| `@deepseek-ai/schemastery` | `~3.19.0` | 允许 `3.19.x` |
| `engines.node` | `>=24.9.0` | DSH 从 `0.2.1-rc.2` 起不再支持 Node 22 |

`devDependencies` 里这四个包用**精确版本**（`0.2.2-rc.3` / `4.1.0` / `3.19.0`），
只为编译时类型对齐；运行时以宿主提供的为准。

---

## 3. 兼容性是怎么判定的（关键，容易踩）

DSH 的规则在 `0.2.1-alpha.3` 变过：

```
若包声明了 engines.dsh  →  只用它判定，**忽略** @deepseek-ai/dsh-* 的 peer
否则                    →  回退到 peerDependencies 判定
```

**本插件不声明 `engines.dsh`**，走 peer 规则。原因有两条：

1. **官方 bundle 全部走 peer**（抽查 `dsh-tool-present`、`dsh-tool-bash`、
   `dsh-experimental-badge-skill-bundle` 均是），跟着官方走风险最低；
2. peer 规则已在 `0.2.2-rc.3` 上**实测通过**（`dsh plugin list` 无警告）。

`^0.2.2-rc.3` 的覆盖范围：

| 目标版本 | 是否满足 | 说明 |
|---|---|---|
| `0.2.2-rc.3` / `0.2.2-rc.4` | ✅ | 同一次版本的后续预发布 |
| `0.2.2`（正式版） | ✅ | `0.2.2 > 0.2.2-rc.3` |
| `0.2.3-rc.1` | ✅ | 同一 minor 内 |
| `0.2.1-rc.4` | ❌ | 低于下界 |
| `0.3.0` | ❌ | 跨 minor，**必须重新核对 API** |

**不匹配时的表现**：插件会被**静默拒绝** —— 看上去"装了但什么都没发生"。
`dsh plugin --profile <p> list` 会打印 `is incompatible with dsh ...`，
并列出被拒绝的 `requirements`。

---

## 4. 升级 DSH 时怎么做

不要只改版本号。按这个顺序：

```bash
# 1) 看当前是否被判不兼容
dsh plugin --profile <profile> list

# 2) 若拒绝，先按官方升级指南核对 API（dsh-doc 技能的 references/upgrading.md）
#    本插件实际用到的表面只有 6 个：
#      SubprocessSpawnSpec / SubprocessHandle（spawn 契约：argv/cwd/stdio/graceMs/signal/env）
#      ToolDefinition / defineTool（output.schema + output.render 是必填）
#      ctx.tools.register
#      ctx.skills.register（可选服务，经 ctx.inject(['skills'])）
#      ctx.inject
#      Config（schemastery）

# 3) 改 peerDependencies 后重装
dsh plugin --profile <profile> add /abs/path/to/stib-planner/dsh-plugin

# 4) 跑自检 + 三类任务（见 §6）
```

**临时豁免**（不推荐）：`dsh plugin --profile <p> allow-version`
—— 有崩溃与数据丢失风险，除非你知道改动不影响上述 6 个表面。

---

## 5. Python 解释器：为什么必须探测

插件调用的是 `python -m stib.cli`，需要 **3.10+**。但常见环境里的默认 Python 往往是 3.9：

| 环境 | 默认 `python` / `python3` | 结果 |
|---|---|---|
| Windows（本机实测） | `python` = **3.9.2** | 低于要求 |
| macOS（Xcode 命令行工具） | `python3` = **3.9.6** | 低于要求 |
| DSH 自带运行时 | **3.12.14** | ✅ |

所以插件按顺序探测，取第一个 ≥3.10 的：

```
python3 → python → DSH 自带运行时（DSH_HOME 与 ~/.dsh 均可）
        → macOS Homebrew（/opt/homebrew/bin/python3、/usr/local/bin/python3）
        → python3.13 / 3.12 / 3.11 / 3.10
```

**注意**：插件运行在 DSH Host 进程里时**看不到 `DSH_HOME`**（它只注入给终端/工具会话），
所以 `~/.dsh` 这个与进程环境无关的兜底是必需的（v0.0.3 修）。

要指定解释器：在插件 config 里设 `python`（绝对路径）。

---

## 6. 自检：怎么判断"这个版本在这台机器上是好的"

```bash
# A. 单元自检（不需要任何外部文件）
python -m stib.cli selftest          # 期望 RESULT: ALL PASS（15/15 + 0 problems + 12/12）

# B. 论文锚点（需要参考文献基因组）
python scripts/verify_against_paper.py --genome data/NC_004663.1.fasta   # 期望 ALL PASS

# C. 脱靶引擎交叉验证（与暴力穷举逐条比对）
python scripts/crosscheck_offtarget.py --seed 42 --length 30000 --n 80 --max-mm 4

# D. 引物特异性 / 任意位点
python scripts/verify_primer_specificity.py     # 期望 13/13
python scripts/verify_insertion_anchors.py      # 期望 39/39

# E. 插件契约（在 dsh-plugin/ 下）
node smoke.mjs && node smoke-path.mjs && node smoke-selfcontained.mjs \
  && node smoke-live.mjs && node smoke-cargo-path.mjs
```

三类插入任务（v0.0.8 起纳入验收）：

| 任务 | 命令 | 期望 |
|---|---|---|
| 保留 CmR | `cargo --site after_selection` | 三片段 6 条引物、CmR 保留 |
| 不保留 CmR | `cargo --site after_le` | 两片段 4 条引物、CmR 移除 |
| 中性位点 | `neutral --genome <Gb>` + `plan --start/--end` | 落点标「区间中央」 |

---

## 7. 数据文件的版本约束

| 文件 | 与谁绑定 | 换掉会怎样 |
|---|---|---|
| `stib/data/construct_psh43.yml` | **必须与载体文件配套** | 坐标对不上 → 切错位置或校验失败 |
| `testdata/psh4300C.dna` | 同上（LE/RE/CmR 坐标都来自它） | — |

换载体时改 spec 里的 `cargo_slot` / `arms` 坐标；`--list-sites` 可先核对解析结果。
**不要**把别的载体的坐标直接套到这份 spec 上。

---

## 8. 版本号语义

- 当前处于 **`0.0.x` 预发布**：API 与输出格式仍可能变，**不承诺兼容**。
- 每次发布都在 `ROADMAP.md` §9 记录改了什么、怎么验证。
- 预发布**不进 `latest`**：若将来发 npm，用 `--tag next`。
