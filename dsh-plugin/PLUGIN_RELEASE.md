# 插件构建与发布要求

> 状态：**源码已通过类型检查与模块冒烟测试**；尚未在真实 DSH 会话中激活。
> 本文列出发布前必须补齐的项与验证方式。

---

## 1. 好消息：插件层第一次通过了真实编译

v0.1 的插件代码是**纯猜测**的。v0.2 依据**官方真实契约**重写后，用 npm 上真实存在的
`0.2.2-rc.3` 类型定义编译，结果：

```
npm install               →  26 packages
npx tsc -p tsconfig.json  →  exit=0（无错误）
```

产物与导出：

```
lib/index.js          16,276 B
lib/types/index.d.ts   3,065 B

node --input-type=module -e "import('./lib/index.js')"
导出: Config, apply, inject, name        ← 符合 DSH 契约
name = "stib-planner"
inject = ["tools","subprocess"]
apply 是函数: true
default 导出存在: false                  ← 符合官方 AGENTS.md 的 bundle 守卫要求
```

### 1.1 v0.1 到底错在哪（契约对照）

| 项 | v0.1（猜的，错） | 真实契约 |
|---|---|---|
| 调用形态 | `spawn({ command, args, cwd, stdio: 'pipe' })` | **`spawn(spec)`**，只有 `spec.argv: string[]`；无 `command`/`args` |
| 必需字段 | 缺 `graceMs` | **`graceMs: number` 必填**（正有限毫秒） |
| `stdio` | 传字符串 `'pipe'` | 每流显式：`{ stdin: 'ignore', stdout: { maxBytes }, stderr: { maxBytes } }` |
| 取输出 | `handle.stdout` 当字节流 `readAll` | collect 模式下**没有** `.stdout`，要用 **`handle.collected.stdout.readFrom(0).text`** |
| 等退出 | `handle.exited` | **`handle.done: Promise<SubprocessOutcome>`** → `.exitCode` |
| 超时 | 未处理 | 契约的 `graceMs` **不是 deadline**，只作用于 terminate 流程 → 须自己起定时器调 `terminate()` |

**依据**：`@deepseek-ai/dsh-subprocess/src/types.ts`（`0.2.2-rc.3`，已复核契约未变）。
`src/index.ts` 顶部完整抄录了该契约，便于日后核对。

---

## 2. 已就绪的部分

| 项 | 状态 |
|---|---|
| `src/index.ts` 按真实契约重写 | ✅ 编译通过 |
| 9 个工具（plan / find_sgrnas / offtarget / batch / open_questions / selftest / find_gene / neutral_sites / cargo_assembly） | ✅ |
| `Config`（schemastery：`python?` / `plannerDir?` / `dataDir?` / `timeoutMs` / `maxOutputBytes`，前三项留空即自动） | ✅ |
| `tsconfig.json` + `npm run build`（tsc） | ✅ |
| `cordis.patch.yml`（bundle 层，含显式 `id`） | ✅ 需填绝对路径 |
| `SKILL.md`（面向无经验课题组） | ✅ |
| `LICENSE`（MIT + 不替代实验判断的声明） | ✅ |
| `README.md`（含版本兼容告警） | ✅ |
| 导出契约冒烟测试 | ✅ `Config, apply, inject, name` |
| peer 版本锁死为**真实存在**的版本 | ✅ 已在 npm 上验证存在 |

---

## 3. 发布前**必须**补齐的项

### 3.1 🔴 生物安全拒绝面（合规底线，当前缺失）

`PLAN.md` §11 要求「在知识库入口即拦截」，**目前完全未实现**。
这是对外发布的**合规必要条件**，不是可选项。

| 需拦截的类别 | 处理方式 |
|---|---|
| 高致病性 / 受管制生物（按国家名录） | 拒绝生成方案，不进入设计流程 |
| 增强致病性、逃逸检测、抗性播散 | 拒绝 |
| 超出机构生物安全等级的操作 | 拒绝并提示所需等级 |

**实现位置**：应在 `stib/knowledge.py` 的入口层做**物种与目标白/黑名单**校验，
而不是写在提示词里（提示词可被绕过）。插件层再加一道 `tools/pre-execute` 守卫。

### 3.2 🟡 在真实 DSH 会话中激活验证

编译通过 ≠ 能装能跑。仍需：

| 步骤 | 验证什么 |
|---|---|
| 1. 填 `cordis.patch.yml` 的绝对路径 | 配置能被解析 |
| 2. 完全退出桌面版 → `dsh plugin --profile desktop add <绝对路径>` | bundle 层能加载；**peer 版本是否被接受** |
| 3. 重开桌面版 | 插件是否出现在列表（若没出现 → 大概率是 `incompatible-version` 静默拒绝） |
| 4. 调 `stib_selftest` | 工具能真正 spawn Python 并拿到输出 |
| 5. 调 `stib_plan` | 端到端：参数传递、collect 模式读输出、中文编码 |

**第 4 步是插件层唯一的真实验证**，也是当前「未验证」的具体所指。

### 3.3 🟡 大输出的 spill 处理

`maxOutputBytes` 默认 1 MiB。`plan` 的 Markdown 输出目前远小于此，
但 `stib_offtarget --panel` 在复杂基因组上可能接近上限。
当前实现会在 `lossy` 时附上 `spillPath`，但**未验证** spill 在桌面版中是否可读。

### 3.4 🟢 发布渠道（你已指示：开源）

| 渠道 | 动作 |
|---|---|
| GitHub | 打 `dsh-plugin` topic（官方 README 指定的发现方式） |
| npm | `npm publish`（`prepublishOnly` 已配置 clean+build） |
| 可选 | 提 PR 到 DSH 仓库的插件列表（若有） |

**开源前请确认**：论文/记录库相关内容已获课题组授权公开；
`knowledge.py` 中引用的是**已发表论文**与**课题组 SOP 的公开部分**，
未包含未发表序列。若含未发表内容，需先做脱敏。

---

## 4. 版本兼容：最容易踩的坑

**规则（自 DSH `0.2.1-alpha.3` 起）**：先读包自身的 `engines.dsh`；
没有才回退到 `@deepseek-ai/dsh-*` 的 `peerDependencies`。本包不声明
`engines.dsh`，走 peer 规则。不匹配时插件行被拒绝，表现为「装了但没反应」。

当前声明（对应 DSH `0.2.2-rc.3`）：

```json
"@deepseek-ai/cordis": "~4.1.0",
"@deepseek-ai/dsh-tools": "^0.2.2-rc.3",
"@deepseek-ai/dsh-subprocess": "^0.2.2-rc.3",
"@deepseek-ai/schemastery": "~3.19.0"
```


### v0.0.3 的两个运行期修复（实机暴露）

**① 插件看不到 `DSH_HOME`，Python 自动探测失效。**
原实现只用 `process.env.DSH_HOME` 定位 DSH 自带运行时。实测插件运行在
Host 进程里时**没有该变量**（它只注入给终端/工具会话），于是探测失败，
只能手工写死解释器路径。
→ 改为多来源：`DSH_HOME` → **`~/.dsh`（与进程环境无关的兜底）** →
`DSH_PROFILE_DIR` 的上两级；并补充 macOS Homebrew 与 `bin/python3` 布局。
`smoke-live.mjs` 现在会**主动删掉 `DSH_HOME`** 来复现 Host 环境。

**② 相对路径只按 `plannerDir` 解析，仓库布局下失效。**
`plannerDir` 默认改成插件包目录后，`testdata/mini.gb` 被解析成
`<插件>/testdata/mini.gb`（不存在），而文件实际在 `<repo>/testdata/`。
→ 改为在**插件目录 → 其上一级 → dataDir → 当前工作目录**依次探测，
命中即用；全部失败时报错列出试过的基准目录。

**③ 顺带：FASTA 读取改用 `utf-8-sig`。**
带 BOM 的文件首行会是 `\ufeff>id`，`startswith(">")` 为假，
整条首记录的头部会被当成序列。已加回归用例。

### 从 `0.2.0-rc.2` 升到 `0.2.2-rc.3` 的实际记录（v0.0.2）

只改版本号是**不够**的，必须重新核对 API。本次核对结果：

| 表面 | 结论 |
|---|---|
| `SubprocessSpawnSpec` | **未变**。`argv` / `cwd` / `stdio` / `graceMs` / `signal` / `env` 全同；`collected.readFrom` 与 `done` 同样未变 |
| `ToolDefinition`（含 `defineTool`） | **未变**，仍要求 `output.schema` + `output.render`；本插件本来就是这个形状 |
| 工具注册 | 未变。`ctx.tools.register` 仍是主入口 |
| 工具呈现 | 只接受 `native` / `ptc`，本插件不涉及 |
| Node 下限 | **变了**：`0.2.1-rc.2` 起要求 Node **24.9+**，`engines.node` 已同步 |
| 兼容性判定 | **变了**：改读 `engines.dsh` 优先（见上） |

**顺带修掉的运行期问题**：旧 `cordis.patch.yml` 把 `python` / `plannerDir`
写成绝对路径。这在别的机器上是死值，且 Windows 上 PATH 里的 `python`
可能是 3.9（CLI 需 3.10+）。现已改为**零配置自定位**：
`plannerDir` 用 `import.meta.url` 指向包内 `stib/`，
`python` 自动探测（`python3` → `python` → DSH 自带运行时，取第一个 3.10+）。

**升级 DSH 后必须重新核对。** 排查命令：

```bash
dsh plugin --profile desktop list
dsh plugin --profile desktop version-exemptions
```

---

## 5. 目录结构

```
dsh-plugin/
├── package.json            peer 精确锁版本；files 声明完整
├── tsconfig.json           tsc 构建配置
├── .gitignore
├── LICENSE                 MIT + 实验判断免责
├── README.md               安装/配置/兼容性/能力边界
├── PLUGIN_RELEASE.md       本文件
├── cordis.patch.yml        bundle 层（需填绝对路径）
├── src/index.ts            插件源码（含真实契约抄录）
├── lib/                    构建产物（tsc 生成，gitignore）
│   ├── index.js
│   └── types/index.d.ts
└── skills/
    └── stib-bacteroides-editing/SKILL.md
```

---

## 6. 发布检查清单

- [x] 源码按真实契约编写
- [x] `tsc` 编译零错误
- [x] 导出契约正确（无 `default`）
- [x] peer 版本为真实存在版本
- [x] LICENSE / README / .gitignore
- [ ] **生物安全拒绝面**（合规底线）
- [ ] **真实 DSH 会话激活 + 6 工具逐一验证**
- [ ] spill 路径在桌面版可读性验证
- [ ] 开源授权确认（论文与 SOP 的公开范围）
- [ ] `cordis.patch.yml` 绝对路径填好（或改为由使用者配置）

**在「生物安全拒绝面」与「真实会话验证」两项完成前，不建议对外发布。**
