/**
 * dsh-stib-planner — 把 stib-planner CLI 包装成 DSH 工具与技能。
 *
 * ═══════════════════════════════════════════════════════════════════════
 * 本文件的 subprocess 调用已按**官方真实契约**重写
 * ═══════════════════════════════════════════════════════════════════════
 *
 * 依据：`@deepseek-ai/dsh-subprocess/src/types.ts`（v0.2.0-rc.2 源码）。
 *
 *   interface SubprocessSpawnSpec {
 *     argv: readonly string[]          // argv[0] 是程序；**不做 shell 解释**
 *     cwd: string
 *     stdio: { stdin, stdout, stderr, control? }   // 全部显式，无默认值
 *     graceMs: number                  // 必填，正有限毫秒
 *     signal?: AbortSignal
 *     env?: NodeJS.ProcessEnv
 *   }
 *
 *   type SubprocessOutputMode = 'pipe' | 'inherit' | SubprocessCollect
 *   interface SubprocessCollect { maxBytes: number; spill?: { maxBytes: number } }
 *
 *   interface SubprocessHandle {
 *     stdin / stdout / stderr : 仅在对应 stdio 为 'pipe' 时存在
 *     collected: { stdout?: SubprocessOutputReader, stderr?: SubprocessOutputReader }
 *     done: Promise<SubprocessOutcome>   // { exitCode: number|null, signal: Signals|null }
 *     terminate(): void
 *     waitForExit(signal?): Promise<boolean>
 *   }
 *   interface SubprocessOutputReader { readFrom(fromByte: number): SubprocessOutputRead }
 *   // readFrom 返回 { text, nextOffset, lossy, spillPath? }；偏移是整流域坐标，多读者互不消费
 *
 * ⚠️ 本文件 v0.1 的写法是**错的**（且是纯猜测）：
 *    ✗ `spawn({ command, args, ... })` —— 真实契约只有 `argv`，且 `graceMs` 必填
 *    ✗ `handle.stdout` 当字节流 readAll —— collect 模式没有 `.stdout`，要用 `.collected`
 *    ✗ `handle.exited` —— 契约里是 `handle.done`
 * 这些错误在真实 DSH 中会**直接抛错或静默拿到空输出**。已在 v0.2 全部修正。
 */
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { homedir } from 'node:os';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import z from '@deepseek-ai/schemastery';
import { defineTool } from '@deepseek-ai/dsh-tools';
export const name = 'stib-planner';
export const inject = ['tools', 'subprocess'];
/**
 * 插件包根目录（`lib/index.js` 的上一级）。
 *
 * 打包时 Python 代码被复制到包内 `stib/`，所以这里同时就是
 * `python -m stib.cli` 的工作目录 —— **不需要用户配置路径**。
 * 用 `import.meta.url` 推导，装到任何机器、任何目录都成立。
 */
const PLUGIN_ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
/** 插件自带的 Python 包所在目录（`python -m stib.cli` 的 cwd）。 */
const BUNDLED_PYTHON_DIR = PLUGIN_ROOT;
/** 最低可用 Python 版本。3.9 会在 CLI 启动时直接失败。 */
const MIN_PYTHON = [3, 10];
export const Config = z.object({
    python: z.string().default(''),
    plannerDir: z.string().default(''),
    dataDir: z.string().default(''),
    timeoutMs: z.number().default(180000),
    maxOutputBytes: z.number().default(1_048_576),
});
/**
 * 运行一次 stib CLI。
 *
 * 用 **collect 模式**（`{ maxBytes }`）取输出，在 `done` 之后用 `readFrom(0)`
 * 一次性读回 —— 这比手搓 pipe 流更省事，且溢出时保留尾部并提供 spill 路径。
 * 偏移是整流域坐标，readFrom 不消费流，因此可重复读。
 */
async function runCli(ctx, config, args, signal) {
    const graceMs = config.timeoutMs ?? 180000;
    const maxBytes = config.maxOutputBytes ?? 1_048_576;
    const plannerDir = resolvePlannerDir(config);
    const python = await resolvePython(ctx, config);
    const handle = ctx.subprocess.spawn({
        argv: [python, '-m', 'stib.cli', ...args],
        cwd: plannerDir,
        stdio: {
            stdin: 'ignore',
            stdout: { maxBytes },
            stderr: { maxBytes },
        },
        graceMs,
        signal,
        // 让子进程按 UTF-8 输出，避免中文报告在 Windows 默认代码页下乱码
        env: { PYTHONIOENCODING: 'utf-8', PYTHONUTF8: '1' },
    });
    // ⚠️ 契约的 graceMs **不是 deadline** —— 它只作用于 terminate 流程与排水。
    //    所以时限要自己管：到点调用 terminate()。
    const timer = setTimeout(() => handle.terminate(), graceMs);
    let exitCode = null;
    try {
        const outcome = await handle.done;
        exitCode = outcome.exitCode;
    }
    finally {
        clearTimeout(timer);
    }
    const readAll = (r) => {
        if (!r)
            return '';
        const out = r.readFrom(0);
        return out.lossy && out.spillPath
            ? `${out.text}\n[输出被截断；完整内容见 ${out.spillPath}]`
            : out.text;
    };
    const stdout = readAll(handle.collected.stdout);
    const stderr = readAll(handle.collected.stderr);
    return { ok: exitCode === 0, stdout, stderr, exitCode };
}
/** 已知的基因组文件扩展名（与 stib/genome.py 的 GENBANK_EXTS/FASTA_EXTS 保持一致）。 */
const GENOME_EXTS = ['.gb', '.gbk', '.gbff', '.genbank', '.fasta', '.fa', '.fna', '.fas'];
/**
 * 把用户给的**输入文件路径**解析成绝对路径。
 *
 * 为什么要这一步：CLI 的工作目录是插件目录，用户写 `testdata/x.fa` 这种
 * 相对路径时，若原样透传就会在插件目录下找不到。而 `accession` 参数
 * 早就做了多基点探测 —— 两者行为必须一致，否则同一个相对路径有的参数
 * 能用、有的报 FileNotFound。
 *
 * 解析顺序：绝对路径 / `~` → 插件目录 → 其上一级（仓库布局）→ 当前工作目录。
 * 全都不存在时返回"插件目录下的那个"，让 Python 报出更具体的错。
 */
function resolveInputPath(config, p) {
    const raw = (p ?? '').trim();
    if (!raw)
        return raw;
    if (/^([A-Za-z]:[\\/]|\/)/.test(raw))
        return raw;
    if (raw.startsWith('~')) {
        try {
            return join(homedir(), raw.slice(1));
        }
        catch {
            // 取不到 home 就按相对路径继续
        }
    }
    const bases = [resolvePlannerDir(config)];
    const parent = dirname(resolvePlannerDir(config));
    if (parent && parent !== resolvePlannerDir(config))
        bases.push(parent);
    try {
        const cwd = process.cwd();
        if (cwd)
            bases.push(cwd);
    }
    catch {
        // 忽略
    }
    for (const b of [...new Set(bases)]) {
        const cand = join(b, raw);
        if (fileExists(cand))
            return cand;
    }
    return join(bases[0], raw);
}
/** 输出目录：相对路径按插件目录解析；留空则用 `<插件目录>/out`。 */
function resolveOutputDir(config, p) {
    const raw = (p ?? '').trim();
    if (!raw)
        return join(resolvePlannerDir(config), 'out');
    if (/^([A-Za-z]:[\\/]|\/)/.test(raw))
        return raw;
    if (raw.startsWith('~')) {
        try {
            return join(homedir(), raw.slice(1));
        }
        catch {
            // fallthrough
        }
    }
    return join(resolvePlannerDir(config), raw);
}
/**
 * 输出**文件**路径：与 [`resolveOutputDir`] 同一套规则。
 *
 * 为什么要单独一个名字：语义不同（文件 vs 目录），但解析逻辑相同 ——
 * 关键是**必须解析**。否则相对路径会被 CLI 按插件目录解释，
 * 而返回给用户的是相对路径字符串，用户按那个路径去找必然找不到文件。
 */
function resolveOutputPath(config, p) {
    return resolveOutputDir(config, p);
}
/**
 * 把用户给的 `accession` 解析成实际基因组文件路径。
 *
 * 支持三种写法 —— 早期版本只会做第 3 种，于是 `.gb` 路径永远进不去：
 *
 *   1. **绝对路径**（含 .gb/.fasta 等扩展名）→ 原样使用
 *   2. **相对路径 / 带目录** → 在插件目录、其上一级（仓库布局）、dataDir、
 *      当前工作目录下依次探测，命中即用
 *   3. **裸 accession**（如 `NC_004663.1`）→ 在 `dataDir` 下按已知扩展名探测：
 *        - 命中唯一文件 → 用它（**这样 .gb / .gbk 缓存也能被发现**）
 *        - 全都没命中 → 报错并列出找过的路径，而不是静默拼一个 `.fasta`
 */
function resolveGenome(config, accession) {
    const dataDir = resolveDataDir(config);
    const a = accession.trim();
    // 0) `~` 展开 —— macOS / Linux 上最常见的写法（Windows 用户一般不会写）
    if (a === '~' || a.startsWith('~/') || a.startsWith('~\\')) {
        try {
            return join(homedir(), a.slice(1));
        }
        catch {
            // 取不到 home 就按原样继续，让后面的分支产生可诊断的报错
        }
    }
    // 1) 绝对路径（Windows 盘符 或 POSIX 根）
    if (/^([A-Za-z]:[\\/]|\/)/.test(a))
        return a;
    const lower = a.toLowerCase();
    // 2) 带目录的相对路径 —— 在多个基点下探测，命中即用
    //
    //    ⚠️ 不能只相对 plannerDir：插件默认的 plannerDir 是**插件包目录**，
    //    而仓库布局里示例数据在 `<repo>/testdata/`（包目录的上一级）。
    //    只用一个基点会让 `testdata/mini.gb` 这种写法直接失效。
    if (a.includes('/') || a.includes('\\')) {
        const pluginDir = resolvePlannerDir(config);
        const bases = [pluginDir];
        const parent = dirname(pluginDir);
        if (parent && parent !== pluginDir)
            bases.push(parent); // 仓库布局
        bases.push(dataDir);
        try {
            const cwd = process.cwd();
            if (cwd)
                bases.push(cwd); // 相对当前工作目录
        }
        catch {
            // 取不到 cwd 就跳过
        }
        for (const base of [...new Set(bases)]) {
            const cand = join(base, a);
            if (fileExists(cand))
                return cand;
        }
        throw new Error(`找不到基因组文件：${a}\n` +
            `已尝试的基准目录：\n` +
            [...new Set(bases)].map((b) => `  - ${b}`).join('\n') +
            `\n\n最省事的做法是直接给出**绝对路径**（.gb / .fasta 均可）。`);
    }
    // 3) 裸文件名：已带已知扩展名 → 去 dataDir 找
    if (GENOME_EXTS.some((e) => lower.endsWith(e))) {
        return join(dataDir, a);
    }
    // 4) 裸 accession → 探测 dataDir 下是否存在任一已知扩展名
    for (const ext of GENOME_EXTS) {
        const cand = join(dataDir, `${a}${ext}`);
        if (fileExists(cand))
            return cand;
    }
    // 探测 `a` 本身就是完整文件名（无扩展名时也可能存在）
    if (fileExists(join(dataDir, a)))
        return join(dataDir, a);
    throw new Error(`找不到基因组文件。已尝试在 ${dataDir} 下查找 ` +
        GENOME_EXTS.map((e) => `${a}${e}`).join(' / ') +
        `。\n请在 dataDir 里放好基因组文件，或直接给出**文件路径**（.gb / .fasta 均可）。`);
}
/** 文件是否存在（用 ESM 导入的 existsSync；插件运行在 Host 进程里）。 */
function fileExists(p) {
    try {
        return existsSync(p);
    }
    catch {
        return false;
    }
}
/** 统一的「工具失败」文本。 */
function failText(r) {
    return [
        `stib-planner 调用失败（exitCode=${r.exitCode ?? 'null'}）`,
        r.stderr.trim() ? `\n--- stderr ---\n${r.stderr.trim()}` : '',
        r.stdout.trim() ? `\n--- stdout ---\n${r.stdout.trim()}` : '',
    ].join('');
}
/**
 * `stib/` 所在目录 —— 运行方式是 `python -m stib.cli`，包名即目录名，
 * 所以该目录**必须直接含 `stib/`**，中隔一层会 ModuleNotFoundError。
 *
 * 默认就是插件包根目录：打包时 Python 代码被复制进包内 `stib/`，
 * 因此**装到任何机器都无需配置**。只有想改用仓库里另一份代码时才需要覆盖。
 */
function resolvePlannerDir(config) {
    const configured = (config.plannerDir ?? '').trim();
    if (configured && fileExists(join(configured, 'stib', 'cli.py')))
        return configured;
    // 1) 包内自带：`npm pack` 时 stage-python 把 Python 代码复制进 dsh-plugin/stib
    if (fileExists(join(BUNDLED_PYTHON_DIR, 'stib', 'cli.py')))
        return BUNDLED_PYTHON_DIR;
    // 2) 仓库根布局：从 git 安装的包，Python 代码在 `<包根>/stib`（与 dsh-plugin/ 同级）
    //    —— 这样 git 安装不必把 Python 代码复制两份。
    const parent = dirname(BUNDLED_PYTHON_DIR);
    if (fileExists(join(parent, 'stib', 'cli.py')))
        return parent;
    return BUNDLED_PYTHON_DIR;
}
/** 基因组缓存目录。默认放在生效的 planner 目录下的 `data/`。 */
function resolveDataDir(config) {
    const configured = (config.dataDir ?? '').trim();
    return configured || join(resolvePlannerDir(config), 'data');
}
/**
 * DSH 主目录（可能多个），用于定位**自带 Python 运行时**。
 *
 * ⚠️ 不能只依赖 `DSH_HOME`：实测插件运行在 Host 进程里时**看不到该变量**
 *    （它只注入给终端/工具会话），所以必须有与进程环境无关的兜底 ——
 *    DSH 的默认 home 就是 `~/.dsh`。
 *
 * 依次来源：
 *   1. `DSH_HOME` 环境变量（若可见）
 *   2. `~/.dsh`（Windows / macOS / Linux 的默认位置）
 *   3. `DSH_PROFILE_DIR` 的上两级（`<home>/profiles/<name>`）
 */
function dshHomes() {
    const out = [];
    const envHome = (process.env.DSH_HOME ?? '').trim();
    if (envHome)
        out.push(envHome);
    try {
        const h = homedir();
        if (h)
            out.push(join(h, '.dsh'));
    }
    catch {
        // 极端情况下取不到 home，忽略
    }
    const profileDir = (process.env.DSH_PROFILE_DIR ?? '').trim();
    if (profileDir)
        out.push(dirname(dirname(profileDir)));
    return [...new Set(out)];
}
/** `dsh-runtimes/<runtime>/dependencies/python/` 下解释器的相对路径（覆盖两平台布局）。 */
const PYTHON_RELATIVE = [
    ['dependencies', 'python', 'python.exe'], // Windows 打包版
    ['dependencies', 'python', 'python'], // 无扩展名的 POSIX 布局
    ['dependencies', 'python', 'bin', 'python3'], // macOS / Linux 打包版
    ['dependencies', 'python', 'bin', 'python'],
];
/**
 * 候选 Python 解释器，按优先级排列。
 *
 * 显式配置时**只用它** —— 用户说了算，失败也要报清楚，不能悄悄换一个。
 *
 * ⚠️ 为什么必须探测而不是直接用 `python`：
 *    - Windows 上 PATH 里的 `python` 常见是 **3.9**（本机实测 3.9.2）
 *    - macOS 上 Xcode 命令行工具自带的 `python3` 是 **3.9.6**
 *    两者都低于 CLI 要求的 3.10，直接用会莫名失败。
 *    DSH 自带运行时是 3.12，所以把它列为主要候选。
 */
function pythonCandidates(config) {
    const explicit = (config.python ?? '').trim();
    if (explicit)
        return [explicit];
    const list = ['python3', 'python'];
    // DSH 自带运行时（枚举 runtime 目录，不写死名字）
    for (const home of dshHomes()) {
        const runtimesDir = join(home, 'dsh-runtimes');
        let runtimes;
        try {
            runtimes = readdirSync(runtimesDir);
        }
        catch {
            continue; // 没有 dsh-runtimes 就跳过
        }
        for (const rt of runtimes) {
            for (const rel of PYTHON_RELATIVE) {
                const p = join(runtimesDir, rt, ...rel);
                if (fileExists(p))
                    list.push(p);
            }
        }
    }
    // macOS 上常见的、比系统自带更新的解释器
    if (process.platform === 'darwin') {
        list.push('/opt/homebrew/bin/python3', // Apple Silicon Homebrew
        '/usr/local/bin/python3');
    }
    // 带小版本的命令名（PATH 上有多个 Python 时往往能命中）
    list.push('python3.13', 'python3.12', 'python3.11', 'python3.10');
    return [...new Set(list)];
}
/** 探测成功的解释器缓存（避免每次调用都探测）。 */
let resolvedPython;
/** 用 `-c` 探一次版本，够新才算通过。 */
async function probePython(ctx, candidate, signal) {
    const probe = `import sys; sys.exit(0 if sys.version_info >= (${MIN_PYTHON[0]}, ${MIN_PYTHON[1]}) else 7)`;
    try {
        const handle = ctx.subprocess.spawn({
            argv: [candidate, '-c', probe],
            cwd: BUNDLED_PYTHON_DIR,
            stdio: { stdin: 'ignore', stdout: 'inherit', stderr: 'inherit' },
            graceMs: 15000,
            signal,
        });
        const timer = setTimeout(() => handle.terminate(), 15000);
        try {
            const outcome = await handle.done;
            return outcome.exitCode === 0;
        }
        finally {
            clearTimeout(timer);
        }
    }
    catch {
        return false;
    }
}
/** 解析出可用的 Python 解释器（带缓存）。 */
async function resolvePython(ctx, config, signal) {
    const explicit = (config.python ?? '').trim();
    if (explicit)
        return explicit; // 用户显式指定：不探测、不替换
    if (resolvedPython)
        return resolvedPython;
    const tried = [];
    for (const candidate of pythonCandidates(config)) {
        tried.push(candidate);
        if (await probePython(ctx, candidate, signal)) {
            resolvedPython = candidate;
            return candidate;
        }
    }
    throw new Error(`找不到 Python ${MIN_PYTHON[0]}.${MIN_PYTHON[1]}+ 解释器。已尝试：\n` +
        tried.map((t) => `  - ${t}`).join('\n') +
        `\n\n已查找的 DSH 主目录：\n` +
        dshHomes().map((h) => `  - ${h}`).join('\n') +
        `\n\n解决办法：安装 Python 3.10+（创建虚拟环境或用系统包管理器也可），\n` +
        `或在插件配置里显式指定解释器路径，例如：\n` +
        `  python: /usr/local/bin/python3        # macOS\n` +
        `  python: C:\\\\Python312\\\\python.exe       # Windows\n` +
        `\n注意：Windows 上 \`py\` 启动器与 macOS 上系统自带的 python3\n` +
        `可能是 3.9，低于本工具要求，需要另装较新版本。`);
}
/**
 * 拆分 SKILL.md 的 YAML frontmatter。
 *
 * 只需要平铺的 `key: value`（本仓库的 SKILL.md 即此形式）；不做通用 YAML 解析，
 * 以免为一处用法引入解析器。找不到 frontmatter 时整个文件都算正文。
 */
function splitFrontmatter(raw) {
    const text = raw.replace(/^\uFEFF/, '');
    const match = /^---\r?\n([\s\S]*?)\r?\n---[ \t]*\r?\n?/.exec(text);
    if (match === null || match[1] === undefined)
        return { meta: {}, body: text };
    const meta = {};
    for (const line of match[1].split(/\r?\n/)) {
        const kv = /^([A-Za-z][A-Za-z0-9_-]*)[ \t]*:[ \t]*(.*)$/.exec(line);
        if (kv !== null && kv[1] !== undefined)
            meta[kv[1]] = (kv[2] ?? '').trim();
    }
    return { meta, body: text.slice(match[0].length) };
}
/** 注册随包技能；任何失败都只少一个技能，不影响工具。 */
function registerBundledSkill(ctx) {
    const registrar = ctx.skills;
    if (registrar === undefined || typeof registrar.register !== 'function')
        return;
    const dir = join(PLUGIN_ROOT, 'skills', 'stib-bacteroides-editing');
    const file = join(dir, 'SKILL.md');
    if (!existsSync(file))
        return;
    try {
        const { meta, body } = splitFrontmatter(readFileSync(file, 'utf8'));
        registrar.register({
            name: meta['name'] ?? 'stib-bacteroides-editing',
            description: meta['description']
                ?? '帮无经验的课题组完成拟杆菌（Bacteroides）的 STIB 基因编辑：sgRNA 设计、脱靶评估、引物、接合流程与排查。',
            ...(meta['whenToUse'] === undefined ? {} : { whenToUse: meta['whenToUse'] }),
            content: body,
            source: 'bundled',
            resourceBase: { kind: 'directory', path: dir },
        });
    }
    catch {
        // 技能注册失败不得影响 9 个工具。
    }
}
export function apply(ctx, config) {
    // ---- tool 1: 生成完整编辑方案 ----
    ctx.tools.register(defineTool({
        name: 'stib_plan',
        description: '为拟杆菌（Bacteroides）生成一份完整的 STIB 基因编辑方案：sgRNA 候选排序（PAM 优选 GTT→GTA）、'
            + '全基因组脱靶评估、可直接下单的寡核苷酸、验证引物、分阶段流程、判定点与风险清单。'
            + '注意：不预测编辑效率；脱靶引擎为 6-mer 块 + 替换枚举（≤4 错配保证不漏检，已交叉验证）；'
            + '未做湿实验验证。',
        parameters: {
            gene: { type: 'string', required: true, description: '目标基因名，如 tdk。**若基因组为 GenBank 格式，此名同时用于自动定位**；'
                    + '若为 FASTA，则仅作标注，需另给 motif 或 start/end' },
            accession: { type: 'string', required: true, description: '参考基因组 accession，如 NC_004663.1（须已取回；.fasta 或 .gb 均可）' },
            motif: { type: 'string', description: '用于定位目标区域的一段序列；与 start/end 二选一' },
            start: { type: 'integer', description: '目标区间起点（0-based），需与 end 同用' },
            end: { type: 'integer', description: '目标区间终点（0-based, exclusive），需与 start 同用' },
            species: { type: 'string', description: '目标物种；影响氯霉素浓度与接合气氛推荐' },
            top: { type: 'integer', description: '返回的 sgRNA 候选数，默认 5' },
            maxMismatch: { type: 'integer', description: '脱靶枚举错配上限，默认 4（≤4 保证不漏检）' },
            community: { type: 'boolean', description: '是否菌群编辑（接合气氛转为厌氧推荐）' },
            detailed: { type: 'boolean', description: 'true = 输出**完整定制实验协议**（动手前核对清单、Golden Gate 与 Gibson 的'
                    + '精确体系体积与程序、接合分步、鉴定流程、时间轴、试剂清单），'
                    + '而不是决策摘要。**给要做实验的人用这个。**' },
            out: { type: 'string', description: '把方案写到这个 .md 路径（留空则只返回文本）。'
                    + '建议同时给 html，方便在浏览器里看图谱。' },
            html: { type: 'string', description: '同时输出**自包含 HTML 可视化报告**到该路径（内联 SVG，双击即可打开，无外部依赖）。'
                    + '含候选位置图、评分条形图、脱靶明细。' },
            jsonPath: { type: 'string', description: '同时输出 JSON 到该路径（供程序消费）' },
        },
        output: {
            schema: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    ok: { type: 'boolean', required: true },
                    content: { type: 'string', required: true },
                },
            },
            render: (_args, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            // ⚠️ 这里**不要**自己做「必须给 motif 或坐标」的前置检查。
            //    Python 侧的 resolve_region 已支持按基因名定位（GenBank），
            //    插件层抢先抛错会让那条路径永远走不到。
            //    缺什么由 Python 报错，且它的提示更具体（会说明该改用 GenBank 还是补坐标）。
            const cli = ['plan', '--gene', args.gene, '--genome', resolveGenome(config, args.accession)];
            if (args.motif)
                cli.push('--motif', args.motif);
            if (args.start !== undefined && args.end !== undefined) {
                cli.push('--start', String(args.start), '--end', String(args.end));
            }
            if (args.species)
                cli.push('--species', args.species);
            if (args.top !== undefined)
                cli.push('--top', String(args.top));
            if (args.maxMismatch !== undefined)
                cli.push('--max-mismatch', String(args.maxMismatch));
            if (args.community)
                cli.push('--community');
            if (args.detailed)
                cli.push('--detailed');
            // 若只要 md 却什么都没给，则连 html 一起出，落在 plannerDir/out/
            // ⚠️ 所有输出路径都解析成绝对路径：CLI 的 cwd 是插件目录，若原样透传
            //    相对路径，文件会落到 `<插件>/out/...`，而返回给用户的却是相对路径 ——
            //    用户按那个路径去找必然找不到。
            const outPath = args.out ? resolveOutputPath(config, args.out) : '';
            const htmlPath = args.html ? resolveOutputPath(config, args.html) : '';
            const jsonPath = args.jsonPath ? resolveOutputPath(config, args.jsonPath) : '';
            if (!outPath && !htmlPath && !jsonPath) {
                cli.push('--out', `${resolvePlannerDir(config)}/out/plan.md`);
                cli.push('--html', `${resolvePlannerDir(config)}/out/plan.html`);
            }
            if (outPath)
                cli.push('--out', outPath);
            if (htmlPath)
                cli.push('--html', htmlPath);
            if (jsonPath)
                cli.push('--json', jsonPath);
            const r = await runCli(ctx, config, cli, exec.signal);
            if (!r.ok)
                return { ok: false, content: failText(r) };
            // 产出文件时，在末尾明确列出路径 —— 用**解析后**的绝对路径
            const wanted = [outPath, htmlPath, jsonPath].filter(Boolean);
            const tail = wanted.length
                ? `\n\n---\n**已写出文件**\n` +
                    wanted.map((f) => `- \`${f}\``).join('\n') +
                    `\n\n> \`.html\` 双击即可打开（无需联网）；\`.md\` 可直接通读或贴进实验记录本。`
                : '';
            return { ok: true, content: r.stdout + tail };
        },
    }));
    // ---- tool 2: 只列 sgRNA 候选 ----
    ctx.tools.register(defineTool({
        name: 'stib_find_sgrnas',
        description: '扫描指定基因组区域的候选 STIB sgRNA（GTN PAM，23 nt spacer），返回候选列表。'
            + '仅做序列层面筛选与块唯一性检查；完整脱靶评估请用 stib_offtarget。',
        parameters: {
            accession: { type: 'string', required: true, description: '基因组：accession 或**文件路径**（.gb / .fasta 均可）' },
            gene: { type: 'string', description: '按基因名定位（需 GenBank 注释）' },
            motif: { type: 'string', description: '目标区域定位序列' },
            start: { type: 'integer' },
            end: { type: 'integer' },
            top: { type: 'integer', description: '显示条数，默认 20' },
        },
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            const cli = ['scan', '--genome', resolveGenome(config, args.accession)];
            if (args.gene)
                cli.push('--gene', args.gene);
            if (args.motif)
                cli.push('--motif', args.motif);
            if (args.start !== undefined && args.end !== undefined) {
                cli.push('--start', String(args.start), '--end', String(args.end));
            }
            if (args.top !== undefined)
                cli.push('--top', String(args.top));
            const r = await runCli(ctx, config, cli, exec.signal);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- tool 3: 全基因组脱靶评估 ----
    ctx.tools.register(defineTool({
        name: 'stib_offtarget',
        description: '对一条或多条 STIB spacer 做全基因组脱靶评估（正负链，≤N 错配枚举），并按错配位置分级：'
            + 'helix（PAM 近端前 17 nt，识别关键区 → 最危险）/ tail（3\' 端外周 6 nt）。'
            + '**任何错配都只增风险，不存在减害**。'
            + '完全相同的额外拷贝会被单独标为串联重复。'
            + '引擎为 6-mer 块 ×3 + 每块 1 处替换枚举：≤4 错配保证不漏检（已用暴力穷举交叉验证）；'
            + '>4 时保证失效并显式警示，不等价于 Cas-OFFinder 穷举。',
        parameters: {
            accession: { type: 'string', required: true, description: '基因组：accession 或**文件路径**（.gb / .fasta 均可）' },
            spacer: { type: 'string', description: '待评估的 23 nt spacer；与 panel 至少提供一个' },
            panel: { type: 'boolean', description: '同时评估内置的已验证 spacer 面板' },
            maxMismatch: { type: 'integer', description: '错配上限，默认 4（≤4 保证不漏检）' },
            top: { type: 'integer', description: '每条最多显示几个命中，默认 8' },
        },
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            if (!args.spacer && !args.panel)
                throw new Error('必须提供 spacer，或设 panel=true');
            const cli = ['offtarget', '--genome', resolveGenome(config, args.accession)];
            if (args.spacer)
                cli.push('--spacer', args.spacer);
            if (args.panel)
                cli.push('--panel');
            if (args.maxMismatch !== undefined)
                cli.push('--max-mismatch', String(args.maxMismatch));
            if (args.top !== undefined)
                cli.push('--top', String(args.top));
            const r = await runCli(ctx, config, cli, exec.signal);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- tool 4: 批量设计 ----
    ctx.tools.register(defineTool({
        name: 'stib_batch',
        description: '批量设计模式：在一个较宽区间（如整个基因簇）内一次给出多条候选 sgRNA，'
            + '含 PAM 优选、落点分级与序列偏好标记。适合「还不知道该改哪个位点」的场景。',
        parameters: {
            accession: { type: 'string', required: true, description: '基因组：accession 或**文件路径**（.gb / .fasta 均可）' },
            gene: { type: 'string', required: true, description: '区域名（用于标注）' },
            motif: { type: 'string', description: '定位序列' },
            start: { type: 'integer' },
            end: { type: 'integer' },
            pad: { type: 'integer', description: 'motif 模式下的两侧扩展，默认 2000' },
            top: { type: 'integer', description: '返回候选数，默认 24' },
            species: { type: 'string' },
        },
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            const cli = ['batch', '--genome', resolveGenome(config, args.accession), '--gene', args.gene];
            if (args.motif)
                cli.push('--motif', args.motif);
            if (args.start !== undefined && args.end !== undefined) {
                cli.push('--start', String(args.start), '--end', String(args.end));
            }
            if (args.pad !== undefined)
                cli.push('--pad', String(args.pad));
            if (args.top !== undefined)
                cli.push('--top', String(args.top));
            if (args.species)
                cli.push('--species', args.species);
            const r = await runCli(ctx, config, cli, exec.signal);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- tool 5: 列出已核实的材料矛盾 ----
    ctx.tools.register(defineTool({
        name: 'stib_open_questions',
        description: '列出 STIB 知识中已核实的材料矛盾（如血液种类、插入窗口两处写法不一致）。'
            + '自 v0.4 起接合气氛不再是未收敛变量 —— 已按菌株耐氧能力确定规则。',
        parameters: {},
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute() {
            const r = await runCli(ctx, config, ['conflicts']);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- tool 6: 自检 ----
    ctx.tools.register(defineTool({
        name: 'stib_selftest',
        description: '运行 stib-planner 自检：用论文已发表的 sgRNA 序列核对序列层，校验知识库的来源锚定与冲突显式约束，'
            + '并用合成底物核对脱靶引擎。首次使用或怀疑环境有问题时调用。',
        parameters: {},
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute() {
            const r = await runCli(ctx, config, ['selftest']);
            return { ok: r.ok, content: r.stdout || r.stderr };
        },
    }));
    // ---- tool 7: 列出/搜索基因（需 GenBank 注释）----
    ctx.tools.register(defineTool({
        name: 'stib_find_gene',
        description: '列出或搜索基因组注释里的基因（gene / locus_tag / product）。'
            + '**当用户说的基因名不确定时，先用这个查准确名字。**'
            + '需要 GenBank 格式基因组；纯 FASTA 没有注释，无法列出基因。',
        parameters: {
            accession: { type: 'string', required: true, description: '基因组：可以是 accession（如 NC_004663.1，在 dataDir 下探测已知扩展名），'
                    + '也可以是**文件路径**（相对 plannerDir 或绝对路径；.gb / .gbk / .fasta 均可）' },
            search: { type: 'string', description: '关键词，匹配 gene / locus_tag / product；留空则列出前 N 条' },
            top: { type: 'integer', description: '显示条数，默认 40' },
        },
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            const cli = ['genes', '--genome', resolveGenome(config, args.accession)];
            if (args.search)
                cli.push('--search', args.search);
            if (args.top !== undefined)
                cli.push('--top', String(args.top));
            const r = await runCli(ctx, config, cli, exec.signal);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- tool 8: 中性位点（插入大片段用）----
    ctx.tools.register(defineTool({
        name: 'stib_neutral_sites',
        description: '查找**中性位点** —— 用于插入大片段 cargo（外源功能基因、通路、报告基因），'
            + '目标是插入不破坏宿主功能。返回按注释与序列特征排序的候选（含评分理由与告警），'
            + '并为每个位点给出保证「预测插入位点落在该区间内」的 sgRNA。'
            + '⚠️ 论文中**已验证**的中性位点会优先列出 —— 它们的价值高于任何序列特征预测。'
            + '⚠️ 需要 GenBank 注释：纯 FASTA 无法区分基因与基因间区。'
            + '⚠️ 这些是候选，不是已验证位点；插入后须确认宿主表型未变。',
        parameters: {
            accession: { type: 'string', required: true, description: '基因组：可以是 accession（如 NC_004663.1，在 dataDir 下探测已知扩展名），'
                    + '也可以是**文件路径**（相对 plannerDir 或绝对路径；.gb / .gbk / .fasta 均可）' },
            minGap: { type: 'integer', description: '最小基因间区长度，默认 120 bp' },
            top: { type: 'integer', description: '候选位点数，默认 8' },
            sgPerSite: { type: 'integer', description: '每个位点显示几条 sgRNA，默认 3' },
            out: { type: 'string', description: '可选：输出 Markdown 路径' },
        },
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            const cli = ['neutral', '--genome', resolveGenome(config, args.accession)];
            if (args.minGap !== undefined)
                cli.push('--min-gap', String(args.minGap));
            if (args.top !== undefined)
                cli.push('--top', String(args.top));
            if (args.sgPerSite !== undefined)
                cli.push('--sg', String(args.sgPerSite));
            if (args.out)
                cli.push('--out', resolveOutputPath(config, args.out));
            const r = await runCli(ctx, config, cli, exec.signal);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- tool 9: cargo 的 Gibson 组装引物 + 可视化图谱 ----
    ctx.tools.register(defineTool({
        name: 'stib_cargo_assembly',
        description: '为插入 cargo 设计 Gibson 组装引物，并输出可在 SnapGene 打开的图谱文件。'
            + '**三片段组装共需 6 条引物**：F1 = GA725 + Cargo1、F2 = cargo、F3 = Cargo2 + GA778'
            + '（GA725 与 GA778 **不是一对**），并逐 junction 检查三处同源区（J1/J2/J3）。'
            + '**硬性校验 cargo 是否落在转座子 LE 与 RE 之间** —— 放错位置会破坏转座子末端识别'
            + '（转座失败，而 PCR、Gibson、菌落 PCR 全都正常，直到接合才发现转化子为零）。'
            + '同时检查 cargo 内的 BsaI 位点（会干扰本载体的 Golden Gate 换 sgRNA）。'
            + '产出 .dna（SnapGene）/ .gb（保底）/ .fa（构建序列）。'
            + '**两个可选插入位点**：after_le（整段替换 LE–RE 之间，cargo 长度无限制，'
            + '但失去 CmR 筛选）与 after_selection（rrnB T1 之后 / NspI 位点附近，保留 CmR 筛选）。'
            + '**cargo 若需表达，应使用 cargoStrand = -1 与 CmR 同向。**'
            + '输出会完整保留骨架的全部注释（sgRNA 盒、酶、启动子、LE/RE、CmR、AmpR、ori）。'
            + '⚠️ Tm 为估算值，首次须梯度 PCR 确认。',
        parameters: {
            cargoFasta: { type: 'string', required: true, description: 'cargo 的 FASTA 文件绝对路径（须只含一条序列）' },
            backbone: { type: 'string', required: true, description: '载体文件绝对路径（.dna 或 .gb）' },
            spec: { type: 'string', description: '构建规范 YAML 路径；留空用内置 psh43 规范' },
            site: { type: 'string', enum: ['after_le', 'after_selection'], description: '插入位点。'
                    + 'after_le = 整段替换 LE 与 RE 之间的内容（8,590–9,784），'
                    + '**cargo 长度无限制**，但会一并失去 PcfxA-CmR（氯霉素筛选）；'
                    + 'after_selection = 插在 rrnB T1 终止子之后 / NspI 位点附近（8,846–8,888），'
                    + '**保留两个 rrnB 终止子与 CmR 筛选**。默认 after_le' },
            cargoStrand: { type: 'integer', enum: [1, -1], description: 'cargo 所在链。-1 = 与 CmR 同向（本载体 CmR 在负链）。'
                    + '**若 cargo 需要表达，应与 CmR 同向，即 -1。** 默认取位点建议值' },
            listSites: { type: 'boolean', description: 'true = 只列出可选插入位点与各自可用空间，不做设计' },
            overlapLen: { type: 'integer', description: 'Gibson overlap 长度，默认 22（NEB 推荐 15–25）' },
            outDir: { type: 'string', description: '输出目录，默认 <plannerDir>/out' },
        },
        output: {
            schema: {
                type: 'object', additionalProperties: false,
                properties: { ok: { type: 'boolean', required: true }, content: { type: 'string', required: true } },
            },
            render: (_a, v) => [{ type: 'text', text: v.content }],
        },
        async execute(args, exec) {
            // ⚠️ 三个路径参数都必须**解析成绝对路径**再交给 CLI ——
            //    否则 CLI 会在插件目录下找相对路径（`testdata/xxx`）而必然失败。
            //    早期版本只解析了 `accession`（基因组），cargoFasta/backbone 原样透传，
            //    于是"同样的相对路径写法，accession 能用、这两个却报 FileNotFound"。
            const cli = [
                'cargo',
                '--cargo', resolveInputPath(config, args.cargoFasta),
                '--backbone', resolveInputPath(config, args.backbone),
            ];
            if (args.spec)
                cli.push('--spec', resolveInputPath(config, args.spec));
            if (args.site)
                cli.push('--site', args.site);
            if (args.cargoStrand !== undefined)
                cli.push('--cargo-strand', String(args.cargoStrand));
            if (args.listSites)
                cli.push('--list-sites');
            if (args.overlapLen !== undefined)
                cli.push('--overlap-len', String(args.overlapLen));
            // 输出目录同理：相对路径按 plannerDir 解析，避免落到意料之外的位置
            cli.push('--out', resolveOutputDir(config, args.outDir));
            const r = await runCli(ctx, config, cli, exec.signal);
            return { ok: r.ok, content: r.ok ? r.stdout : failText(r) };
        },
    }));
    // ---- 技能：让随包的 SKILL.md 进入技能目录 ----
    // 用 ctx.inject 而不是加进模块级 inject：skills 缺失时本插件仍应正常提供 9 个工具。
    //
    // ⚠️ 防御：`ctx.inject` 是 cordis 的可选服务注入 API。极少数精简宿主
    //    （以及冒烟测试的 mock ctx）没有它 —— 此时**只少一个技能**，
    //    9 个工具必须照常注册，不能让整个 apply() 抛错。
    const injectFn = ctx.inject;
    if (typeof injectFn === 'function') {
        try {
            injectFn.call(ctx, ['skills'], (scoped) => {
                registerBundledSkill(scoped);
            });
        }
        catch {
            // 技能通道不可用时静默降级
        }
    }
}
