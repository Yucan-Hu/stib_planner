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
import type { Context } from '@deepseek-ai/cordis';
import z from '@deepseek-ai/schemastery';
export declare const name = "stib-planner";
export declare const inject: string[];
/** 插件配置。 */
export interface Config {
    /**
     * Python 解释器（绝对路径，或 PATH 上的名字）。**留空即自动探测**：
     * 依次尝试 `python3`、`python`、DSH 自带的运行时，取第一个 3.10+ 的。
     */
    python?: string;
    /**
     * `stib/` 所在目录。**留空即使用插件自带的包**（`<插件目录>/stib`）。
     * 仅当你想改用仓库里的另一份 Python 代码时才需要设置。
     */
    plannerDir?: string;
    /** 基因组缓存目录；缺省为 `<plannerDir>/data`。 */
    dataDir?: string;
    /** 单次调用的超时（毫秒）。同时用作 spec.graceMs。 */
    timeoutMs?: number;
    /** 收集输出的内存上限（字节）。超出保留**尾部**。 */
    maxOutputBytes?: number;
}
export declare const Config: z<Config>;
export declare function apply(ctx: Context, config: Config): void;
