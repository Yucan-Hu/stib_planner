/**
 * 端到端验证「零配置」：不提供 python / plannerDir，直接调用工具。
 *
 * 与 smoke.mjs 的区别：这里 `ctx.subprocess.spawn` 用**真实** child_process
 * 实现（按 dsh-subprocess 的 collect 契约），因此能验证：
 *   1. Python 自动探测是否真的找到 3.10+ 解释器
 *   2. 插件是否用 `import.meta.url` 定位到自带的 `stib/`
 *   3. CLI 是否真的跑起来并返回 ALL PASS
 *
 * 用法：node smoke-live.mjs
 */
import { spawn } from 'node:child_process'
import { existsSync } from 'node:fs'
import { dirname, join, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
const PKG = resolve(here)

/** 按 dsh-subprocess 的 collect 契约实现 spawn。 */
function realSpawn(spec) {
  const [cmd, ...args] = spec.argv
  const maxOut = spec.stdio?.stdout?.maxBytes ?? 1 << 20
  const maxErr = spec.stdio?.stderr?.maxBytes ?? 1 << 20

  const child = spawn(cmd, args, {
    cwd: spec.cwd,
    env: spec.env ? { ...process.env, ...spec.env } : process.env,
    stdio: ['ignore', 'pipe', 'pipe'],
    // Windows 上不能直接 spawn .cmd/.bat；用 shell 解析 PATH 上的名字，
    // 与 DSH 自己的 provider 行为接近（绝对路径时 shell 也无害）。
    shell: false,
  })

  const bufs = { out: [], err: [] }
  let outLen = 0
  let errLen = 0
  child.stdout.on('data', (d) => {
    if (outLen < maxOut) { bufs.out.push(d); outLen += d.length }
  })
  child.stderr.on('data', (d) => {
    if (errLen < maxErr) { bufs.err.push(d); errLen += d.length }
  })

  const done = new Promise((res, rej) => {
    child.on('error', rej)
    child.on('close', (code, signal) => res({ exitCode: code, signal }))
  })

  const reader = (which) => ({
    readFrom() {
      const b = Buffer.concat(bufs[which])
      return { text: b.toString('utf8'), nextOffset: b.length, lossy: false }
    },
  })

  return {
    stdin: undefined,
    stdout: undefined,
    stderr: undefined,
    control: undefined,
    collected: { stdout: reader('out'), stderr: reader('err') },
    done,
    terminate() { try { child.kill() } catch {} },
    async waitForExit() { await done; return true },
  }
}

const registered = []
const ctx = {
  tools: { register(d) { registered.push(d); return () => {} } },
  subprocess: { spawn: realSpawn },
}

const mod = await import('./lib/index.js')

// ★ 关键：config 里**不写** python / plannerDir / dataDir
const config = { timeoutMs: 300000, maxOutputBytes: 2 << 20 }
mod.apply(ctx, config)

console.log('=== 1. 工具注册 ===')
console.log('  注册数:', registered.length, '（期望 9）')

const self = registered.find((t) => t.name === 'stib_selftest')
if (!self) { console.log('  ✗ 找不到 stib_selftest'); process.exit(1) }

console.log()
console.log('=== 2. 无配置调用 stib_selftest（真实子进程）===')
let result
try {
  result = await self.execute({}, { signal: undefined })
} catch (e) {
  console.log('  ✗ 抛错:', e.message)
  process.exit(1)
}

const blocks = self.output.render({}, result)
const text = blocks.map((b) => b.text ?? '').join('')

console.log('  返回字节:', text.length)
const pass = /ALL PASS/.test(text)
console.log(`  ${pass ? 'PASS' : 'FAIL'}  CLI 在零配置下跑通`)

if (!pass) {
  console.log('--- 输出片段 ---')
  console.log(text.split('\n').slice(0, 25).join('\n'))
}

console.log()
console.log('=== 3. 连带验证一个真正依赖 Python 的工具 ===')
const genes = registered.find((t) => t.name === 'stib_find_gene')
if (genes && existsSync(join(PKG, '..', 'testdata', 'mini.gb'))) {
  const mini = join(PKG, '..', 'testdata', 'mini.gb')
  let r2
  try {
    r2 = await genes.execute({ accession: mini, search: 'thr' }, { signal: undefined })
    const t2 = genes.output.render({}, r2).map((b) => b.text ?? '').join('')
    const ok2 = /thrA/.test(t2)
    console.log(`  ${ok2 ? 'PASS' : 'FAIL'}  stib_find_gene 按基因名搜到 thrA`)
    if (!ok2) console.log(t2.slice(0, 600))
  } catch (e) {
    console.log('  ✗ 抛错:', e.message)
  }
} else {
  console.log('  [SKIP] 缺少 testdata/mini.gb')
}

console.log()
console.log(pass ? '===== 零配置端到端通过 =====' : '===== 失败 =====')
process.exit(pass ? 0 : 1)
