/**
 * 端到端验证「npm 包自包含」：
 *   1. 用打包后的插件目录作为 plannerDir（模拟用户安装后的目录结构）
 *   2. 用 mock ctx 捕获 argv，再用**真实 Python** 执行该 argv
 *   3. 确认工具在只有包内文件的情况下能跑通
 *
 * ⚠️ plannerDir 不硬编码本机路径 —— 从脚本自身位置推导，
 *    这样任何机器上 clone 后都能直接跑。
 */

import { spawnSync } from 'node:child_process'
import { existsSync } from 'node:fs'
import { dirname, join, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
const PKG = resolve(here)                       // dsh-plugin/
const DIST = join(PKG, 'stib')                  // 包内自带的 Python 包（目录名即包名）

console.log('=== 1. 检查包内是否自带 Python 代码 ===')
const need = ['cli.py', 'planner.py', 'genome.py', join('data', 'construct_psh43.yml')]
let ok = true
for (const f of need) {
  const p = join(DIST, f)
  const e = existsSync(p)
  if (!e) ok = false
  console.log(`  ${e ? 'OK ' : 'MISS'}  stib/${f}`)
}
if (!ok) {
  console.log()
  console.log('  提示：先跑 `npm run stage:python` 把 ../stib 复制进包内再测。')
}

console.log()
console.log('=== 2. mock ctx 捕获 argv（plannerDir 指向插件目录本身）===')
const cap = []
const reg = []
const mod = await import('./lib/index.js')
mod.apply({
  tools: { register(d) { reg.push(d); return () => {} } },
  subprocess: {
    spawn(spec) {
      cap.push(spec)
      return {
        collected: {
          stdout: { readFrom: () => ({ text: '', nextOffset: 0, lossy: false }) },
          stderr: { readFrom: () => ({ text: '', nextOffset: 0, lossy: false }) },
        },
        done: Promise.resolve({ exitCode: 0, signal: null }),
        terminate() {}, waitForExit: async () => true,
      }
    },
  },
}, {
  python: process.env.PYTHON || 'python',
  plannerDir: PKG,
  timeoutMs: 30000, maxOutputBytes: 65536,
})

const self = reg.find((t) => t.name === 'stib_selftest')
await self.execute({}, { signal: undefined })
const spec = cap[0]

console.log('  argv :', JSON.stringify(spec.argv))
console.log('  cwd  :', spec.cwd)
const cwdOk = spec.cwd === PKG
console.log(`  ${cwdOk ? 'PASS' : 'FAIL'}  cwd 即 plannerDir（其下有 stib/，故 -m stib.cli 可用）`)

console.log()
console.log('=== 3. 用真实 Python 执行该 argv（只有包内文件可用）===')
if (!process.env.PYTHON) {
  console.log('  [SKIP] 未设 PYTHON 环境变量；请用 PYTHON=<解释器路径> 运行本脚本')
} else {
  const r = spawnSync(spec.argv[0], spec.argv.slice(1), {
    cwd: spec.cwd, encoding: 'utf8', timeout: 300000,
    env: { ...process.env, PYTHONIOENCODING: 'utf-8', PYTHONUTF8: '1' },
  })
  const out = (r.stdout || '') + (r.stderr || '')
  const pass = /RESULT: ALL PASS/.test(out)
  console.log(`  exit=${r.status}  ${pass ? 'PASS —— 打包副本可独立运行' : 'FAIL'}`)
  if (!pass) console.log(out.split('\n').slice(-15).join('\n'))
  ok = ok && pass && cwdOk
}

console.log()
console.log(ok ? '===== 自包含验证通过 =====' : '===== 仍有问题 =====')
process.exit(ok ? 0 : 1)
