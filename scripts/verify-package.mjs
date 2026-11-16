/**
 * 验证「从 git 安装」这份包是否真的能用。
 *
 * 背景：`dsh plugin add <git-url>` 只是转发给 `pnpm add`，装的是**仓库根**。
 * 所以仓库根必须自己是一个可用的组合包：
 *   - package.json 声明 dsh.bundle.patch
 *   - 带**编译好的** dsh-plugin/lib（pnpm 默认不跑依赖的构建脚本）
 *   - 带 Python 包（根目录的 stib/，不再重复放一份）
 *
 * 本脚本做四件事：
 *   1. 在仓库根 npm pack
 *   2. 解压到临时目录
 *   3. 检查关键文件是否都在包里
 *   4. **从解压后的布局**加载插件、用真实 Python 跑 stib_selftest
 *
 * 用法：PYTHON=<3.10+ 解释器> node scripts/verify-package.mjs
 */
import { spawn, spawnSync } from 'node:child_process'
import { existsSync, mkdtempSync, readFileSync, rmSync, symlinkSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { dirname, join, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
const REPO = resolve(here, '..')

const rows = []
const check = (ok, name, detail = '') => rows.push([ok, name, detail])

/** 按 DSH collect 契约跑一个子进程，供 mock ctx.subprocess 使用。 */
function realSpawn(spec) {
  const [cmd, ...args] = spec.argv
  const child = spawn(cmd, args, {
    cwd: spec.cwd,
    env: spec.env ? { ...process.env, ...spec.env } : process.env,
    stdio: ['ignore', 'pipe', 'pipe'],
  })
  const out = []
  const err = []
  child.stdout.on('data', (d) => out.push(d))
  child.stderr.on('data', (d) => err.push(d))
  const done = new Promise((res, rej) => {
    child.on('error', rej)
    child.on('close', (code, signal) => res({ exitCode: code, signal }))
  })
  const stream = (buf) => ({
    readFrom: () => ({ text: Buffer.concat(buf).toString('utf8'), nextOffset: 0, lossy: false }),
  })
  return {
    collected: { stdout: stream(out), stderr: stream(err) },
    done,
    terminate() {
      try {
        child.kill()
      } catch {
        /* 已退出 */
      }
    },
    async waitForExit() {
      await done
      return true
    },
  }
}

// ── 1. pack ──
// ⚠️ Windows 上 `npm` 是 `npm.cmd`，spawnSync 不加 shell 会直接 ENOENT
//    （而且 stderr 为 null，报错信息看起来像 "undefined"）。
const isWin = process.platform === 'win32'
const pack = spawnSync('npm', ['pack', '--silent'], {
  cwd: REPO,
  encoding: 'utf8',
  shell: isWin,
})
if (pack.status !== 0) {
  console.error('npm pack 失败:', pack.stderr || pack.stdout || pack.error)
  process.exit(1)
}
const name = (pack.stdout || '').trim().split(/\r?\n/).pop()
const tgz = join(REPO, name)
check(existsSync(tgz), 'npm pack 成功', name)

// ── 2. 解压 ──
const tmp = mkdtempSync(join(tmpdir(), 'dsh-stib-pkg-'))
const untar = spawnSync('tar', ['-xzf', tgz, '-C', tmp], { encoding: 'utf8' })
if (untar.status !== 0) {
  console.error('解压失败:', untar.stderr)
  rmSync(tmp, { recursive: true, force: true })
  process.exit(1)
}
const PKG = join(tmp, 'package')
check(existsSync(PKG), '解压得到 package/')

// ── 3. 关键文件与清单 ──
const need = [
  'package.json',
  'dsh-plugin/lib/index.js',
  'dsh-plugin/cordis.patch.yml',
  'dsh-plugin/skills/stib-bacteroides-editing/SKILL.md',
  'stib/cli.py',
  'stib/data/construct_psh43.yml',
  // 载体与示例：没有它们，装完连一个可跑的 cargo 组装都没有
  'testdata/psh4300C.dna',
  'testdata/mini.gb',
  'testdata/demo_cargo.fa',
  'USER_GUIDE.md',
  'VERSIONS.md',
]
for (const f of need) check(existsSync(join(PKG, f)), `包内含 ${f}`)

const manifest = JSON.parse(readFileSync(join(PKG, 'package.json'), 'utf8'))
check(manifest.name === 'dsh-stib-planner', '包名是 dsh-stib-planner', manifest.name)
const patchRel = manifest.dsh?.bundle?.patch
check(Boolean(patchRel), '声明了 dsh.bundle.patch', String(patchRel))
check(Boolean(patchRel) && existsSync(join(PKG, patchRel)), 'patch 指向的文件存在', String(patchRel))

// ── 4. 从解压布局真实执行 ──
// ⚠️ peerDependencies 由宿主（DSH profile）提供，临时目录里没有。
//    这里把开发目录的 node_modules 链进来，模拟"宿主已装好 peer"的真实情形。
//    若跳过这步，import 会在解析 @deepseek-ai/schemastery 时失败 ——
//    那不是包的问题，是验证环境缺少宿主依赖。
const devModules = join(REPO, 'dsh-plugin', 'node_modules')
if (existsSync(devModules)) {
  try {
    symlinkSync(devModules, join(PKG, 'node_modules'), isWin ? 'junction' : 'dir')
    check(true, '已链入 peerDependencies（模拟宿主环境）')
  } catch (e) {
    check(false, '链入 peerDependencies 失败', e.message)
  }
} else {
  check(false, '开发目录缺少 node_modules，无法模拟宿主 peer')
}

const reg = []
const cap = []
const mod = await import(
  `file://${join(PKG, 'dsh-plugin/lib/index.js').replace(/\\/g, '/')}`
)
mod.apply(
  {
    tools: {
      register(d) {
        reg.push(d)
        return () => {}
      },
    },
    subprocess: {
      spawn(spec) {
        cap.push(spec)
        return realSpawn(spec)
      },
    },
  },
  { python: process.env.PYTHON || 'python' },
)

check(reg.length === 9, '从包内布局注册了 9 个工具', String(reg.length))

const self = reg.find((t) => t.name === 'stib_selftest')
if (!self) {
  check(false, '找到 stib_selftest 工具')
} else {
  let text = ''
  try {
    const v = await self.execute({}, { signal: undefined })
    text = self.output.render({}, v).map((b) => b.text ?? '').join('')
  } catch (e) {
    text = `EXEC_ERROR: ${e.message}`
  }
  const ok = /ALL PASS/.test(text)
  check(ok, '包内布局下 stib_selftest 返回 ALL PASS', ok ? '' : text.slice(0, 300))
  check(cap[0]?.cwd === PKG, 'cwd 解析到包根（其下有 stib/）', String(cap[0]?.cwd))
}

// ── 汇总 ──
rmSync(tmp, { recursive: true, force: true })
try {
  rmSync(tgz, { force: true })
} catch {
  /* ignore */
}

const width = Math.max(...rows.map(([, n]) => n.length))
let fails = 0
console.log('\n=== 打包布局验证（模拟 pnpm add github:...） ===')
for (const [ok, n, d] of rows) {
  if (!ok) fails++
  console.log(`  [${ok ? 'PASS' : 'FAIL'}] ${n.padEnd(width)}  ${d}`)
}
console.log(`\n  ${rows.length - fails}/${rows.length} passed`)
console.log('RESULT:', fails === 0 ? 'ALL PASS' : `${fails} FAILED`)
process.exit(fails === 0 ? 0 : 1)
