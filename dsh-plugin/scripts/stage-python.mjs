/**
 * prepublish 脚本：把 Python 侧的 stib/ 包复制进插件目录，
 * 使 npm 包**自包含** —— 别人 `npm install` 或 `dsh plugin add` 后
 * 无需再单独获取 Python 代码即可使用 9 个工具。
 *
 * ⚠️ 为什么需要这一步：
 *   npm 的 `files` **无法包含父目录**（实测 `"../stib"` 被静默忽略）。
 *   插件在 ../dsh-plugin，Python 包在 ../stib，所以必须在打包前复制进来。
 *
 * ⚠️ 复制进来的是**分发副本**，不是源码权威：
 *   主副本始终是 ../../stib（仓库根下的 stib/），开发只改那一份。
 *   本脚本每次打包都会覆盖 stib/，避免两份漂移。
 */

import { cpSync, rmSync, existsSync, readdirSync, statSync } from 'node:fs'
import { dirname, join } from 'node:path'
import { fileURLToPath } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
// ⚠️ 本脚本在 dsh-plugin/scripts/ 下，所以仓库根的 stib/ 要上**两**级：
//      <repo>/dsh-plugin/scripts/stage-python.mjs  →  <repo>/stib
//    分发副本放在 dsh-plugin/stib/ —— **目录名必须是 stib**，
//    因为运行方式是 `python -m stib.cli`，包名即目录名。
const src = join(here, '..', '..', 'stib')
const dst = join(here, '..', 'stib')

if (!existsSync(src)) {
  console.error(`[stage-python] 找不到 Python 包: ${src}`)
  console.error('  期望仓库结构：<repo>/stib/ 与 <repo>/dsh-plugin/ 同级。')
  process.exit(1)
}

rmSync(dst, { recursive: true, force: true })
cpSync(src, dst, { recursive: true })

// ⚠️ 不依赖 cpSync 的 filter —— 实测它会漏掉 __pycache__（**pyc 被打进了 npm 包**）。
//    改为复制后**显式清理**，语义确定。
function pruneCache(dir) {
  let removed = 0
  for (const e of readdirSync(dir)) {
    const p = join(dir, e)
    const st = statSync(p)
    if (st.isDirectory()) {
      if (e === '__pycache__') {
        rmSync(p, { recursive: true, force: true })
        removed += 1
      } else {
        removed += pruneCache(p)
      }
    } else if (e.endsWith('.pyc')) {
      rmSync(p, { force: true })
      removed += 1
    }
  }
  return removed
}
const pruned = pruneCache(dst)

// 统计并报告，便于在 CI 日志里确认真的打包了
function walk(dir) {
  let n = 0
  let bytes = 0
  for (const e of readdirSync(dir)) {
    const p = join(dir, e)
    const st = statSync(p)
    if (st.isDirectory()) {
      const sub = walk(p)
      n += sub.n
      bytes += sub.bytes
    } else if (e.endsWith('.py')) {
      n += 1
      bytes += st.size
    }
  }
  return { n, bytes }
}

const { n, bytes } = walk(dst)
console.log(`[stage-python] 已复制 ${n} 个 .py（${(bytes / 1024).toFixed(1)} KB），清理缓存 ${pruned} 项`)
console.log('[stage-python] 主副本仍是 ../stib —— 开发只改那一份，本副本每次打包重建')
