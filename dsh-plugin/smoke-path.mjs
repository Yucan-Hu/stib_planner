/**
 * 验证 resolveGenome 的路径解析：用 mock ctx 捕获传给 spawn 的 argv，
 * 检查基因组路径是否被正确解析（这是 agent 报的 bug 的直接复核）。
 */

import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const captured = []
const registered = []

function fakeCtx() {
  return {
    tools: { register(def) { registered.push(def); return () => {} } },
    subprocess: {
      spawn(spec) {
        captured.push(spec)
        return {
          collected: {
            stdout: { readFrom: () => ({ text: '(mock)', nextOffset: 0, lossy: false }) },
            stderr: { readFrom: () => ({ text: '', nextOffset: 0, lossy: false }) },
          },
          done: Promise.resolve({ exitCode: 0, signal: null }),
          terminate() {},
          waitForExit: async () => true,
        }
      },
    },
  }
}

const mod = await import('./lib/index.js')
const ctx = fakeCtx()
// 从脚本自身位置推导仓库根，不硬编码本机路径
const PLANNER = resolve(dirname(fileURLToPath(import.meta.url)), '..')
mod.apply(ctx, {
  python: process.env.PYTHON || 'python',
  plannerDir: PLANNER,
  dataDir: `${PLANNER}/data`,
  timeoutMs: 30000, maxOutputBytes: 65536,
})

const tools = registered
const findGene = tools.find((t) => t.name === 'stib_find_gene')
const plan = tools.find((t) => t.name === 'stib_plan')
const neutral = tools.find((t) => t.name === 'stib_neutral_sites')

console.log('=== 工具数 ===')
console.log('  注册:', tools.length)

/** 跑一次并返回 argv 里的 --genome 值。 */
async function genomeArg(tool, args) {
  captured.length = 0
  try {
    await tool.execute(args, { signal: undefined })
  } catch (e) {
    return { error: e.message }
  }
  const argv = captured[0]?.argv || []
  const i = argv.indexOf('--genome')
  return { genome: i >= 0 ? argv[i + 1] : null }
}

console.log()
console.log('=== resolveGenome 各种输入（这是 bug 的直接复核）===')
const cases = [
  ['testdata/mini.gb', 'agent 报的那条：相对路径 + .gb'],
  ['testdata/psh4300C.dna', '相对路径 + .dna'],
  ['testdata/demo_cargo.fa', '相对路径 + .fa'],
  ['E:/abs/path/genome.gb', '绝对路径（正斜杠）'],
  ['E:\\\\abs\\\\path\\\\genome.gbk', '绝对路径（反斜杠）'],
  ['NC_004663.1', '裸 accession（data/ 下有 .fasta）'],
  ['mini.gb', '裸文件名（带扩展名）'],
]
for (const [acc, desc] of cases) {
  const r = await genomeArg(findGene, { accession: acc, search: 'x' })
  console.log(`  ${desc}`)
  console.log(`    入参: ${acc}`)
  console.log(`    解析: ${r.genome || '(报错) ' + (r.error || '').split('\\n')[0]}`)
}

console.log()
console.log('=== 不存在的文件应给出可操作报错 ===')
const bad = await genomeArg(findGene, { accession: 'NC_999999.9', search: 'x' })
console.log('  ' + (bad.error ? bad.error.split('\\n').slice(0, 3).join('\\n    ') : `未报错，解析为 ${bad.genome}`))

console.log()
console.log('=== 抽查：不是裸 accession 时不应再拼 .fasta ===')
const r1 = await genomeArg(findGene, { accession: 'testdata/mini.gb', search: 'tdk' })
const ok1 = r1.genome && !r1.genome.endsWith('.fasta') && r1.genome.endsWith('mini.gb')
console.log(`  testdata/mini.gb → ${r1.genome}`)
console.log(`  ${ok1 ? 'PASS' : 'FAIL'}  不再被强加 .fasta`)

const r2 = await genomeArg(plan, { gene: 'thrA', accession: 'testdata/mini.gb' })
console.log(`  stib_plan + mini.gb → ${r2.genome}`)
const ok2 = r2.genome && r2.genome.endsWith('mini.gb')
console.log(`  ${ok2 ? 'PASS' : 'FAIL'}`)

const r3 = await genomeArg(neutral, { accession: 'testdata/mini.gb' })
console.log(`  stib_neutral_sites + mini.gb → ${r3.genome}`)
const ok3 = r3.genome && r3.genome.endsWith('mini.gb')
console.log(`  ${ok3 ? 'PASS' : 'FAIL'}`)

console.log()
const allOk = ok1 && ok2 && ok3
console.log(allOk ? '===== 路径解析全部通过 =====' : '===== 仍有问题 =====')
process.exit(allOk ? 0 : 1)
