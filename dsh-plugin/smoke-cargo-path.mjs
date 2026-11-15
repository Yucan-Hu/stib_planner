/** 验证：插件的 cargo 工具会把相对路径解析成绝对路径。 */
import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
const PKG = resolve(here)

const cap = []
const reg = []
const mod = await import('./lib/index.js')
mod.apply(
  {
    tools: { register(d) { reg.push(d); return () => {} } },
    subprocess: {
      spawn(spec) {
        cap.push(spec)
        return {
          collected: {
            stdout: { readFrom: () => ({ text: '(mock)', nextOffset: 0, lossy: false }) },
            stderr: { readFrom: () => ({ text: '', nextOffset: 0, lossy: false }) },
          },
          done: Promise.resolve({ exitCode: 0, signal: null }),
          terminate() {}, waitForExit: async () => true,
        }
      },
    },
  },
  { python: 'python', plannerDir: PKG, timeoutMs: 30000, maxOutputBytes: 65536 },
)

const tool = reg.find((t) => t.name === 'stib_cargo_assembly')
console.log('=== 用相对路径调用 stib_cargo_assembly ===')
await tool.execute(
  { cargoFasta: 'testdata/cargo600.fa', backbone: 'testdata/psh4300C.dna',
    site: 'after_selection', cargoStrand: -1, outDir: 'out/rel_test' },
  { signal: undefined },
)
const argv = cap[0].argv
const iC = argv.indexOf('--cargo')
const iB = argv.indexOf('--backbone')
const iO = argv.indexOf('--out')
console.log('  --cargo    =', argv[iC + 1])
console.log('  --backbone =', argv[iB + 1])
console.log('  --out      =', argv[iO + 1])

const abs = (p) => p.startsWith('/') || /^[A-Za-z]:[\\/]/.test(p)
const ok = abs(argv[iC + 1]) && abs(argv[iB + 1]) && abs(argv[iO + 1])
console.log()
console.log(ok ? '  PASS  三个路径都已解析成绝对路径' : '  FAIL  仍有相对路径')
process.exit(ok ? 0 : 1)
