/**
 * 插件冒烟测试：用 mock Context 调用 apply()，检查工具注册是否成功。
 *
 * 这比 dump-config 更接近真实：
 *   - 真的执行了插件的 apply()
 *   - 真的走了 ctx.tools.register() 与 ctx.subprocess.spawn() 的调用形态
 *   - defineTool 会对参数 DSL 与输出 schema 做**真实校验**（格式错会抛错）
 *
 * 不做的：不启动完整 Cordis 运行时、不真正 spawn Python（那是下一步）。
 */

const registered = []
const spawned = []

function fakeCtx() {
  return {
    tools: {
      register(def) {
        registered.push(def)
        return () => {}
      },
    },
    subprocess: {
      spawn(spec) {
        spawned.push(spec)
        // 返回一个符合契约的最小 handle
        return {
          stdin: undefined,
          stdout: undefined,
          stderr: undefined,
          control: undefined,
          collected: {
            stdout: { readFrom: () => ({ text: '(mock stdout)', nextOffset: 0, lossy: false }) },
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
console.log('=== 1. 导出契约 ===')
console.log('导出:', Object.keys(mod).sort().join(', '))
console.log('无 default 导出:', !('default' in mod))
console.log('name =', mod.name)
console.log('inject =', JSON.stringify(mod.inject))

console.log()
console.log('=== 2. 调用 apply()（真实注册全部工具）===')
const config = {
  python: 'python',
  plannerDir: 'E:/tmp/stib_planner',
  dataDir: 'E:/tmp/stib_planner/data',
  timeoutMs: 30000,
  maxOutputBytes: 65536,
}
try {
  mod.apply(fakeCtx(), config)
  console.log('apply() 未抛错 ✅')
} catch (e) {
  console.log('apply() 抛错 ❌', e && e.message)
  process.exit(1)
}

console.log()
console.log('=== 3. 注册结果 ===')
console.log('注册工具数:', registered.length, '（期望 9）')
for (const d of registered) {
  const params = d.parameters ? Object.keys(d.parameters) : []
  console.log(`  - ${d.name}  (参数: ${params.length ? params.join(', ') : '无'})`)
}

console.log()
console.log('=== 4. 工具名与预期比对 ===')
const want = ['stib_plan', 'stib_find_sgrnas', 'stib_offtarget', 'stib_batch',
  'stib_open_questions', 'stib_selftest',
  'stib_find_gene', 'stib_neutral_sites', 'stib_cargo_assembly']
const got = registered.map(d => d.name)
const missing = want.filter(w => !got.includes(w))
const extra = got.filter(g => !want.includes(g))
console.log('期望:', want.join(', '))
console.log('实际:', got.join(', '))
console.log('缺失:', missing.length ? missing.join(', ') : '无')
console.log('多余:', extra.length ? extra.join(', ') : '无')

console.log()
console.log('=== 5. 执行一个工具，验证 spawn 调用形态 ===')
const selftest = registered.find(d => d.name === 'stib_selftest')
try {
  const res = await selftest.execute({}, { signal: undefined })
  console.log('execute() 返回:', JSON.stringify(res))
} catch (e) {
  console.log('execute() 抛错 ❌', e && e.message)
  process.exit(1)
}

console.log()
console.log('=== 6. 实际传给 ctx.subprocess.spawn 的 spec（契约核对）===')
console.log(JSON.stringify(spawned[0], null, 2))

console.log()
console.log('=== 7. 契约符合性断言 ===')
const s = spawned[0]
const checks = [
  ['argv 是数组且 argv[0] 是程序', Array.isArray(s.argv) && typeof s.argv[0] === 'string'],
  ['无 command/args 字段（契约只要 argv）', !('command' in s) && !('args' in s)],
  ['cwd 存在', typeof s.cwd === 'string' && s.cwd.length > 0],
  ['graceMs 存在且为正', typeof s.graceMs === 'number' && s.graceMs > 0],
  ['stdio 三流均为显式对象/字面量', !!s.stdio && 'stdin' in s.stdio && 'stdout' in s.stdio && 'stderr' in s.stdio],
  ['stdout 为 collect 对象（含 maxBytes）', typeof s.stdio.stdout === 'object' && typeof s.stdio.stdout.maxBytes === 'number'],
  ['stdin 为 ignore', s.stdio.stdin === 'ignore'],
  ['env 为对象', s.env === undefined || typeof s.env === 'object'],
]
let ok = true
for (const [name, pass] of checks) {
  console.log(`  ${pass ? 'PASS' : 'FAIL'}  ${name}`)
  if (!pass) ok = false
}

console.log()
console.log(ok && missing.length === 0 ? '===== 全部通过 =====' : '===== 存在问题 =====')
process.exit(ok && missing.length === 0 ? 0 : 1)
