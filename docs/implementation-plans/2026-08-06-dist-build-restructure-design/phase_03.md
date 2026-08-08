# Dist Build Restructure Implementation Plan

## Phase 3: The build pipeline

**Goal:** Replace `scripts/rewriteBuilds.ts`, `scripts/buildUMD.ts`, and
both `tsc` emit passes with a single `scripts/build.ts` that produces
`dist/browser/`, `dist/node/`, and `dist/umd/`.

**Architecture:** esbuild bundles, terser minifies. An esbuild resolve
plugin performs the `./node-worker` to `./worker` swap for the browser
target, which is the reason the pipeline must use the JS API rather than
CLI flags. All minified output goes through one terser configuration, the
same one that has been proven against `wcln` in the existing UMD build.

**Tech Stack:** esbuild 0.25.12 (JS API, pinned `^0.25.0` in
`package.json`), terser 5.49.2, TypeScript 6.0.3, tsx. No task upgrades
esbuild; 0.28.1 is the current release but is not what is installed.

**Scope:** Phase 3 of 7.

**Codebase verified:** 2026-08-06

**Depends on:** Phase 2. The declaration pass reads `tsconfig.build.json`
with its new `scripts` exclusion **and its new `rootDir: "src"`**, and
the node bundle assumes `src/node-worker.ts` has a static import. Task 3
fails outright if that `rootDir` is missing.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC1: The build produces the specified dist/ layout

- **dist-build-restructure-design.AC1.1 Success:** `npm run build` writes
  `dist/browser/index.js`, `index.js.map`, `index.min.js`,
  `index.min.js.map`, and `index.d.ts`.
- **dist-build-restructure-design.AC1.2 Success:** `npm run build` writes
  the same five files under `dist/node/`.
- **dist-build-restructure-design.AC1.3 Success:** `npm run build` writes
  a minified `dist/umd/fflate.js`.
- **dist-build-restructure-design.AC1.4 Success:** The build creates no
  `lib/`, `esm/`, or top level `umd/` directory.
- **dist-build-restructure-design.AC1.5 Success:** The build removes a
  stale `dist/` before emitting, so no artifact from a previous run
  survives.

### dist-build-restructure-design.AC2: The worker swap is correct per build target

- **dist-build-restructure-design.AC2.1 Success:**
  `dist/browser/index.js` contains the `src/worker.ts` implementation,
  identifiable by `URL.createObjectURL` and `Blob`.
- **dist-build-restructure-design.AC2.3 Failure:**
  `dist/browser/index.js` contains no reference to `worker_threads` at
  module scope.

---

## Investigation findings

Verified against the working tree on 2026-08-06, plus two behaviours
confirmed by running esbuild and terser directly.

Confirmed as the design describes:

- `src/index.ts:15` is `import wk from './node-worker';`, the sole worker
  import.
- `src/index.ts:1037` carries the comment "The reason we can't just use
  the original variable names is minifiers mangling the toplevel scope."
- `wcln` is defined at `src/index.ts:1040-1062` and called at
  `src/index.ts:1088` and `:1089`.
- The thunks whose array literal text is parsed are `bInflt`
  (`src/index.ts:1096`), `bDflt` (`:1097`), `gze` (`:1100`), `guze`
  (`:1102`), `zle` (`:1104`), and `zule` (`:1106`).
- `scripts/buildUMD.ts:14-23` uses exactly
  `{ mangle: { toplevel: true }, compress: { passes: 5, unsafe: true,
  pure_getters: true }, sourceMap: false }`.
- `scripts/rewriteBuilds.ts:39` prepends
  `"import { createRequire } from 'module';\nvar require = createRequire('/');\n"`
  to `esm/index.mjs`.
- Nothing in the repository today exercises the async APIs against a
  minified bundle. The design is correct that this is new exposure.

Corrections and additions the design did not cover:

1. **`write: false` does return source maps as separate output files.**
   Verified by running the installed esbuild 0.25.12 with
   `{ write: false, sourcemap: true }`: `result.outputFiles` had length 2,
   containing `index.js.map` and `index.js`. The build script can
   therefore hold bundles in memory, feed them to terser, and write both
   the bundle and its map itself. No temporary files are needed.

2. **The resolve plugin must join against `args.resolveDir`, not a
   hard coded `src/` prefix.** `args.resolveDir` is already the importing
   file's directory, so the correct expression is
   `path.resolve(args.resolveDir, 'worker.ts')`. Verified end to end: the
   emitted browser bundle contained the `worker.ts` implementation and
   none of the `node-worker.ts` implementation.

3. **The terser configuration survives the `wcln` pattern on ESM
   output.** Verified against a reduced model of the `bInflt` pattern
   minified with the proven options plus `module: true`. The array
   literal emerged as `[n,t,r,e,s,c,i,l]`, every entry a bare mangled
   identifier, with mangling consistent between the definitions and the
   literal. This raises confidence but does not settle the question: the
   model is eight entries, whereas `bDflt` has thirty and the real file
   is 3,885 lines. Phase 6 is the actual gate.

4. **The UMD artifact narrows from dual-target to browser-only.**
   `scripts/buildUMD.ts:34-37` currently inlines both worker
   implementations and picks between them at runtime with a
   `typeof module` test, so today's UMD works under Node as well as in a
   browser. The design specifies UMD as "esbuild bundle of the browser
   variant", which drops the Node path. This is consistent with the
   design casting UMD as "a fallback for script-tag and legacy-bundler
   consumers", but it is a deliberate behaviour reduction and is recorded
   here so it is not mistaken for a regression.

5. **The UMD output path and filename both change.**
   `scripts/buildUMD.ts:39` writes `umd/index.js`. The design specifies
   `dist/umd/fflate.js`.

6. **The declaration emit produces more than `index.d.ts`.** With
   `declaration: true` and `declarationMap: true` from Phase 2,
   `tsc --emitDeclarationOnly -p tsconfig.build.json` writes
   `index.d.ts`, `index.d.ts.map`, `worker.d.ts`, `worker.d.ts.map`,
   `node-worker.d.ts`, and `node-worker.d.ts.map` into `dist/`. Task 3
   keeps only `index.d.ts` and removes the rest, since both published
   builds are bundled and neither worker module is separately importable.

7. **`keepNames` must be `false`, not `true`.** Corrected during
   execution on 2026-08-06; the code blocks in Tasks 4, 5 and 6 have been
   amended. esbuild's `keepNames` wraps every function as
   `__name(fn, 'name')`. The async APIs stringify functions and eval them
   inside a worker, where the `__name` helper is not in scope, so
   `deflate()` dies with `ReferenceError: __name is not defined` while
   every sync API keeps working. Confirmed by building the identical
   source both ways. Name recovery for the worker payload is `wcln`'s
   job, not esbuild's, and nothing in `src/` reads
   `Function.prototype.name`.

8. **Verify the async round trip from a file, not from
   `node --input-type=module -e`.** The `-e` form puts
   `--input-type=module` into `process.execArgv`, and worker threads
   inherit `execArgv`. The inherited flag makes the worker parse its
   eval'd payload as ESM, so the payload's sloppy-mode implicit globals
   (`u8=Uint8Array` and friends) throw
   `ReferenceError: u8 is not defined`. Isolated directly: the same
   worker eval succeeds with `execArgv: []` and fails with the inherited
   value, and a plain `.mjs` parent works because its `execArgv` is
   empty. This is an artifact of the verification command only. The
   shipped bundle is unaffected, and any harness that checks the async
   path must attach a `.catch` and assert the exit code, since an
   unhandled rejection is easy to misread as a pass.

---

<!-- START_SUBCOMPONENT_A (tasks 1-2) -->

<!-- START_TASK_1 -->
### Task 1: Create the build script skeleton and clean step

**Verifies:** dist-build-restructure-design.AC1.5

**Files:**
- Create: `scripts/build.ts`

**Implementation:**

Create `scripts/build.ts` with the module scaffolding and the clean step.
Later tasks in this phase append to it. `import.meta.dirname` is used
throughout because the package is `"type": "module"`.

```ts
import * as esbuild from 'esbuild'
import { minify, type MinifyOptions } from 'terser'
import { execFileSync } from 'child_process'
import { rmSync, mkdirSync, writeFileSync, readFileSync } from 'fs'
import { resolve, join } from 'path'

const root = resolve(import.meta.dirname, '..')
function p (...parts:Array<string>):string {
    return join(root, ...parts)
}

// Step 1: remove any output from a previous run.
rmSync(p('dist'), { recursive: true, force: true })
mkdirSync(p('dist/browser'), { recursive: true })
mkdirSync(p('dist/node'), { recursive: true })
mkdirSync(p('dist/umd'), { recursive: true })

console.log('cleaned dist/')
```

**Verification:**

Run: `npx tsx scripts/build.ts`
Expected: prints `cleaned dist/` and exits 0.

Run: `ls dist`
Expected: `browser`, `node`, `umd`.

**Commit:** `build: add build script skeleton with clean step`
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Add the shared terser configuration and a write helper

**Verifies:** None (infrastructure)

**Files:**
- Modify: `scripts/build.ts` (append)

**Implementation:**

Append the terser options and a helper that minifies a bundle and writes
both the minified file and its chained source map.

These options are copied from `scripts/buildUMD.ts:14-23`, which is the
configuration already proven safe against `wcln` in this codebase. The
one addition is `module: true` for the two ESM targets, which the design
requires. esbuild's own minifier is deliberately not used: mixing two
manglers over code whose correctness depends on consistent mangling is
the risk the design calls out.

```ts
// Proven against wcln by the existing UMD build. Do not add compress
// options without re-running the Phase 6 browser suite against the
// minified bundle: any option that substitutes a constant's VALUE for
// its NAME inside the bInflt/bDflt array literals breaks name recovery
// silently, leaving sync APIs working while async APIs corrupt output.
function terserOpts (esm:boolean):MinifyOptions {
    return {
        mangle: { toplevel: true },
        compress: { passes: 5, unsafe: true, pure_getters: true },
        module: esm
    }
}

type Bundle = { code:string; map:string }

// esbuild with write:false returns the bundle and its map as separate
// entries in outputFiles. Confirmed against esbuild 0.25.12.
function collect (out:esbuild.BuildResult):Bundle {
    const js = out.outputFiles!.find(f => !f.path.endsWith('.map'))!
    const map = out.outputFiles!.find(f => f.path.endsWith('.map'))!
    return { code: js.text, map: map.text }
}

async function writeBundle (
    dir:string,
    name:string,
    bundle:Bundle,
    esm:boolean
):Promise<void> {
    writeFileSync(p(dir, name + '.js'), bundle.code)
    writeFileSync(p(dir, name + '.js.map'), bundle.map)

    const min = await minify(bundle.code, {
        ...terserOpts(esm),
        sourceMap: {
            content: bundle.map,
            url: name + '.min.js.map'
        }
    })
    if (!min.code) throw new Error('terser produced no output for ' + name)

    writeFileSync(p(dir, name + '.min.js'), min.code)
    // terser types map as string|RawSourceMap, but it returns a string
    // whenever sourceMap.content is supplied, which it always is here.
    writeFileSync(p(dir, name + '.min.js.map'), min.map as string)

    console.log('wrote ' + dir + '/' + name + '{.js,.min.js} and maps')
}
```

**Verification:**

Run: `npx tsx scripts/build.ts`
Expected: still prints `cleaned dist/` and exits 0. Nothing calls
`writeBundle` yet, so no additional output.

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep "error TS" | grep -c "scripts/build.ts" || true`
Expected: `0`. The `grep "error TS"` filter is required: the root config
sets `listFiles: true`, so `tsc` prints every file it loads and an
unfiltered grep matches the filename itself.

**Commit:** `build: add shared terser config and write helper`
<!-- END_TASK_2 -->

<!-- END_SUBCOMPONENT_A -->

<!-- START_TASK_3 -->
### Task 3: Emit declarations

**Verifies:** dist-build-restructure-design.AC1.1,
dist-build-restructure-design.AC1.2

**Files:**
- Modify: `scripts/build.ts` (append)

**Implementation:**

One `tsc` pass emits declarations for `src/` only, because
`tsconfig.build.json` excludes `example`, `test`, and `scripts`. Phase 2
Task 2 sets `"rootDir": "src"` on that config, which is what makes the
declaration land at `dist/index.d.ts`. Without it TypeScript raises
`TS5011` and emits to `dist/src/index.d.ts`, and this task fails twice
over: the `tsc` call exits non-zero and the file it reads is absent.

The public API is identical between the browser and node builds, so the
single declaration is copied to both. `worker.d.ts` and
`node-worker.d.ts` are discarded because both builds are bundled and
neither worker module is separately importable.

```ts
// Step 2: declarations. Invoke the local tsc through node rather than
// npx: npx resolves to npx.cmd on Windows, which execFileSync will not
// find without a shell.
execFileSync(
    process.execPath,
    [
        p('node_modules/typescript/bin/tsc'),
        '--emitDeclarationOnly',
        '--project',
        'tsconfig.build.json'
    ],
    { cwd: root, stdio: 'inherit' }
)

const dts = readFileSync(p('dist/index.d.ts'), 'utf8')

// The bundled builds expose no separate worker entry point, so only
// index.d.ts is published. Declaration maps are dropped with it: they
// would point at src/ paths that are not shipped.
for (const f of [
    'index.d.ts', 'index.d.ts.map',
    'worker.d.ts', 'worker.d.ts.map',
    'node-worker.d.ts', 'node-worker.d.ts.map'
]) {
    rmSync(p('dist', f), { force: true })
}

writeFileSync(p('dist/browser/index.d.ts'), dts)
writeFileSync(p('dist/node/index.d.ts'), dts)

console.log('wrote declarations')
```

**Verification:**

Run: `npx tsx scripts/build.ts`
Expected: exits 0 and prints `wrote declarations`.

Run: `ls dist/browser/index.d.ts dist/node/index.d.ts`
Expected: both exist.

Run: `ls dist/*.d.ts 2>&1`
Expected: "No such file or directory". The intermediate declarations must
have been removed.

Run: `grep -nE "from '\./(node-)?worker'" dist/browser/index.d.ts`
Expected: no output. If this matches, `index.d.ts` is not self-contained
and references a worker module that is not published. Stop and report it
rather than working around it; it means the public API surface leaks an
internal type.

**Commit:** `build: emit declarations to both dist targets`
<!-- END_TASK_3 -->

<!-- START_SUBCOMPONENT_B (tasks 4-5) -->

<!-- START_TASK_4 -->
### Task 4: Bundle the browser target with the worker swap

**Verifies:** dist-build-restructure-design.AC1.1,
dist-build-restructure-design.AC2.1,
dist-build-restructure-design.AC2.3

**Files:**
- Modify: `scripts/build.ts` (append)

**Implementation:**

This is the reason the pipeline uses the esbuild JS API. `src/index.ts:15`
imports `./node-worker`, and the browser build must resolve that to
`src/worker.ts`. esbuild's `--alias:` flag rewrites bare package
specifiers only, never relative paths, so a resolve plugin is required.

`args.resolveDir` is already the importing file's directory, so the
target is `worker.ts` relative to it. Verified against esbuild 0.25.12:
the emitted bundle contained the `worker.ts` implementation and none of
`node-worker.ts`.

```ts
// Step 3: browser bundle.
// esbuild's --alias only rewrites bare specifiers, so the relative
// ./node-worker to ./worker swap needs a resolve plugin.
const workerSwap:esbuild.Plugin = {
    name: 'browser-worker-swap',
    setup (build) {
        build.onResolve({ filter: /^\.\/node-worker$/ }, args => ({
            path: resolve(args.resolveDir, 'worker.ts')
        }))
    }
}

const browser = collect(await esbuild.build({
    entryPoints: [p('src/index.ts')],
    outfile: p('dist/browser/index.js'),
    bundle: true,
    write: false,
    format: 'esm',
    platform: 'browser',
    target: 'es2022',
    keepNames: false,
    sourcemap: true,
    plugins: [workerSwap]
}))

await writeBundle('dist/browser', 'index', browser, true)
```

**Verification:**

Run: `npx tsx scripts/build.ts`
Expected: exits 0 and reports the browser bundle was written.

Run: `ls dist/browser`
Expected: `index.d.ts`, `index.js`, `index.js.map`, `index.min.js`,
`index.min.js.map`.

Run: `grep -c "createObjectURL" dist/browser/index.js`
Expected: at least `1`. This proves the browser worker was substituted.

Run: `grep -c "worker_threads" dist/browser/index.js || true`
Expected: `0`. Any match means the swap did not take effect and the
browser bundle carries the Node implementation.

Run: `node --input-type=module -e "import('./dist/browser/index.js').then(m => console.log(typeof m.deflateSync, typeof m.deflate))"`
Expected: prints `function function`.

**Commit:** `build: bundle browser target with worker swap`
<!-- END_TASK_4 -->

<!-- START_TASK_5 -->
### Task 5: Bundle the node target

**Verifies:** dist-build-restructure-design.AC1.2

**Files:**
- Modify: `scripts/build.ts` (append)

**Implementation:**

The node bundle keeps `src/node-worker.ts`, so no plugin is involved.
`node:worker_threads` is listed explicitly in `external` rather than
relying on `platform: 'node'` alone, so the intent is visible in the
config and does not depend on esbuild's builtin detection.

```ts
// Step 4: node bundle.
const node = collect(await esbuild.build({
    entryPoints: [p('src/index.ts')],
    outfile: p('dist/node/index.js'),
    bundle: true,
    write: false,
    format: 'esm',
    platform: 'node',
    target: 'es2022',
    keepNames: false,
    sourcemap: true,
    external: ['node:worker_threads']
}))

await writeBundle('dist/node', 'index', node, true)
```

**Verification:**

Run: `npx tsx scripts/build.ts`
Expected: exits 0 and reports the node bundle was written.

Run: `ls dist/node`
Expected: `index.d.ts`, `index.js`, `index.js.map`, `index.min.js`,
`index.min.js.map`.

Run: `grep -c "node:worker_threads" dist/node/index.js`
Expected: at least `1`, and it must appear as an `import` rather than an
inlined module body.

Run: `grep -c "createObjectURL" dist/node/index.js || true`
Expected: `0`. The node bundle must not carry the browser worker.

Run: `node --input-type=module -e "import('./dist/node/index.js').then(m => { const d = m.deflateSync(new Uint8Array([1,2,3])); console.log(m.inflateSync(d).length) })"`
Expected: prints `3`.

Run: `node --input-type=module -e "import('./dist/node/index.js').then(m => new Promise((res, rej) => m.deflate(new Uint8Array([1,2,3]), (e, d) => e ? rej(e) : m.inflate(d, (e2, r) => e2 ? rej(e2) : res(r.length))))).then(console.log)"`
Expected: prints `3`. This exercises the async path, which is what
actually spawns a worker.

**Commit:** `build: bundle node target`
<!-- END_TASK_5 -->

<!-- END_SUBCOMPONENT_B -->

<!-- START_TASK_6 -->
### Task 6: Build the UMD artifact

**Verifies:** dist-build-restructure-design.AC1.3

**Files:**
- Modify: `scripts/build.ts` (append)

**Implementation:**

The UMD build bundles the browser variant, so it reuses the same resolve
plugin. esbuild emits `var fflate = (() => { ... })()` for
`format: 'iife'` with `globalName: 'fflate'`; the wrapper below declares
that inside a factory function and returns it.

Note that this narrows UMD to browser semantics. The previous
`scripts/buildUMD.ts` inlined both worker implementations and chose
between them at runtime. The design specifies the browser variant, which
matches UMD's role as a script-tag and legacy-bundler fallback.

UMD is CommonJS-shaped, so `module: false` is passed to terser here,
unlike the two ESM targets.

```ts
// Step 5: UMD. Bundles the browser variant, so the swap plugin applies.
const umd = await esbuild.build({
    entryPoints: [p('src/index.ts')],
    outfile: p('dist/umd/fflate.js'),
    bundle: true,
    write: false,
    format: 'iife',
    globalName: 'fflate',
    platform: 'browser',
    target: 'es2022',
    keepNames: false,
    sourcemap: false,
    plugins: [workerSwap]
})

const iife = umd.outputFiles!.find(f => !f.path.endsWith('.map'))!.text

const wrapped = [
    '(function(root, factory){',
    "  if (typeof module === 'object' && typeof exports === 'object')",
    '    module.exports = factory();',
    "  else if (typeof define === 'function' && define.amd)",
    '    define([], factory);',
    '  else root.fflate = factory();',
    "})(typeof self !== 'undefined' ? self : this, function(){",
    iife,
    '  return fflate;',
    '});'
].join('\n')

// Step 6: minify. UMD is CommonJS shaped, so module:false here.
const umdMin = await minify(wrapped, terserOpts(false))
if (!umdMin.code) throw new Error('terser produced no output for umd')

writeFileSync(p('dist/umd/fflate.js'), umdMin.code)

// The root package.json is "type": "module", so node would parse this
// .js as ESM, where the wrapper's `this` is undefined and loading throws
// "Cannot set properties of undefined". Scoping the directory back to
// commonjs makes the artifact loadable by both require() and import()
// while keeping the filename that unpkg points at.
writeFileSync(p('dist/umd/package.json'), '{ "type": "commonjs" }\n')

console.log('wrote dist/umd/fflate.js')
```

**Verification:**

Run: `npx tsx scripts/build.ts`
Expected: exits 0 and prints `wrote dist/umd/fflate.js`.

Run: `ls dist/umd`
Expected: `fflate.js` and `package.json`. The latter scopes the directory
to commonjs so the artifact is loadable; see the comment in the code
block above.

The sibling `dist/umd/package.json` written above is what makes this
artifact loadable. Without it the root `"type": "module"` makes Node
parse the `.js` as ESM, where the wrapper's `this` is `undefined` and
loading throws `TypeError: Cannot set properties of undefined (setting
'fflate')`. `require()` fails too, and there is no Node version or flag
on which it quietly returns an empty namespace. This passage previously
claimed there was; it was corrected after measurement, by moving the
file aside on node v25.8.2:

- `node t.cjs` -- `TypeError: Cannot set properties of undefined
  (setting 'fflate')`, identical to the `import()` path. Under
  `require(esm)`, default since 22.12, the file is parsed as ESM, so
  the wrapper's `this` is undefined and it throws on assignment.
- `node --no-experimental-require-module t.cjs` --
  `Error [ERR_REQUIRE_ESM]`. This is the pre-22.12 behaviour class:
  a different error, raised earlier, before the wrapper is evaluated
  at all, so `this` never enters the picture.

Both generations throw, and the remediation is the same either way:
keep the sibling `package.json`.
Scoping the directory to commonjs fixes both, and keeps the filename the
design names and `unpkg` points at, so do not rename to `.cjs`.

Measured with the sibling present: `require('.../dist/umd/fflate.js')`
and `import('.../dist/umd/fflate.js')` both yield the fflate object and
round trip correctly. The artifact remains reachable from a script tag,
from AMD loaders and from legacy bundlers, which is the role the design
assigns it.

Verify the global branch, which is the branch script-tag consumers take:

```bash
node -e "
const s = require('fs').readFileSync('dist/umd/fflate.js', 'utf8');
const g = {};
// module, exports and define must be shadowed as undefined. Under
// 'node -e' they are real objects on globalThis, so without these
// parameters the UMD wrapper takes its CommonJS branch and never
// assigns the global, leaving g.fflate undefined.
new Function('self', 'module', 'exports', 'define', s).call(g, g);
const f = g.fflate;
const d = f.deflateSync(new Uint8Array([1, 2, 3]));
console.log('umd roundtrip:', f.inflateSync(d).length);
"
```
Expected: prints `umd roundtrip: 3`. Note that `this` must not appear in
the parameter list: `new Function('self', 'this', s)` is a `SyntaxError`,
because `this` is not a valid parameter name.

Run: `grep -c "worker_threads" dist/umd/fflate.js || true`
Expected: `0`.

**Commit:** `build: produce umd artifact`
<!-- END_TASK_6 -->

<!-- START_TASK_7 -->
### Task 7: Delete the superseded build scripts

**Verifies:** dist-build-restructure-design.AC1.4

**Files:**
- Delete: `scripts/rewriteBuilds.ts`
- Delete: `scripts/buildUMD.ts`

**Implementation:**

Both scripts are fully superseded. `scripts/rewriteBuilds.ts` existed to
perform the worker swap with regular expressions over emitted JavaScript,
which the esbuild resolve plugin now does at resolution time. Its
`createRequire` banner is unnecessary because Phase 2 removed the only
module-scope `require`. `scripts/buildUMD.ts` is replaced by Task 6.

```bash
git rm scripts/rewriteBuilds.ts scripts/buildUMD.ts
```

Leave `package.json` alone. Its `build:rewrite` and `build:umd` scripts
still reference these files and are removed in Phase 4, which rewrites
the whole script block.

**Verification:**

Run: `ls scripts`
Expected: `build.ts` and `cpGHPages.ts` only.

Run: `grep -rn "rewriteBuilds\|buildUMD" --include='*.ts' --include='*.js' . --exclude-dir=node_modules --exclude-dir=docs --exclude-dir=dist`
Expected: no output.

**Commit:** `build: remove superseded build scripts`
<!-- END_TASK_7 -->

<!-- START_TASK_8 -->
### Task 8: Verify the complete build output

**Verifies:** dist-build-restructure-design.AC1.1,
dist-build-restructure-design.AC1.2,
dist-build-restructure-design.AC1.3,
dist-build-restructure-design.AC1.4,
dist-build-restructure-design.AC1.5

**Files:**
- No files changed. This task is a gate.

**Implementation:**

Run the full build from a clean tree and confirm the complete layout.

**Verification:**

Run: `rm -rf dist && npx tsx scripts/build.ts && find dist -type f | sort`
Expected exactly these twelve paths:

```
dist/browser/index.d.ts
dist/browser/index.js
dist/browser/index.js.map
dist/browser/index.min.js
dist/browser/index.min.js.map
dist/node/index.d.ts
dist/node/index.js
dist/node/index.js.map
dist/node/index.min.js
dist/node/index.min.js.map
dist/umd/fflate.js
dist/umd/package.json
```

Run: `ls lib esm umd 2>&1`
Expected: three "No such file or directory" errors. Stale directories
from the old pipeline are removed in Phase 7; if they still exist here,
confirm the build did not recreate them by checking their mtimes are
older than this run.

Run: `node -e "
const {statSync} = require('fs');
for (const f of ['dist/browser/index.js','dist/browser/index.min.js']) {
  console.log(f, statSync(f).size);
}"`
Expected: the `.min.js` size is smaller than the `.js` size. If they are
equal, terser did not run.

Run: `head -c 200 dist/browser/index.min.js`
Expected: minified output on a single line. Confirm it is not identical
to the head of `dist/browser/index.js`.

Run: `tail -1 dist/browser/index.min.js | grep -c "sourceMappingURL" || true`
Expected: `0` or `1`. Terser does not append the comment automatically;
either is acceptable because the map file is written alongside.

**Commit:** `build: verify complete dist output`
<!-- END_TASK_8 -->

---

## Phase 3 completion criteria

- `npx tsx scripts/build.ts` exits 0 from a clean tree and produces
  exactly the twelve files listed in Task 8.
- `dist/browser/index.js` contains `createObjectURL` and no
  `worker_threads`.
- `dist/node/index.js` imports `node:worker_threads` and contains no
  `createObjectURL`.
- `dist/umd/fflate.js` loads via the global branch, and the sibling
  `dist/umd/package.json` makes it loadable by `require()` and `import()`
  as well; see Task 6.
- Round-trip compression works through the node bundle, synchronously and
  asynchronously.
- `scripts/rewriteBuilds.ts` and `scripts/buildUMD.ts` are deleted.

The `wcln` hazard is **not** settled by this phase. The minified bundles
are produced here but nothing yet exercises their async paths. Phase 6 is
the gate for that, and it is the reason the browser suite is part of this
work rather than a follow-up.

`npm run build` still points at the old script chain until Phase 4.
