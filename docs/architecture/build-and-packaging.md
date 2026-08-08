# Build and packaging

Every shipped artifact, how it is produced, and the invariants that fail
silently.

Sources:

- [`scripts/build.ts`](../../scripts/build.ts) -- the whole pipeline.
- [`package.json`](../../package.json) -- `exports`, `files`, scripts.
- [`tsconfig.build.json`](../../tsconfig.build.json) -- declaration emit.
- [`test/browser/harness.ts`](../../test/browser/harness.ts) -- guards the
  timeout invariant.
- [`specs/2026-08-06-dist-build-restructure-design.md`](../../specs/2026-08-06-dist-build-restructure-design.md)
  -- rationale for the three-output layout.

## Pipeline

`npm run build` runs `tsx scripts/build.ts`, which is the single source of
truth for `dist/`. It wipes the directory first, so partial rebuilds are not
a thing.

1. Clean `dist/` and recreate `browser/`, `node/`, `umd/`.
2. Emit declarations: the local `tsc --emitDeclarationOnly` against
   `tsconfig.build.json`, invoked through `process.execPath` rather than
   `npx` so Windows does not need a shell. The single `dist/index.d.ts` is
   read, then copied to `dist/browser/` and `dist/node/`. Worker declarations
   and all `.d.ts.map` files are deleted -- the maps would point at `src/`
   paths that are not published.
3. Browser bundle: esbuild, ESM, `platform:browser`, with the
   `browser-worker-swap` plugin.
4. Node bundle: esbuild, ESM, `platform:node`, `node:worker_threads` marked
   external.
5. UMD: esbuild IIFE with `globalName:'fflate'`, wrapped by hand in a
   CommonJS/AMD/global factory.
6. Minify. Each of steps 3 and 4 is minified in `writeBundle` with
   `module:true`; the UMD wrapper is minified with `module:false`.

Terser options are fixed at `mangle:{toplevel:true}` and
`compress:{passes:5, unsafe:true, pure_getters:true}`.

## Artifacts

| Path | Format | Declarations | Reachable via |
| --- | --- | --- | --- |
| `dist/browser/index.js` | ESM | `dist/browser/index.d.ts` | `.` (default), `./browser`, `module` |
| `dist/browser/index.min.js` | ESM | reuses the plain `.d.ts` | `./min` |
| `dist/node/index.js` | ESM | `dist/node/index.d.ts` | `.` (node), `./node`, `main` |
| `dist/node/index.min.js` | ESM | reuses the plain node `.d.ts` | `./node/min` |
| `dist/umd/fflate.js` | UMD | none, deliberately | `./umd`, `unpkg` |
| `dist/umd/package.json` | -- | -- | scopes the directory to CommonJS |

Every bundle also emits a `.js.map`, except UMD, which is built with
`sourcemap:false`.

Only `.` is platform-conditional. `./node`, `./node/min`, `./browser`, and
`./min` are all explicit picks, so `./min` means the minified **browser**
bundle on every platform. A node consumer importing `./min` gets working
sync APIs and async APIs that throw on `URL.createObjectURL`. `./node/min`
is the node-safe equivalent, added as its own subpath rather than by making
`./min` conditional, which would have silently redirected existing
consumers.

`test/3-node-min.ts` is the only coverage of `dist/node/index.min.js`. It
reaches the file by relative path rather than through the subpath, so it
exercises the bundle but not the `exports` entry.

Both minified bundles reuse the plain declarations for their tier, so a
drift between a plain and a minified bundle is invisible to the type
checker. Runtime coverage is what catches it: the browser suites run under
both `plain` and `min` labels, and `test/3-node-min.ts` compares export
names across the two node bundles.

`dist/umd/package.json` containing `{ "type": "commonjs" }` is load-bearing,
not build residue. The root package is `"type": "module"`, so without that
file node parses the UMD wrapper as ESM, where the wrapper's `this` is
undefined and loading throws `Cannot set properties of undefined (setting
'fflate')`. Verified on node v25.8.2 by removing and restoring it: present,
both `require()` and `import()` succeed; absent, both throw.

`./umd` carries no `types` entry, unlike every other subpath. Adding one
would type check code that fails at runtime, because the UMD runtime shapes
do not match the ESM ones.

## Invariants

These break silently. Sync APIs keep passing while async output corrupts, so
a green test run is not on its own evidence that they hold.

### `keepNames` must stay `false`

In all three esbuild configs. It wraps each function as
`__name(fn, 'name')`. The async APIs stringify functions and eval them inside
a worker where that helper is not in scope, giving
`ReferenceError: __name is not defined`. Name recovery for the worker payload
is `wcln`'s job, not esbuild's.

### Terser `compress` options are load-bearing

`wcln` recovers identifier names by parsing the **text** of the `bInflt` and
`bDflt` array literals (see [components.md](components.md)). Any option that
substitutes a constant's value for its name inside those literals breaks
recovery. Do not add compress options without re-running the browser suite
against the minified bundle.

### The two timeout numbers are a matched pair

`--timeout` in the `test:browser` script and `withTimeout`'s default in
`test/browser/util.ts` (`DEFAULT_TIMEOUT_MS`, currently 2000, plus
`TIMEOUT_MARGIN_MS` of 500).

tapout ends a run after `max(500, min(3000, floor(timeout * 0.2)))` ms of
console silence, and these suites are silent while polling. A stalled test
goes quiet, auto-finish fires, and every later test is dropped -- including
the entire minified pass -- while the run still reports PASS.
`test/browser/harness.ts` asserts the relationship rather than either
number, and fails if one drifts. Change them together.

### Two type assertions hide missing members

`Unzip.push` builds its result with a literal ending `} as UnzipFile`
(`src/index.ts:4033`). The assertion suppresses both the excess- and the
missing-property check a plain annotation would apply, so a new required
member on `UnzipFile` leaves the build green while the property is simply
absent at runtime. Adding `mtime` hit exactly this: `npm run build` and
`tsc --noEmit` both passed, and only the browser suite caught it. Every
required member added to `UnzipFile` must also be added to that literal.

The parallel `} as ZIFE` at `src/index.ts:3494` is deliberate -- `c`, `b`,
and `h` are assigned after construction -- but it hides the same class of
mistake.

Conversely, adding a required member to `UnzipFileInfo` **does** break the
build, in a place that looks unrelated. Declaration emit runs against
`tsconfig.build.json`, and the filter arguments in `unzip` and `unzipSync`
are contextually typed, so both call sites must change in the same commit.
Do not plan those as separate steps.

## Test surfaces

| Suite | Command | Covers |
| --- | --- | --- |
| `test/index.ts` -> `0-valid`, `1-size`, `2-perf`, `3-node-min` | `npm run test:node` | Correctness, size, perf against source; both node bundles. Numeric prefixes set run order -- tapzero has no CLI runner and flushes in registration order. |
| `test/browser/index.ts` -> `harness`, then `async`, `streams`, `zip` | `npm run test:browser` | Runs every suite twice, `plain` and `min`, against `dist/browser/`. The minified pass is the only thing that catches name-recovery breakage. |

`test/browser/index.ts` imports from `../../dist/browser/index.js`, and
`dist/` is gitignored, so `test:browser` fails on a clean checkout unless
`npm run build` ran first. `npm test` chains them; the standalone script
does not.

Run only one npm invocation at a time against this repo. They share `dist/`,
and a concurrent run produces an `npm test` exit 1 that looks like a real
failure and is not.

## Scripts that are destructive

`npm run gh-pages` (`scripts/cpGHPages.ts`) checks out the `gh-pages` branch
and `unlinkSync`s every top-level file before copying `public/` over it and
committing. `npm run build-docs` writes typedoc output into `docs/`, which
also holds hand-maintained plans and this inventory; it is safe only because
of `--cleanOutputDir false`.
