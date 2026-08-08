# fflate (@substrate-system/fflate)

Last verified: 2026-08-08

Fork of fflate. High performance (de)compression. ESM-only package,
bundled with esbuild at `es2022` into three outputs under `dist/`.

## Tech Stack

- Language: TypeScript, `target: es2022`
- Bundler: esbuild; minifier: terser
- Runner: tsx
- Tests: tapzero (assertions), tapout (browser runner)
- Example app: vite + preact

## Commands

- `npm run build` -- `tsx scripts/build.ts`. Single source of truth for
  every artifact. Wipes and rebuilds `dist/`.
- `npm test` -- build, then node suite, then browser suite.
- `npm run test:node` -- `tsx test/index.ts`.
- `npm run test:browser` -- bundles `test/browser/index.ts` and pipes it
  to tapout. Requires a prior `npm run build`; see Gotchas.
- `npm start` -- vite dev server for `example/`.

Run only ONE npm invocation at a time against this repo. They share
`dist/`, and a concurrent run yields an `npm test` exit 1 that looks
like a real failure and is not.

## Never Run These

- `npm run gh-pages` / `scripts/cpGHPages.ts` -- checks out the
  `gh-pages` branch and `unlinkSync`s every top-level file before
  copying and committing. Destructive to that branch.
- `npm run build-docs` -- writes typedoc output into `docs/`, which
  also holds implementation plans.

## Project Structure

- `src/` -- `index.ts` (all public exports), `worker.ts` (browser
  worker), `node-worker.ts` (node worker). The build swaps
  `./node-worker` for `./worker` in browser and UMD bundles.
- `scripts/` -- exactly `build.ts` and `cpGHPages.ts`.
- `test/` -- node suite (`index.ts` imports `0-valid`, `1-size`,
  `2-perf` in that order; numeric prefixes set run order).
- `test/browser/` -- zip, streams, and async coverage. Runs every suite
  twice, plain and minified.
- `example/` -- vite demo app.
- `docs/` -- mixed. `README.md`, `classes/`, `functions/`,
  `interfaces/`, `type-aliases/`, and `variables/` are GENERATED
  typedoc output. `architecture/`, `ARCHITECTURE.md`,
  `implementation-plans/`, and `test-plans/` are hand-maintained.
- `specs/` -- design documents.
- `dist/` -- GENERATED, gitignored.

Do not hand-edit anything in `dist/`, or the generated parts of
`docs/` listed above. `build-docs` passes `--cleanOutputDir false`,
which is the only reason the hand-maintained files under `docs/`
survive a regeneration. Keep that flag.

`docs/architecture/` is the runtime inventory: what exists, where it
lives, and its current contract. Start at `docs/architecture/INDEX.md`.
Update it when you change the build pipeline, the shipped artifacts,
the worker execution path, or the caches in `src/index.ts` and
`src/worker.ts`.

## Packaging

`"type": "module"`, `files: ["dist"]`. The `exports` map advertises no
`require` condition; subpaths are `.`, `./node`, `./node/min`,
`./browser`, `./min`, `./umd`, `./package.json`.

`./min` is the minified BROWSER build, not a platform-conditional one.
Only `.` switches on the `node` condition; every other subpath is an
explicit pick. A node consumer importing `./min` therefore gets the
browser bundle, whose async APIs call `URL.createObjectURL` and throw
under node while the sync APIs keep working. That is why `./node/min`
exists as a separate subpath rather than `./min` being made
conditional -- changing `./min` would silently redirect existing
consumers.

`./umd` is the exception. It is a bare string pointing at
`./dist/umd/fflate.js`, and `build.ts` writes
`dist/umd/package.json` containing `{ "type": "commonjs" }`. That file
is load-bearing, not build residue -- someone will eventually try to
delete it as a stray. Without it, the root `"type": "module"` makes
node parse the UMD wrapper as ESM, where the wrapper's `this` is
undefined and loading throws `Cannot set properties of undefined
(setting 'fflate')`. Measured on node v25.8.2, by removing the file
and restoring it: with it present both `require()` and `import()` of
`dist/umd/fflate.js` succeed; with it removed BOTH throw that same
`TypeError`. Node's `require(esm)` support does not rescue this --
the wrapper fails on `this` before the module format matters. Keep
the file.

## Invariants

These are silent when broken. Sync APIs keep passing while async
output corrupts, so tests staying green is not evidence.

- `keepNames` must stay `false` in every esbuild config. It wraps each
  function as `__name(fn, 'name')`, and the async APIs stringify
  functions and eval them inside a worker where that helper is not in
  scope, giving `ReferenceError: __name is not defined`. Name recovery
  for the worker payload is `wcln`'s job, not esbuild's.
- `DeflateState.b` is the index of the first legal back-reference target
  in `dat`, and `dflt`'s `maxd` bound is `Math.min(32767, i - (st.b || 0))`.
  It is 0 for `deflateSync` and the one-shot paths, where `dat` is the
  stream. The streaming classes deflate out of a 98304-byte scratch
  buffer whose leading 32768 bytes are reserved lookback space, so their
  base starts at 32768 (or at the dictionary start) and drops to 0 at the
  wrap copy in `Deflate.push`. Any new code that moves data within that
  buffer, or that reserves space in front of it, must set `s.b` in the
  same edit. Getting it wrong emits a distance past the start of the
  stream, which every inflater rejects -- but only for inputs whose tail
  happens to match the zero pad, so almost all tests stay green. See
  `specs/2026-08-08-streaming-deflate-invalid-distance.md` and
  `test/10-stream-window.ts`. A base that is too LARGE is silent in the
  other direction: it just costs ratio, so keep that test's long-range
  match assertion.
- Do not add terser `compress` options without re-running the browser
  suite against the minified bundle. `wcln` recovers identifier names
  by parsing the TEXT of the `bInflt`/`bDflt` array literals. Any
  option that substitutes a constant's value for its name inside those
  literals breaks recovery silently.
- The `--timeout` in the `test:browser` script and `withTimeout`'s
  default in `test/browser/util.ts` are a matched pair. tapout ends a
  run after `max(500, min(3000, floor(timeout * 0.2)))` ms of console
  silence, and these suites are silent while polling. A stalled test
  goes quiet, auto-finish fires, and every later test is dropped --
  including the entire minified pass -- while the run still reports
  PASS. `test/browser/harness.ts` asserts the relationship and fails if
  either number drifts. Change them together.
- Every required member you add to `UnzipFile` must also be added to
  the object literal in `Unzip.push`. That literal ends `} as UnzipFile`
  (`src/index.ts:4033`), and the assertion suppresses the excess- and
  missing-property checks a plain annotation would apply, so the build
  stays green while the property is simply absent at runtime. Adding
  `mtime` hit exactly this: `npm run build` and `tsc --noEmit` both
  passed with the interface member declared and the literal never
  updated, and only the browser suite caught it. The parallel
  `} as ZIFE` at `src/index.ts:3494` is deliberate -- `c`, `b`, and `h`
  are assigned after construction -- but it hides the same class of
  mistake.
- Conversely, adding a required member to `UnzipFileInfo` DOES break the
  build, and breaks it in a place that looks unrelated. `npm run build`
  runs `tsc --emitDeclarationOnly` against `tsconfig.build.json`, and the
  filter arguments in `unzip` and `unzipSync` are contextually typed, so
  both call sites must be updated in the same commit or declaration emit
  fails. Do not plan those two as separate steps.

## Gotchas

- `dist/` is gitignored, and `test/browser/index.ts` imports from
  `../../dist/browser/index.js`. `npm run test:browser` therefore fails
  on a clean checkout unless `npm run build` ran first. `npm test`
  chains them; the standalone script does not.
- Never verify anything in this repo with `node --input-type=module
  -e`. That form puts `--input-type=module` into `process.execArgv`,
  worker threads inherit `execArgv`, and the worker then parses its
  eval'd payload as ESM, where fflate's sloppy-mode implicit globals
  throw `ReferenceError: u8 is not defined`. It is an artifact of the
  command, not a defect in the code. Use a real `.mjs` or `.cjs` file
  and assert the exit code.
- `./min` has no declarations of its own. Its `types` entry reuses
  `dist/browser/index.d.ts`, so a drift between the plain and minified
  bundles would not be caught by the type checker. The runtime coverage
  is the full test suites in `test/browser/` running under both `plain`
  and `min` labels; `test/browser/min.ts` holds a type-level cast
  documenting that the two export names are expected to match.
- `./umd` deliberately carries NO `types` entry, unlike every other
  subpath. Adding one type checks code that fails at runtime: the UMD
  wrapper's runtime shapes do not match the ESM ones.

## Key Files

- `scripts/build.ts` -- the whole pipeline, heavily commented with the
  reasoning behind each non-obvious flag.
- `test/browser/harness.ts` -- guards the timeout invariant above.
- `test/3-node-min.ts` -- the only coverage of `dist/node/index.min.js`,
  which ships but which the browser suite cannot reach. Its async case
  goes through `worker_threads`, so it is the node-side counterpart to
  the browser `min` pass.
- `test/browser/util.ts` -- `withTimeout`, `autoFinishWindow`, and the
  timeout constants.
