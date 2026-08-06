# Build restructure: single `dist/`, ESM-only, template-aligned

Date: 2026-08-06
Status: approved, not yet implemented

## Goal

Restructure `@substrate-system/fflate` to match the layout and npm script
conventions of `../template-ts-browser`:

- All compiled output goes to `dist/`, organised into subdirectories.
- CommonJS is dropped as a published format.
- UMD, node, and browser artifacts are all still built.
- The template's npm lifecycle hooks are adopted.

## Why this is not a mechanical change

Three properties of the current codebase constrain the design. Each is
load-bearing and any implementation must preserve it.

### 1. The worker swap

`src/index.ts` imports a worker factory:

```ts
import wk from './node-worker';
```

The browser build must substitute `src/worker.ts` for `src/node-worker.ts`.
This single relative-path swap is the entire reason node and browser are
separate artifacts, and it is what `scripts/rewriteBuilds.ts` exists to do
via regular expressions over emitted JavaScript.

esbuild's `--alias:` flag only rewrites bare package specifiers, not
relative paths, so the CLI cannot express this. A resolve plugin using the
esbuild JS API can, which is why the build moves from npm-script CLI
invocations to a single script.

### 2. `wcln` recovers variable names from function source text

`src/index.ts:1040` contains the most fragile code in the project:

```ts
const wcln = (fn: () => unknown[], fnStr: string, td: Record<string, unknown>) => {
  const dt = fn();
  const st = fn.toString();
  const ks = st.slice(st.indexOf('[') + 1, st.lastIndexOf(']'))
    .replace(/\s+/g, '').split(',');
  ...
}
```

The async APIs construct worker source at runtime by calling `.toString()`
on thunks such as:

```ts
const bInflt = () => [u8, u16, i32, fleb, fdeb, clim, fl, fd, /* ... */];
```

It then parses the identifier names back out of the array-literal *text*
and pairs them positionally with the evaluated values, reassigning each
under the recovered name inside the worker. The comment at `src/index.ts:1037`
states the contract directly: "The reason we can't just use the original
variable names is minifiers mangling the toplevel scope."

Two consequences for the build:

- Any minifier in the pipeline must rename the top-level scope
  **consistently**, so that the names in the array-literal text still match
  the names referenced inside the serialised function bodies.
- No minifier may **inline a constant into that array literal**. If `fleb`
  were substituted for its value, name recovery would silently produce
  wrong keys.

The current build sidesteps this: `lib/` and `esm/` ship unminified `tsc`
output, and terser (with `mangle.toplevel`) is applied only when building
UMD. Adding minified ESM outputs is therefore **new exposure**, and is the
single highest-risk part of this work.

### 3. The node test harness depends on a CJS build existing

`test/util.ts:167`:

```ts
const fflate = resolve(here, '..', 'lib', 'index.cjs');
```

The benchmark harness generates `worker_threads` workers whose source
`require()`s that absolute path (`test/util.ts:113-114`), deliberately
bypassing the `exports` map. `worker_threads` with `{ eval: true }`
evaluates as CommonJS, so dropping CJS breaks the harness at its
foundation, independently of any test-runner migration.

## Design

### Output layout

```
dist/
  browser/
    index.js        index.js.map
    index.min.js    index.min.js.map
    index.d.ts
  node/
    index.js        index.js.map
    index.min.js    index.min.js.map
    index.d.ts
  umd/
    fflate.js       (minified)
```

### tsconfig

`tsconfig.json` currently targets the old CommonJS build and must change:

| option | from | to |
| --- | --- | --- |
| `target` | `es5` | `es2022` |
| `module` | `commonjs` | `es2022` |
| `outDir` | `lib/` | `dist/` |
| `lib` | (implicit) | `["ES2022", "DOM", "WebWorker"]` |
| `types` | `["node"]` | `["node"]` (kept) |

`lib` must be set explicitly. `src/worker.ts` references `Worker`,
`URL.createObjectURL`, and `Blob`, which the current configuration only
picks up via the default lib for an ES5 target.

`tsconfig.esm.json` and `tsconfig.demo.json` are deleted.
`tsconfig.build.json` is retained; its `exclude` list gains `example` in
place of the removed `demo`.

### Build pipeline

A single `scripts/build.ts`, run under `tsx`, replaces `scripts/rewriteBuilds.ts`,
`scripts/buildUMD.ts`, and both `tsc` emit passes.

Steps:

1. Remove `dist/`.
2. **Declarations.** One `tsc --emitDeclarationOnly` pass. `index.d.ts` is
   copied to both `dist/browser/` and `dist/node/`; the public API is
   identical between the two. Emitted `worker.d.ts` and `node-worker.d.ts`
   are discarded, since both builds are bundled.
3. **Browser.** esbuild bundle of `src/index.ts`, `format: 'esm'`,
   `platform: 'browser'`, `target: 'es2022'`, `keepNames: true`, with the
   resolve plugin mapping `./node-worker` to `./worker`.
4. **Node.** esbuild bundle of `src/index.ts`, `format: 'esm'`,
   `platform: 'node'`, `target: 'es2022'`, `keepNames: true`, with
   `node:worker_threads` marked external.
5. **UMD.** esbuild bundle of the browser variant, `format: 'iife'` with a
   `globalName`, wrapped in a UMD preamble exposing the global `fflate`.
6. **Minification.** terser, applied to the esbuild output of steps 3-5.

#### esbuild bundles, terser minifies

esbuild's own `--minify` is **not** used. All minified artifacts go through
terser with the settings `scripts/buildUMD.ts` has already proven safe
against `wcln`:

```js
{ mangle: { toplevel: true },
  compress: { passes: 5, unsafe: true, pure_getters: true } }
```

ESM outputs additionally set `module: true`.

Rationale: mixing esbuild's minifier for `.min.js` with terser for UMD
would apply two different mangling strategies to code whose correctness
depends on mangling being consistent. One mangler, already validated
against this codebase, is the lower-risk choice.

#### Target

ES2022, up from the current ES5. This is a deliberate compatibility
reduction, consistent with dropping CommonJS. It also stops `worker.ts`'s
`||=` from being downleveled.

### `src/node-worker.ts`

Replace the `require('worker_threads')` in a `try`/`catch` with a static
import:

```ts
import { Worker, isMarkedAsUntransferable } from 'node:worker_threads';
```

The fallback branch that returns the "async operations unsupported - update
to Node 12+" error is removed along with it; the import cannot fail on any
Node version this package now targets. The `createRequire` banner that
`rewriteBuilds.ts` prepends to `esm/index.mjs` is no longer needed.

The browser build never sees this file, because the resolve plugin swaps in
`src/worker.ts`.

### package.json

```json
"main":   "./dist/node/index.js",
"module": "./dist/browser/index.js",
"types":  "./dist/browser/index.d.ts",
"unpkg":  "./dist/umd/fflate.js",
"files":  ["dist"],
"exports": {
  ".": {
    "node":    { "types": "./dist/node/index.d.ts",
                 "default": "./dist/node/index.js" },
    "default": { "types": "./dist/browser/index.d.ts",
                 "default": "./dist/browser/index.js" }
  },
  "./node":    { "types": "./dist/node/index.d.ts",
                 "default": "./dist/node/index.js" },
  "./browser": { "types": "./dist/browser/index.d.ts",
                 "default": "./dist/browser/index.js" },
  "./min":     { "types": "./dist/browser/index.d.ts",
                 "default": "./dist/browser/index.min.js" },
  "./umd":     "./dist/umd/fflate.js",
  "./package.json": "./package.json"
}
```

Scripts:

| script | purpose |
| --- | --- |
| `build` | `tsx scripts/build.ts` |
| `build-example` | vite build of `example/` into `public/` |
| `gh-pages` | `tsx scripts/cpGHPages.ts` |
| `build-docs` | typedoc, output `docs/` (unchanged) |
| `start` | vite dev server |
| `test` | `build`, then `test:node` and `test:browser` |
| `test:node` | `tsx test/index.ts` |
| `test:browser` | `esbuild test/browser/index.ts --bundle \| tapout` |
| `toc` | markdown-toc |
| `version` | toc, auto-changelog, stage |
| `postversion` | push with tags, publish |
| `prepublishOnly` | `build` |

The half-merged template scripts currently in `package.json`
(`build-cjs`, `build-esm`, `build-esm:min`, `build-cjs:min`, `build:lib`,
`build:umd`, `build:rewrite`, `build:demo`, `script`, `//build`) are all
removed. The generic `script` runner (`tsx scripts/$SC.ts`) goes away with
them; its only remaining caller, `cpGHPages.ts`, gets the named `gh-pages`
script above.

`.gitignore` and `.npmignore` lose their `lib/`, `esm/`, and `umd/` entries.

### Linting

No change. fflate has no `eslint.config.js` and none is added, so there is
no `lint` script and no `preversion` hook (the template's `preversion`
exists solely to run lint). `newneostandard`, `eslint`, and
`typescript-eslint` are not added as devDependencies.

This is a deliberate deviation from the template. `src/index.ts` is 3,885
lines of vendored upstream code: 192 lines exceed 80 columns and 1,289 end
in a semicolon. The template config sets 4-space indent, no-space
`key-spacing`, and standard's no-semicolon style, so autofixing would
rewrite substantially every line. The repository tracks
`git remote upstream` at `101arrowz/fflate` and merges from it; a
wholesale reformat would make every future merge a conflict.

### Tests

`npm test` builds, then runs both suites.

**Node suite** (`test:node`, unchanged runner `tsx test/index.ts`).
Preserves the existing correctness and benchmark coverage against zlib,
pako, uzip, tiny-inflate, and jszip, none of which has a browser
equivalent. Changes:

- `test/util.ts:167` points at `dist/node/index.js` instead of
  `lib/index.cjs`.
- The `wc()` worker generator emits ES module source and spawns workers via
  a `data:text/javascript` URL rather than `{ eval: true }`, which is
  CommonJS-only. Comparison libraries that ship CommonJS are reached
  through Node's ESM-CommonJS interop.
- The fixture downloader, `perf_hooks` timing, and
  `test/results/*.json` writing are untouched.

**Browser suite** (`test:browser`, new). `test/browser/index.ts`, bundled
by esbuild and piped to `tapout`, importing from `dist/browser/`. This is
new coverage: the browser build is currently untested.

The three files that are presently empty TODO stubs **move** out of the
node suite and into `test/browser/`, where they are written for the first
time:

| from | to | current contents |
| --- | --- | --- |
| `test/3-zip.ts` | `test/browser/zip.ts` | `// TODO: test ZIP` |
| `test/4-streams.ts` | `test/browser/streams.ts` | `// TODO: test all streams (including ZIP)` |
| `test/5-async.ts` | `test/browser/async.ts` | `// TODO: test all async operations` |

`test/index.ts` therefore drops its `./3-zip.js`, `./4-streams.js`, and
`./5-async.js` imports, retaining only `./0-valid.js`, `./1-size.js`, and
`./2-perf.js`.

The browser suite must exercise the async and streaming APIs **against
`dist/browser/index.min.js`**, not only the unminified bundle. That is the
only real check on the `wcln` hazard described above, and the reason the
browser suite is worth writing as part of this change rather than after it.

### example/

`demo/` moves to `example/` and builds with vite plus `@preact/preset-vite`,
following the template's `vite.config.js` (root `example`, `publicDir`
`_public`, output `public/`). The application itself is kept: `App.tsx`,
`index.tsx`, `index.css`, `components/code-box` (prism-based),
`components/file-picker`, and `util/workers.ts`.

`demo/sw.ts` imports `manifest` and `version` from `@parcel/service-worker`.
It switches to `vite-plugin-pwa` in `injectManifest` mode; the service
worker body already exists and only the manifest source changes.

`parcel` and `@parcel/service-worker` are removed from devDependencies.
`build-example` runs vite with `--base="/fflate"` (the template ships a
placeholder `/repo-name`).

`scripts/cpGHPages.ts` **must** be repointed from `dist/` to `public/`.
Under parcel the demo was emitted to `dist/`, which cpGHPages then copied
to the `gh-pages` branch. After this change `dist/` holds the library and
the vite config emits the demo to `public/`; leaving cpGHPages untouched
would publish the compiled library to `gh-pages` instead of the demo.

`public/` is added to `.gitignore`, which does not currently list it.

## Risks

1. **Terser-minified ESM and `wcln`.** The highest-risk item. Mitigated by
   running the browser async suite against the minified bundle. If terser's
   `compress` inlines a constant into a `bInflt`/`bDflt` array literal,
   name recovery breaks silently -- the sync APIs keep working while the
   async ones corrupt output. If this occurs, the fallback is to disable
   the offending `compress` option for ESM outputs, or ship
   `dist/*/index.min.js` built with `mangle` only.
2. **ES2022 target** narrows the supported runtime range relative to the
   current ES5 output.
3. **Module workers in the node test harness.** Reworking `wc()` is the
   least mechanical part of the test changes; CommonJS comparison libraries
   under ESM interop is where problems would surface.
4. **Dropping CommonJS is breaking** for any consumer using `require()`.
   The UMD artifact remains as a fallback for script-tag and legacy-bundler
   consumers.

## Out of scope

- Reformatting or linting `src/`.
- Changing typedoc configuration or the tracked `docs/` output.
- Any change to the compression algorithms themselves.
