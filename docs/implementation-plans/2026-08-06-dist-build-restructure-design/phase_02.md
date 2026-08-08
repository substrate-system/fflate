# Dist Build Restructure Implementation Plan

## Phase 2: tsconfig restructure and node-worker static import

**Goal:** Replace the four TypeScript configurations with a single root
`tsconfig.json` targeting ES2022 ESM, and convert `src/node-worker.ts`
from a guarded `require` to a static `node:worker_threads` import.

**Architecture:** One root config covers `example`, `src`, `test`, and
`scripts`. `tsconfig.build.json` narrows that to `src` for the
declaration emit that Phase 3 drives. The worker shim loses its
unreachable fallback branch, which in turn removes the need for the
`createRequire` banner that `scripts/rewriteBuilds.ts` prepends today.

**Tech Stack:** TypeScript 6.0.3, Node types, Vite client types.

**Scope:** Phase 2 of 7.

**Codebase verified:** 2026-08-06

**Depends on:** Phase 1. The new `types` array includes `vite/client`,
which does not resolve until Phase 1 installs Vite, and the new `include`
array names `example`, which does not exist until Phase 1 renames
`demo/`.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC8: Type checking passes and obsolete configuration is gone

- **dist-build-restructure-design.AC8.1 Success:** `npx tsc --noEmit`
  against the root `tsconfig.json` reports zero errors across `example`,
  `src`, `test`, and `scripts`.
- **dist-build-restructure-design.AC8.3 Success:** `tsconfig.esm.json`,
  `tsconfig.demo.json`, and `test/tsconfig.json` no longer exist, and no
  file references them.

### dist-build-restructure-design.AC2: The worker swap is correct per build target

- **dist-build-restructure-design.AC2.2 Success:** `src/node-worker.ts`
  obtains `Worker` and `isMarkedAsUntransferable` from a static
  `node:worker_threads` import, with no `require` at module scope and no
  unsupported-runtime fallback branch.

---

## Investigation findings

Verified against the working tree on 2026-08-06.

Confirmed as the design describes:

- `tsconfig.esm.json`, `tsconfig.demo.json`, and `test/tsconfig.json` all
  exist. `test/tsconfig.json:4-5` sets `module` and `moduleResolution` to
  `nodenext`, which does conflict with `Bundler` resolution.
- `tsconfig.build.json` exists, extends `./tsconfig.json`, and excludes
  `["example", "test"]`. The design adds `scripts` to that list.
- `src/node-worker.ts:6-9` wraps `require('worker_threads')` in a
  `try`/`catch` with an empty catch, and lines 24 to 31 are the
  "async operations unsupported" fallback.
- `src/index.ts:15` is `import wk from './node-worker';`.
- `src/worker.ts:4` uses `||=` and `URL.createObjectURL` with `Blob`.
- `src/` contains exactly three files: `index.ts`, `node-worker.ts`,
  `worker.ts`.
- TypeScript is 6.0.3 and supports `moduleResolution: Bundler`.
- The current configuration type checks with zero errors.

Corrections and additions the design did not cover:

1. **The design's tsconfig omits JSX settings, which breaks the
   `example` include.** The design's `include` array lists `example`, and
   `example/` contains four `.tsx` files (`App.tsx`, `index.tsx`,
   `components/code-box/index.tsx`,
   `components/file-picker/index.tsx`). The design's `compilerOptions`
   block has no `jsx` key. Today those files are only ever checked by
   `tsconfig.demo.json`, which sets `"jsx": "preserve"` and which this
   phase deletes. Without a `jsx` setting the root config fails on every
   `.tsx` file. Task 1 adds `"jsx": "react-jsx"` with
   `"jsxImportSource": "preact"`, which matches the `@preact/preset-vite`
   runtime configured in Phase 1.

2. **The design's tsconfig omits `paths`, which the React specifiers
   need.** The example imports from `react` and `react-dom/client` while
   only `preact` is installed. Vite resolves these through aliases at
   build time, but `tsc` resolves them through `paths`. Task 1 adds the
   mapping. Without it, `tsc --noEmit` fails with `TS2307: Cannot find
   module 'react'`.

3. **`tsconfig.demo.json:11` excludes `demo/sw.ts` from type checking.**
   The service worker was never type checked under the old setup. The new
   root config has no such exclusion, so `example/sw.ts` is now checked.
   Phase 1 Task 6 already gives it a self-contained
   `ServiceWorkerGlobalScope` declaration and the design's `lib` array
   includes `WebWorker`, so this should pass. Task 4 verifies it.

4. **The current root `tsconfig.json:4` already sets
   `moduleResolution: "bundler"`** (lowercase). The design writes
   `"Bundler"`. These are equivalent; casing is not significant.

5. **`rootDir` must be set explicitly on `tsconfig.build.json`, and the
   design's replacement drops it.** TypeScript does not silently infer a
   usable root here. Measured on TypeScript 6.0.3 with `include:
   ["src/**/*"]` and `declarationDir` set:

   | configuration | result |
   | --- | --- |
   | no `rootDir` | `error TS5011: The common source directory ... is './src'. The 'rootDir' setting must be explicitly set`, and declarations emitted to `<dir>/src/index.d.ts` |
   | `"rootDir": "src"` | zero errors, declarations emitted to `<dir>/index.d.ts` |

   Phase 3 Task 3 reads `dist/index.d.ts` and would fail twice over
   without this: the `tsc` invocation exits non-zero on TS5011, and the
   file it wants is at `dist/src/index.d.ts`.

   `rootDir` belongs on `tsconfig.build.json` only, never on the root
   config. The root config includes `example`, `test`, and `scripts`
   alongside `src`, so a `src` root there is wrong for every other
   directory.

6. **The current config sets `"esModuleInterop": false` and
   `"ignoreDeprecations": "6.0"`**, neither of which appears in the
   design's replacement. Dropping `ignoreDeprecations` is safe once
   `target` and `module` move off `es5`/`commonjs`, which is what
   triggered the deprecation warning. Task 4 verifies this.

7. **`src/node-worker.ts:4` contains a `require` inside a string
   literal**, not at module scope:

   ```
   const workerAdd = ";var __w=require('worker_threads');..."
   ```

   That string is appended to worker source and evaluated by
   `new Worker(..., { eval: true })`, which Node evaluates as CommonJS
   regardless of the host package's `type` field. It must be left exactly
   as it is. Only the module-scope `require` on line 7 is replaced. This
   is why the design is correct that the `createRequire` banner becomes
   unnecessary: the banner existed solely for the line 7 `require`.

8. **`src/worker.ts` and `src/node-worker.ts` export identical
   signatures**, confirmed field by field:
   `<T>(c:string, id:number, msg:unknown, transfer:ArrayBuffer[], cb:(err:Error, msg:T) => void) => Worker`.
   The Phase 3 resolve plugin swap is therefore type safe.

---

<!-- START_SUBCOMPONENT_A (tasks 1-2) -->

<!-- START_TASK_1 -->
### Task 1: Replace the root tsconfig.json

**Verifies:** dist-build-restructure-design.AC8.1

**Files:**
- Modify: `tsconfig.json` (full replacement)

**Implementation:**

Replace the entire contents of `tsconfig.json`. This is the design's
configuration plus the two additions recorded in investigation findings 1
and 2, which the design's version is missing.

```json
{
  "compilerOptions": {
    "listFiles": true,
    "module": "ES2022",
    "target": "ES2022",
    "moduleResolution": "Bundler",
    "lib": ["ES2022", "DOM", "WebWorker"],
    "types": ["node", "vite/client"],
    "allowJs": false,
    "skipLibCheck": true,
    "outDir": "dist",
    "allowSyntheticDefaultImports": true,
    "experimentalDecorators": true,
    "emitDecoratorMetadata": true,
    "strict": false,
    "noImplicitAny": true,
    "sourceMap": true,
    "forceConsistentCasingInFileNames": true,
    "resolveJsonModule": true,
    "isolatedModules": true,
    "declaration": true,
    "declarationDir": "dist",
    "declarationMap": true,
    "jsx": "react-jsx",
    "jsxImportSource": "preact",
    "paths": {
      "react": ["./node_modules/preact/compat/"],
      "react/jsx-runtime": ["./node_modules/preact/jsx-runtime"],
      "react-dom": ["./node_modules/preact/compat/"],
      "react-dom/client": ["./node_modules/preact/compat/client"],
      "react-dom/*": ["./node_modules/preact/compat/*"]
    }
  },
  "include": [
    "example",
    "src/**/*",
    "test",
    "scripts"
  ]
}
```

Two keys here are not in the design document and are required. `jsx` and
`jsxImportSource` are needed because `include` names `example`, which
holds four `.tsx` files that were previously only checked by the
`tsconfig.demo.json` this phase deletes. `paths` is needed because the
example imports `react` and `react-dom/client`, which are not installed;
Vite aliases them at build time but `tsc` needs the mapping.

Note that `react-dom/client` maps to `preact/compat/client`, not
`preact/compat`. `createRoot` is exported only from the client entry, and
mapping to `preact/compat` gives
`TS2305: Module '"react"' has no exported member 'createRoot'`.

Note that `listFiles: true` is carried over from the template and makes
`tsc` print every file it loads. Build output is verbose by design.

**Verification:**

Do not run a bare `tsc` yet; Task 2 and Task 3 must land first, because
`test/tsconfig.json` still shadows resolution for the test directory and
`src/node-worker.ts` still has a module-scope `require`. Task 4 runs the
full check.

Run: `node -e "JSON.parse(require('fs').readFileSync('tsconfig.json','utf8').replace(/\/\/.*$/gm,''))" && echo "valid json"`
Expected: prints `valid json`.

**Commit:** `build: replace tsconfig with es2022 esm configuration`
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Delete the obsolete TypeScript configurations

**Verifies:** dist-build-restructure-design.AC8.3

**Files:**
- Delete: `tsconfig.esm.json`
- Delete: `tsconfig.demo.json`
- Delete: `test/tsconfig.json`
- Modify: `tsconfig.build.json`

**Implementation:**

`tsconfig.esm.json` drove the second `tsc` emit pass that Phase 3
replaces. `tsconfig.demo.json` type checked the Parcel demo that Phase 1
replaced. `test/tsconfig.json` sets `nodenext` module resolution, which
conflicts with the root `Bundler` setting; the node suite runs under
`tsx`, which does not need it.

```bash
git rm tsconfig.esm.json tsconfig.demo.json test/tsconfig.json
```

Then replace the contents of `tsconfig.build.json` with:

```json
{
    "extends": "./tsconfig.json",
    "compilerOptions": {
        "rootDir": "src"
    },
    "exclude": [
        "example",
        "test",
        "scripts"
    ]
}
```

Two changes. `scripts` joins the exclude list, per the design, which
keeps the Phase 3 declaration emit scoped to `src`.

`rootDir` is not in the design and is required. Without it TypeScript
6.0.3 raises `TS5011` and emits to `dist/src/index.d.ts`; with it the
emit is clean and lands at `dist/index.d.ts`, which is the path Phase 3
Task 3 reads. It is set here rather than on the root config because the
root config also includes `example`, `test`, and `scripts`, for which a
`src` root would be wrong.

**Verification:**

Run: `ls tsconfig.esm.json tsconfig.demo.json test/tsconfig.json 2>&1`
Expected: three "No such file or directory" errors.

Run: `grep -rn "tsconfig.esm\|tsconfig.demo\|test/tsconfig" --include='*.json' --include='*.ts' --include='*.js' . --exclude-dir=node_modules --exclude-dir=docs`
Expected: no output. If `package.json` still references `tsconfig.esm.json`
through the `build:lib` script, that reference is removed in Phase 4;
record it and continue.

**Commit:** `build: delete obsolete tsconfig files`
<!-- END_TASK_2 -->

<!-- END_SUBCOMPONENT_A -->

<!-- START_TASK_3 -->
### Task 3: Convert node-worker.ts to a static import

**Verifies:** dist-build-restructure-design.AC2.2

**Files:**
- Modify: `src/node-worker.ts` (full replacement)

**Implementation:**

Replace the entire contents of `src/node-worker.ts` with the following.

The guarded `require` and the fallback branch both go. The fallback
returned an "async operations unsupported, update to Node 12+" error,
which cannot be reached on any Node version this package now targets. Its
removal is what makes the `createRequire` banner in
`scripts/rewriteBuilds.ts` unnecessary.

The `workerAdd` string on the first line is unchanged and must stay byte
for byte identical. Its `require('worker_threads')` runs inside the
worker, which Node evaluates as CommonJS because of `{ eval: true }`, so
it is unaffected by this package being `"type": "module"`.

```ts
import { Worker, isMarkedAsUntransferable } from 'node:worker_threads'

const workerAdd = ";var __w=require('worker_threads');__w.parentPort.on('message',function(m){onmessage({data:m})}),postMessage=function(m,t){__w.parentPort.postMessage(m,t)},close=process.exit;self=global";

export default <T>(
    c:string,
    _:number,
    msg:unknown,
    transfer:ArrayBuffer[],
    cb:(err:Error, msg:T) => void
) => {
    let done = false
    const w = new Worker(c + workerAdd, { eval: true })
        .on('error', e => cb(e as Error, null))
        .on('message', m => cb(null, m))
        .on('exit', c => {
            if (c && !done) cb(new Error('exited with code ' + c), null)
        })
    if (isMarkedAsUntransferable) {
        transfer = transfer.filter(t => !isMarkedAsUntransferable(t))
    }
    w.postMessage(msg, transfer)
    w.terminate = () => {
        done = true
        return Worker.prototype.terminate.call(w)
    }
    return w
}
```

The exported signature is unchanged, which matters because
`src/worker.ts` must remain interchangeable with this module for the
Phase 3 resolve plugin swap.

**Verification:**

Run: `grep -n "require(" src/node-worker.ts`
Expected: exactly one match, on the `workerAdd` line, inside the string
literal. No match at module scope.

Run: `grep -n "async operations unsupported" src/node-worker.ts`
Expected: no output.

Run: `grep -n "^import { Worker, isMarkedAsUntransferable } from 'node:worker_threads'" src/node-worker.ts`
Expected: one match on line 1.

**Commit:** `refactor: use static node:worker_threads import`
<!-- END_TASK_3 -->

<!-- START_TASK_3B -->
### Task 3B: Resolve the newly exposed example type errors

**Verifies:** dist-build-restructure-design.AC8.1

**Files:**
- Create: `example/components/code-box/prism.d.ts`
- Modify: `example/components/code-box/index.tsx:314`
- Modify: `example/components/file-picker/index.tsx:106`

**Implementation:**

Bringing `example` under the root config exposes three errors that
`tsconfig.demo.json` never surfaced, because it did not set
`noImplicitAny` and did not use preact's types. All three are measured;
none is optional, since Task 4 gates on zero errors.

**Error 1, `TS7016` on `prism.js`.** `example/components/code-box/`
contains a vendored `prism.js`, and `allowJs` is `false` while
`noImplicitAny` is `true`, so importing it has no declaration to bind to.
Add `example/components/code-box/prism.d.ts` declaring the module's
surface. Determine the actual shape from the import site in `index.tsx`
rather than guessing; if only `highlight` and `languages` are used, the
declaration needs only those. Do not set `allowJs: true` to dodge this:
that would pull `prism.js` into the program and produce a much larger
error surface in vendored code.

**Error 2, `TS2322` at `code-box/index.tsx:314`.** `spellCheck` is not
present on preact/compat's `TextareaHTMLAttributes`. Preact's DOM
attribute casing differs from React's. Use the lowercase DOM attribute
name that preact expects, verifying against
`node_modules/preact/src/jsx.d.ts` rather than assuming the spelling.

**Error 3, `TS2698` at `file-picker/index.tsx:106`.** A spread of a value
TypeScript cannot prove is an object type. Narrow the value at the spread
site. Do not cast to `any`.

These are all in `example/`, which is this repository's own demo code,
not vendored upstream. Editing it is in scope, unlike `src/`.

**Verification:**

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep -E "error TS" | grep "example/"`
Expected: no output.

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep -cE "error TS7016|error TS2322|error TS2698" || true`
Expected: `0`.

**Commit:** `fix: resolve example type errors under strict resolution`
<!-- END_TASK_3B -->

<!-- START_TASK_4 -->
### Task 4: Verify the full type check

**Verifies:** dist-build-restructure-design.AC8.1

**Files:**
- No files changed. This task is a gate.

**Implementation:**

Run the type check across all four included directories and drive the
error count to zero.

**Verification:**

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep -E "error TS" | head -40`
Expected: no output.

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep -cE "error TS" || true`
Expected: `0`.

If errors appear, resolve them by category rather than by adding
suppressions:

- `TS2688: Cannot find type definition file for 'vite/client'` means
  Phase 1 Task 1 did not complete. Install Vite.
- `TS17004` or `TS6142` on `.tsx` files means the `jsx` option is missing
  or misspelled. Re-check Task 1.
- `TS2307: Cannot find module 'react'` means the `paths` block is missing
  or `baseUrl` semantics differ in this TypeScript version. Verify the
  `paths` entries resolve by running
  `ls node_modules/preact/compat/index.d.ts`.
- `TS2305: Module '"react"' has no exported member 'createRoot'` means
  `react-dom/client` is mapped to `preact/compat` instead of
  `preact/compat/client`. Re-check Task 1.
- `TS7016`, `TS2322`, or `TS2698` in `example/` means Task 3B did not
  complete.
- Errors in `src/index.ts` are unexpected. The design measured zero
  errors at `strict: false` with `noImplicitAny: true`, and this was
  reconfirmed on 2026-08-06. Do not edit `src/index.ts` to silence them;
  it is vendored upstream code and this repository merges from
  `git remote upstream`. Report the errors instead.

Run: `npx tsc --noEmit --project tsconfig.build.json 2>&1 | grep -cE "error TS" || true`
Expected: `0`. A `TS5011` here means `rootDir` is missing from
`tsconfig.build.json`; see Task 2.

Run: `npx tsc --emitDeclarationOnly --project tsconfig.build.json && ls dist/index.d.ts && rm -rf dist`
Expected: `dist/index.d.ts` exists. This is a dry run of the emit that
Phase 3 Task 3 depends on; if the file appears at `dist/src/index.d.ts`
instead, `rootDir` is not taking effect.

**Commit:** `build: verify typecheck across all sources`
<!-- END_TASK_4 -->

---

## Phase 2 completion criteria

- `npx tsc --noEmit --project tsconfig.json` reports zero errors.
- `npx tsc --noEmit --project tsconfig.build.json` reports zero errors.
- `tsconfig.esm.json`, `tsconfig.demo.json`, and `test/tsconfig.json` are
  deleted.
- `tsconfig.build.json` excludes `example`, `test`, and `scripts`.
- `src/node-worker.ts` has no module-scope `require` and no fallback
  branch, and its `workerAdd` string is byte for byte unchanged.

`npm run build` still runs the old pipeline and will now fail, because
`build:lib` invokes `tsc --project tsconfig.esm.json` against a deleted
file. That is expected and is resolved in Phase 3, which replaces the
pipeline, and Phase 4, which rewrites the scripts. No published artifact
is produced between here and Phase 4.
