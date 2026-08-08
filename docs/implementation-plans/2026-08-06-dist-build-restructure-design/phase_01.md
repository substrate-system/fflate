# Dist Build Restructure Implementation Plan

## Phase 1: Dependencies and example/ migration

**Goal:** Replace the Parcel demo toolchain with Vite, move `demo/` to
`example/`, and install the dependency set the later phases require.

**Architecture:** The example app is a Preact application that still uses
React import specifiers. `@preact/preset-vite` aliases those to
`preact/compat`. Vite builds from `example/` into `public/` at the repo
root. The service worker keeps its hand-written body and takes its
precache list from the manifest that `vite-plugin-pwa` injects.

**Tech Stack:** vite 8, @preact/preset-vite 2, vite-plugin-pwa 1, preact
10, TypeScript 6.

**Scope:** Phase 1 of 7.

**Codebase verified:** 2026-08-06

---

## Verification conventions

These apply to every phase in this plan.

`grep -c` exits non-zero when the count is zero, so a check whose
expected value is `0` is written with `|| true` appended. Under a runner
that stops on first error, omitting it turns a passing check into a
reported failure.

`git grep -E` does not honour `\b`. Where a word boundary is needed, use
`git grep -P`. A `-E` pattern containing `\b` silently matches nothing
and reports a false pass.

Invoking `tsc` on a named file while `tsconfig.json` exists fails with
`TS5112`. Add `--ignoreConfig` for single-file checks.

Scratch files go in the session scratchpad
(`/private/tmp/claude-501/-Users-nick-code-fflate/e8c45c98-8ea8-4895-8a65-48cc175c6ce5/scratchpad/`),
not `/tmp`.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC7: The example app builds and deploys from `public/`

- **dist-build-restructure-design.AC7.1 Success:** `npm run build-example`
  exits 0 and writes `index.html` plus hashed assets into `public/` at the
  repository root.
- **dist-build-restructure-design.AC7.2 Success:** The example's React
  import specifiers, including the `react-dom/client` subpath used by
  `example/index.tsx`, resolve to `preact/compat` at build time with
  neither `react` nor `react-dom` installed.
- **dist-build-restructure-design.AC7.3 Success:** The built service worker
  precaches the application assets using the manifest injected by
  `vite-plugin-pwa` rather than `@parcel/service-worker`.
- **dist-build-restructure-design.AC7.4 Success:** `scripts/cpGHPages.ts`
  reads the directory that the example build writes (`public/`), not
  `dist/`.

---

## Investigation findings

These were verified against the working tree on 2026-08-06. They correct
or extend the design document.

Confirmed as the design describes:

- `demo/` exists with `App.tsx`, `index.tsx`, `index.css`,
  `components/code-box`, `components/file-picker`, `util/workers.ts`,
  `sw.ts`, `index.html`, `favicon.ico`, `global.d.ts`, `augment.d.ts`.
- `demo/sw.ts:2` imports `{ manifest, version }` from
  `@parcel/service-worker`.
- `parcel@^2.9.3` and `@parcel/service-worker@^2.9.3` are both in
  devDependencies.
- `scripts/cpGHPages.ts:18` and `:20` read and copy from `dist`.
- `.gitignore` does not list `public/`.
- No `example/`, no `vite.config.*`, and no `public/` exist yet.

Corrections and additions the design did not cover:

1. **The example is written against React, and React is not installed.**
   `demo/index.tsx:1`, `demo/App.tsx:1`,
   `demo/components/code-box/index.tsx:1`, and
   `demo/components/file-picker/index.tsx:1` import from `react`;
   `demo/index.tsx:3` imports `createRoot` from `react-dom/client`.
   Neither `react` nor `react-dom` resolves (`require.resolve` fails for
   both). Only `preact@^10.25.4` is installed. The demo therefore cannot
   build today. `@preact/preset-vite` resolves this, but see the next
   point.

2. **`@preact/preset-vite` does not alias the `react-dom/client`
   subpath, and the alias target is `preact/compat/client`.** The preset
   aliases `react`, `react-dom`, `react/jsx-runtime`, and
   `react-dom/test-utils` only. `example/index.tsx` imports `createRoot`
   from `react-dom/client`, so an explicit `resolve.alias` entry is
   required. That entry must point at `preact/compat/client`, not
   `preact/compat`: verified that `createRoot` appears in no emitted
   `.js` or `.d.ts` under `node_modules/preact/compat/dist/` (only in the
   embedded sources of the `.map` files), that
   `node_modules/preact/compat/client.mjs` exists, and that preact's
   `exports` map declares `./compat/client` as a separate entry.
   Targeting `preact/compat` yields
   `TS2305: Module '"react"' has no exported member 'createRoot'`.

3. **`vite-plugin-pwa` declares `workbox-build` and `workbox-window` as
   peer dependencies** (`^7.4.1` each). They must be installed even though
   the service worker body here does not import any `workbox-*` runtime
   package.

4. **`vite-plugin-pwa` provides no equivalent of Parcel's `version`
   export.** It injects an array of `{ url, revision }` entries at
   `self.__WB_MANIFEST`. The existing `sw.ts` body uses `version` as the
   cache name. Task 6 keeps the hand-written body and supplies the cache
   name from a build-time constant, which is the smallest change
   consistent with the design's statement that "the service worker body
   already exists and only the manifest source changes".

5. **`package.json` already sets `"type": "module"`.** No change is needed
   for the ESM `.js` output that later phases emit, but it does mean
   `vite.config.js` and `scripts/cpGHPages.ts` cannot use `__dirname`.
   `scripts/cpGHPages.ts:5` currently uses `__dirname` and must be
   converted.

6. **`scripts/cpGHPages.ts` has two defects beyond the output path.**
   Line 15 calls `statSync(f)` with a bare entry name rather than
   `statSync(to(f))`, so it stats relative to the process working
   directory. Line 24 checks out `master`, which is wrong on any other
   branch. Both are fixed in Task 7.

7. **`demo/index.tsx:5` reads `process.env.NODE_ENV`**, which Vite does not
   define. It becomes `import.meta.env.PROD`.

8. **`demo/index.tsx:7` registers the service worker** with
   `new URL('sw.ts', import.meta.url)`, a Parcel idiom. Vite emits `sw.js`
   at the site root, so the registration path changes.

9. **The example imports the library by relative root path.**
   `demo/components/code-box/sandbox.ts:1` uses `import * as fflate from
   '../../..'` and `stream-adapter.ts:1` uses `import { AsyncDeflate } from
   '../../..'`. The directory depth is unchanged by the move, but resolving
   through the root `package.json` would couple the example build to the
   library build. Task 4 adds an alias to `src/index.ts` so the example
   builds from source and this phase does not depend on Phases 3 or 4.

10. **Stale build output is present.** `dist/`, `lib/`, `esm/`, and `umd/`
    all exist in the working tree from the previous half-merged build.
    They are removed in Phase 7, not here.

11. **No `CLAUDE.md` or `AGENTS.md` exists anywhere in the repository.**
    Project code style comes from the user's global instructions: 80
    column limit, no space between colon and type annotation, ternary
    operator with the `?` ending the line, no em dashes or arrow
    characters in comments.

---

<!-- START_SUBCOMPONENT_A (tasks 1-3) -->

<!-- START_TASK_1 -->
### Task 1: Install the Vite toolchain and remove Parcel

**Verifies:** None (infrastructure)

**Files:**
- Modify: `package.json` (devDependencies)

**Step 1: Remove the Parcel dependencies**

```bash
npm uninstall parcel @parcel/service-worker
```

**Step 2: Install the Vite toolchain**

`workbox-build` and `workbox-window` are peer dependencies of
`vite-plugin-pwa` and must be present even though the service worker body
does not import them.

```bash
npm install --save-dev \
  vite@^8.2.1 \
  @preact/preset-vite@^2.10.6 \
  vite-plugin-pwa@^1.3.0 \
  workbox-build@^7.4.1 \
  workbox-window@^7.4.1
```

**Step 3: Remove the unused React type packages**

`preact/compat` ships its own types, so `@types/react` and
`@types/react-dom` are no longer referenced once the alias is in place.

```bash
npm uninstall @types/react @types/react-dom
```

**Step 4: Verify**

Run: `npm ls vite @preact/preset-vite vite-plugin-pwa preact`
Expected: all four resolve with no `UNMET DEPENDENCY` lines.

Run: `node -e "require.resolve('vite')"`
Expected: exits 0 with no output.

Run: `node -e "try{require.resolve('parcel');console.log('STILL PRESENT')}catch(e){console.log('removed')}"`
Expected: prints `removed`.

**Step 5: Commit**

```bash
git add package.json package-lock.json
git commit -m "build: replace parcel toolchain with vite"
```
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Move demo/ to example/ and add the static asset directory

**Verifies:** None (infrastructure)

**Files:**
- Rename: `demo/` to `example/`
- Create: `example/_public/`
- Move: `example/favicon.ico` to `example/_public/favicon.ico`
- Delete: `example/package.json`

**Step 1: Move the directory with git so history follows**

```bash
git mv demo example
```

**Step 2: Create the static asset directory and move the favicon**

The design sets Vite's `publicDir` to `_public`. `example/index.html:8`
references `favicon.ico` at the site root, which is where Vite copies
`publicDir` contents.

```bash
mkdir -p example/_public
git mv example/favicon.ico example/_public/favicon.ico
```

**Step 3: Delete the nested package.json**

`example/package.json` contains only `{"sideEffects": true}`, which was a
Parcel bundling hint. Vite does not use it, and leaving a nested
`package.json` inside the Vite root changes module resolution.

```bash
git rm example/package.json
```

**Step 4: Verify the layout**

Run: `ls example example/_public`
Expected: `example/` lists `App.tsx`, `augment.d.ts`, `components`,
`global.d.ts`, `index.css`, `index.html`, `index.tsx`, `sw.ts`, `util`,
`_public`. `example/_public/` lists `favicon.ico`. Neither
`example/package.json` nor `example/favicon.ico` is present.

Run: `test ! -d demo && echo "demo removed"`
Expected: prints `demo removed`.

**Step 5: Commit**

```bash
git add -A example
git commit -m "refactor: move demo to example and add _public"
```
<!-- END_TASK_2 -->

<!-- START_TASK_3 -->
### Task 3: Point the example entry at Vite idioms

**Verifies:** None (infrastructure; verified operationally in Task 8)

**Files:**
- Modify: `example/index.tsx:5` and `example/index.tsx:7`

**Implementation:**

Two Parcel-specific idioms must change. `process.env.NODE_ENV` is not
defined by Vite; the equivalent is `import.meta.env.PROD`. The service
worker is registered by URL relative to the module under Parcel, whereas
`vite-plugin-pwa` emits `sw.js` at the site root.

Replace lines 5 to 9 of `example/index.tsx` with:

```tsx
if (import.meta.env.PROD) {
    if ('serviceWorker' in navigator) {
        navigator.serviceWorker.register('/sw.js', { type: 'module' })
    }
}
```

Leave the `react` and `react-dom/client` imports on lines 1 and 3 alone.
They are resolved by the aliases configured in Task 4.

**Verification:**

Run: `grep -n "process.env" example/index.tsx`
Expected: no output.

Run: `grep -n "import.meta.env.PROD" example/index.tsx`
Expected: one match.

**Commit:** `refactor: use vite env and sw registration in example`
<!-- END_TASK_3 -->

<!-- END_SUBCOMPONENT_A -->

<!-- START_SUBCOMPONENT_B (tasks 4-6) -->

<!-- START_TASK_4 -->
### Task 4: Create vite.config.js

**Verifies:** dist-build-restructure-design.AC7.1,
dist-build-restructure-design.AC7.2

**Files:**
- Create: `vite.config.js`

**Implementation:**

Create `vite.config.js` at the repository root with exactly this content.

Four things here are load bearing and are explained in the comments: the
`react-dom/client` alias that the Preact preset does not provide, the
`fflate` alias that decouples the example from the library build, the
`emptyOutDir` flag that Vite requires when `outDir` sits outside `root`,
and `import.meta.dirname` in place of `__dirname` because the package is
`"type": "module"`.

```js
import { defineConfig } from 'vite'
import preact from '@preact/preset-vite'
import { VitePWA } from 'vite-plugin-pwa'
import { resolve } from 'path'

const root = import.meta.dirname

export default defineConfig({
    root: 'example',
    publicDir: '_public',

    resolve: {
        // Array form, because the worker swap below matches on a regex.
        alias: [
            // preset-vite aliases react, react-dom, react/jsx-runtime
            // and react-dom/test-utils, but not the react-dom/client
            // subpath that example/index.tsx imports createRoot from.
            // The target is preact/compat/client, not preact/compat:
            // createRoot lives in the client entry only.
            { find: 'react-dom/client', replacement: 'preact/compat/client' },

            // The example imports the library as '../../..', which would
            // resolve through the root package.json and couple this
            // build to the library build. Resolve it to source instead.
            { find: 'fflate', replacement: resolve(root, 'src/index.ts') },

            // Because the line above pulls in src/index.ts, this build
            // sees its `import wk from './node-worker'`. From Phase 2
            // onward that module statically imports node:worker_threads,
            // which must never reach a browser bundle. Apply the same
            // swap the library build performs.
            {
                find: /^\.\/node-worker$/,
                replacement: resolve(root, 'src/worker.ts')
            }
        ]
    },

    plugins: [
        preact(),
        VitePWA({
            strategies: 'injectManifest',
            srcDir: '.',
            filename: 'sw.ts',
            injectRegister: false,
            injectManifest: {
                injectionPoint: 'self.__WB_MANIFEST'
            }
        })
    ],

    build: {
        // outDir is outside root, so Vite refuses to clear it unless
        // emptyOutDir is set explicitly.
        outDir: resolve(root, 'public'),
        emptyOutDir: true
    }
})
```

**Verification:**

Run: `node --input-type=module -e "import('./vite.config.js').then(m => console.log(typeof m.default))"`
Expected: prints `object`.

**Commit:** `build: add vite config for example app`
<!-- END_TASK_4 -->

<!-- START_TASK_5 -->
### Task 5: Redirect the example's library imports to the alias

**Verifies:** dist-build-restructure-design.AC7.1

**Files:**
- Modify: `example/components/code-box/sandbox.ts:1`
- Modify: `example/components/code-box/stream-adapter.ts:1`

**Implementation:**

Both files import the library by walking up to the repository root. Change
them to use the `fflate` alias registered in Task 4, so the example builds
from `src/index.ts` and does not depend on Phase 3 or Phase 4 having run.

In `example/components/code-box/sandbox.ts`, change line 1 from
`import * as fflate from '../../..';` to:

```ts
import * as fflate from 'fflate'
```

In `example/components/code-box/stream-adapter.ts`, change line 1 from
`import { AsyncDeflate } from '../../..';` to:

```ts
import { AsyncDeflate } from 'fflate'
```

**Verification:**

Run: `grep -rn "from '\.\./\.\./\.\.'" example/`
Expected: no output.

Run: `grep -rn "from 'fflate'" example/`
Expected: two matches, one per file above.

**Commit:** `refactor: import library via alias in example`
<!-- END_TASK_5 -->

<!-- START_TASK_6 -->
### Task 6: Convert the service worker to the injected manifest

**Verifies:** dist-build-restructure-design.AC7.3

**Files:**
- Modify: `example/sw.ts:1-9`

**Implementation:**

`vite-plugin-pwa` in `injectManifest` mode replaces the
`self.__WB_MANIFEST` token with an array of `{ url, revision }` entries.
It provides no counterpart to Parcel's `version` export, so the cache name
is derived from the manifest itself: hashing the concatenated revisions
gives a value that changes exactly when the precached set changes, which
is the property the old `version` string provided.

Keep lines 10 onward of `example/sw.ts` unchanged. Replace lines 1 to 9
with:

```ts
/// <reference lib="webworker" />

// Required. Removing the @parcel/service-worker import leaves this file
// with no import or export, which makes it a global script rather than a
// module. In a global script `declare const self` collides with
// lib.dom's `declare var self`, giving TS2451 "Cannot redeclare
// block-scoped variable 'self'" followed by TS2339 on __WB_MANIFEST.
export {}

declare const self:ServiceWorkerGlobalScope & {
    __WB_MANIFEST:{ url:string; revision:string|null }[]
}

const manifest = self.__WB_MANIFEST

// vite-plugin-pwa injects per asset revisions but no build version.
// Derive a stable cache name from the revisions so it changes when,
// and only when, the precached set changes.
const precacheVersion = 'fflate-' + manifest
    .map(e => e.revision ?? e.url)
    .join('|')
    .split('')
    .reduce((h, c) => (((h << 5) - h) + c.charCodeAt(0)) | 0, 0)
    .toString(36)

const precacheFiles = manifest
    .map(e => e.url)
    .filter(u => !/\.(ico)$/.test(u))

const ch = () => caches.open(precacheVersion)
```

Note that the original line 7 declared `const sw = self as unknown as
ServiceWorkerGlobalScope`. The block above declares `self` directly with
the correct type, so every later reference to `sw` in the file must become
`self`. There are four such references, on the original lines 11, 20, 26,
and 30.

**Verification:**

Run: `grep -n "@parcel/service-worker" example/sw.ts`
Expected: no output.

Run: `grep -c "sw\." example/sw.ts || true`
Expected: `0`.

Run: `npx tsc --noEmit --ignoreConfig --lib ES2022,DOM,WebWorker --skipLibCheck example/sw.ts`
Expected: exits 0 with no diagnostics. `--ignoreConfig` is required:
naming a file on the command line while `tsconfig.json` exists otherwise
fails with `TS5112`.

**Commit:** `refactor: precache from injected manifest in service worker`
<!-- END_TASK_6 -->

<!-- END_SUBCOMPONENT_B -->

<!-- START_TASK_7 -->
### Task 7: Repoint cpGHPages at public/ and fix its path handling

**Verifies:** dist-build-restructure-design.AC7.4

**Files:**
- Modify: `scripts/cpGHPages.ts` (whole file)

**Implementation:**

Three defects are fixed together, because the script cannot be verified
without all three. The output directory moves from `dist` to `public`
(the design requires this, since `dist/` now holds the library). Line 15
stats a bare entry name instead of a repository-relative path. Line 24
returns to `master` unconditionally, which is wrong on any other branch.
`__dirname` is also unavailable because the package is `"type": "module"`.

Replace the entire contents of `scripts/cpGHPages.ts` with:

```ts
import { simpleGit } from 'simple-git'
import { resolve, join } from 'path'
import { copyFileSync, readdirSync, statSync, unlinkSync } from 'fs'

const baseDir = resolve(import.meta.dirname, '..')
const to = (...paths:string[]) => join(baseDir, ...paths)
const git = simpleGit()

const branch = (await git.revparse(['--abbrev-ref', 'HEAD'])).trim()
const log = await git.log({ from: 'HEAD~1', to: 'HEAD' })
const hash = log.latest!.hash.slice(0, 7)

await git.checkout('gh-pages')

for (const f of readdirSync(to('.'))) {
    // statSync needs the repository relative path, not the bare entry.
    if (statSync(to(f)).isFile()) unlinkSync(to(f))
}

const files = readdirSync(to('public'))
for (const f of files) {
    copyFileSync(to('public', f), to(f))
}

await git.add(files)
await git.commit('Build demo from ' + hash)

// Return to the branch we started on, not an assumed master.
await git.checkout(branch)
```

**Verification:**

This script mutates git branches, so do not execute it as a verification
step. Verify statically.

Run: `grep -n "'dist'" scripts/cpGHPages.ts`
Expected: no output.

Run: `grep -n "__dirname" scripts/cpGHPages.ts`
Expected: no output.

Run: `npx tsc --noEmit --ignoreConfig --module ES2022 --target ES2022 --moduleResolution Bundler --types node --skipLibCheck scripts/cpGHPages.ts`
Expected: exits 0 with no diagnostics.

**Commit:** `fix: copy example build from public in cpGHPages`
<!-- END_TASK_7 -->

<!-- START_TASK_8 -->
### Task 8: Wire the example scripts and ignore public/

**Verifies:** dist-build-restructure-design.AC7.1,
dist-build-restructure-design.AC7.2,
dist-build-restructure-design.AC7.3

**Files:**
- Modify: `package.json` (scripts)
- Modify: `.gitignore`

**Implementation:**

Set these three script entries in `package.json`. Leave every other script
alone; Phase 4 rewrites the full script block.

```json
"start": "vite",
"build-example": "vite build --base=\"/fflate\"",
"gh-pages": "tsx scripts/cpGHPages.ts"
```

The design replaces the template's `/repo-name` placeholder with
`/fflate`. Remove the old `build:demo` script, whose Parcel invocation no
longer has a toolchain, and the `SC=cpGHPages npm run script` indirection
that `gh-pages` now supersedes.

Append `public/` to `.gitignore`, which does not currently list it:

```
public/
```

**Verification:**

Run: `npm run build-example`
Expected: exits 0. Vite reports the example root and writes to the
`public` directory. No unresolved-import errors for `react`,
`react-dom/client`, or `fflate`.

Run: `ls public`
Expected: contains `index.html`, an `assets/` directory, `favicon.ico`,
and `sw.js`.

Run: `grep -c "__WB_MANIFEST" public/sw.js || true`
Expected: `0`. The token must have been replaced by the real manifest
array at build time; a remaining token means injection did not run.

Run: `grep -o "precacheVersion" public/sw.js | head -1`
Expected: either a match or no output. Both are acceptable, since the
identifier may be minified. This is a smoke check only.

Run: `git status --porcelain public`
Expected: no output, confirming `public/` is ignored.

**Commit:** `build: add vite example scripts and ignore public`
<!-- END_TASK_8 -->

---

## Phase 1 completion criteria

- `npm run build-example` exits 0 and populates `public/`.
- `npm run start` serves the example without unresolved imports.
- No file under `example/` references `@parcel/service-worker`,
  `process.env`, or `'../../..'`.
- `parcel`, `@parcel/service-worker`, `@types/react`, and
  `@types/react-dom` are absent from `package.json`.
- `scripts/cpGHPages.ts` reads `public/` and contains no `__dirname`.
- `public/` is git-ignored.

The repository build (`npm run build`) is still the old half-merged
pipeline at the end of this phase. Phase 2 replaces the TypeScript
configuration and Phase 3 replaces the build pipeline itself.
