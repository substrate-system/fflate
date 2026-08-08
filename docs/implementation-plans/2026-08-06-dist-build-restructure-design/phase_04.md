# Dist Build Restructure Implementation Plan

## Phase 4: package.json entry points and scripts

**Goal:** Point every published entry at `dist/`, drop CommonJS from the
`exports` map, and replace the half-merged script block with the
design's script table.

**Architecture:** A conditional `exports` map selects the node bundle
under the `node` condition and the browser bundle otherwise, with
explicit `./node`, `./browser`, `./min`, and `./umd` subpaths. `files` is
narrowed to `dist`.

**Tech Stack:** npm conditional exports, Node 25 resolution.

**Scope:** Phase 4 of 7.

**Codebase verified:** 2026-08-06

**Depends on:** Phase 3, which produces the files this map points at.
Verification here loads the real artifacts.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC4: Published entry points resolve

- **dist-build-restructure-design.AC4.1 Success:** Importing the package
  under the `node` condition resolves to `dist/node/index.js`, and its
  types to `dist/node/index.d.ts`.
- **dist-build-restructure-design.AC4.2 Success:** Importing under any
  other condition resolves to `dist/browser/index.js`, and its types to
  `dist/browser/index.d.ts`.
- **dist-build-restructure-design.AC4.3 Success:** The `./node`,
  `./browser`, `./min`, `./umd`, and `./package.json` subpaths all
  resolve to existing files.
- **dist-build-restructure-design.AC4.4 Success:** `npm pack` includes
  only `dist/` and the always-included metadata files.
- **dist-build-restructure-design.AC4.5 Failure:** No `require`
  condition is advertised anywhere in the `exports` map.

---

## Investigation findings

Verified against the working tree on 2026-08-06.

Confirmed as the design describes:

- The current `exports` map advertises `require` conditions pointing at
  `./lib/node.cjs`, `./lib/browser.cjs`, and `./lib/index.cjs`, and
  `import` conditions pointing at `./esm/index.mjs` and
  `./esm/browser.js`. None of these paths is produced by the Phase 3
  build.
- `main` is `dist/index.js`, which the Phase 3 build does not produce
  either; the design changes it to `./dist/node/index.js`.
- `files` is `["./dist/*", "./esm/*", "./lib/*", "./umd/*"]`.
- Every half-merged script named in the design is present: `build-cjs`,
  `build-esm`, `build-esm:min`, `build-cjs:min`, `build:lib`,
  `build:umd`, `build:rewrite`, `build:demo`, `script`, and `//build`.

Corrections and additions the design did not cover:

1. **The design's script table renames `build:docs` to `build-docs`
   without saying so.** The current script is `build:docs`. The design's
   table lists `build-docs` and marks it "(unchanged)", referring to the
   typedoc invocation rather than the script name. Task 2 renames the key
   and preserves the command verbatim.

2. **`test:valid` is not in the design's table, and `prepublishOnly`
   already exists.** `test:valid` (`tsx test/0-valid.ts`) is a
   convenience entry the design mentions neither way; it is kept, since
   the design's list of removals is explicit and does not include it.
   There is no `prepack` script to remove: `package.json` already sets
   `"prepublishOnly": "npm run build"`, exactly what the design's table
   specifies, so that row is satisfied by an existing field rather than a
   change.

3. **`type: "module"` is already set** and must stay. It is what makes
   the ESM `.js` output in `dist/` resolvable.

4. **`sideEffects` is absent.** The design does not ask for it and it is
   not added. Adding `"sideEffects": false` would be wrong here anyway:
   `src/worker.ts` and the async APIs rely on module-scope state.

5. **`./umd` needs a sibling `dist/umd/package.json`, and no `types`
   entry.** Added during execution on 2026-08-06. The root package is
   `"type": "module"`, so node parses `dist/umd/fflate.js` as ESM, where
   the UMD wrapper's `this` is `undefined`: importing it throws
   `TypeError: Cannot set properties of undefined (setting 'fflate')` and
   requiring it returns an empty object. Phase 3 now emits
   `dist/umd/package.json` containing `{ "type": "commonjs" }`, which
   makes the artifact load correctly under both `require()` and
   `import()` while keeping the filename `unpkg` points at.

   `./umd` deliberately carries NO `types` entry, unlike every other
   subpath. Adding one was tried and reverted. The only declaration
   available is `dist/browser/index.d.ts`, which sits under the root
   `"type": "module"` and is therefore ESM-shaped, while the artifact it
   would describe is CommonJS. That combination type checks code that
   fails at runtime: under `moduleResolution: nodenext`,
   `import { zlibSync } from '.../umd'` compiles but throws
   `SyntaxError: Named export 'zlibSync' not found`, and
   `import * as umd` compiles but yields an object whose members are all
   `undefined` -- turning a loud failure into a silent one. Only the
   default import works. Omitting `types` makes TypeScript emit `TS7016`
   (`Could not find a declaration file`) under `noImplicitAny`; with that
   flag off the import degrades to `any` rather than to a wrong type,
   which is still the better of the two failure modes. Bundlers are
   unaffected either way, since their CJS interop resolves the members at
   runtime.

6. **`build-docs` must not clean its output directory.** Added during
   execution on 2026-08-06. TypeDoc defaults `cleanOutputDir` to true and
   the command writes `--out docs/`, which is also where
   `docs/implementation-plans/` lives, so running the script would delete
   this plan. `--cleanOutputDir false` is added to the command.

7. **`build-example` keeps the trailing slash in its base.** Corrected
   during execution on 2026-08-06; the Task 2 script block has been
   amended from `--base="/fflate"` to `--base="/fflate/"`. Vite does not
   normalise `import.meta.env.BASE_URL` to a trailing slash, so with
   `/fflate` the service worker registration in `example/index.tsx`
   concatenates to `/fflatesw.js` and 404s on the deployed site, leaving
   nothing precached. Asset URLs in `index.html` stay correct either way,
   because Vite path-joins those separately, which is what makes the bug
   easy to miss. Phase 1 fixed this; restating the script block without
   the slash would silently revert it.

---

<!-- START_TASK_1 -->
### Task 1: Rewrite the entry point fields and exports map

**Verifies:** dist-build-restructure-design.AC4.1,
dist-build-restructure-design.AC4.2,
dist-build-restructure-design.AC4.3,
dist-build-restructure-design.AC4.5

**Files:**
- Modify: `package.json` (`main`, `module`, `types`, `unpkg`, `files`,
  `exports`)

**Implementation:**

Replace the `main` field and the entire `exports` object, and add
`module`, `types`, and `unpkg`. Set `files` to `["dist"]`.

Within a condition object, `types` must come first so TypeScript sees it
before falling through to `default`. There is no `require` condition
anywhere: dropping CommonJS is the point, and an explicit `require` entry
pointing at ESM would only restate what `default` already covers.

Note that omitting `require` does not make CommonJS consumers fail.
`default` matches the `require` condition too, so on Node 20.19+ or
22.12+ a CJS caller gets the ESM node build through require-ESM, which
was measured to work; older Node throws a clear `ERR_REQUIRE_ESM`.
AC4.5 is about not advertising a `require` condition, not about blocking
CommonJS.

```json
"main": "./dist/node/index.js",
"module": "./dist/browser/index.js",
"types": "./dist/browser/index.d.ts",
"unpkg": "./dist/umd/fflate.js",
"files": ["dist"],
"exports": {
  ".": {
    "node": {
      "types": "./dist/node/index.d.ts",
      "default": "./dist/node/index.js"
    },
    "default": {
      "types": "./dist/browser/index.d.ts",
      "default": "./dist/browser/index.js"
    }
  },
  "./node": {
    "types": "./dist/node/index.d.ts",
    "default": "./dist/node/index.js"
  },
  "./browser": {
    "types": "./dist/browser/index.d.ts",
    "default": "./dist/browser/index.js"
  },
  "./min": {
    "types": "./dist/browser/index.d.ts",
    "default": "./dist/browser/index.min.js"
  },
  "./umd": "./dist/umd/fflate.js",
  "./package.json": "./package.json"
}
```

Leave `"type": "module"` in place.

**Verification:**

The build from Phase 3 must be present. Run `npx tsx scripts/build.ts`
first if `dist/` is absent.

Run: `node -e "console.log(require.resolve ? 'ok' : '')" && node --input-type=module -e "import('@substrate-system/fflate').then(m => console.log('node cond:', typeof m.deflateSync))"`
Expected: the import resolves. If it fails with `ERR_MODULE_NOT_FOUND`,
the package is not linked into its own `node_modules`; use the explicit
resolution check below instead.

Run: `node -e "console.log(require('module').createRequire(process.cwd()+'/x.js').resolve('./dist/node/index.js'))"`
Expected: prints an absolute path ending in `dist/node/index.js`.

Run these five and expect every one to print an existing path:
```bash
for s in dist/node/index.js dist/node/index.d.ts dist/browser/index.js \
         dist/browser/index.d.ts dist/browser/index.min.js dist/umd/fflate.js; do
  test -f "$s" && echo "ok $s" || echo "MISSING $s"
done
```
Expected: six `ok` lines, no `MISSING`.

Run: `node -e "const e=require('./package.json').exports; console.log(JSON.stringify(e).includes('require') ? 'HAS REQUIRE' : 'no require condition')"`
Expected: prints `no require condition`.

Run: `node -e "const e=require('./package.json').exports; for (const k of ['.','./node','./browser','./min','./umd','./package.json']) if(!(k in e)) console.log('MISSING KEY',k)"`
Expected: no output.

**Commit:** `feat: point exports map at dist and drop commonjs`
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Replace the script block

**Verifies:** None (infrastructure; the scripts are exercised by Phases
5, 6, and 7)

**Files:**
- Modify: `package.json` (`scripts`)

**Implementation:**

Replace the entire `scripts` object with the following. `start`,
`build-example`, and `gh-pages` were already set in Phase 1 Task 8 and
are restated here so the final block is unambiguous.

```json
"scripts": {
  "build": "tsx scripts/build.ts",
  "build-example": "vite build --base=\"/fflate/\"",
  "build-docs": "typedoc --plugin typedoc-plugin-markdown --hideBreadcrumbs --readme none --disableSources --excludePrivate --excludeProtected --expandParameters --githubPages false --cleanOutputDir false --out docs/ src/index.ts",
  "gh-pages": "tsx scripts/cpGHPages.ts",
  "start": "vite",
  "test": "npm run build && npm run test:node && npm run test:browser",
  "test:node": "tsx test/index.ts",
  "test:valid": "tsx test/0-valid.ts",
  "test:browser": "esbuild test/browser/index.ts --bundle | tapout --timeout 30000",
  "toc": "markdown-toc --maxdepth 3 -i README.md",
  "version": "npm run toc && auto-changelog -p --template keepachangelog --breaking-pattern 'BREAKING CHANGE:' && git add CHANGELOG.md README.md",
  "postversion": "git push --follow-tags && npm publish",
  "prepublishOnly": "npm run build"
}
```

Every one of `build-cjs`, `build-esm`, `build-esm:min`, `build-cjs:min`,
`build:lib`, `build:umd`, `build:rewrite`, `build:demo`, `build:docs`,
`script`, and `//build` is gone. (`prepack` is in the removal check below
only as a guard; it does not currently exist.) The generic `script` runner
(`tsx scripts/$SC.ts`) goes with them; its only remaining caller was
`cpGHPages.ts`, which now has the named `gh-pages` entry.

`--timeout 30000` on `test:browser` is load-bearing, not a slow-machine
allowance. tapout ends a run after
`max(500, min(3000, floor(timeout * 0.2)))` ms with no console output and
resets that timer on every line printed. At tapout's default timeout of
5000 that window is 1000 ms; the browser suites are silent while
polling, so a stalled test goes quiet, auto-finish fires, and every
later test is dropped -- including the entire minified pass, which is
the one thing the suite exists to run. Raising the timeout puts the
window at its 3000 ms cap, above `withTimeout`'s 2000 ms default in
`test/browser/util.ts`, so a stall rejects first: the catch prints
`not ok`, that output resets the window, and the run continues.

The two values are a pair; changing either alone reintroduces the
truncation, and the failure presents as a PASS, so nothing else would
catch it. Phase 6 Task 6 adds `test/browser/harness.ts`, which reads
the `--timeout` back out of this script and asserts it leaves room for
`withTimeout`'s default. Remove the timeout flag, or lower it below
12500, and that test goes `not ok` -- 12500 is where
`floor(t * 0.2)` stops clearing `withTimeout`'s default plus its
margin. Raising it is safe: the window is capped at 3000 ms, so a
larger `--timeout` only buys whole-run headroom.

Note that `test:browser` does not work yet. It requires
`test/browser/index.ts` and the `tapout` binary, both of which arrive in
Phase 6. `npm test` therefore fails at its third step until Phase 6
completes; `npm run test:node` is the working gate in the meantime.

**Verification:**

Run: `node -e "const s=require('./package.json').scripts; const gone=['build-cjs','build-esm','build-esm:min','build-cjs:min','build:lib','build:umd','build:rewrite','build:demo','build:docs','script','//build','prepack']; const left=gone.filter(k=>k in s); console.log(left.length?'STILL PRESENT: '+left.join(', '):'all removed')"`
Expected: prints `all removed`.

Run: `npm run build`
Expected: exits 0 and produces the twelve files from Phase 3 Task 8.

Run: `npm run test:node`
Expected: fails. `test/util.ts` still points at `lib/index.cjs`, which no
longer exists. This is expected and is fixed in Phase 5. Record the
failure and continue.

**Commit:** `build: replace script block with dist pipeline scripts`
<!-- END_TASK_2 -->

<!-- START_TASK_3 -->
### Task 3: Trim the ignore files and verify the package contents

**Verifies:** dist-build-restructure-design.AC4.4

**Files:**
- Modify: `.gitignore`
- Modify: `.npmignore`

**Implementation:**

Remove the `lib/`, `esm/`, and `umd/` entries from `.gitignore`, keeping
`dist/`, `public/`, `node_modules/`, `.DS_STORE`, and the Rust entries.
The resulting file:

```
node_modules/
dist/
public/
.DS_STORE
# Rust version - available when ready
rs/*/target/
rs/*/Cargo.lock
```

The duplicate `.DS_STORE` line and the `.parcel-cache` entry go too; the
former was a duplicate and the latter belonged to the toolchain Phase 1
removed.

Replace `.npmignore` with:

```
node_modules
.env
.DS_Store
```

The `!dist/`, `!esm/`, `!umd/`, and `!lib/` negations are removed. They
are inert now: `files: ["dist"]` in `package.json` takes precedence over
`.npmignore` for deciding what is published.

**Verification:**

Run: `npm pack --dry-run 2>&1 | grep -E "^npm notice" | grep -vE "package size|unpacked size|shasum|integrity|total files|filename|=== |name:|version:"`
Expected: every listed file is under `dist/`, plus `package.json`,
`README.md`, and `LICENSE`. No `src/`, `test/`, `example/`, `scripts/`,
`lib/`, `esm/`, or `umd/` entries.

Run:
```bash
npm pack --dry-run --json 2>/dev/null | node -e "
let s = '';
process.stdin.on('data', d => s += d).on('end', () => {
  const files = JSON.parse(s)[0].files.map(f => f.path);
  const allowed = /^(dist\/|package\.json$|README\.md$|LICENSE$)/;
  const bad = files.filter(f => !allowed.test(f));
  console.log(bad.length ? 'UNEXPECTED: ' + bad.join(', ') : 'contents ok');
});
"
```
Expected: prints `contents ok`.

Run: `git check-ignore -q public/ && git check-ignore -q dist/ && echo "both ignored"`
Expected: prints `both ignored`.

**Commit:** `chore: trim ignore files for dist only publishing`
<!-- END_TASK_3 -->

---

## Phase 4 completion criteria

- `npm run build` exits 0.
- Every path named in the `exports` map exists on disk after a build.
- The `exports` map contains no `require` condition.
- `npm pack --dry-run` lists only `dist/` plus package metadata.
- All twelve superseded scripts are gone from `package.json`.
- `.gitignore` no longer lists `lib/`, `esm/`, or `umd/`, and does list
  `dist/` and `public/`.

`npm test` fails at `test:node` until Phase 5 and at `test:browser` until
Phase 6. That is the expected intermediate state.
