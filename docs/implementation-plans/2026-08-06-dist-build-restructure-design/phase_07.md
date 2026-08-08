# Dist Build Restructure Implementation Plan

## Phase 7: Legacy removal and final verification

**Goal:** Delete the stale build output left by the old pipeline, confirm
nothing in the repository still references it, and verify the published
package end to end.

**Architecture:** No new code. This phase removes directories, sweeps for
dangling references, and runs the acceptance checks that span more than
one phase.

**Tech Stack:** npm pack, node resolution.

**Scope:** Phase 7 of 7.

**Codebase verified:** 2026-08-06

**Depends on:** Phases 1 through 6, all complete and passing.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC8: Type checking passes and obsolete configuration is gone

- **dist-build-restructure-design.AC8.2 Success:** No file in the
  repository references `scripts/rewriteBuilds.ts`,
  `scripts/buildUMD.ts`, the `lib/`, `esm/`, or top level `umd/`
  directories, or the `SC=` script runner.

This phase also re-verifies, as a combined gate, criteria already
implemented earlier:
`dist-build-restructure-design.AC1.1` through
`dist-build-restructure-design.AC1.5`,
`dist-build-restructure-design.AC4.1` through
`dist-build-restructure-design.AC4.4`,
`dist-build-restructure-design.AC5.1`,
`dist-build-restructure-design.AC6.1`, and
`dist-build-restructure-design.AC6.2`.

---

## Investigation findings

Verified on 2026-08-06.

1. **`lib/`, `esm/`, `umd/`, and `dist/` all exist in the working tree**
   as output from the previous half-merged build. All four are
   git-ignored today, so they are untracked local artifacts rather than
   committed files. Removing them is a filesystem operation, not a
   commit.

2. **`umd/` is git-ignored but was a published directory.** The old
   `files` array listed `./umd/*`. Phase 4 narrowed `files` to `["dist"]`,
   so nothing outside `dist/` is published regardless of what remains on
   disk.

3. **The `SC=` indirection had exactly one remaining caller**,
   `cpGHPages.ts`, which Phase 1 replaced with the named `gh-pages`
   script. Phase 4 removed the `script` entry itself.

4. **`docs/` is tracked typedoc output** and is explicitly out of scope
   per the design. The `build-docs` script and its typedoc invocation are
   unchanged. Do not regenerate `docs/` as part of this phase; a
   regeneration would produce a large unrelated diff.

5. **The README documents seven entry points that this work removes**,
   and no phase updated it. Measured: `README.md:69` shows
   `require('fflate')`; `:81` and `:82` link CDN paths for the old
   package name, `:82` specifically to `umd/index.js`; `:87` and `:97`
   use skypack with the old name; `:96` references `lib/index.d.ts`; and
   `:106` and `:108` import `fflate/esm/browser.js` and `fflate/esm`.
   The design's Risk 4 calls dropping CommonJS "breaking for any consumer
   using `require()`" but assigns no documentation work. Task 2B covers
   both.

---

<!-- START_TASK_1 -->
### Task 1: Remove stale build output

**Verifies:** dist-build-restructure-design.AC1.4

**Files:**
- Delete: `lib/`, `esm/`, `umd/` (untracked build output)

**Implementation:**

These are the outputs of the pipeline that Phase 3 replaced. All three
are git-ignored, so this is a working tree cleanup with no commit.

Confirm they are untracked before deleting, so nothing tracked is lost:

```bash
git ls-files lib esm umd
```

Expected: no output. If any path is listed, stop: it is tracked, and
deleting it is a committed change that must be reviewed rather than
swept.

Then remove them:

```bash
rm -rf lib esm umd
```

**Verification:**

Run: `ls lib esm umd 2>&1`
Expected: three "No such file or directory" errors.

Run: `npm run build && find dist -type d -maxdepth 1 | sort`
Expected: `dist`, `dist/browser`, `dist/node`, `dist/umd`. The build must
not recreate the top level `lib`, `esm`, or `umd`.

Run: `ls lib esm umd 2>&1`
Expected: still three errors, after the build.

**Commit:** No commit. These paths are git-ignored; `git status` should
be unchanged by this task.
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Sweep for dangling references

**Verifies:** dist-build-restructure-design.AC8.2

**Files:**
- Modify: any file the sweep finds. Expected: none.

**Implementation:**

Search the tracked sources for references to anything the restructure
removed. Every command below should produce no output. Where one does,
fix the reference.

`docs/` is excluded throughout: it is generated typedoc output, out of
scope per the design, and will contain historical references.

`specs/` is excluded as well. The design document is tracked and
describes the very things being removed, so five of the six commands hit
it. Following the "fix the reference" instruction there would mean
editing the source of truth this plan is measured against.

`README.md` is excluded too, and Task 2B rewrites it immediately after.
Without that exclusion this sweep fails on `README.md:106`, a reference
the plan already knows about and has scheduled a fix for. Run Task 2B
before treating any README hit as outstanding.

**Verification:**

Run each and expect no output:

```bash
git grep -n "rewriteBuilds\|buildUMD" -- . ':!docs' ':!specs'
git grep -n "SC=" -- . ':!docs' ':!specs'
git grep -nP "\blib/index\.cjs\b|\besm/index\.mjs\b|\besm/browser\.js\b" -- . ':!docs' ':!specs' ':!README.md'
git grep -nP "tsconfig\.(esm|demo)\.json" -- . ':!docs' ':!specs'
git grep -n "@parcel/service-worker\|parcel " -- . ':!docs' ':!specs' ':!package-lock.json'
git grep -n "build:lib\|build:umd\|build:rewrite\|build:demo\|build-cjs\|build-esm" -- . ':!docs' ':!specs'
```

Run: `git grep -nP "'\.\./\.\./\.\.'" -- example`
Expected: no output. The example imports via the `fflate` alias.

Run: `git grep -n "from 'react'" -- example | wc -l`
Expected: greater than `0`. The React specifiers are intentionally kept
and resolved through `preact/compat`; this confirms Phase 1 did not
rewrite them by mistake.

Run: `ls scripts`
Expected: exactly `build.ts` and `cpGHPages.ts`.

**Commit:** `chore: remove dangling references to old build` (only if the
sweep found something; otherwise no commit)
<!-- END_TASK_2 -->

<!-- START_TASK_2B -->
### Task 2B: Update the README for the new entry points

**Verifies:** dist-build-restructure-design.AC8.2

**Files:**
- Modify: `README.md` lines 69, 81, 82, 87, 96, 97, 106, 108

**Implementation:**

The README documents entry points that this restructure removes. Every
one of these is measured against the current file and is wrong after
Phase 4:

| line | current text | problem |
| --- | --- | --- |
| 69 | `const fflate = require('fflate');` | CommonJS is dropped. There is no `require` condition. |
| 82 | `<script src="https://cdn.jsdelivr.net/npm/fflate@0.8.3/umd/index.js">` | UMD moved to `dist/umd/fflate.js`, and the package is `@substrate-system/fflate`. |
| 96 | `@deno-types=".../fflate@0.8.3/lib/index.d.ts"` | `lib/` no longer exists; types are at `dist/browser/index.d.ts`. |
| 106 | `import * as fflate from 'fflate/esm/browser.js';` | `esm/` no longer exists; the subpath is `@substrate-system/fflate/browser`. |
| 108 | `import * as fflate from 'fflate/esm';` | `esm/` no longer exists. |
| 81 | `<script src="https://unpkg.com/fflate@0.8.3">` | Wrong package name; `unpkg` now resolves to `dist/umd/fflate.js`. |
| 87, 97 | skypack `fflate@0.8.3` | Wrong package name. |

Update each to the new layout. Replace the line 69 CommonJS example with
the ESM equivalent rather than deleting it, so the README still shows a
Node usage example.

Add a short note to the README recording that dropping CommonJS is a
breaking change for `require()` consumers, and that the UMD artifact at
`dist/umd/fflate.js` remains the fallback for script-tag and
legacy-bundler use. This is the design's Risk 4, which otherwise has no
home in the plan.

Do not regenerate `docs/`. It is typedoc output and out of scope.

**Verification:**

Run: `grep -nE "require\('fflate'\)|umd/index\.js|lib/index\.d\.ts|fflate/esm|esm/browser\.js|esm/index\.mjs" README.md || true`
Expected: no output. Note `fflate/esm` rather than only `esm/browser.js`:
line 108 imports the bare `fflate/esm` subpath and the narrower pattern
misses it.

Run: `grep -c "fflate@0.8.3" README.md || true`
Expected: `0`. The pinned old-package CDN references on lines 81, 82, 87,
96, and 97 must all be updated.

Run: `grep -c "dist/umd/fflate.js" README.md || true`
Expected: at least `1`.

Run: `npm run toc`
Expected: exits 0. If the edits changed any heading, this keeps the table
of contents in sync.

**Commit:** `docs: update readme for dist layout and esm only`
<!-- END_TASK_2B -->

<!-- START_TASK_3 -->
### Task 3: Verify the published package

**Verifies:** dist-build-restructure-design.AC4.1,
dist-build-restructure-design.AC4.2,
dist-build-restructure-design.AC4.3,
dist-build-restructure-design.AC4.4

**Files:**
- No files changed. This task is a gate.

**Implementation:**

Pack the package and install it into a scratch directory, so resolution
is tested the way a consumer experiences it rather than through the
working tree.

```bash
npm run build
npm pack
# note the emitted tarball name, e.g. substrate-system-fflate-0.8.17.tgz
mkdir -p /private/tmp/claude-501/-Users-nick-code-fflate/e8c45c98-8ea8-4895-8a65-48cc175c6ce5/scratchpad/fflate-consume && cd /private/tmp/claude-501/-Users-nick-code-fflate/e8c45c98-8ea8-4895-8a65-48cc175c6ce5/scratchpad/fflate-consume
npm init -y >/dev/null
npm install /Users/nick/code/fflate/substrate-system-fflate-*.tgz
```

**Verification:**

From `/private/tmp/claude-501/-Users-nick-code-fflate/e8c45c98-8ea8-4895-8a65-48cc175c6ce5/scratchpad/fflate-consume`, run each check.

**Write each of these to a real file and run it. Do NOT use
`node --input-type=module -e`.** That form puts `--input-type=module`
into `process.execArgv`; worker threads inherit `execArgv` and then
parse their eval'd payload as ESM, so fflate's sloppy-mode implicit
globals throw `ReferenceError: u8 is not defined`. The async check
below spawns a worker and would fail for that reason alone, reporting
a defect that is not there.

Node condition resolves to the node build and round-trips. `sync.mjs`:
```js
import * as f from '@substrate-system/fflate'
const d = f.deflateSync(new TextEncoder().encode('hello hello hello'))
console.log('sync roundtrip:', new TextDecoder().decode(f.inflateSync(d)))
```
Run `node sync.mjs`. Expected: `sync roundtrip: hello hello hello`,
exit 0.

The async path works, which is what spawns a worker. `async.mjs`:
```js
import * as f from '@substrate-system/fflate'
const src = new TextEncoder().encode('async round trip payload')
f.deflate(src, (e, d) => {
  if (e) throw e
  f.inflate(d, (e2, r) => {
    if (e2) throw e2
    console.log('async roundtrip:', new TextDecoder().decode(r))
  })
})
```
Run `node async.mjs`. Expected: `async roundtrip: async round trip
payload`, exit 0.

Subpath exports resolve. `subpaths.mjs`:
```js
const subs = ['', '/node', '/browser', '/min', '/umd']
for (const s of subs) {
  console.log('ok', await import.meta.resolve('@substrate-system/fflate' + s))
}
```
Run `node subpaths.mjs`. Expected: five `ok` lines -- the bare
specifier plus the four subpaths. The bare one and `/node` both resolve
to `dist/node/index.js` under the `node` condition. Any rejection means
an `exports` entry points at a file that was not packed.

CommonJS is gone from the `exports` map. Note carefully what this does
and does not mean.

Do **not** assert that `require()` throws. Node has supported
`require()` of ESM unflagged since 22.12, and this tree runs 25.8.2.
Measured against the exact `exports` map from Phase 4: `require()`
resolves conditions `["require", "node", ...]`, matches the `node` key,
falls through to `default`, and loads `dist/node/index.js` as ESM
successfully. The bundle has no top-level await, so nothing prevents it.

What the design actually requires is that no `require` condition is
advertised, which Phase 4 Task 1 asserts directly. Re-assert it here
against the packed artifact rather than against the working tree:

`norequire.cjs`:
```js
const e = require('@substrate-system/fflate/package.json').exports
console.log(JSON.stringify(e).includes('"require"')
  ? 'UNEXPECTED: require condition present'
  : 'no require condition')
```
Run `node norequire.cjs`. Expected: `no require condition`. A `.cjs`
file is safe here: it only reads the JSON, and spawns no worker.

That `require()` happens to work anyway on modern Node is a property of
the runtime, not of this package, and it is a compatibility bonus rather
than something to rely on or advertise. The design's Risk 4 stands:
consumers on Node below 22.12 lose `require()` support.

Only `dist/` was published:
```bash
find node_modules/@substrate-system/fflate -type d -maxdepth 1 | sort
```
Expected: the package root and `dist` only. No `src`, `test`, `example`,
`scripts`, `lib`, `esm`, or top level `umd`. At the file level the
tarball also carries `LICENSE`, `README.md` and `package.json`, which
npm includes regardless of the `files` array.

Types resolve:
```bash
ls node_modules/@substrate-system/fflate/dist/browser/index.d.ts \
   node_modules/@substrate-system/fflate/dist/node/index.d.ts
```
Expected: both exist.

Clean up:
```bash
cd /Users/nick/code/fflate && rm -rf /private/tmp/claude-501/-Users-nick-code-fflate/e8c45c98-8ea8-4895-8a65-48cc175c6ce5/scratchpad/fflate-consume substrate-system-fflate-*.tgz
```

**Commit:** No commit. This task only verifies.
<!-- END_TASK_3 -->

<!-- START_TASK_4 -->
### Task 4: Full acceptance run

**Verifies:** dist-build-restructure-design.AC1.1 through
dist-build-restructure-design.AC1.5,
dist-build-restructure-design.AC5.1,
dist-build-restructure-design.AC6.1,
dist-build-restructure-design.AC6.2,
dist-build-restructure-design.AC8.1

**Files:**
- No files changed. This task is the final gate.

**Implementation:**

Run everything from a clean state.

**Verification:**

Run: `rm -rf dist public && npm run build && npm test; echo "exit=$?"`
Expected: `exit=0`. Build, node suite, and browser suite all pass.

Run: `find dist -type f | sort`
Expected exactly:
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

`dist/umd/package.json` is `{ "type": "commonjs" }` and is
load-bearing, not a stray. The root package.json is `"type": "module"`,
so without it node parses `fflate.js` as ESM, where the UMD wrapper's
`this` is undefined and loading throws "Cannot set properties of
undefined". Scoping the directory back to commonjs makes the artifact
loadable by both `require()` and `import()` while keeping the filename
`unpkg` points at. Added in `292f07f`; this expected list predated it.

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep -cE "error TS" || true`
Expected: `0`.

Run: `npm run build-example && ls public/index.html && echo "example ok"`
Expected: prints the path and `example ok`.

Run: `git status --porcelain`
Expected: no entries for `dist/`, `public/`, `lib/`, `esm/`, or `umd/`.
Only intentional source changes appear.

Run: `npm run test:browser 2>&1 | grep -c "^not ok"`
Expected: `0`.

Run: `npm run test:browser 2>&1 | grep -c "^# min >"` compared against
`npm run test:browser 2>&1 | grep -c "^# plain >"`
Expected: both equal (both `103`). This is the final confirmation that
both the plain and minified passes completed and the `wcln` hazard did
not fire against the minified bundle.

**Commit:** `chore: complete dist build restructure`
<!-- END_TASK_4 -->

---

## Phase 7 completion criteria

- `lib/`, `esm/`, and top level `umd/` do not exist and are not
  recreated by a build.
- No tracked file outside `docs/` references the removed scripts,
  directories, or the `SC=` runner.
- A packed and installed tarball resolves under the `node` and default
  conditions and all four subpaths, round-trips data synchronously and
  asynchronously, and advertises no `require` condition in its `exports`
  map. It does not need to reject `require()`; see Task 3.
- `npm pack` publishes only `dist/`.
- `npm test` exits 0 from a clean tree.
- `npm run build-example` produces `public/index.html`.
- `README.md` documents the `dist/` entry points and records the
  CommonJS removal as breaking.

## Release path -- do NOT use `npm version`

Added during the Phase 7 review. The `version` field was bumped
`0.8.17` to `0.9.0` by editing `package.json` directly, because this
repo's `version` script chains
`toc && auto-changelog && git add CHANGELOG.md README.md` and firing it
mid-review would have staged files behind the reviewer's back.

The consequence is that `npm version 0.9.0` now refuses: the field
already reads `0.9.0`. Both halves of what that script would have done
have been performed by hand and committed, so the release is:

```bash
git tag -a v0.9.0 -m 'v0.9.0'
git push --follow-tags
npm publish
```

**The `-a` is load-bearing.** `git push --follow-tags` pushes only
ANNOTATED tags. A lightweight `git tag v0.9.0` is skipped silently:
measured against a real bare remote, the push printed
`Everything up-to-date` and `git ls-remote --tags origin` returned
nothing, with no error and no non-zero exit. The operator would then
publish 0.9.0 to the registry with no tag on GitHub, which turns the
`compare/v0.8.17...v0.9.0` link in `CHANGELOG.md` into a 404 and
leaves the release uncuttable. `npm version` creates annotated tags,
which is why this fork's own releases have them: `git cat-file -t
v0.8.17` and `v0.8.16` both return `tag`. Do not read that as true of
every tag here -- 13 of the 23 predate the fork and are lightweight.
The hand-rolled path has to reproduce the annotation deliberately.
`git push origin v0.9.0` as a separate command works too, and does
push a lightweight tag.

The bump was required, not cosmetic. `@substrate-system/fflate@0.8.17`
is already on the registry carrying the pre-restructure layout, so
without it `npm publish` is rejected for an existing version and every
`dist/` path the README documents stays unreachable. Measured: the
published 0.8.17 contains `esm/`, `lib/` and `umd/index.js`, has no
`unpkg` field, and its `exports` map still has `require` conditions.

Two related facts, both measured:

- `CHANGELOG.md` was regenerated with
  `npx auto-changelog -p --template keepachangelog --breaking-pattern
  'BREAKING CHANGE:'`. The `-p` flag reads the version from
  `package.json`, so no tag is needed first.
- No commit between `v0.8.17` and the review carried the
  `BREAKING CHANGE:` trailer that `--breaking-pattern` keys on, so the
  changelog would have presented 0.9.0 as an ordinary release -- with
  no mention of the CommonJS removal in the one place a `require()`
  consumer would look. Empty commit `770c5c4` carries the trailer.
- The keepachangelog template renders the commit SUBJECT, not the
  trailer body. A subject describing the bookkeeping ("record the
  ESM-only migration as a breaking change") therefore produced a
  changelog line that told a `require()` consumer nothing. The subject
  states the break directly instead: `chore!: drop CommonJS support;
  entry points move under dist/`. Put the break in the subject line of
  any future breaking commit for the same reason.

`auto-changelog`'s default `commitLimit` of 3 is kept deliberately.
It selects by diff size rather than significance, so the two entries
below the breaking one are `setup` and a test commit, which
under-describes the largest release in this fork's history. Raising
it with `--commit-limit false` is not the fix: regeneration rewrites
the whole file, so every one of the 23 older sections would expand
too. The breaking entry leads its section, and the README's
`## Usage` section carries the migration detail.

## Out of scope, per the design

- Reformatting or linting `src/`. No `eslint.config.js` is added and no
  `lint` or `preversion` script exists.
- Enabling TypeScript `strict`, which would require roughly 180 edits to
  vendored upstream code and would make every future merge from
  `git remote upstream` a conflict.
- Changing typedoc configuration or regenerating the tracked `docs/`
  output.
- Any change to the compression algorithms.
