# Dist Build Restructure Implementation Plan

## Phase 6: Browser test suite

**Goal:** Write the ZIP, streaming, and async coverage that has never
existed, and run it in a real browser against both `dist/browser/index.js`
and `dist/browser/index.min.js`.

**Architecture:** Assertions live in shared suite modules that take the
fflate namespace as a parameter. A single entry point imports both the
plain and minified bundles and runs every suite twice, once per bundle,
so a failure names which build broke. esbuild bundles the entry and pipes
it to `tapout`, which drives Playwright.

**Tech Stack:** @substrate-system/tapzero, @substrate-system/tapout,
Playwright, esbuild.

**Scope:** Phase 6 of 7.

**Codebase verified:** 2026-08-06

**Depends on:** Phases 3, 4, and 5. It consumes the built bundles and the
`test/browser/` directory that Phase 5 created.

**Note on AC8.1 after this phase.** `test/browser/index.ts` value-imports
`../../dist/browser/index.js`, so from here on a bare `tsc --noEmit`
fails with `TS2307` unless a build has run. Every gate that checks
AC8.1 from this point runs `npm run build` first, including Phase 7
Task 4. Deriving the `Fflate` *type* from `src/` (Task 2) removes the
build dependency for type resolution but not for the value import.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC3: Minified output preserves wcln name recovery

- **dist-build-restructure-design.AC3.1 Success:** The async callback
  APIs (`deflate`, `inflate`, `gzip`, `gunzip`, `zlib`, `unzlib`)
  round-trip data correctly when called on
  `dist/browser/index.min.js`.
- **dist-build-restructure-design.AC3.2 Success:** The async streaming
  classes (`AsyncDeflate`, `AsyncInflate`, `AsyncGzip`, `AsyncGunzip`)
  round-trip chunked input correctly on the minified bundle.
- **dist-build-restructure-design.AC3.3 Success:** The async ZIP APIs
  (`zip`, `unzip`, `AsyncZipDeflate`, `AsyncUnzipInflate`) round-trip a
  multi-file archive correctly on the minified bundle.
- **dist-build-restructure-design.AC3.4 Failure:** No async operation on
  the minified bundle produces output that differs from the same
  operation on the unminified bundle.

### dist-build-restructure-design.AC6: The browser suite runs and passes

- **dist-build-restructure-design.AC6.1 Success:** `npm run test:browser`
  exits 0, running the ZIP, stream, and async suites against
  `dist/browser/index.js`.
- **dist-build-restructure-design.AC6.2 Success:** The same suites pass
  against `dist/browser/index.min.js` in the same run.

---

## Investigation findings

Verified on 2026-08-06.

1. **The `tapout` on npm is the wrong package.** The design's script is
   `esbuild test/browser/index.ts --bundle | tapout`. The bare `tapout`
   package on npm is version 2.0.0, last modified in 2022, described as
   "Tap output stream". It is not a browser runner. The browser runner is
   **`@substrate-system/tapout`**, currently 0.0.41, described as "Run
   tests in a browser from the command line", which installs a `tapout`
   binary and depends on Playwright. The design's script text is correct
   only once the scoped package is the thing installed. This pairs with
   `@substrate-system/tapzero`, already a devDependency, from the same
   org.

2. **`@substrate-system/tapout` pulls in Playwright** (`playwright` and
   `@playwright/test`, both `^1.58.2`). Playwright browsers are a
   separate download, so `npx playwright install chromium` is a required
   setup step and CI must account for it.

3. **It reads the bundle from stdin** and exits 1 on any failure,
   including uncaught exceptions and console errors, which makes the
   pipeline in the design's script correct as written.

4. **`@substrate-system/tapzero` declares no browser condition** in its
   `exports` map, only `import` and `require`. esbuild resolves `import`
   and inlines it, which works, but the package is not explicitly built
   for a browser target. If it turns out to reference node builtins at
   module scope, the bundle will fail; Task 1 verifies this early rather
   than discovering it after the suites are written.

5. **The public API surface was enumerated from `src/index.ts`.** The
   names the suites exercise are: sync `deflateSync`, `inflateSync`,
   `gzipSync`, `gunzipSync`, `zlibSync`, `unzlibSync`, `zipSync`,
   `unzipSync`; async callback `deflate`, `inflate`, `gzip`, `gunzip`,
   `zlib`, `unzlib`, `decompress`, `zip`, `unzip`; streaming `Deflate`,
   `Inflate`, `Gzip`, `Gunzip`, `Zlib`, `Unzlib`, `Decompress` and their
   `Async` counterparts; ZIP streaming `Zip`, `Unzip`,
   `ZipPassThrough`, `ZipDeflate`, `AsyncZipDeflate`,
   `UnzipPassThrough`, `UnzipInflate`, `AsyncUnzipInflate`; helpers
   `strToU8`, `strFromU8`, `EncodeUTF8`, `DecodeUTF8`.

6. **Bundling the minified file is still a valid test of the hazard, with
   one caveat worth recording.** esbuild renames identifiers when it
   bundles two modules that share names, so the code the browser runs is
   not byte-identical to the shipped `index.min.js`. That does not
   weaken the test: the `wcln` hazard is terser substituting a
   constant's *value* for its *name* inside the `bInflt`/`bDflt` array
   literals, which would already have happened in the file esbuild reads.
   esbuild's renaming is consistent across the whole module, which is the
   property `wcln` depends on. What this setup cannot catch is a defect
   introduced by esbuild's own renaming, which is not a risk this design
   raises.

7. **A reduced model of the hazard was measured and passed.** Minifying a
   `bInflt`-shaped thunk with the proven terser options plus
   `module: true` produced `[n,t,r,e,s,c,i,l]`, all bare identifiers. The
   real `bDflt` has thirty entries against the model's eight, so this
   phase remains the actual gate.

---

<!-- START_SUBCOMPONENT_A (tasks 1-2) -->

<!-- START_TASK_1 -->
### Task 1: Install the browser runner and verify the pipeline

**Verifies:** None (infrastructure)

**Files:**
- Modify: `package.json` (devDependencies)

**Implementation:**

Install the scoped package. Installing bare `tapout` gets an unrelated
2022 library and the suite will not run.

```bash
npm install --save-dev @substrate-system/tapout
npx playwright install chromium
```

Before writing any suites, prove the toolchain end to end with a
throwaway test, so a tapzero-in-browser problem surfaces now rather than
after three suites exist.

The probe files must live **inside the repository**. esbuild resolves
bare specifiers by walking up from the entry file, so a probe in the
session scratchpad cannot find `node_modules` and fails with
`Could not resolve "@substrate-system/tapzero"`. This is the one case
where the scratchpad convention does not apply.

```bash
mkdir -p test/browser
cat > test/browser/_smoke.ts <<'EOF'
import { test } from '@substrate-system/tapzero'
test('smoke', t => { t.ok(true, 'runner works') })
EOF
npx esbuild test/browser/_smoke.ts --bundle | npx tapout
```

**Verification:**

Expected: the smoke command prints TAP output containing `ok 1` and exits
0.

Run: `npx esbuild test/browser/_smoke.ts --bundle | npx tapout; echo "exit=$?"`
Expected: `exit=0`.

Then confirm the failure path, because a runner that cannot fail is worse
than no runner:

```bash
cat > test/browser/_fail.ts <<'EOF'
import { test } from '@substrate-system/tapzero'
test('deliberate failure', t => { t.ok(false, 'should fail') })
EOF
npx esbuild test/browser/_fail.ts --bundle | npx tapout; echo "exit=$?"
```
Expected: `exit=1` and a `not ok` line.

Clean up: `rm -f test/browser/_smoke.ts test/browser/_fail.ts`
Expected: both probe files are gone before Task 2 begins. They must not
be committed.

If the smoke test fails to bundle because tapzero references node
builtins, stop and report it. The fallback is to bundle with
`--platform=browser --define:process.env.NODE_ENV='"test"'`, but confirm
the cause before adding flags.

**Commit:** `test: add browser test runner`
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Create the shared fixture and suite scaffolding

**Verifies:** None (infrastructure)

**Files:**
- Create: `test/browser/util.ts`

**Implementation:**

Create `test/browser/util.ts` holding the fixture data and helpers every
suite shares. Fixtures are generated in the browser rather than
downloaded: the node suite's multi-megabyte network fixtures are wrong
here, and deterministic generated data makes failures reproducible.

The module must export:

- `type Fflate = typeof import('../../src/index.js')`, the namespace type
  every suite takes as a parameter. Note it is derived from **source**,
  not from `dist/`. Deriving it from the built bundle would make
  `tsc --noEmit` fail whenever a build has not run, and no phase
  establishes that ordering for the typecheck.
- `fixtures()`, returning a small set of `Uint8Array` inputs that
  between them cover the compressor's branches: highly compressible
  (a long run of one byte), incompressible (a deterministic
  pseudo-random sequence, seeded so runs are comparable), text (UTF-8
  including multi-byte characters), and empty.
- `largeCompressible()` and `largeIncompressible()`, large fixtures
  for testing the async worker paths. The two paths want fixtures with
  OPPOSITE properties, so do not share one:

  | consumer | gate | fixture |
  |---|---|---|
  | `zip()` | `originalSize >= 160000` | `largeCompressible()` |
  | `unzip()` | `originalSize >= 524288` AND `compressed <= 0.8 * original` | `largeCompressible()` |
  | `AsyncUnzipInflate` | the COMPRESSED size, passed as `sz`, `>= 320000` | `largeIncompressible()` |

  `largeCompressible()` must therefore be at least 524288 bytes, not
  just 160000.

  `largeIncompressible()` must be measurably incompressible. Do not
  generate it with an LCG of the form
  `state = (state * 1103515245 + 12345) & 0x7fffffff`: in JS that
  product exceeds `2**53`, so the low bits are lost to float rounding
  before the mask applies and the sequence degenerates. That version
  deflated 600000 bytes to 19536 and silently sent `AsyncUnzipInflate`
  down its synchronous branch. Use a generator that stays in exact
  32-bit integer arithmetic (xorshift128), and verify the result by
  measuring `deflateSync(fixture).length`, not by assuming.
- `eq(a:Uint8Array, b:Uint8Array):boolean`, a length-then-bytes
  comparison.
- `chunks(data:Uint8Array, n:number):Uint8Array[]`, splitting input into
  `n` roughly equal slices for the streaming suites.
- `concat(parts:Uint8Array[]):Uint8Array`, joining chunks from streaming
  classes for comparison.
- `withTimeout(fn, ms?)`, racing a promise-returning function against a
  timeout. One definition, imported by all three suites. Its default
  must stay below tapout's auto-finish window -- see the `--timeout`
  note in phase_04.md, and Task 6 below.
- `firstDiff(a, b):string | null`, returning null when two arrays are
  equal and otherwise the first differing offset with both byte values,
  or a length mismatch. Byte-identity assertions use this instead of a
  bare boolean, so a two-byte header skew is distinguishable from the
  wholesale corruption a `wcln` break produces.

Use a seeded generator, not `Math.random()`, so a failure can be
reproduced.

**Verification:**

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep "error TS" | grep -c "test/browser/util.ts" || true`
Expected: `0`.

Run: `npx esbuild test/browser/util.ts --bundle --outfile=/dev/null`
Expected: exits 0.

**Commit:** `test: add browser suite fixtures and helpers`
<!-- END_TASK_2 -->

<!-- END_SUBCOMPONENT_A -->

<!-- START_SUBCOMPONENT_B (tasks 3-5) -->

<!-- START_TASK_3 -->
### Task 3: Write the async suite

**Verifies:** dist-build-restructure-design.AC3.1,
dist-build-restructure-design.AC3.4

**Files:**
- Modify: `test/browser/async.ts` (currently the moved TODO stub)

**Implementation:**

Replace the stub. Export a single function
`asyncSuite(f:Fflate, label:string)` that registers tapzero tests named
with `label` as a prefix, so the same assertions run against both bundles
and a failure identifies which one.

**Naming convention, used by all three suites and asserted in Task 6:**
every test name is `` `${label} > ${description}` ``. The `# min > ` and
`# plain > ` comment-line prefixes are what Task 6 greps to count how many
TESTS per bundle actually ran -- 103 each, across three suites -- and so
to prove both passes completed rather than one being truncated.

This suite is the primary gate on the `wcln` hazard. Every async
callback API spawns a worker whose source is assembled by `wcln` from
`fn.toString()` text, so a name-recovery break shows up here as corrupt
output rather than as an exception.

**Testing:**

Tests must verify these AC cases:

- `dist-build-restructure-design.AC3.1`: for each of `deflate`/`inflate`,
  `gzip`/`gunzip`, and `zlib`/`unzlib`, and for each fixture, compressing
  then decompressing returns bytes equal to the input. Assert on the
  bytes, not on the absence of an error, because a `wcln` break yields
  wrong output rather than a thrown error.
- `dist-build-restructure-design.AC3.1`: `decompress` correctly
  auto-detects and reverses all three container formats.
- `dist-build-restructure-design.AC3.4`: for each fixture, output from
  the async API equals output from the corresponding sync API on the same
  bundle. This is the assertion that catches a silent name-recovery
  failure, since the sync path does not go through `wcln` at all.
- The callback error slot is `null` on success. The codebase passes
  `null` into callback error slots throughout, and a worker that failed
  to build would surface here.

Each async call must be wrapped in a promise with a timeout, via the
shared `withTimeout` from Task 2. A `wcln` break can hang a worker
rather than failing it, and an un-timed-out hang never reports at all.

Do NOT give it a 10-second timeout, even though that matches the node
suite's `cProm.timeout(10000)`. The browser runner has a constraint the
node suite does not: tapout ends a run after a window of console
silence, and a timeout longer than that window means auto-finish wins
and the run is truncated rather than reported. `withTimeout` defaults to
2000 ms against a 3000 ms window; see the `--timeout` note in
phase_04.md for why the two values are a pair, and Task 6 for the test
that holds them together.

Follow the tapzero idiom already used in `test/0-valid.ts`:
`test(name, async t => { ... t.ok(value, message) })`.

**Verification:**

Run: `npx esbuild test/browser/async.ts --bundle --outfile=/dev/null`
Expected: exits 0.

Full-suite verification happens in Task 6.

**Commit:** `test: add browser async suite`
<!-- END_TASK_3 -->

<!-- START_TASK_4 -->
### Task 4: Write the streaming suite

**Verifies:** dist-build-restructure-design.AC3.2,
dist-build-restructure-design.AC3.4

**Files:**
- Modify: `test/browser/streams.ts` (currently the moved TODO stub)

**Implementation:**

Replace the stub. Export `streamSuite(f:Fflate, label:string)` following
the same parameterised shape as Task 3.

**Testing:**

Tests must verify these AC cases:

- `dist-build-restructure-design.AC3.2`: each of `AsyncDeflate`,
  `AsyncInflate`, `AsyncGzip`, `AsyncGunzip`, `AsyncZlib`, and
  `AsyncUnzlib` accumulates chunked input pushed via `push(chunk, final)`
  and produces, through its `ondata` handler, output that round-trips to
  the original. Use `chunks()` from Task 2 with more than one chunk, and
  include a case where the final chunk is empty, since `push(empty, true)`
  is a real usage pattern that exercises the terminating branch.
- `dist-build-restructure-design.AC3.2`: the synchronous streaming
  classes (`Deflate`, `Inflate`, `Gzip`, `Gunzip`, `Decompress`) do the
  same. These do not use workers, so a divergence between them and their
  `Async` counterparts localises a fault to the worker path.
- `dist-build-restructure-design.AC3.4`: async streamed output equals
  SYNC STREAMED output over the same input at the same chunking.

  Do not compare a stream against one-shot `gzipSync`/`zlibSync`. That
  is not a real invariant and the first implementation failed on it:
  chunking changes DEFLATE block boundaries, so on incompressible input
  even the sync stream differs from the one-shot output (561 vs 535
  bytes, measured 2026-08-07). Async-equals-sync at identical chunking
  is the comparison that isolates the worker path, which is what this
  AC is for.
- `EncodeUTF8` and `DecodeUTF8` round-trip text containing multi-byte
  characters, including a character split across a chunk boundary. That
  boundary case is the reason those classes hold state.

**Verification:**

Run: `npx esbuild test/browser/streams.ts --bundle --outfile=/dev/null`
Expected: exits 0.

**Commit:** `test: add browser streaming suite`
<!-- END_TASK_4 -->

<!-- START_TASK_5 -->
### Task 5: Write the ZIP suite

**Verifies:** dist-build-restructure-design.AC3.3,
dist-build-restructure-design.AC3.4

**Files:**
- Modify: `test/browser/zip.ts` (currently the moved TODO stub)

**Implementation:**

Replace the stub. Export `zipSuite(f:Fflate, label:string)`.

**Testing:**

Tests must verify these AC cases:

- `dist-build-restructure-design.AC3.3`: `zipSync` then `unzipSync` on a
  multi-file archive, including a nested directory path and an empty
  file, returns every original file byte for byte with its path intact.
- `dist-build-restructure-design.AC3.3`: the async `zip` and `unzip`
  produce the same result as their sync counterparts. This is the ZIP
  arm of the `wcln` gate, since async ZIP compression runs in a worker.
- `dist-build-restructure-design.AC3.3`: the streaming `Zip` and `Unzip`
  classes, driven with `ZipDeflate`, `AsyncZipDeflate`,
  `ZipPassThrough`, `UnzipInflate`, and `AsyncUnzipInflate`, round-trip
  a multi-file archive.
- `dist-build-restructure-design.AC3.4`: an archive produced by the
  async path is byte-identical to one produced by the sync path at the
  same compression level, and each is readable by the other's decoder.

  Pass a fixed `mtime` to BOTH calls. fflate stamps wall-clock time
  into every ZIP header (`src/index.ts:2917`), and the DOS time field
  has 2-second resolution, so unpinned the two archives differ
  whenever the calls straddle a tick.

  The same applies to every gzip byte comparison, in Tasks 3 and 4 as
  well as here: `gzh` writes `Date.now() / 1000` into header bytes 4-7
  (`src/index.ts:1199`) at 1-second resolution. Measured on the empty
  fixture: 4 mismatches in 200 at offset 4 unpinned, 0 in 200 pinned.
  This bit `gzip sync/async match` and `AsyncGzip vs sync`, one-shot
  and streaming respectively. zlib and raw deflate carry no timestamp
  and need no pinning. Measured: at 600000 bytes the
  archives diverge at byte offsets 10 and 650 -- the mod-time field of
  the local file header and of the central directory record --
  1 time in 60 unaligned, 29 in 60 when aligned to a tick. Pinned:
  0 in 60. Note the async overload is `zip(files, opts, cb)`; passing
  `undefined` as `opts` throws.

  Report the first differing offset on mismatch rather than asserting
  a bare boolean. A `wcln` name-recovery break corrupts wholesale; a
  header skew differs in two bytes. Without a diff the next maintainer
  cannot tell them apart, and phase_06.md sends them to the terser
  options in `scripts/build.ts` for either one.

  Cover both directions of "readable by the other's decoder": the
  async decoder over a sync-produced archive, and `unzipSync` over an
  async-produced one.
- A `ZipPassThrough` entry (stored, not deflated) survives a round trip,
  since stored entries take a different code path from deflated ones.

Every test that claims to exercise an async branch must prove it took
that branch, with an assertion that fails when it did not. Round-trip
equality alone proves nothing here: the sync fallback returns the same
bytes, so the test passes either way. Two specific traps:

- A flag set inside the `Unzip` `onfile` callback is vacuous. That
  callback fires for the synchronous decoder too, so the assertion can
  never fail.
- An assertion placed after a poll that only resolves once the asserted
  value is truthy is vacuous for the same reason. Poll on a separate
  `settled` flag that BOTH the success and the error branch set, then
  assert. Polling on the success value alone also means an error spins
  until the timeout instead of failing fast.
- `Unzip` keeps the decoder instance private, so wrap the registered
  constructor to observe it. `AsyncUnzipInflate` assigns `terminate`
  only on its async branch, which makes it the one externally visible
  difference from the `Inflate` fallback.

Also assert the fixture property that gates the branch (the sizes in
Task 2's table), so a fixture that drifts fails loudly instead of
silently downgrading to the sync path. Confirm each new guard bites by
temporarily breaking the fixture and watching the test go `not ok`.

**Verification:**

Run: `npx esbuild test/browser/zip.ts --bundle --outfile=/dev/null`
Expected: exits 0.

**Commit:** `test: add browser zip suite`
<!-- END_TASK_5 -->

<!-- END_SUBCOMPONENT_B -->

<!-- START_TASK_6 -->
### Task 6: Wire the entry point and run against both bundles

**Verifies:** dist-build-restructure-design.AC6.1,
dist-build-restructure-design.AC6.2,
dist-build-restructure-design.AC3.1,
dist-build-restructure-design.AC3.2,
dist-build-restructure-design.AC3.3

**Files:**
- Create: `test/browser/min.ts`
- Create: `test/browser/index.ts`

**Implementation:**

Two files. First create `test/browser/min.ts`, which is the only place
the minified bundle is referenced.

`dist/browser/index.min.js` has no declaration file, and `allowJs` is
`false`, so importing it directly from `index.ts` fails the typecheck
with `TS2307` and would break Phase 7 Task 4's re-verification of AC8.1.
Isolating the untyped import in one small module keeps the suppression
narrow and gives the value a real type.

```ts
// The minified bundle ships without declarations, and generating them
// would be meaningless: its public API is identical to the plain
// bundle's by construction. Assert that here, in one place.
// @ts-expect-error no declarations are emitted for the minified bundle
import * as minBundle from '../../dist/browser/index.min.js'
import type { Fflate } from './util.js'

// The double cast (as unknown as Fflate) is necessary because the
// dist/ bundle's types are structurally different from src/ due to
// private members on classes. The minified bundle's runtime value is
// identical to the plain bundle's by construction, so the cast is safe.
export const min = minBundle as unknown as Fflate
```

Then create the entry point. It imports both bundles and runs all three
suites against each.

```ts
// Runs every suite twice: once against the plain bundle and once
// against the minified one.
//
// The minified pass is the point of this file. The async APIs build
// worker source by parsing identifier names out of array literal TEXT
// via Function.prototype.toString, so a minifier that inlines a
// constant into one of those literals breaks name recovery silently:
// the sync APIs keep working while the async ones corrupt output.
// Nothing else in this repository would catch that.
import * as plain from '../../dist/browser/index.js'
import { min } from './min.js'
import { asyncSuite } from './async.js'
import { streamSuite } from './streams.js'
import { zipSuite } from './zip.js'
import { harnessSuite } from './harness.js'

// Runs once, not per bundle: it checks the runner's own configuration,
// which is what decides whether the minified pass below runs at all.
harnessSuite()

for (const [label, f] of [
    ['plain', plain],
    ['min', min]
] as const) {
    asyncSuite(f, label)
    streamSuite(f, label)
    zipSuite(f, label)
}
```

Also add `test/browser/harness.ts`, exporting `harnessSuite()`. It
guards the one invariant that, when broken, makes this suite lie
rather than fail: `withTimeout`'s default must stay below tapout's
auto-finish window, or a stalled test truncates the run and the
minified pass never registers. Measured at tapout's default timeout
with one deliberately stalled test: zero `# min >` headers, no plan
line, and the run reported as a PASS.

The invariant spans two files, so assert the relationship rather than
either number: read `--timeout` back out of package.json's
`test:browser` script (esbuild inlines the JSON import), recompute
tapout's window as `max(500, min(3000, floor(timeout * 0.2)))`, and
assert `DEFAULT_TIMEOUT_MS` fits under it with a margin -- the two
timers do not start together, so a bare `<` would pass a default that
still loses the race. Accept `-t` as well as `--timeout`, since tapout
treats them as aliases, and take the LAST occurrence, since tapout
does. Fail loudly if neither flag is present, because tapout's default
of 5000 puts the window at 1000 ms.

The `test:browser` script set in Phase 4 now reads
`esbuild test/browser/index.ts --bundle | tapout --timeout 30000`.
Confirm it is unchanged; see the note in phase_04.md for why the
`--timeout` is load-bearing.

**Verification:**

The build must exist. Run `npm run build` first.

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep "error TS" | grep -c "test/browser" || true`
Expected: `0`. This is the check that the isolated `min.ts` exists to
satisfy; a `TS2307` here means the minified bundle was imported directly.

Run: `npm run build && npm run test:browser; echo "exit=$?"`
Expected: `exit=0`, with TAP output containing both `plain` and `min`
prefixed test names and no `not ok` lines.

Run: `npm run test:browser 2>&1 | grep -c "^not ok"`
Expected: `0`.

Run: `npm run test:browser 2>&1 | grep -c "^# min >"` compared against
`npm run test:browser 2>&1 | grep -c "^# plain >"`
Expected: equal counts (both `103` on healthy tree). If min is less than
plain, the minified pass was truncated before completion, indicating a
timeout or crash. **If either check fails while the `plain` equivalents
pass, the `wcln` hazard has fired.** Do not work around it in the tests.
Per the design's risk register, the response is to narrow the terser
`compress` options for the ESM outputs in `scripts/build.ts`, or fall
back to shipping `dist/*/index.min.js` built with `mangle` only, then
re-run this suite. Record which option was taken and why.

Run: `npm test; echo "exit=$?"`
Expected: `exit=0`. Build, node suite, and browser suite all pass. This
is the first point in the plan where the full `npm test` chain works.

**Commit:** `test: run browser suites against plain and minified bundles`
<!-- END_TASK_6 -->

---

## Phase 6 completion criteria

- `npm run test:browser` exits 0.
- Both `plain` and `min` prefixed tests appear in the output.
- `npm test` exits 0 end to end.
- The three files that were single-line TODO comments now contain real
  ZIP, streaming, and async coverage.
- The async, streaming, and ZIP APIs are proven to produce output
  identical to their sync counterparts on the minified bundle, which is
  the only real check on the `wcln` hazard.
