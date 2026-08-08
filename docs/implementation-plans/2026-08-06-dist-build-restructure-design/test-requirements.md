# Test Requirements: dist build restructure

Date: 2026-08-06

Source design: `specs/2026-08-06-dist-build-restructure-design.md`
Source phases: `phase_01.md` through `phase_07.md` in this directory.

All acceptance criteria are scoped with the prefix
`dist-build-restructure-design.`. This document uses the short form
(`AC1.1`) in headings and the full scoped form in the mapping tables.

## How this plan is tested

This is a build restructure, not a feature. The overwhelming majority of
its acceptance criteria are verified by **build-verification**: shell
commands that run the build, then assert on the shape and content of the
emitted artifacts, on `package.json` fields, or on the absence of
removed files. The plan creates **no unit tests** and no test framework
harness for the build script itself. Inventing unit tests here would
misrepresent what the phases actually produce.

There is exactly one body of real behavioural test code in the whole
plan: the browser suite written in Phase 6, at `test/browser/async.ts`,
`test/browser/streams.ts`, and `test/browser/zip.ts`, driven by
`test/browser/index.ts`. That suite exists to gate the `wcln` hazard and
covers AC3.1 through AC3.4 and AC6.1 and AC6.2. It is documented in the
most detail below.

The node suite (AC5) is **pre-existing** coverage. Phase 5 repoints it at
the new bundle; it writes no new assertions.

Test-type vocabulary used below:

| type | meaning in this plan |
| --- | --- |
| build-verification | shell command asserting on build output or config |
| behavioural | real assertions in test code, run by a test runner |
| e2e | packed tarball installed and consumed as a real dependency |
| human | manual step, justified in place |

---

## AC1: The build produces the specified dist/ layout

Produced by Phase 3. Re-verified as a combined gate by Phase 7 Task 4.
All five criteria are build-verification. There is no test file; the
verification commands are the tests.

### AC1.1

> **dist-build-restructure-design.AC1.1 Success:** `npm run build` writes
> `dist/browser/index.js`, `index.js.map`, `index.min.js`,
> `index.min.js.map`, and `index.d.ts`.

Automated, build-verification. Produced by Phase 3 Task 3 (declarations)
and Phase 3 Task 4 (browser bundle), gated by Phase 3 Task 8.

Commands, verbatim from the phase:

- `ls dist/browser` (Task 4) -- expected `index.d.ts`, `index.js`,
  `index.js.map`, `index.min.js`, `index.min.js.map`.
- `ls dist/browser/index.d.ts dist/node/index.d.ts` (Task 3) -- both
  exist.
- `rm -rf dist && npx tsx scripts/build.ts && find dist -type f | sort`
  (Task 8) -- expected exactly the twelve listed paths.
- `find dist -type f | sort` (Phase 7 Task 4) -- same twelve paths.

Task 8 additionally asserts terser actually ran, which is what makes
`index.min.js` distinct from `index.js`:

- `node -e "const {statSync} = require('fs'); ..."` -- `.min.js` size
  must be smaller than `.js` size; equal sizes mean terser did not run.
- `head -c 200 dist/browser/index.min.js` -- minified single line, not
  identical to the head of `index.js`.

### AC1.2

> **dist-build-restructure-design.AC1.2 Success:** `npm run build` writes
> the same five files under `dist/node/`.

Automated, build-verification. Produced by Phase 3 Task 5, gated by Phase
3 Task 8.

- `ls dist/node` (Task 5) -- the same five names.
- `find dist -type f | sort` (Task 8, Phase 7 Task 4) -- twelve paths.

Task 5 also runs two functional smoke checks against the emitted bundle,
which is the closest thing to a unit test in Phase 3:

- `node --input-type=module -e "import('./dist/node/index.js').then(m =>
  { const d = m.deflateSync(new Uint8Array([1,2,3]));
  console.log(m.inflateSync(d).length) })"` -- prints `3`.
- the async equivalent, chaining `m.deflate` into `m.inflate` -- prints
  `3`. This one matters because it is the first point in the plan where a
  worker is actually spawned from a built artifact.

### AC1.3

> **dist-build-restructure-design.AC1.3 Success:** `npm run build` writes
> a minified `dist/umd/fflate.js`.

Automated, build-verification. Produced by Phase 3 Task 6.

- `ls dist/umd` -- `fflate.js` and `package.json`. The second is
  load-bearing, not residue: it contains `{ "type": "commonjs" }`, and
  without it the root `"type": "module"` makes node parse the UMD
  wrapper as ESM, where `this` is undefined and BOTH `require()` and
  `import()` throw `Cannot set properties of undefined (setting
  'fflate')`. Measured on node v25.8.2.
- the `new Function('self', 'module', 'exports', 'define', s)` harness in
  Task 6 -- prints `umd roundtrip: 3`.

Note the deliberate constraint recorded in Phase 3 Task 6: the artifact
must **not** be verified with `require()`. The package is
`"type": "module"`, so Node parses `dist/umd/fflate.js` as ESM and a
`require`-based check returns a namespace of `undefined` exports rather
than throwing, which would pass a naive assertion while proving nothing.
The global-branch harness exists specifically to avoid that trap, and
`module`, `exports`, and `define` must be shadowed as parameters or the
UMD wrapper takes its CommonJS branch and never assigns the global.

### AC1.4

> **dist-build-restructure-design.AC1.4 Success:** The build creates no
> `lib/`, `esm/`, or top level `umd/` directory.

Automated, build-verification. Negative criterion, verified by absence.
Produced by Phase 3 Task 7 (deletes the superseded scripts) and Phase 7
Task 1 (removes the stale directories).

- `ls lib esm umd 2>&1` (Phase 3 Task 8, Phase 7 Task 1) -- three "No
  such file or directory" errors.
- `npm run build && find dist -type d -maxdepth 1 | sort` (Phase 7 Task
  1) -- `dist`, `dist/browser`, `dist/node`, `dist/umd`, then `ls lib esm
  umd 2>&1` again, still three errors after the build. This ordering is
  the actual test: it proves the build does not recreate them, which
  `ls` alone before a build would not.
- `git ls-files lib esm umd` (Phase 7 Task 1) -- no output, a
  precondition check confirming nothing tracked is being deleted.

### AC1.5

> **dist-build-restructure-design.AC1.5 Success:** The build removes a
> stale `dist/` before emitting, so no artifact from a previous run
> survives.

Automated, build-verification. Produced by Phase 3 Task 1 (the
`rmSync(p('dist'), { recursive: true, force: true })` clean step).

- `npx tsx scripts/build.ts` then `ls dist` (Task 1) -- prints
  `cleaned dist/`, then `browser`, `node`, `umd`.
- `find dist -type f | sort` after a full run (Task 8) -- exactly twelve
  paths. Any surviving stale artifact from a previous run appears as a
  thirteenth path, which is how this criterion is actually caught.

---

## AC2: The worker swap is correct per build target

AC2.2 is produced by Phase 2; AC2.1 and AC2.3 by Phase 3. All
build-verification or static inspection.

### AC2.1

> **dist-build-restructure-design.AC2.1 Success:**
> `dist/browser/index.js` contains the `src/worker.ts` implementation,
> identifiable by `URL.createObjectURL` and `Blob`.

Automated, build-verification. Produced by Phase 3 Task 4 (the
`browser-worker-swap` esbuild resolve plugin).

- `grep -c "createObjectURL" dist/browser/index.js` -- at least `1`.
- `node --input-type=module -e "import('./dist/browser/index.js').then(m
  => console.log(typeof m.deflateSync, typeof m.deflate))"` -- prints
  `function function`.

The complementary check lives on the node side, in Phase 3 Task 5:
`grep -c "createObjectURL" dist/node/index.js || true` expects `0`. The
pair together is what proves the swap is per-target rather than
unconditional.

### AC2.2

> **dist-build-restructure-design.AC2.2 Success:** `src/node-worker.ts`
> obtains `Worker` and `isMarkedAsUntransferable` from a static
> `node:worker_threads` import, with no `require` at module scope and no
> unsupported-runtime fallback branch.

Automated, static source verification. Produced by Phase 2 Task 3.

- `grep -n "^import { Worker, isMarkedAsUntransferable } from
  'node:worker_threads'" src/node-worker.ts` -- one match on line 1.
- `grep -n "require(" src/node-worker.ts` -- exactly one match, on the
  `workerAdd` line, inside the string literal. This check carries real
  content: the `workerAdd` string contains
  `require('worker_threads')` which **must stay byte for byte
  unchanged**, because it is evaluated inside a `{ eval: true }` worker
  that Node treats as CommonJS regardless of the host package's `type`
  field. A verification that simply expected zero `require` matches would
  drive an implementer to break the harness.
- `grep -n "async operations unsupported" src/node-worker.ts` -- no
  output, confirming the fallback branch is gone.

Behavioural confirmation is indirect but real: Phase 3 Task 5's async
round-trip through `dist/node/index.js` spawns a worker through this
module, and the whole node suite (AC5) depends on it.

### AC2.3

> **dist-build-restructure-design.AC2.3 Failure:**
> `dist/browser/index.js` contains no reference to `worker_threads` at
> module scope.

Automated, build-verification. **Negative criterion, verified by grep for
absence.** Produced by Phase 3 Task 4.

- `grep -c "worker_threads" dist/browser/index.js || true` -- expected
  `0`. The `|| true` is required by the plan's verification conventions:
  `grep -c` exits non-zero on a zero count, so without it a passing check
  reports as a failure under a fail-fast runner.
- `grep -c "worker_threads" dist/umd/fflate.js || true` (Phase 3 Task 6)
  -- expected `0`. The UMD artifact bundles the browser variant and
  reuses the same resolve plugin, so it inherits the same requirement.

Limitation worth recording: the grep is over the whole file, not
scoped to module scope, so it is strictly stronger than the criterion
states. That is fine -- a stronger check cannot produce a false pass.
It could produce a false failure if the string `worker_threads` ever
appeared legitimately in browser output, which it does not, because the
`workerAdd` string containing it lives only in `src/node-worker.ts` and
the plugin ensures that file never enters the browser graph.

---

## AC3: Minified output preserves wcln name recovery

**This is the only AC group with real behavioural test code, and it is
the highest-risk part of the whole change.** Produced by Phase 6, Tasks 3
through 6.

### Why this group is different

`src/index.ts:1040` defines `wcln`, which reconstructs variable names by
calling `.toString()` on a thunk such as `bInflt`, slicing the text
between `[` and `]`, splitting on commas, and pairing the recovered
identifier strings positionally with the evaluated values. Two build
properties must hold or this breaks:

1. the minifier must mangle the top-level scope **consistently**, so the
   names in the array-literal text still match the names referenced
   inside the serialised function bodies; and
2. the minifier must not **inline a constant into that array literal** --
   substituting a value for a name silently produces wrong keys.

The current repository ships unminified `tsc` output for `lib/` and
`esm/` and applies terser only to UMD, so minified ESM is new exposure.
A break is silent: the sync APIs keep working while the async ones
corrupt output. No grep, size check, or type check can detect it. Only
executing the async paths against the minified bundle can.

Phase 3 investigation finding 3 and Phase 6 investigation finding 7 both
record a reduced model of the `bInflt` pattern that was minified with the
production terser options plus `module: true` and emerged as
`[n,t,r,e,s,c,i,l]`, all bare identifiers. That raises confidence but
does not settle it: the model has eight entries, `bDflt` has thirty, and
`src/index.ts` is 3,885 lines. **Phase 6 is the actual gate.**

### Test architecture

| file | role |
| --- | --- |
| `test/browser/util.ts` | `Fflate` type, `fixtures()`, `eq()`, `chunks()` |
| `test/browser/async.ts` | exports `asyncSuite(f, label)` |
| `test/browser/streams.ts` | exports `streamSuite(f, label)` |
| `test/browser/zip.ts` | exports `zipSuite(f, label)` |
| `test/browser/min.ts` | isolates the untyped minified import |
| `test/browser/harness.ts` | guards the run's own completeness: asserts the `--timeout` in the `test:browser` script and `withTimeout`'s default stay a compatible pair |
| `test/browser/index.ts` | runs all three suites against both bundles |

Each suite is a function taking the fflate namespace and a label, not a
module with hard-coded imports. That parameterisation is what lets the
identical assertions run twice.

`test/browser/index.ts` iterates
`[['plain', plain], ['min', min]] as const` and calls `asyncSuite`,
`streamSuite`, and `zipSuite` for each. `plain` is
`import * as plain from '../../dist/browser/index.js'`; `min` comes from
`test/browser/min.ts`, which carries a `@ts-expect-error` on
`import * as minBundle from '../../dist/browser/index.min.js'` and
re-exports it as `Fflate`. The minified bundle ships no declarations and
`allowJs` is `false`, so isolating that one import keeps the suppression
narrow and preserves AC8.1.

`Fflate` in `test/browser/util.ts` is
`typeof import('../../src/index.js')` -- derived from **source**, not
from `dist/`. Deriving it from the built bundle would make
`tsc --noEmit` fail whenever a build had not run.

**Test naming convention and header anchors:** every test name is
`` `${label} > ${description}` ``, with `label` being `plain` or `min`.
The `# plain > ` and `# min > ` comment-line **prefixes** are
load-bearing: Task 6 greps these headers to count how many TESTS per
bundle actually ran -- 103 each, across three suites. The gates rely on
the header counts being equal between plain and min to prove neither
was truncated. A naming change that removes the `> ` breaks the gate
silently.

Read the two gates TOGETHER with the adjacent `exit=$?` check, and
against the absolute `103` rather than only against each other. A run
that produces no output at all -- a failed esbuild bundle, say -- gives
`^not ok` of 0 and header counts of 0 == 0, which satisfies both gates
in isolation. The exit code is what catches that case. Measured by
deleting `dist/browser/index.min.js` and re-running.

Fixtures are generated in the browser from a **seeded** generator, not
`Math.random()`, and cover highly compressible (a long run of one byte),
incompressible (deterministic pseudo-random), UTF-8 text including
multi-byte characters, and empty input.

Every async call is wrapped in a promise with a 10 second timeout,
matching the node suite's `cProm.timeout(10000)`. A `wcln` break can hang
a worker rather than failing it, and an un-timed-out hang surfaces as a
generic tapout timeout with no useful detail.

### AC3.1

> **dist-build-restructure-design.AC3.1 Success:** The async callback
> APIs (`deflate`, `inflate`, `gzip`, `gunzip`, `zlib`, `unzlib`)
> round-trip data correctly when called on `dist/browser/index.min.js`.

Automated, behavioural. Test file `test/browser/async.ts`, produced by
Phase 6 Task 3, executed by Phase 6 Task 6.

Required assertions:

- for each of `deflate`/`inflate`, `gzip`/`gunzip`, and `zlib`/`unzlib`,
  and for each fixture, compressing then decompressing returns bytes
  equal to the input. **Assert on the bytes, not on the absence of an
  error** -- a `wcln` break yields wrong output rather than a thrown
  error.
- `decompress` correctly auto-detects and reverses all three container
  formats.
- the callback error slot is `null` on success.

Run command: `npm run test:browser`, which is
`esbuild test/browser/index.ts --bundle | tapout --timeout 30000`. The runner is
`@substrate-system/tapout` (**not** the bare `tapout` package on npm,
which is an unrelated 2022 stream library); it reads the bundle from
stdin, drives Playwright, and exits 1 on any failure including uncaught
exceptions and console errors.

Per-file bundle check during development:
`npx esbuild test/browser/async.ts --bundle --outfile=/dev/null` exits 0.

### AC3.2

> **dist-build-restructure-design.AC3.2 Success:** The async streaming
> classes (`AsyncDeflate`, `AsyncInflate`, `AsyncGzip`, `AsyncGunzip`)
> round-trip chunked input correctly on the minified bundle.

Automated, behavioural. Test file `test/browser/streams.ts`, produced by
Phase 6 Task 4, executed by Phase 6 Task 6.

Required assertions:

- `AsyncDeflate`, `AsyncInflate`, `AsyncGzip`, `AsyncGunzip`,
  `AsyncZlib`, and `AsyncUnzlib` each accumulate chunked input pushed via
  `push(chunk, final)` and produce, through `ondata`, output that
  round-trips to the original. Use `chunks()` with more than one chunk,
  and include a case where the final chunk is empty, since
  `push(empty, true)` exercises the terminating branch.
- the synchronous streaming classes (`Deflate`, `Inflate`, `Gzip`,
  `Gunzip`, `Decompress`) do the same. These do not use workers, so a
  divergence between them and their `Async` counterparts **localises a
  fault to the worker path** -- that is the diagnostic value of including
  them.
- `EncodeUTF8` and `DecodeUTF8` round-trip text containing multi-byte
  characters, including a character split across a chunk boundary.

Note the AC text names four classes; the phase requires six, adding
`AsyncZlib` and `AsyncUnzlib`. The implementation is broader than the
criterion, which is correct.

### AC3.3

> **dist-build-restructure-design.AC3.3 Success:** The async ZIP APIs
> (`zip`, `unzip`, `AsyncZipDeflate`, `AsyncUnzipInflate`) round-trip a
> multi-file archive correctly on the minified bundle.

Automated, behavioural. Test file `test/browser/zip.ts`, produced by
Phase 6 Task 5, executed by Phase 6 Task 6.

Required assertions:

- `zipSync` then `unzipSync` on a multi-file archive, including a nested
  directory path and an empty file, returns every original file byte for
  byte with its path intact.
- async `zip` and `unzip` produce the same result as their sync
  counterparts. This is the ZIP arm of the `wcln` gate, since async ZIP
  compression runs in a worker.
- streaming `Zip` and `Unzip`, driven with `ZipDeflate`,
  `AsyncZipDeflate`, `ZipPassThrough`, `UnzipInflate`, and
  `AsyncUnzipInflate`, round-trip a multi-file archive.
- a `ZipPassThrough` entry (stored, not deflated) survives a round trip,
  since stored entries take a different code path from deflated ones.

### AC3.4

> **dist-build-restructure-design.AC3.4 Failure:** No async operation on
> the minified bundle produces output that differs from the same
> operation on the unminified bundle.

Automated, behavioural. Negative criterion, but unlike AC2.3 and AC4.5 it
is **not** verified by grep -- it is verified by differential assertion
inside all three suites (Phase 6 Tasks 3, 4, and 5).

The differential assertions required:

- async suite: for each fixture, output from the async API equals output
  from the corresponding sync API **on the same bundle**. This is the
  assertion that catches a silent name-recovery failure, because the sync
  path does not go through `wcln` at all.
- stream suite: streamed output equals one-shot sync output over the same
  input.
- zip suite: an archive produced by the async path is byte-identical to
  one produced by the sync path at the same compression level, and each
  is readable by the other's decoder.

Because all three suites run under both the `plain` and `min` labels, the
cross-bundle comparison the AC states is obtained transitively: async
equals sync on `plain`, async equals sync on `min`, and the sync path is
unaffected by the hazard.

The run-level gate, Phase 6 Task 6:

```
npm run test:browser 2>&1 | grep -c "^not ok"
```

Expected `0`.

```
npm run test:browser 2>&1 | grep -c "^# min >" == npm run test:browser 2>&1 | grep -c "^# plain >"
```

Expected both equal. If the min count is less than the plain count,
the minified pass was truncated, indicating a timeout or crash.

**If either check fails while the `plain` equivalents pass, the `wcln`
hazard has fired.** The plan's prescribed response is not to work around
it in the tests: narrow the terser `compress` options for the ESM
outputs in `scripts/build.ts`, or fall back to shipping `dist/*/index.min.js`
built with `mangle` only, then re-run. The comment block above
`terserOpts` in `scripts/build.ts` (Phase 3 Task 2) records the same
constraint at the point of change.

Phase 7 Task 4 re-runs both checks as the final acceptance gates.

### Known limitation of the AC3 gate

Phase 6 investigation finding 6 records it: esbuild renames identifiers
when bundling `dist/browser/index.min.js` into the test bundle, so the
code the browser runs is not byte-identical to the shipped file. This
does not weaken the test, because the hazard is terser substituting a
constant's *value* for its *name* inside the `bInflt`/`bDflt` array
literals, which would already have happened in the file esbuild reads,
and esbuild's renaming is consistent across the module -- the property
`wcln` depends on. What this setup cannot catch is a defect introduced by
esbuild's own renaming, which is not a risk the design raises. No
additional test is planned for it.

---

## AC4: Published entry points resolve

Produced by Phase 4. Re-verified against a packed and installed tarball
by Phase 7 Task 3, which is the stronger of the two passes because it
tests resolution the way a consumer experiences it.

### AC4.1

> **dist-build-restructure-design.AC4.1 Success:** Importing the package
> under the `node` condition resolves to `dist/node/index.js`, and its
> types to `dist/node/index.d.ts`.

Automated, build-verification (Phase 4 Task 1) plus e2e (Phase 7 Task 3).

Phase 4 Task 1:

- `node --input-type=module -e
  "import('@substrate-system/fflate').then(m => console.log('node cond:',
  typeof m.deflateSync))"` -- the import resolves.
- the six-path `for` loop asserting each of `dist/node/index.js`,
  `dist/node/index.d.ts`, `dist/browser/index.js`,
  `dist/browser/index.d.ts`, `dist/browser/index.min.js`, and
  `dist/umd/fflate.js` exists -- six `ok` lines, no `MISSING`.

Phase 7 Task 3, from an installed tarball in a scratch directory:

- `import * as f from '@substrate-system/fflate'` then `deflateSync` /
  `inflateSync` -- prints `sync roundtrip: hello hello hello`.
- the same, exercising `f.deflate` into `f.inflate` -- prints
  `async roundtrip: async round trip payload`. This is the e2e proof that
  the node worker path survives packaging.
- `ls node_modules/@substrate-system/fflate/dist/browser/index.d.ts
  node_modules/@substrate-system/fflate/dist/node/index.d.ts` -- both
  exist.

### AC4.2

> **dist-build-restructure-design.AC4.2 Success:** Importing under any
> other condition resolves to `dist/browser/index.js`, and its types to
> `dist/browser/index.d.ts`.

Automated but **partially**: build-verification plus static inspection of
the `exports` map. This is the weakest of the AC4 criteria.

What is actually run:

- the file-existence loop in Phase 4 Task 1 confirms both target files
  exist.
- `node -e "const e=require('./package.json').exports; for (const k of
  ['.','./node','./browser','./min','./umd','./package.json']) if(!(k in
  e)) console.log('MISSING KEY',k)"` -- no output.
- Phase 7 Task 3's `import.meta.resolve('@substrate-system/fflate' +
  '/browser')` returns a path.

What is **not** run: no command resolves the bare `.` entry under a
non-node condition. Node's own resolver always supplies `node`, so the
`default` branch of the `.` key cannot be exercised from a Node process.
The effective verification is inspection of the `exports` map plus proof
that `./browser` resolves to the same file the `default` branch names.

That gap is acceptable and needs no human step, because the browser
default branch is exercised in substance by AC6: `test/browser/index.ts`
imports `../../dist/browser/index.js` directly and the whole browser
suite runs against it in a real browser. What is unverified is the
`exports` map wiring, not the artifact.

### AC4.3

> **dist-build-restructure-design.AC4.3 Success:** The `./node`,
> `./browser`, `./min`, `./umd`, and `./package.json` subpaths all
> resolve to existing files.

Automated, build-verification (Phase 4 Task 1) plus e2e (Phase 7 Task 3).

- Phase 4 Task 1: the `MISSING KEY` loop over all six keys, plus the
  six-path file-existence loop.
- Phase 7 Task 3, against the installed tarball:

```
node --input-type=module -e "
const subs = ['/node', '/browser', '/min', '/umd']
Promise.all(subs.map(s =>
  import.meta.resolve('@substrate-system/fflate' + s)
)).then(rs => rs.forEach(r => console.log('ok', r)))
"
```

Expected four `ok` lines. Any rejection means a subpath in the `exports`
map points at a file that was not packed. This is the check that catches
a mismatch between `exports` and `files`, which the working-tree checks
structurally cannot.

### AC4.4

> **dist-build-restructure-design.AC4.4 Success:** `npm pack` includes
> only `dist/` and the always-included metadata files.

Automated, build-verification (Phase 4 Task 3) plus e2e (Phase 7 Task 3).

Phase 4 Task 3:

- `npm pack --dry-run 2>&1 | grep -E "^npm notice" | grep -vE "package
  size|unpacked size|shasum|integrity|total files|filename|=== |name:
  |version:"` -- every listed file under `dist/`, plus `package.json`,
  `README.md`, `LICENSE`.
- the `npm pack --dry-run --json` pipeline testing every path against
  `/^(dist\/|package\.json$|README\.md$|LICENSE$)/` -- prints
  `contents ok`. This is the machine-checkable form and the one to rely
  on; the grep above is for human reading.
- `git check-ignore -q public/ && git check-ignore -q dist/ && echo
  "both ignored"` -- prints `both ignored`. Both details matter. `-q`
  accepts only ONE pathname (`fatal: --quiet is only valid with a
  single pathname`), so the two-argument form always exits 128 and the
  `&&` never fires. And the trailing slashes are required: `.gitignore`
  writes the patterns as `public/` and `dist/`, so `check-ignore -q
  public` exits 1 while `check-ignore -q public/` exits 0.

Phase 7 Task 3, against the installed tarball:

- `find node_modules/@substrate-system/fflate -type d -maxdepth 1 | sort`
  -- the package root and `dist` only. No `src`, `test`, `example`,
  `scripts`, `lib`, `esm`, or `umd`.

### AC4.5

> **dist-build-restructure-design.AC4.5 Failure:** No `require`
> condition is advertised anywhere in the `exports` map.

Automated, static inspection of `package.json`. **Negative criterion,
verified by string absence.**

Phase 4 Task 1:

```
node -e "const e=require('./package.json').exports;
console.log(JSON.stringify(e).includes('require')
  ? 'HAS REQUIRE' : 'no require condition')"
```

Expected `no require condition`.

Phase 7 Task 3 re-asserts it against the **packed artifact** rather than
the working tree, reading
`require('@substrate-system/fflate/package.json').exports`.

Phase 7 Task 3 records an important non-test: **do not assert that
`require()` throws.** Node has supported `require()` of ESM unflagged
since 22.12, and this tree runs 25.8.2. Measured against the exact
`exports` map from Phase 4, `require()` resolves conditions
`["require", "node", ...]`, matches the `node` key, falls through to
`default`, and loads `dist/node/index.js` as ESM successfully. A test
expecting a throw would fail on modern Node while the package is
correct. The criterion is about what is advertised, and that is what is
tested.

---

## AC5: The node suite passes against the new build

**Covered by the existing node suite, not by new tests.** Phase 5 changes
how `test/util.ts` loads fflate and which suites `test/index.ts` imports.
It writes no new assertions. The behavioural coverage comes from
`test/0-valid.ts`, `test/1-size.ts`, and `test/2-perf.ts`, which already
existed and are unchanged.

### AC5.1

> **dist-build-restructure-design.AC5.1 Success:** `npm run test:node`
> exits 0 with every assertion passing, loading fflate from
> `dist/node/index.js`.

Automated, existing behavioural suite plus build-verification of the
wiring. Produced by Phase 5 Task 1, gated by Phase 5 Task 3 and Phase 7
Task 4.

Run command: `npm run test:node` (`tsx test/index.ts`).

- `npm run build && npm run test:node` (Phase 5 Task 3) -- exits 0, TAP
  output ends with a plan line and no `not ok` entries.
- `npm run test:node 2>&1 | grep -c "^not ok" || true` -- `0`.
- `npm run test:node 2>&1 | grep -cE
  "ERR_UNSUPPORTED_RESOLVE_REQUEST|Cannot find module|
  ERR_MODULE_NOT_FOUND" || true` -- `0`. Any match means a worker failed
  to load its module, which is the specific failure mode the loader
  rework risks.

Wiring checks (Phase 5 Task 1):

- `grep -n "lib', 'index.cjs\|lib/index.cjs" test/util.ts` -- no output.
- `grep -n "pathToFileURL" test/util.ts` -- two matches, the import and
  the fflate constant.
- `npx tsx test/0-valid.ts` -- exits 0. This is the cheap early gate: it
  exercises the fflate worker (dynamic `import()` of the ESM bundle)
  against the zlib worker (`require` of a builtin), proving both loader
  paths in one run.

Recorded correction: the design proposed rewriting `wc()` to spawn
`data:text/javascript` module workers. Phase 5 investigation finding 1
measured that against Node 25.8.2 and it fails --
`ERR_UNSUPPORTED_RESOLVE_REQUEST`, because Node's ESM resolver does not
resolve bare specifiers relative to a `data:` URL, and every comparison
library is loaded by bare specifier. The implemented approach keeps
`{ eval: true }` CommonJS workers and changes only how fflate is loaded.

### AC5.2

> **dist-build-restructure-design.AC5.2 Success:** The benchmark workers
> spawn for all eight fflate methods defined in `test/util.ts`
> (`deflate`, `inflate`, `gzip`, `gunzip`, `zlib`, `unzlib`, `zip`,
> `unzip`) and return timings for the six that `test/2-perf.ts` runs,
> plus the pako, uzip, tiny-inflate, and zlib comparisons.

Reworded during the test-coverage analysis. The original said "spawn and
return results for all eight", which the tree falsifies:
`test/2-perf.ts:26` is `if (l == 'zip' || l == 'unzip') continue`, so
`timings.json` carries six fflate keys -- `deflate`, `inflate`, `gzip`,
`gunzip`, `zlib`, `unzlib`. Measured.

The skip is pre-existing upstream behaviour that this branch does not
touch, and for `unzip` it is structural rather than an oversight:
`unzip` returns an `Unzipped` plain object with no `.buffer`, so the
harness's fixed `postMessage(buf, [buf.buffer])` transfer list rejects
it with `DataCloneError`. `test/util.ts:210-219` documents exactly
that. Building a benchmark worker that could carry it is out of scope
-- the design's Tests section preserves existing coverage rather than
extending it -- and the ZIP behaviour this criterion gestures at is
covered far more strongly by `test/browser/zip.ts`, which runs 2 x 15
ZIP tests including worker-branch assertions.

Automated, existing behavioural suite. Gated by Phase 5 Task 3.

- `test -f test/results/timings.json && test -f
  test/results/longTimings.json && echo "results written"` -- prints
  `results written`, confirming the benchmark suites ran to completion
  rather than being skipped.
- `node -e "const t=require('./test/results/timings.json'); const
  libs=new Set(Object.keys(t).map(k=>k.split('.')[0]));
  console.log([...libs].sort().join(','))"` -- expected
  `fflate,pako,tinyInflate,uzip,zlib`. A missing library means its worker
  never produced a result.

**Partial coverage, flagged.** The library-set assertion is exact, but no
command enumerates the fflate methods individually; the criterion's
method list is covered only insofar as `0-valid.ts` and the benchmark
suites happen to call them. Strengthening this would mean asserting on
the second segment of the `timings.json` keys as well as the first --
which, given the amended wording above, means asserting the six that
run rather than the eight that are defined. The plan does not do so,
and this is noted rather than invented.

Also recorded: the design lists jszip as a comparison library, but the
harness has **no jszip worker** -- `test/util.ts` defines workers for
fflate, pako, uzip, tinyInflate, and zlib only. Phase 5 deliberately does
not add one, since the design's Tests section describes preserving
existing coverage rather than extending it. The expected library list
above reflects the harness as it is.

### AC5.3

> **dist-build-restructure-design.AC5.3 Success:** `test/index.ts`
> imports no ZIP, stream, or one-shot async suite. Those moved to
> `test/browser/`, where they run against the real browser build
> including the minified bundle.

Automated, static source verification. Produced by Phase 5 Task 2.

The criterion is about WHICH suites live on the node side, not about
the count of numbered files. It was originally worded as "imports only
`./0-valid.js`, `./1-size.js`, and `./2-perf.js`", which the final
review falsified: `test/3-node-min.ts` was added afterwards to cover
`dist/node/index.min.js`, the one shipped artifact in the `wcln`
silent-failure class that the browser-only suite cannot reach. That is
a legitimate node-side test -- it exercises a node bundle through
`worker_threads` -- so the AC is widened rather than the test removed.

- `grep -cE "3-zip|4-streams|5-async" test/index.ts || true` -- `0`.
  This is the load-bearing check: the three suites that moved must not
  come back.
- `ls test` -- `0-valid.ts`, `1-size.ts`, `2-perf.ts`, `3-node-min.ts`,
  `browser`, `data`, `index.ts`, `results`, `util.ts`.
- `ls test/browser` -- `async.ts`, `harness.ts`, `index.ts`, `min.ts`,
  `streams.ts`, `util.ts`, `zip.ts`.

---

## AC6: The browser suite runs and passes

Produced by Phase 6 Task 6. Re-verified by Phase 7 Task 4. Shares its
test files with AC3; see the AC3 test architecture section above.

### AC6.1

> **dist-build-restructure-design.AC6.1 Success:** `npm run test:browser`
> exits 0, running the ZIP, stream, and async suites against
> `dist/browser/index.js`.

Automated, behavioural. Run command: `npm run test:browser`
(`esbuild test/browser/index.ts --bundle | tapout --timeout 30000`).

- `npm run build && npm run test:browser; echo "exit=$?"` -- `exit=0`,
  with TAP output containing both `plain` and `min` prefixed test names
  and no `not ok` lines.
- `npm run test:browser 2>&1 | grep -c "^not ok" || true` -- `0`.

Setup dependency: `@substrate-system/tapout` pulls in Playwright, whose
browsers are a separate download, so `npx playwright install chromium`
is a required setup step and CI must account for it (Phase 6 Task 1).

Phase 6 Task 1 also verifies the **failure path** before any suite is
written, by bundling a throwaway `test/browser/_fail.ts` containing
`t.ok(false, ...)` and asserting `exit=1` and a `not ok` line. A runner
that cannot fail is worse than no runner. Both probe files
(`_smoke.ts`, `_fail.ts`) are deleted before Task 2 and must not be
committed. They must live inside the repository, not the scratchpad,
because esbuild resolves bare specifiers by walking up from the entry
file.

### AC6.2

> **dist-build-restructure-design.AC6.2 Success:** The same suites pass
> against `dist/browser/index.min.js` in the same run.

Automated, behavioural. The critical assertion is that the minified pass
actually **ran to completion**, not merely that nothing failed:

- `npm run test:browser 2>&1 | grep -c "^not ok"` -- expected `0`. Any
  failing assertion, whether in `plain` or `min`, indicates a defect.
- `npm run test:browser 2>&1 | grep -c "^# min >"` -- expected to equal
  `npm run test:browser 2>&1 | grep -c "^# plain >"`. Both should be
  `103` on the healthy tree. If the min count is less, the run was
  truncated before the min pass completed, indicating a timeout or crash
  in the worker bootstrap. This check closes both the assertion hole
  (a label appears only on comment lines, not on `not ok` lines) and the
  truncation hole documented in phase_04.md.
- `npm test; echo "exit=$?"` -- `exit=0`. This is the first point in the
  plan where the full chain (build, node suite, browser suite) works.

---

## AC7: The example app builds and deploys from public/

Produced by Phase 1. AC7.1 and AC7.2 are build-verification via
`npm run build-example`. AC7.4 is static inspection. **AC7.3 has the
weakest automated coverage in the plan and requires a human step.**

### AC7.1

> **dist-build-restructure-design.AC7.1 Success:** `npm run build-example`
> exits 0 and writes `index.html` plus hashed assets into `public/` at the
> repository root.

Automated, build-verification. Produced by Phase 1 Tasks 4, 5, and 8.

- `npm run build-example` (Phase 1 Task 8) -- exits 0, Vite reports the
  example root and writes to `public`, with no unresolved-import errors
  for `react`, `react-dom/client`, or `fflate`.
- `ls public` -- contains `index.html`, an `assets/` directory,
  `favicon.ico`, and `sw.js`.
- `git status --porcelain public` -- no output, confirming `public/` is
  git-ignored.
- `npm run build-example && ls public/index.html && echo "example ok"`
  (Phase 7 Task 4) -- prints the path and `example ok`.

Supporting static checks: `grep -rn "from '\.\./\.\./\.\.'" example/`
gives no output and `grep -rn "from 'fflate'" example/` gives two matches
(Phase 1 Task 5), confirming the example builds through the `fflate`
alias to `src/index.ts` rather than through the root `package.json`,
which is what decouples this phase from Phases 3 and 4.

### AC7.2

> **dist-build-restructure-design.AC7.2 Success:** The example's React
> import specifiers, including the `react-dom/client` subpath used by
> `example/index.tsx`, resolve to `preact/compat` at build time with
> neither `react` nor `react-dom` installed.

Automated, build-verification. Produced by Phase 1 Tasks 1 and 4.

This is verified by the **success** of `npm run build-example`: an
unresolved `react` or `react-dom/client` is a hard Vite build error, so a
zero exit is proof of resolution. There is no separate assertion, and
none is needed.

Two supporting checks:

- `node -e "try{require.resolve('parcel');console.log('STILL
  PRESENT')}catch(e){console.log('removed')}"` (Phase 1 Task 1) -- prints
  `removed`.
- `git grep -n "from 'react'" -- example | wc -l` (Phase 7 Task 2) --
  greater than `0`. This confirms the React specifiers were
  intentionally kept and resolved through `preact/compat`, rather than
  rewritten to `preact` imports by mistake. Without it, a phase that
  "fixed" the imports would pass every other check while defeating the
  criterion.

Type-level coverage comes from AC8.1: the root `tsconfig.json` `paths`
block maps `react-dom/client` to `preact/compat/client`, and a
`tsc --noEmit` with zero errors proves that mapping resolves. Mapping to
`preact/compat` instead produces `TS2305: Module '"react"' has no
exported member 'createRoot'`, so this is a real check, not a formality.

### AC7.3

> **dist-build-restructure-design.AC7.3 Success:** The built service
> worker precaches the application assets using the manifest injected by
> `vite-plugin-pwa` rather than `@parcel/service-worker`.

**Split coverage: partly automated, partly human. This is the weakest
automated coverage in the plan.**

Automated portion (Phase 1 Tasks 6 and 8):

- `grep -n "@parcel/service-worker" example/sw.ts` -- no output.
- `grep -c "sw\." example/sw.ts || true` -- `0`, confirming the four
  `sw.` references were converted to `self.`.
- `npx tsc --noEmit --ignoreConfig --lib ES2022,DOM,WebWorker
  --skipLibCheck example/sw.ts` -- exits 0.
- `grep -c "__WB_MANIFEST" public/sw.js || true` -- **`0`**. A remaining
  token means injection did not run.
- `grep -o "precacheVersion" public/sw.js | head -1` -- either a match or
  no output; the phase explicitly calls this a smoke check only, since
  the identifier may be minified.

**What the automation does not establish.** The only real assertion is
that the `self.__WB_MANIFEST` token is absent from the built `sw.js`. An
absent token proves the plugin rewrote something. It does not prove the
manifest contains the right entries, that the derived `precacheVersion`
cache name is stable across builds and changes when the asset set
changes, that the service worker installs without throwing, or that a
subsequent navigation is served from the cache. Nothing in the plan
registers the service worker in a browser.

**Human verification required.**

Justification: verifying service worker install and cache-serving needs a
real browser with a secure context, a full navigation, and an offline
transition. The repository has a browser runner (`@substrate-system/
tapout` over Playwright), but it drives a bundled test module against
`about:blank`-style pages, not a served origin with a registered service
worker. Wiring a served-origin Playwright fixture is a larger piece of
work than the phase scopes, and the design does not ask for it.

Steps:

1. `npm run build-example`
2. Serve `public/` over `http://localhost` with any static server (a
   secure context; `localhost` qualifies without TLS).
3. Open the site in Chrome. In DevTools, Application -> Service Workers,
   confirm one worker is registered, its status is `activated`, and its
   script URL ends in `/sw.js`.
4. In Application -> Cache Storage, confirm a cache exists whose name
   begins `fflate-`, and that it contains `index.html` and the hashed
   `assets/` entries. Confirm `favicon.ico` is absent. Note the reason
   changed during the migration: the old `demo/sw.ts` filtered `.ico`
   out of `precacheFiles` explicitly, but `example/sw.ts:26` is now a
   bare `manifest.map(e => e.url)`, and the exclusion comes from
   vite-plugin-pwa's default `injectManifest` glob instead. So this
   step confirms the plugin default is in effect, not that a filter in
   this repo's own code applied. If `favicon.ico` ever appears here,
   the fix is a `globIgnores` entry or a filter in `example/sw.ts`, not
   a hunt for a filter that no longer exists.
5. Reload with DevTools Network throttling set to Offline. The page must
   still render. This is the assertion that the precache is real rather
   than merely populated.
6. Note the cache name. Rebuild with `npm run build-example` and reload.
   The cache name must be unchanged if no source changed. Then edit any
   file under `example/`, rebuild, reload, and confirm the cache name
   changed. This exercises the property the removed Parcel `version`
   export used to provide, and is the part most likely to be wrong.

Record the result in the PR description. If step 5 or step 6 fails, the
defect is in the `precacheVersion` derivation in `example/sw.ts`, not in
the Vite configuration.

### AC7.4

> **dist-build-restructure-design.AC7.4 Success:** `scripts/cpGHPages.ts`
> reads the directory that the example build writes (`public/`), not
> `dist/`.

Automated, static source verification. Produced by Phase 1 Task 7.

The script is deliberately **not executed** as a verification step: it
checks out the `gh-pages` branch and commits, so running it mutates
repository state. Verification is static:

- `grep -n "'dist'" scripts/cpGHPages.ts` -- no output.
- `grep -n "__dirname" scripts/cpGHPages.ts` -- no output (the package is
  `"type": "module"`).
- `npx tsc --noEmit --ignoreConfig --module ES2022 --target ES2022
  --moduleResolution Bundler --types node --skipLibCheck
  scripts/cpGHPages.ts` -- exits 0 with no diagnostics.

Not covered by automation, and accepted as such: that the script actually
produces a correct `gh-pages` commit. Task 7 also fixes two defects
beyond the output path -- `statSync(f)` becoming `statSync(to(f))`, and
the unconditional `git checkout('master')` becoming a return to the
branch the run started on. Neither fix is exercised by any command. The
first time the script runs for real, confirm by hand that the `gh-pages`
branch received the example build and that the working branch was
restored afterwards.

---

## AC8: Type checking passes and obsolete configuration is gone

AC8.1 and AC8.3 are produced by Phase 2; AC8.2 by Phase 7.

### AC8.1

> **dist-build-restructure-design.AC8.1 Success:** `npx tsc --noEmit`
> against the root `tsconfig.json` reports zero errors across `example`,
> `src`, `test`, and `scripts`.

Automated, build-verification. Produced by Phase 2 Tasks 1 and 3B, gated
by Phase 2 Task 4 and re-verified by Phase 7 Task 4.

- `npx tsc --noEmit --project tsconfig.json 2>&1 | grep -cE "error TS" ||
  true` -- `0`.
- `npx tsc --noEmit --project tsconfig.build.json 2>&1 | grep -cE "error
  TS" || true` -- `0`. A `TS5011` here means `rootDir` is missing from
  `tsconfig.build.json`.
- `npx tsc --emitDeclarationOnly --project tsconfig.build.json && ls
  dist/index.d.ts && rm -rf dist` -- `dist/index.d.ts` exists. This is a
  dry run of the emit Phase 3 Task 3 depends on; a file at
  `dist/src/index.d.ts` instead means `rootDir` is not taking effect.

Two conventions matter for anyone re-running this:

1. The root config sets `listFiles: true`, so `tsc` prints every file it
   loads. An unfiltered grep for a filename matches the file listing
   itself. Every check must filter on `error TS` first. Phase 3 Task 2
   calls this out explicitly.
2. **From Phase 6 onward, a bare `tsc --noEmit` fails with `TS2307`
   unless a build has run**, because `test/browser/index.ts`
   value-imports `../../dist/browser/index.js`. Every gate checking AC8.1
   after Phase 6 must run `npm run build` first, including Phase 7 Task
   4. Deriving the `Fflate` type from `src/` removes the build dependency
   for type resolution but not for the value import.

Phase 2 Task 4 enumerates the expected error categories and their causes
(`TS2688` for a missing Vite install, `TS17004`/`TS6142` for a missing
`jsx` option, `TS2307` for a missing `paths` block, `TS2305` for
`react-dom/client` mapped to the wrong preact entry, `TS7016`/`TS2322`/
`TS2698` for incomplete Task 3B). It also states the rule that errors in
`src/index.ts` are unexpected and must be reported rather than silenced:
`src/` is vendored upstream code and this repository merges from
`git remote upstream`.

### AC8.2

> **dist-build-restructure-design.AC8.2 Success:** No file in the
> repository references `scripts/rewriteBuilds.ts`,
> `scripts/buildUMD.ts`, the `lib/`, `esm/`, or top level `umd/`
> directories, or the `SC=` script runner.

Automated, static repository sweep. **Negative criterion, verified by
grep for absence.** Produced by Phase 7 Tasks 2 and 2B.

Six `git grep` commands, each expected to produce no output:

```
git grep -n "rewriteBuilds\|buildUMD" -- . ':!docs' ':!specs'
git grep -n "SC=" -- . ':!docs' ':!specs'
git grep -nP "\blib/index\.cjs\b|\besm/index\.mjs\b|\besm/browser\.js\b" \
  -- . ':!docs' ':!specs' ':!README.md'
git grep -nP "tsconfig\.(esm|demo)\.json" -- . ':!docs' ':!specs'
git grep -n "@parcel/service-worker\|parcel " \
  -- . ':!docs' ':!specs' ':!package-lock.json'
git grep -n "build:lib\|build:umd\|build:rewrite\|build:demo\|\
build-cjs\|build-esm" -- . ':!docs' ':!specs'
```

Plus `ls scripts` -- exactly `build.ts` and `cpGHPages.ts`.

Two conventions are load-bearing here and are stated in the Phase 1
verification-conventions section:

- **`git grep -E` does not honour `\b`.** Where a word boundary is
  needed, use `git grep -P`. A `-E` pattern containing `\b` silently
  matches nothing and reports a false pass. Two of the six commands above
  use `-P` for exactly this reason.
- Path exclusions are deliberate, not laziness. `docs/` is generated
  typedoc output and out of scope. `specs/` holds the design document,
  which describes the very things being removed, so five of the six
  commands would hit it and "fixing the reference" would mean editing the
  source of truth. `README.md` is excluded because Task 2B rewrites it
  immediately after.

Phase 7 Task 2B covers the README, which documents seven entry points
this work removes:

- `grep -nE "require\('fflate'\)|umd/index\.js|lib/index\.d\.ts|
  fflate/esm|esm/browser\.js|esm/index\.mjs" README.md || true` -- no
  output. Note `fflate/esm` and not only `esm/browser.js`: line 108
  imports the bare `fflate/esm` subpath and the narrower pattern misses
  it.
- `grep -c "fflate@0.8.3" README.md || true` -- `0`.
- `grep -c "dist/umd/fflate.js" README.md || true` -- at least `1`.
- `npm run toc` -- exits 0.

Task 2B also requires a prose note recording that dropping CommonJS is
breaking for `require()` consumers. That note's presence is not asserted
by any grep; reviewer inspection covers it.

Phase 3 Task 7 contributes the complementary check at deletion time:
`grep -rn "rewriteBuilds\|buildUMD" --include='*.ts' --include='*.js' .
--exclude-dir=node_modules --exclude-dir=docs --exclude-dir=dist` -- no
output. Phase 4 Task 2 contributes the script-block check, asserting all
twelve superseded npm scripts are gone and printing `all removed`.

Final state check, Phase 7 Task 4: `git status --porcelain` -- no entries
for `dist/`, `public/`, `lib/`, `esm/`, or `umd/`.

### AC8.3

> **dist-build-restructure-design.AC8.3 Success:** `tsconfig.esm.json`,
> `tsconfig.demo.json`, and `test/tsconfig.json` no longer exist, and no
> file references them.

Automated, static verification. **Negative criterion.** Produced by Phase
2 Task 2.

- `ls tsconfig.esm.json tsconfig.demo.json test/tsconfig.json 2>&1` --
  three "No such file or directory" errors.
- `grep -rn "tsconfig.esm\|tsconfig.demo\|test/tsconfig" --include='*
  .json' --include='*.ts' --include='*.js' . --exclude-dir=node_modules
  --exclude-dir=docs` -- no output.

Phase 2 Task 2 records one expected transient: if `package.json` still
references `tsconfig.esm.json` through the `build:lib` script at that
point, the reference is removed in Phase 4 Task 2, which replaces the
whole script block. Record it and continue rather than treating it as a
failure. Phase 7 Task 2's `tsconfig\.(esm|demo)\.json` sweep is the
final assertion after Phase 4 has landed.

---

## Coverage summary

| AC | type | where |
| --- | --- | --- |
| AC1.1 | build-verification | Phase 3 T3, T4, T8 |
| AC1.2 | build-verification | Phase 3 T5, T8 |
| AC1.3 | build-verification | Phase 3 T6 |
| AC1.4 | build-verification (absence) | Phase 3 T7, T8; Phase 7 T1 |
| AC1.5 | build-verification | Phase 3 T1, T8 |
| AC2.1 | build-verification | Phase 3 T4 |
| AC2.2 | static source | Phase 2 T3 |
| AC2.3 | grep for absence | Phase 3 T4, T6 |
| AC3.1 | behavioural | `test/browser/async.ts`, Phase 6 T3, T6 |
| AC3.2 | behavioural | `test/browser/streams.ts`, Phase 6 T4, T6 |
| AC3.3 | behavioural | `test/browser/zip.ts`, Phase 6 T5, T6 |
| AC3.4 | behavioural (differential) | all three suites, Phase 6 T6 |
| AC4.1 | build-verification + e2e | Phase 4 T1; Phase 7 T3 |
| AC4.2 | inspection + build-verification | Phase 4 T1; Phase 7 T3 |
| AC4.3 | build-verification + e2e | Phase 4 T1; Phase 7 T3 |
| AC4.4 | build-verification + e2e | Phase 4 T3; Phase 7 T3 |
| AC4.5 | inspection for absence | Phase 4 T1; Phase 7 T3 |
| AC5.1 | existing node suite | Phase 5 T1, T3; Phase 7 T4 |
| AC5.2 | existing node suite (partial) | Phase 5 T3 |
| AC5.3 | static source | Phase 5 T2 |
| AC6.1 | behavioural | Phase 6 T6; Phase 7 T4 |
| AC6.2 | behavioural | Phase 6 T6; Phase 7 T4 |
| AC7.1 | build-verification | Phase 1 T8; Phase 7 T4 |
| AC7.2 | build-verification | Phase 1 T1, T4, T8; Phase 7 T2 |
| AC7.3 | grep only, **plus human** | Phase 1 T6, T8 |
| AC7.4 | static source | Phase 1 T7 |
| AC8.1 | build-verification | Phase 2 T1, T3B, T4; Phase 7 T4 |
| AC8.2 | grep sweep for absence | Phase 7 T2, T2B |
| AC8.3 | static, absence | Phase 2 T2; Phase 7 T2 |

Human verification is required for exactly one criterion, AC7.3, and for
the first real run of `scripts/cpGHPages.ts` under AC7.4.

The negative criteria -- AC1.4, AC2.3, AC4.5, AC8.2, AC8.3 -- are all
verified by grep or file-absence rather than by executing anything. Every
one of them needs the `|| true` and `-P` conventions from the Phase 1
verification-conventions section, or it will report a false result.

---

## Deliberately not tested

Drawn from the design's "Out of scope" section and reaffirmed in Phase 7.

**No linting of `src/`.** fflate has no `eslint.config.js` and none is
added, so there is no `lint` script and no `preversion` hook. This is a
deliberate deviation from `template-ts-browser`. `src/index.ts` is 3,885
lines of vendored upstream code: 192 lines exceed 80 columns and 1,289
end in a semicolon, against a template config that sets 4-space indent,
no-space `key-spacing`, and standard's no-semicolon style. Autofixing
would rewrite substantially every line, and this repository merges from
`git remote upstream` at `101arrowz/fflate`, so a wholesale reformat
would make every future merge a conflict. `newneostandard`, `eslint`, and
`typescript-eslint` are not added as devDependencies, so there is no
lint check to include in any gate.

**No TypeScript `strict` mode, and no test asserting it.** The design
measured 180 errors under `strict: true` against `src/**/*`: 54 from
`strictPropertyInitialization`, roughly 101 from `strictNullChecks` (the
codebase passes `null` into callback error slots throughout), and 25 from
`strictFunctionTypes`, `noImplicitThis`, and
`useUnknownInCatchVariables`. Every one is in vendored upstream code. The
chosen configuration is `strict: false, noImplicitAny: true`, which
measures zero errors, and AC8.1 is written against that. Raising
strictness is out of scope for the same upstream-merge reason as linting.

**No typedoc changes and no assertion about `docs/`.** `docs/` is tracked
typedoc output. The `build-docs` script and its typedoc invocation are
unchanged apart from the key rename from `build:docs`. Phase 7 explicitly
instructs against regenerating `docs/`, since that would produce a large
unrelated diff, and `docs/` is excluded from every reference sweep
because it will legitimately contain historical references to the removed
paths.

**No change to the compression algorithms, and no unit tests for them.**
The plan adds no unit tests for deflate, inflate, gzip, zlib, or the ZIP
container logic as algorithms. They are vendored upstream and unchanged
by this work, and the existing `test/0-valid.ts` already validates
correctness against zlib, pako, uzip, and tiny-inflate. The new browser
suite tests the algorithms only incidentally: its assertions are about
whether the *build* preserved them, which is why every AC3 assertion is
framed as async-equals-sync or round-trip-equals-input rather than as a
check against known compressed output. A test asserting specific byte
output for a given input would fail on any legitimate upstream
compression change and would be testing the wrong thing.

**No test that the UMD artifact loads under Node's `require()`.** Phase 3
Task 6 records this as a real limitation of the design's chosen filename
rather than a defect: `dist/umd/fflate.js` sits inside a
`"type": "module"` package, so Node parses it as ESM. It is reachable
from a script tag, from AMD loaders, and from legacy bundlers, which is
the role the design assigns it. Node consumers use the `node` export
condition. Renaming to `.cjs` to make a `require` test pass is explicitly
ruled out -- the design names the file and `unpkg` points at it.

**No test that `require('@substrate-system/fflate')` fails.** See AC4.5.
Node has supported `require()` of ESM unflagged since 22.12 and it
succeeds against this package's `exports` map. The criterion is about
what is advertised, not what the runtime happens to permit.

**No jszip benchmark worker.** The design's Tests section lists jszip
among the comparison libraries, but `test/util.ts` has never defined a
jszip worker. Phase 5 does not add one, since the design describes
preserving existing coverage rather than extending it.

**No served-origin service worker test.** See AC7.3. This is the one
place where the absence of automation leaves a real gap rather than a
justified omission, and it is why AC7.3 carries a human step.
