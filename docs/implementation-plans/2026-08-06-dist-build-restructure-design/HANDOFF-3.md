# Handoff 3: dist build restructure, all phases executed

Written 2026-08-07. Supersedes `HANDOFF.md` and `HANDOFF-2.md` for
anything they disagree on.

**All seven phases are executed and committed. Phase 6 is approved
after six review cycles. What remains is Phase 7's review and the
four closing steps.**

## Current state

- Branch `dist-build-restructure`, HEAD `6031716`. No upstream; the
  branch is local only. Working tree clean.
- `npm test` exits 0: node 28/28, browser 306/306, zero `not ok`,
  103 `# plain >` and 103 `# min >` headers, plans `1..28` and
  `1..306`.
- `npm run build` exits 0. `npx tsc --noEmit -p tsconfig.json` reports
  0 errors, **but only after `npm run build`** -- see gotcha 1.
- `npm run build-example` exits 0 and produces `public/index.html`.
- Phase 7's own gates all pass; see the `6031716` commit message for
  the measurements.

## How to resume

Invoke `ed3d-plan-and-execute:executing-an-implementation-plan` with:

- plan directory:
  `/Users/nick/code/fflate/docs/implementation-plans/2026-08-06-dist-build-restructure-design/`
- working directory: `/Users/nick/code/fflate/`

Skip phase discovery and task creation for phases 1 to 7; every task
is done. Resume at Phase 7's code review.

There is no `.ed3d/implementation-plan-guidance.md`, so do not pass an
IMPLEMENTATION_GUIDANCE path. `test-requirements.md` DOES exist in the
plan directory and is used by the final test-analyst step.

## Remaining work

### 1. Phase 7 code review

Phase 7's execution is committed but not reviewed. Review
`e394490..6031716`, which is the README rewrite (Task 2B), the
`example/sw.ts` comment from the Task 2 sweep, and the Task 1/3/4
verification plus two `phase_07.md` amendments.

Nothing in Phase 7 changed library code. The risk surface is the
README's accuracy about the new entry points and whether the
verification actually exercised the packed artifact rather than the
working tree.

### 2. Invoke `ed3d-extending-claude:project-claude-librarian`

No `CLAUDE.md` or `AGENTS.md` exists anywhere in this repo, so it may
report nothing.

### 3. Final code review over `0146a83..HEAD`

With an AC coverage check: verify every
`dist-build-restructure-design.AC*` in the design is covered by some
phase. Loop to zero issues.

### 4. `ed3d-plan-and-execute:test-analyst`

Only after step 3 reaches zero issues. Dispatch with
`TEST_REQUIREMENTS_PATH` set to `test-requirements.md` in the plan
directory. On PASS, write the returned human test plan to
`docs/test-plans/2026-08-06-dist-build-restructure-design.md` and
commit it.

### 5. `finishing-a-development-branch`

Not before step 4.

## Gotchas that cost real time

All measured, not assumed. Items 1 to 11 carry forward from
`HANDOFF-2.md`; 12 to 17 are new this session.

1. **`tsc --noEmit` REQUIRES a prior `npm run build`.**
   `test/browser/index.ts` value-imports `../../dist/browser/index.js`,
   so a bare typecheck fails with `TS2307`. This affects every AC8.1
   gate.

2. **Never verify with `node --input-type=module -e`.** That form puts
   `--input-type=module` into `process.execArgv`; worker threads
   inherit `execArgv` and then parse their eval'd payload as ESM, so
   fflate's sloppy-mode implicit globals throw `ReferenceError: u8 is
   not defined`. It is an artifact of the command, not a real defect.
   Use a real `.mjs`/`.cjs`/`.mts` file and assert the exit code.
   `phase_07.md` Task 3 originally prescribed this form for four
   checks, one of which spawns a worker; it was rewritten in
   `6031716`.

3. **`keepNames` must stay `false` in `scripts/build.ts`.** esbuild's
   `keepNames` wraps functions as `__name(fn, 'name')`; the async APIs
   stringify functions and eval them in a worker where that helper is
   out of scope, so every async API dies with `ReferenceError: __name
   is not defined` while sync APIs keep working.

4. **`src/index.ts` is vendored from upstream but WAS deliberately
   edited**, on the user's explicit instruction, to remove the
   `Deflate.call(this, ...)` idiom. That is sanctioned. It WILL
   conflict on future upstream merges. `src/worker.ts` remains
   untouched; `src/node-worker.ts` is editable.

5. **Do not run `scripts/cpGHPages.ts`, `npm run gh-pages`, or
   `npm run build-docs`.** The first checks out `gh-pages`, deletes
   every top-level file and commits. The last writes into `docs/`,
   where this plan lives.

6. **`build-example` must keep the trailing slash: `--base="/fflate/"`.**
   Vite does not normalise `import.meta.env.BASE_URL`, so `/fflate`
   makes the service worker register at `/fflatesw.js`, which 404s.
   Verified in `6031716`: the built bundle registers `/fflate/sw.js`.

7. **`./umd` deliberately carries NO `types` entry.** Adding one type
   checks code that fails at runtime.

8. **Async ZIP APIs transfer their input buffer to the worker, which
   detaches it.** `test/browser/zip.ts` regenerates fixtures per
   access (`const fx = () => fixtures()`) and pushes `.slice()`
   copies. A shared fixture array gets emptied for every later test.

9. **`Zip.add(file)` installs the entry's `ondata` handler.** Pushing
   before adding fails with `no stream handler`.

10. **`UnzipFile.start()` instantiates the REGISTERED decoder.**
    Calling it without `unzipper.register(...)` throws `ctr is not a
    constructor`; not calling it at all means no data flows and the
    test hangs. See gotcha 13 for what a hang now does.

11. **Streamed output does not equal one-shot output.** Chunking
    changes DEFLATE block boundaries; on incompressible input even the
    SYNC stream differs from one-shot (561 vs 535 bytes, measured).
    Compare async stream against SYNC STREAM at identical chunking.

12. **`dist/umd/package.json` is load-bearing.** It contains
    `{ "type": "commonjs" }`. The root `package.json` is
    `"type": "module"`, so without it node parses `dist/umd/fflate.js`
    as ESM, where the UMD wrapper's `this` is undefined and loading
    throws "Cannot set properties of undefined". Do not delete it as
    build residue.

13. **`withTimeout`'s default and tapout's `--timeout` are a pair.**
    tapout ends a run after `max(500, min(3000, floor(timeout * 0.2)))`
    ms of console silence and resets that timer on every printed line.
    The browser suites are silent while polling, so a stalled test goes
    quiet, auto-finish fires, and every later test is dropped --
    including the entire minified pass. Measured at tapout's default
    timeout: zero `# min >` headers, no plan line, **and the run
    reported as a PASS**. `test:browser` therefore passes
    `--timeout 30000` and `withTimeout` defaults to 2000 ms with a
    500 ms margin. `test/browser/harness.ts` asserts the relationship
    and fails if either value drifts; it accepts `-t` as well as
    `--timeout` and takes the last occurrence, matching tapout.

14. **Anything comparing compressed bytes for equality must pin
    `mtime`.** `src/index.ts` has exactly two wall-clock reads: line
    1199 writes `Date.now() / 1000` into gzip header bytes 4-7 at
    1-second resolution, and line 2917 writes the DOS time into every
    ZIP header at 2-second resolution. Measured unpinned: gzip 4
    mismatches in 200 at offset 4; zip 1 in 60 at offsets 10 and 650,
    rising to 29 in 60 when the calls are aligned to a tick. Pinned:
    0 in 200 and 0 in 60. zlib and raw deflate carry no timestamp and
    need no pinning.

15. **Do not generate incompressible test data with an LCG.**
    `state = (state * 1103515245 + 12345) & 0x7fffffff` degenerates in
    JS: the product exceeds `2**53`, so the low bits are lost to float
    rounding before the mask applies. It deflated 600000 bytes to
    19536, which silently sent `AsyncUnzipInflate` down its
    synchronous branch. `test/browser/util.ts` uses xorshift128 in
    exact 32-bit arithmetic. Re-measure with `deflateSync` if you
    change a fixture.

16. **The async paths want fixtures with opposite properties.**

    | consumer | gate | fixture |
    |---|---|---|
    | `zip()` | `originalSize >= 160000` | `largeCompressible()` |
    | `unzip()` | `originalSize >= 524288` AND `compressed <= 0.8 * original` | `largeCompressible()` |
    | `AsyncUnzipInflate` | COMPRESSED size, passed as `sz`, `>= 320000` | `largeIncompressible()` |

    Each test asserts its own gate, so a fixture that drifts fails
    loudly instead of silently downgrading to the sync path.

17. **Run only ONE npm invocation at a time against this repo.** They
    share `dist/`. A concurrent run produced an `npm test` exit 1 that
    looked like a real failure and was not.

## Watch for vacuous assertions

This is the defect this plan produced most. Seven instances across six
review cycles, one of them introduced by a fix for the same class.

The shapes seen so far:

- `t.ok(x, ...)` placed after a poll that only resolves when `x` is
  truthy. It can never fail.
- A flag set inside the `Unzip` `onfile` callback used to prove the
  async decoder ran. That callback fires for the sync decoder too.
- An error captured into a variable but asserted after an `await` that
  only resolves on success, so the assertion is skipped in exactly the
  case it was written for. TypeScript will narrow the variable to
  `null` and type the truthy branch `never`, which is the tell.
- A test that names an async path but feeds it a fixture below the
  gate, so it runs entirely on the main thread and would pass
  identically on the sync fallback.

The working pattern, used in thirteen places in `test/browser/`: a
`settled` flag that BOTH the success and the error branch set, polled
instead of the payload, then assert the error is null and the content
matches. **When you add a guard, prove it bites** by breaking the thing
it guards and confirming `not ok` in both the plain and min passes.

## Process notes

The plan documents have been amended repeatedly wherever the spec was
wrong. **Keep doing that** -- if you change prescribed code, change the
phase file's block to match, or the next run reintroduces the defect.
Amendments so far: `phase_04.md` (the `test:browser` timeout pairing),
`phase_05.md` (AC5.2, the error-marshalling rationale, Task 1's
line-number anchors, the carried-forward note), `phase_06.md` (AC3.4,
Task 6's entry-point block, Task 2's exports list, the fixture table,
the LCG trap, the gzip mtime rule, the anti-vacuous-assertion rules,
Task 3's timeout) and `phase_07.md` (Task 4's expected dist listing,
Task 3's verification commands).

Subagent output on this plan has needed heavy verification.
**Verify claims by measurement before accepting them**, especially any
statement of the form "X cannot work" or "X is now covered". Concrete
failures seen: an agent added regex post-processing to
`scripts/build.ts` on an unverified premise and committed it with a
message claiming a fix that never worked; another dropped all async
streaming coverage and reported it as a "design decision"; a third
skipped the one measurement it was asked for; a fourth replaced one
vacuous assertion with a different vacuous assertion and deleted the
error check while reporting the issue fixed.

The reviewers, by contrast, were consistently good once given the
measured record and told not to re-derive it. Two of the three most
valuable findings in this phase came from a reviewer probing rather
than reading.

## Phase history

| Phase | Tasks | Review cycles | Outcome |
|---|---|---|---|
| 1 Dependencies and example/ migration | 8 | 4 | approved |
| 2 tsconfig and node-worker static import | 4 | 2 | approved |
| 3 The build pipeline | 8 | 4 | approved |
| 4 package.json entry points and scripts | 3 | 4 | approved |
| 5 Node test suite migration | 3 | 5 | approved |
| 6 Browser test suite | 6 | 6 | approved |
| 7 Legacy removal and final verification | 5 | 0 | executed, NOT reviewed |

## What this plan actually bought

The build restructure was the stated goal, but the tests it forced
were worth more.

Phase 6 caught that eight exported classes -- `Gzip`, `Gunzip`,
`Zlib`, `Unzlib`, `Compress`, `AsyncDecompress`, `ZipDeflate`,
`AsyncZipDeflate` -- threw on construction in the shipped package,
because Phase 3 moved from `tsc` at `target: es5` to esbuild at
`es2022` and esbuild cannot lower classes to ES5 by any route. That
bug would have shipped.

The review cycles on Phase 6 then found that the suite meant to guard
the minified bundle could report a clean pass while never running the
minified pass at all, and that several tests claiming to exercise
worker paths never left the main thread. Both were silent. Neither
would have been found by reading the code.
