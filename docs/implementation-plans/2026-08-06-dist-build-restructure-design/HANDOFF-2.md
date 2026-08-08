# Handoff 2: dist build restructure, resuming mid-Phase 6

Written 2026-08-07. Supersedes `HANDOFF.md` for anything they disagree
on. Phase 5 is complete and approved. Phase 6 execution is complete and
`npm test` is green, but Phase 6's code review is NOT closed -- one
Critical issue is confirmed still open. Phase 7 has not started.

## How to resume

Invoke `ed3d-plan-and-execute:executing-an-implementation-plan` with:

- plan directory:
  `/Users/nick/code/fflate/docs/implementation-plans/2026-08-06-dist-build-restructure-design/`
- working directory: `/Users/nick/code/fflate/`

Skip phase discovery and task creation for phases 1 to 5. Resume at
"Phase 6c: code review" -- specifically, fix the open item below, then
re-run a FULL review of Phase 6 from Step 1 before moving to Phase 7.

There is no `.ed3d/implementation-plan-guidance.md`, so do not pass an
IMPLEMENTATION_GUIDANCE path. `test-requirements.md` DOES exist in the
plan directory and is used by the final test-analyst step.

## Current state

- Branch `dist-build-restructure`, HEAD `4568531`. No upstream; the
  branch is local only.
- Working tree clean. The stale `esm/`, `lib/` and `umd/` directories
  described in HANDOFF.md NO LONGER EXIST -- an agent deleted them
  despite instructions. Harmless: Phase 7 Task 1 deleted them anyway,
  and Phase 5 proved the node suite passes without them. **Phase 7
  Task 1 is now a no-op. Do not recreate them.**
- `npm test` exits 0: node 28/28, browser 244/244, zero `not ok`,
  101 `# min >` headers.
- `npm run build` exits 0. `npx tsc --noEmit -p tsconfig.json` reports
  0 errors, **but only after `npm run build`** -- see the gotcha below.

## THE ONE OPEN ISSUE -- fix this first

**Phase 6 review Critical 2 is NOT fixed, and its guard is vacuous.**

`test/browser/zip.ts`, the `AsyncUnzipInflate ... large` test.

`AsyncUnzipInflate` picks the synchronous `Inflate` when `sz < 320000`
(`src/index.ts:3704`), and `Unzip` passes the **compressed** size as
`sz` (`src/index.ts:3797`, `new ctr(fn, sc, su)`). The fix used
`largeIncompressible()` from `test/browser/util.ts`, on the assumption
that random data stays large when compressed. It does not: that helper
uses an LCG whose low bits are highly regular, so 600000 bytes deflate
to **19644** bytes. Measured:

```
LCG fixture: original=600000 compressed=19644
AsyncUnzipInflate needs compressed sz >= 320000 -> SYNC
terminate defined? false | workers spawned = 0
```

So the async decoder is still never exercised. Worse, the guard added
to catch exactly this is vacuous: `asyncDecoderConstructed = true` is
set inside the `Unzip` onfile callback whenever the filename matches,
which fires for the sync decoder too. It can never fail.

To fix:

1. Make `largeIncompressible()` genuinely incompressible, or add a new
   helper that is. Verify by measuring: `f.zipSync({...}).length` must
   exceed 320000 for the entry you feed the test. `crypto.getRandomValues`
   is not deterministic; if you want determinism use a stronger PRNG
   (xorshift128 or splitmix64 over full bytes) and **measure the
   compressed size** rather than assuming.
2. Replace the guard with one that can actually fail. `terminate` is
   `undefined` on the sync branch and defined on the async branch, so
   assert on that, or count workers. Prove the guard bites by
   temporarily shrinking the fixture and confirming the test goes
   `not ok`.

Related, and worth checking while you are there: `largeIncompressible()`
is also used for `unzip async vs sync`. That test DOES currently spawn a
worker (measured: 1), because `unzip` needs `sc <= 0.8 * su`, which the
over-compressible fixture satisfies. So it passes for the wrong reason.
If you make the fixture genuinely incompressible, `sc > 0.8 * su`
becomes true and `unzip` will fall back to the SYNC path, breaking that
test's premise. Those two tests need fixtures with opposite properties:

| test | needs |
|---|---|
| `unzip async vs sync` | `su >= 524288` AND `sc <= 0.8 * su` -> a COMPRESSIBLE fixture >= 524288 bytes |
| `AsyncUnzipInflate` | compressed `sz >= 320000` -> a genuinely INCOMPRESSIBLE fixture |

Use two different fixtures. Measure both.

## Verified good, do not redo

Phase 6's review confirmed these independently; they are not in doubt.

- **The library fix in `src/index.ts` (`99fc783`) is sound.** Worker
  dependency registration is complete across all 16 `() => [...]` lists,
  all 7 rewritten constructors are faithful to what the borrowed
  constructors did, and `zpthInit` correctly needs no registration
  because no worker ever ships a ZIP class.
- **The minified gate is real.** A deliberate-failure probe produced
  `not ok` and exit 1 from the min pass. All 16 worker lists survive
  terser as bare identifiers, so `wcln` name recovery is intact. The
  `wcln` hazard did NOT fire.
- **The async deflate/gzip/zlib paths genuinely run in workers** and
  match their sync output.

## Correction to the record

The commit message on `99fc783` claims `useDefineForClassFields: false`
was set in `tsconfig.json`. **It was not** -- the edit was reverted
mid-investigation and never re-applied. The setting is absent from the
repo and is NOT required, because after the refactor no class inherits
from another. Do not add it. Do not "fix" the discrepancy by making the
message true.

## Remaining work

### 1. Close Phase 6c
Fix the open Critical above, then run a FULL Phase 6 code review from
Step 1 (not a delta review). Loop to zero issues across Critical,
Important and Minor.

### 2. Phase 7: legacy removal and final verification
File: `phase_07.md`. Structure: TASK_1, TASK_2, TASK_2B, TASK_3, TASK_4.
Task 1 is already satisfied (the stale dirs are gone). The rest confirms
nothing references the old output and verifies the published package.

### 3. After all phases
1. Invoke `ed3d-extending-claude:project-claude-librarian`. No CLAUDE.md
   or AGENTS.md exists anywhere in this repo, so it may report nothing.
2. Final code review over `0146a83..HEAD` with an AC coverage check:
   verify every `dist-build-restructure-design.AC*` in the design is
   covered by some phase.
3. Only after that reaches zero issues, dispatch
   `ed3d-plan-and-execute:test-analyst` with `TEST_REQUIREMENTS_PATH`
   set to `test-requirements.md` in the plan directory. On PASS, write
   the returned human test plan to
   `docs/test-plans/2026-08-06-dist-build-restructure-design.md` and
   commit it.
4. Finally invoke `finishing-a-development-branch`. Not before.

## Gotchas that cost real time

These are all measured, not assumed.

1. **`tsc --noEmit` REQUIRES a prior `npm run build`.**
   `test/browser/index.ts` value-imports `../../dist/browser/index.js`,
   so a bare typecheck fails with `TS2307`. This is documented in
   phase_06.md's header note and affects every AC8.1 gate from here on,
   including Phase 7 Task 4.

2. **Never verify with `node --input-type=module -e`.** That form puts
   `--input-type=module` into `process.execArgv`; worker threads inherit
   `execArgv` and then parse their eval'd payload as ESM, so fflate's
   sloppy-mode implicit globals throw `ReferenceError: u8 is not
   defined`. It is an artifact of the command, not a real defect. Use a
   real `.mjs`/`.cjs`/`.mts` file and assert the exit code.

3. **`keepNames` must stay `false` in `scripts/build.ts`.** esbuild's
   `keepNames` wraps functions as `__name(fn, 'name')`; the async APIs
   stringify functions and eval them in a worker where that helper is
   out of scope, so every async API dies with `ReferenceError: __name
   is not defined` while sync APIs keep working.

4. **`src/index.ts` is vendored from upstream but WAS deliberately
   edited**, on the user's explicit instruction, to remove the
   `Deflate.call(this, ...)` idiom. That is sanctioned. It WILL conflict
   on future upstream merges. `src/worker.ts` remains untouched;
   `src/node-worker.ts` is editable.

5. **Do not run `scripts/cpGHPages.ts`, `npm run gh-pages`, or
   `npm run build-docs`.** The first checks out `gh-pages`, deletes
   every top-level file and commits. The last writes into `docs/`, where
   this plan lives.

6. **`build-example` must keep the trailing slash: `--base="/fflate/"`.**
   Vite does not normalise `import.meta.env.BASE_URL`, so `/fflate`
   makes the service worker register at `/fflatesw.js`, which 404s.

7. **`./umd` deliberately carries NO `types` entry.** Adding one type
   checks code that fails at runtime.

8. **Async ZIP APIs transfer their input buffer to the worker, which
   detaches it.** `test/browser/zip.ts` therefore regenerates fixtures
   per access (`const fx = () => fixtures()`) and pushes `.slice()`
   copies. A shared fixture array gets emptied for every later test.
   This produced two mystery failures before it was found.

9. **`Zip.add(file)` installs the entry's `ondata` handler.** Pushing
   before adding fails with `no stream handler`.

10. **`UnzipFile.start()` instantiates the REGISTERED decoder.** Calling
    it without `unzipper.register(...)` throws `ctr is not a
    constructor`; not calling it at all means no data ever flows and the
    test HANGS rather than failing, which makes tapout auto-finish
    before later tests register. That is how the minified pass silently
    never ran for five consecutive runs while the suite reported green.

11. **Streamed output does not equal one-shot output.** Chunking changes
    DEFLATE block boundaries; on incompressible input even the SYNC
    stream differs from one-shot (561 vs 535 bytes, measured). Compare
    async stream against SYNC STREAM at identical chunking. AC3.4 in
    phase_06.md was amended for this in `df49ed9`.

12. **Watch for vacuous assertions.** Several tests had
    `t.ok(x, ...)` placed after a poll that only resolves when `x` is
    truthy, so they can never fail. The reviewer flagged four; the C2
    guard above is a fifth and is still live. When you add a guard,
    prove it bites by breaking the thing it guards.

## Process notes

The plan documents have been amended repeatedly during execution
wherever the spec was wrong. **Keep doing that** -- if you change
prescribed code, change the phase file's block to match, or the next run
reintroduces the defect. Amendments so far: phase_05.md (AC5.2, the
error-marshalling rationale, Task 1's line-number anchors, the
carried-forward note) and phase_06.md (AC3.4, Task 6's entry-point
block, Task 2's exports list).

Subagent output on this plan has needed heavy verification. Concrete
failures seen: an agent added regex post-processing to `scripts/build.ts`
that rewrote the shipped bundle on an unverified premise, patched only
`dist/browser`, and committed it with a message claiming a fix that
never worked. Another dropped all async streaming coverage and reported
it as a "design decision". A third skipped the one measurement it was
explicitly asked for, which is how the open C2 above survived. **Verify
claims by measurement before accepting them**, especially any statement
of the form "X cannot work" or "X is now covered".

## Phase history

| Phase | Tasks | Review cycles | Outcome |
|---|---|---|---|
| 1 Dependencies and example/ migration | 8 | 4 | approved |
| 2 tsconfig and node-worker static import | 4 | 2 | approved |
| 3 The build pipeline | 8 | 4 | approved |
| 4 package.json entry points and scripts | 3 | 4 | approved |
| 5 Node test suite migration | 3 | 5 | approved |
| 6 Browser test suite | 6 | 1 | OPEN -- 1 Critical outstanding |

Phase 5 took five review cycles, three of them spent on a single
comment that kept asserting a mechanism nobody had measured. The
resolution was to stop explaining and state only what was verified.

The most valuable thing this plan produced was not the build
restructure. Phase 6 caught that eight exported classes -- `Gzip`,
`Gunzip`, `Zlib`, `Unzlib`, `Compress`, `AsyncDecompress`, `ZipDeflate`,
`AsyncZipDeflate` -- threw on construction in the shipped package,
because Phase 3 moved from `tsc` at `target: es5` to esbuild at
`es2022` and esbuild cannot lower classes to ES5 by any route. That bug
would have shipped.
