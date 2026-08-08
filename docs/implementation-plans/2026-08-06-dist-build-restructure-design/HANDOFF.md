# Handoff: dist build restructure execution

Written 2026-08-06 at a clean stopping point. Phases 1 through 4 are
complete and code-review approved. Phases 5, 6 and 7 remain, followed by
the plan's final review sequence.

## How to resume

Invoke the `ed3d-plan-and-execute:executing-an-implementation-plan`
skill with:

- plan directory:
  `/Users/nick/code/fflate/docs/implementation-plans/2026-08-06-dist-build-restructure-design/`
- working directory: `/Users/nick/code/fflate/`

Skip its phase-discovery and task-creation steps for phases 1 to 4;
resume at "Phase 5a: read the phase file". The skill's rules still
apply: read one phase at a time, execute every task in it, then run ONE
code review for the phase and loop until zero issues. Fix every issue
including Minor ones.

There is no `.ed3d/implementation-plan-guidance.md` in this repo, so do
not pass an IMPLEMENTATION_GUIDANCE path to the reviewer.
`test-requirements.md` DOES exist in the plan directory and is used by
the final test-analyst step.

Session scratchpad used so far (any new session should use its own):
`/tmp/exec-2026-08-06-dist-build-restructure-design-42af18d4`

## Current state

- Branch `dist-build-restructure`, HEAD `f49c06a`, 36 commits ahead of
  the starting commit `0146a83`.
- Working tree clean apart from three untracked stale directories,
  `esm/`, `lib/` and `umd/`. These are output of the OLD pipeline. Phase
  7 deletes them. Do not delete them early -- see the warning below.
- `npm run build` exits 0 and emits exactly TWELVE files:
  `dist/browser/{index.d.ts,index.js,index.js.map,index.min.js,index.min.js.map}`,
  the same five under `dist/node/`, plus `dist/umd/fflate.js` and
  `dist/umd/package.json`.
- `npx tsc --noEmit -p tsconfig.json` reports 0 errors.
- `npm run test:node` passes 25/25 (but see the warning below).
- `npm pack --dry-run` lists 15 files: `dist/` plus `package.json`,
  `README.md`, `LICENSE`.
- `npm test` as a whole FAILS until Phase 6, because `test:browser`
  needs `test/browser/index.ts` and the `tapout` binary. Expected.

## Read this before touching anything

These cost real debugging time. Several are defects in the plan itself,
not in the code.

1. **`npm run test:node` is green only because the stale `lib/` still
   exists.** `test/util.ts:167` builds a CommonJS path from separate
   segments, `resolve(here, '..', 'lib', 'index.cjs')`, so a grep for
   `lib/` does NOT find it. Move `lib/` aside and the suite aborts with
   `Cannot find module '.../lib/index.cjs' at [worker eval]:2:29`. It
   would fail on a clean checkout and in CI. **This is Phase 5's job.**
   The replacement must be CommonJS loadable by `require()` from a
   `worker_threads` eval string, so `dist/node/index.js` will NOT serve
   -- it is ESM inside a `"type": "module"` package. This is recorded in
   phase_05.md under "Carried forward from the Phase 4 review". Verify
   any fix with `lib/` moved aside.

2. **Never verify the async path with `node --input-type=module -e`.**
   That form puts `--input-type=module` into `process.execArgv`, worker
   threads inherit `execArgv`, and the worker then parses its eval'd
   payload as ESM, so fflate's sloppy-mode implicit globals throw
   `ReferenceError: u8 is not defined`. It is an artifact of the command,
   not the bundle. Verify from a real `.mjs` file, attach a `.catch`, and
   assert the exit code -- an unhandled rejection reads like a pass if
   you only look at stdout. Recorded as phase_03.md investigation
   finding 8.

3. **`keepNames` must stay `false` in `scripts/build.ts`.** esbuild's
   `keepNames` wraps functions as `__name(fn, 'name')`; the async APIs
   stringify functions and eval them in a worker where that helper is
   out of scope, so `deflate()` dies with
   `ReferenceError: __name is not defined` while every sync API keeps
   working. The plan originally specified `true`; all three blocks have
   been amended. Recorded as phase_03.md investigation finding 7.

4. **`build-example` must keep the trailing slash:
   `--base="/fflate/"`.** Vite does not normalise
   `import.meta.env.BASE_URL`, so `/fflate` makes the service worker
   register at `/fflatesw.js`, which 404s and precaches nothing. Asset
   URLs in `index.html` look correct either way, which is what makes it
   easy to miss. Phase 4's task block originally restated it without the
   slash; amended, recorded as phase_04.md finding 7.

5. **`./umd` deliberately carries NO `types` entry.** Adding one type
   checks code that fails at runtime. Recorded as phase_04.md finding 5.

6. **`src/index.ts` and `src/worker.ts` are vendored upstream** and this
   repo merges from `git remote upstream`. Do not edit them.
   `src/node-worker.ts` IS editable and was legitimately changed in
   Phase 2.

7. **Do not run `scripts/cpGHPages.ts` or `npm run gh-pages`.** It
   checks out `gh-pages`, deletes every top-level file, and commits. No
   `gh-pages` branch exists in this clone, so it cannot succeed here
   anyway. Verify it statically only.

8. **`npm run build-docs` writes `--out docs/`**, which is where this
   plan lives. `--cleanOutputDir false` has been added to stop it
   deleting the plan, but still avoid running it casually.

9. The plan documents have been amended during execution wherever the
   spec was wrong. Keep doing that -- if you change prescribed code,
   change the phase file's block to match, or the next run reintroduces
   the defect.

## Remaining work

### Phase 5: Node test suite migration
File: `phase_05.md`. Structure: SUBCOMPONENT_A (tasks 1-2), TASK_3.
Goal: repoint the node benchmark harness at `dist/node/index.js` and move
three empty stub suites out of the node runner. Start by reading the
"Carried forward from the Phase 4 review" section at the top.

### Phase 6: Browser test suite
File: `phase_06.md`. Structure: SUBCOMPONENT_A (tasks 1-2),
SUBCOMPONENT_B (tasks 3-5), TASK_6.
Goal: write ZIP, streaming and async coverage that has never existed and
run it in a real browser against both `dist/browser/index.js` and
`dist/browser/index.min.js`.

**This phase is the real gate for the `wcln` minification hazard.** The
async APIs recover mangled names by parsing array-literal text; any
terser option that substitutes a constant's value for its name breaks
this silently, leaving sync APIs working while async APIs corrupt output.
Encouraging early signal: the minified NODE bundle already round trips
sync and async correctly. The browser path is still unproven.

Phase 6 also introduces `test/browser/index.ts` and the `tapout` binary,
which is what finally makes `npm test` pass end to end.

### Phase 7: Legacy removal and final verification
File: `phase_07.md`. Structure: TASK_1, TASK_2, TASK_2B, TASK_3, TASK_4.
Goal: delete the stale `lib/`, `esm/`, `umd/` output, confirm nothing
references it, and verify the published package end to end. Phase 5 must
have fixed `test/util.ts` first or the node suite breaks here.

### After all phases

1. Invoke the `ed3d-extending-claude:project-claude-librarian` subagent
   to update CLAUDE.md files. No CLAUDE.md or AGENTS.md exists anywhere
   in this repo today, so it may report nothing to do.
2. Final code review over the whole range, base `0146a83` to HEAD, with
   an AC coverage check: verify every acceptance criterion in the design
   (scoped `dist-build-restructure-design.AC*`) is covered by some phase.
3. Then, and only after that review reaches zero issues, dispatch
   `ed3d-plan-and-execute:test-analyst` with
   `TEST_REQUIREMENTS_PATH` set to `test-requirements.md` in the plan
   directory. If it returns PASS it also returns a human test plan; write
   that to `docs/test-plans/2026-08-06-dist-build-restructure-design.md`
   and commit it.
4. Finally invoke `finishing-a-development-branch`. Not before.

## Phases already done

| Phase | Tasks | Review cycles | Outcome |
|---|---|---|---|
| 1 Dependencies and example/ migration | 8 | 4 | approved |
| 2 tsconfig and node-worker static import | 4 | 2 | approved |
| 3 The build pipeline | 8 | 4 | approved |
| 4 package.json entry points and scripts | 3 | 4 | approved |

No compromises were made and no review issue was left outstanding in any
of them. Every phase ended at zero issues across Critical, Important and
Minor.

Notable defects caught, all fixed and verified:

- `cpGHPages.ts` would have checked out `gh-pages`, deleted every root
  file, then thrown on the first `copyFileSync` because Vite's output
  nests an `assets/` directory -- leaving the user stranded on the wrong
  branch with the repo emptied.
- The service worker could never register on the deployed site, twice,
  for two different reasons.
- `git.add(['.'])` would have committed `node_modules/` to `gh-pages`,
  because the branch switch removes `.gitignore` from the working tree.
- The build's `keepNames: true` broke every async API while leaving sync
  APIs working.
- The `./umd` export was unloadable in both module systems.
- `build-docs` would have deleted this plan.

Two of these came from the plan document rather than from an
implementation mistake, which is why amending the plan as you go matters.
