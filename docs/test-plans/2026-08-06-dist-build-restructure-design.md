# Human Test Plan: dist build restructure

Generated 2026-08-07 from
`docs/implementation-plans/2026-08-06-dist-build-restructure-design/test-requirements.md`,
after the test-coverage analysis returned PASS on all 29 automated
acceptance criteria.

Almost everything in this plan is automated. What follows is the
remainder: two criteria that cannot be automated without either a served
origin or a destructive branch operation, plus two end-to-end scenarios
that close documented gaps in the automation.

## Prerequisites

1. Node 25.x. Everything below was measured on 25.8.2.
2. Clean working tree on `dist-build-restructure`.
3. `npm ci`, then `npx playwright install chromium`.
4. `npm test` passing: exit 0, node 34/34, browser 306/306.
5. Chrome with DevTools, and a static server such as `npx serve public`
   or `python3 -m http.server`.

**Do not run `npm run gh-pages`, `scripts/cpGHPages.ts`, or
`npm run build-docs` during sections 1, 3, or 4.** The first two check
out `gh-pages`, delete every top-level file, and commit. The last writes
into `docs/`, where this plan and the implementation plan both live.

Run only one npm invocation at a time. They share `dist/`, and a
concurrent run produces a failure that looks real and is not.

## 1. Service worker and offline precache (AC7.3)

The genuinely outstanding item. The browser suite drives a bundled
module rather than a served origin with a registered worker, so no
automated check covers registration, the cache contents, or the offline
transition.

Run `npm run build-example` first. It should report `precache 4 entries`
and write `../public/sw.js`.

| Step | Action | Expected |
| --- | --- | --- |
| 1 | Serve `public/` on `http://localhost:3000` and open it in Chrome | Page renders. `localhost` is a secure context, so no TLS is needed |
| 2 | DevTools -> Application -> Service Workers | Exactly one worker, status `activated`, script URL ending `/sw.js` |
| 3 | Application -> Cache Storage | A cache whose name begins `fflate-`, containing `index.html` plus the three hashed `assets/` entries (`index-*.css`, `index-*.js`, `workers-*.js`) |
| 4 | Look for `favicon.ico` in that cache | Absent. See the note below -- the reason is not what an earlier version of this plan said it was |
| 5 | Note the cache name. Set Network throttling to Offline and reload | Page still renders fully. This is the assertion that the precache is real rather than merely populated |
| 6 | Back online. Re-run `npm run build-example` with no source change, reload | Cache name unchanged |
| 7 | Edit any file under `example/`, re-run `npm run build-example`, reload | Cache name changed |

**On step 4.** `favicon.ico` is excluded by vite-plugin-pwa's default
`injectManifest` glob, not by anything in this repository.
`example/sw.ts:26` is a bare `manifest.map(e => e.url)`. The old
`demo/sw.ts` did filter `.ico` explicitly, and the rationale outlived
the filter. If `favicon.ico` ever appears here, the fix is a
`globIgnores` entry or a new filter in `example/sw.ts`, not a hunt for
a filter that no longer exists.

**On steps 6 and 7.** These are the part most likely to be wrong. They
exercise the property that Parcel's removed `version` export used to
provide. If step 5 or step 7 fails, look at the `precacheVersion`
derivation in `example/sw.ts`, not at the Vite config.

Revert the step 7 edit afterwards and confirm `git status --porcelain`
is empty.

## 2. First gh-pages deploy (AC7.4)

Not runnable before merge. `scripts/cpGHPages.ts` checks out `gh-pages`,
deletes every top-level file, and commits. Verify on the first real run.

| Step | Action | Expected |
| --- | --- | --- |
| 1 | Note the current branch name | Needed for step 4 |
| 2 | `npm run build-example`, then `npm run gh-pages` | Exits 0 |
| 3 | `git log -1 --stat gh-pages` | The commit carries `index.html`, `assets/`, `favicon.ico`, and `sw.js` from `public/` -- not from `dist/` |
| 4 | `git branch --show-current` | The branch from step 1, **not** `master` |
| 5 | Load the published GitHub Pages URL | The example app loads under the `/fflate/` base path |

Step 4 is the one to watch. Phase 1 Task 7 replaced an unconditional
`git checkout('master')` with a return to the starting branch, and no
automated check exercises it. That task also changed `statSync(f)` to
`statSync(to(f))`; neither fix is covered by any command.

The trailing slash in `--base="/fflate/"` is load-bearing for step 5.
Vite does not normalise `import.meta.env.BASE_URL`, so `/fflate` would
register the service worker at `/fflatesw.js`, which 404s.

## 3. UMD artifact in a real browser

The build verifies `dist/umd/fflate.js` through a
`new Function('self', 'module', 'exports', 'define', s)` harness, never
in a browser. The design assigns this file the script-tag and `unpkg`
role, so confirm it in the medium it ships for.

Create a scratch directory containing a copy of `dist/umd/fflate.js` and
an HTML file:

```html
<script src="fflate.js"></script>
<script>
  const d = fflate.deflateSync(new TextEncoder().encode('umd works'))
  console.log(new TextDecoder().decode(fflate.inflateSync(d)))
</script>
```

Serve it over `http://localhost` and open the console.

Expected: `umd works`, and `typeof fflate === 'object'` on `window`.

**Do not test this file with Node `require()`.** The package is
`"type": "module"`, so Node parses it as ESM and the check proves
nothing. The sibling `dist/umd/package.json` containing
`{ "type": "commonjs" }` is what makes the shipped artifact loadable;
with it removed, both `require()` and `import()` throw
`TypeError: Cannot set properties of undefined (setting 'fflate')`.

## 4. Bundler resolution of the `default` condition

This closes the plan's one acknowledged automation gap, in AC4.2. No
Node process can resolve the bare `.` entry under a non-node condition,
because Node's resolver always supplies `node`. The browser artifact
itself is well covered by the browser suite; the `exports` wiring for
that branch is not.

In a scratch directory:

1. `npm init -y`
2. `npm pack` in the repo, then install the resulting tarball here
3. Add `src/index.js` with
   `import { deflateSync } from '@substrate-system/fflate'` and a
   trivial call
4. `npx vite build`
5. Inspect the output bundle

Expected: the build succeeds and the bundled code contains
`URL.createObjectURL`, proving the browser variant was selected. A
bundle containing `worker_threads` means the `default` branch is
mis-wired.

The consumer-facing half of AC4 is already verified automatically
against a packed and installed tarball: sync and async round-trips from
a real dependency, all four subpaths resolving, both `.d.ts` files
shipping, only `dist/` in the installed root, and no `require`
condition in the packed `exports`. This step is narrowly about the
`default` branch under a non-node condition.

## What is already covered automatically

Recorded so nobody re-does it by hand. All 29 sub-criteria passed the
coverage analysis; these are the ones a human might otherwise be
tempted to spot-check.

| Area | Automated gate |
| --- | --- |
| dist layout | `npm run build` then `find dist -type f` -- exactly 12 paths; terser confirmed by size (browser 66753 -> 31221, node 67030 -> 31409) |
| Worker swap per target | `createObjectURL` count 1 in browser and 0 in node; `worker_threads` count 0 in browser and umd |
| `wcln` hazard | `test/browser/` runs all three suites twice, plain and minified, 306 assertions. `test/3-node-min.ts` covers `dist/node/index.min.js`, which the browser suite cannot reach |
| Published entry points | packed tarball installed into a scratch project: round-trips, subpath resolution, file list, no `require` condition |
| Node suite | `npm run test:node`, 34/34 |
| Type checking | `npx tsc --noEmit -p tsconfig.json`, 0 errors -- **after** a build, since `test/browser/index.ts` value-imports `dist/browser/index.js` |
| Legacy removal | six `git grep` sweeps, all silent; `ls scripts` is exactly `build.ts` and `cpGHPages.ts` |

## Reporting

Record the section 1 result in the PR description. Sections 2 through 4
can be reported wherever the branch is being tracked; section 2 cannot
run until after merge.
